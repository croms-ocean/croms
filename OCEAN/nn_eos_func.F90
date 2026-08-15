!======================================================================
! CROCO Neural Network Equation of State
!======================================================================
!
!  Replaces the nonlinear polynomial EOS (Shchepetkin 2003) with a
!  per-dt trained MLP.  (BASIC))
!======================================================================

#include "cppdefs.h"

MODULE nn_eos_func_mod

  IMPLICIT NONE
  PRIVATE

#ifdef MPI
  include 'mpif.h'
#endif

  !--------------------------------------------------------------------
  ! Tuneable parameters
  !--------------------------------------------------------------------
  INTEGER,  PARAMETER :: NN_H      = 32
  INTEGER,  PARAMETER :: N_BATCH   = 2000
  REAL(8),  PARAMETER :: LR        = 8.D-3
  REAL(8),  PARAMETER :: W_RHO     = 0.3D0   ! weight on rho loss ***
  REAL(8),  PARAMETER :: W_GRAD    = 0.7D0   ! weight on gradient loss **
  REAL(8),  PARAMETER :: ADAM_B1   = 0.9D0
  REAL(8),  PARAMETER :: ADAM_B2   = 0.999D0
  REAL(8),  PARAMETER :: ADAM_EPS  = 1.D-8

  !--------------------------------------------------------------------
  ! Reference EOS polynomial coefficients (Shchepetkin 2003)
  !--------------------------------------------------------------------
  REAL(8), PARAMETER :: r00=999.842594D0,   r01=6.793952D-2, r02=-9.095290D-3
  REAL(8), PARAMETER :: r03=1.001685D-4,    r04=-1.120083D-6, r05=6.536332D-9
  REAL(8), PARAMETER :: r10=0.824493D0,     r11=-4.08990D-3,  r12=7.64380D-5
  REAL(8), PARAMETER :: r13=-8.24670D-7,    r14=5.38750D-9
  REAL(8), PARAMETER :: rS0=-5.72466D-3,    rS1=1.02270D-4,   rS2=-1.65460D-6
  REAL(8), PARAMETER :: r20=4.8314D-4
  REAL(8), PARAMETER :: rho0_ref=1025.D0    ! Boussinesq Ref.

  !--------------------------------------------------------------------
  ! Network weights: 2 -> NN_H -> NN_H -> 3
  !--------------------------------------------------------------------
  REAL(8), SAVE :: W1(NN_H,2),    b1(NN_H)
  REAL(8), SAVE :: W2(NN_H,NN_H), b2(NN_H)
  REAL(8), SAVE :: W3(3,NN_H),    b3(3)
  REAL(8), SAVE :: x_mean(2)=0.D0, x_std(2)=1.D0
  REAL(8), SAVE :: y_mean(3)=0.D0, y_std(3)=1.D0

  ! Adam moments
  REAL(8), SAVE :: mW1(NN_H,2),    vW1(NN_H,2)
  REAL(8), SAVE :: mW2(NN_H,NN_H), vW2(NN_H,NN_H)
  REAL(8), SAVE :: mW3(3,NN_H),    vW3(3,NN_H)
  REAL(8), SAVE :: mb1(NN_H), vb1(NN_H)
  REAL(8), SAVE :: mb2(NN_H), vb2(NN_H)
  REAL(8), SAVE :: mb3(3),    vb3(3)
  INTEGER, SAVE :: adam_t = 0

  ! Buffer
  REAL(4), ALLOCATABLE, SAVE :: buf_T(:), buf_S(:), buf_z(:)
  REAL(4), ALLOCATABLE, SAVE :: buf_rho(:), buf_dT(:), buf_dS(:)
  INTEGER, SAVE :: buf_n     = 0
  INTEGER, SAVE :: buf_max   = 0
  INTEGER, SAVE :: last_iic  = -9999
  INTEGER, SAVE :: nn_ready  = 0
  LOGICAL, SAVE :: buf_alloc = .FALSE.

  CHARACTER(LEN=256), SAVE :: wfile = ''

  PUBLIC :: nn_eos_init
  PUBLIC :: nn_eos_func
  PUBLIC :: nn_eos_finalize

CONTAINS

!=======================================================================
  SUBROUTINE nn_eos_init(run_title, my_rank)
!=======================================================================
    CHARACTER(LEN=*), INTENT(IN) :: run_title
    INTEGER,          INTENT(IN) :: my_rank

    INTEGER    :: ios
    INTEGER(8) :: rng
    LOGICAL    :: fexist
    CHARACTER(LEN=256) :: base

    ! Gen nn_weights_eos_<TITLE>.txt
    CALL make_fname(run_title, base)
    wfile = 'nn_eos_' // TRIM(base)

    IF (my_rank == 0) THEN
      WRITE(*,'(A,A)') ' [nn_eos] weight file: ', TRIM(wfile)
      INQUIRE(FILE=TRIM(wfile), EXIST=fexist)
      IF (fexist) THEN
        CALL load_weights(TRIM(wfile), ios)
        IF (ios == 0) THEN
          nn_ready = 1
          WRITE(*,'(A)') ' [nn_eos] weights loaded -> NN EOS active from step 1'
        ELSE
          rng = 99_8; CALL init_weights(rng)
          WRITE(*,'(A)') ' [nn_eos] bad file -> polynomial EOS used, training in bg'
        END IF
      ELSE
        rng = 99_8; CALL init_weights(rng)
        WRITE(*,'(A)') ' [nn_eos] no file -> polynomial EOS this run, NN next run'
      END IF
    END IF

!-------------------BLOCK-MPI----------------------
#ifdef MPI
    BLOCK
      INTEGER :: mpi_err
      CALL MPI_Bcast(nn_ready, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
      IF (nn_ready == 1) CALL bcast_weights()
    END BLOCK
#endif
!-------------------BLOCK-MPI----------------------

  END SUBROUTINE nn_eos_init


!=======================================================================
  SUBROUTINE nn_eos_func(Tt, Ts, depth, rho1_out, drho_dT, drho_dS, &
                          my_rank, timestep, Istr, Iend, Jstr, Jend, Nz)
!=======================================================================
!  Call per cell in rho_eos_tile.
!=======================================================================
    REAL,    INTENT(IN)  :: Tt, Ts, depth
    REAL,    INTENT(OUT) :: rho1_out, drho_dT, drho_dS
    INTEGER, INTENT(IN)  :: my_rank, timestep
    INTEGER, INTENT(IN)  :: Istr, Iend, Jstr, Jend, Nz

    REAL(8) :: h1(NN_H), h2(NN_H), h3(3), acc
    REAL(8) :: rho1_ref, dT_ref, dS_ref
    REAL(8) :: sqTs, Td, Sd
    INTEGER :: r, c

    ! Auto-allocate buffer
    IF (my_rank == 0 .AND. .NOT. buf_alloc) THEN
      buf_max = (Iend-Istr+1) * (Jend-Jstr+1) * Nz + 1000
      ALLOCATE(buf_T(buf_max), buf_S(buf_max), buf_z(buf_max))
      ALLOCATE(buf_rho(buf_max), buf_dT(buf_max), buf_dS(buf_max))
      buf_n = 0; buf_alloc = .TRUE.
      WRITE(*,'(A,I8)') ' [nn_eos] buffer auto-sized: ', buf_max
    END IF

    ! Timestep boundary (<prev)
    IF (my_rank == 0 .AND. timestep /= last_iic) THEN
      IF (buf_n >= 10) CALL do_train(buf_n)
      buf_n = 0; last_iic = timestep
    END IF

    ! Refpolynomial EOS
    Td   = DBLE(Tt)
    Sd   = MAX(DBLE(Ts), 0.D0)
    sqTs = SQRT(Sd)
    rho1_ref = (r00 - rho0_ref) &
             + Td*(r01+Td*(r02+Td*(r03+Td*(r04+Td*r05)))) &
             + Sd*(r10+Td*(r11+Td*(r12+Td*(r13+Td*r14))) &
                  + sqTs*(rS0+Td*(rS1+Td*rS2)) + Sd*r20)

    ! Ana gradients of polynomial (ext)
    dT_ref = r01 + Td*(2.D0*r02 + Td*(3.D0*r03 + Td*(4.D0*r04 + Td*5.D0*r05))) &
           + Sd*(r11 + Td*(2.D0*r12 + Td*(3.D0*r13 + Td*4.D0*r14)) &
                + sqTs*(rS1 + Td*2.D0*rS2))
    dS_ref = r10 + Td*(r11+Td*(r12+Td*(r13+Td*r14))) &
           + 1.5D0*sqTs*(rS0+Td*(rS1+Td*rS2)) + 2.D0*Sd*r20

    ! Buffer
    IF (my_rank == 0 .AND. buf_alloc .AND. buf_n < buf_max) THEN
      buf_n       = buf_n + 1
      buf_T(buf_n) = Tt
      buf_S(buf_n) = Ts
      buf_z(buf_n) = depth
      buf_rho(buf_n) = REAL(rho1_ref)
      buf_dT (buf_n) = REAL(dT_ref)
      buf_dS (buf_n) = REAL(dS_ref)
    END IF

    ! Return
    IF (nn_ready == 1) THEN
      h1(1) = (DBLE(Tt) - x_mean(1)) / (x_std(1) + 1.D-30)
      h1(2) = (DBLE(Ts) - x_mean(2)) / (x_std(2) + 1.D-30)

      !!!----------------------------------TODO(sm-loop)----------------------------!!!
      !!!------ depth not used: surface EOS labels have no pressure dependence -----!!!
      !!!---------------------------------------------------------------------------!!!

      DO r = 1, NN_H
        acc = b1(r)
        DO c = 1, 2; acc = acc + W1(r,c)*h1(c); END DO
        h2(r) = TANH(acc)
      END DO
      DO r = 1, NN_H
        acc = b2(r)
        DO c = 1, NN_H; acc = acc + W2(r,c)*h2(c); END DO
        h1(r) = TANH(acc)
      END DO
      DO r = 1, 3
        acc = b3(r)
        DO c = 1, NN_H; acc = acc + W3(r,c)*h1(c); END DO
        h3(r) = acc
      END DO

      rho1_out = REAL(h3(1)*y_std(1) + y_mean(1))
      drho_dT  = REAL(h3(2)*y_std(2) + y_mean(2))
      drho_dS  = REAL(h3(3)*y_std(3) + y_mean(3))
    ELSE
      rho1_out = REAL(rho1_ref)
      drho_dT  = REAL(dT_ref)
      drho_dS  = REAL(dS_ref)
    END IF

  END SUBROUTINE nn_eos_func


!=======================================================================
  SUBROUTINE nn_eos_finalize(my_rank)
!=======================================================================
    INTEGER, INTENT(IN) :: my_rank
    IF (my_rank /= 0) RETURN
    IF (buf_alloc .AND. buf_n >= 10) CALL do_train(buf_n)
    CALL save_weights(TRIM(wfile))
    WRITE(*,'(A,A)') ' [nn_eos] weights saved -> ', TRIM(wfile)
    WRITE(*,'(A,I8)') ' [nn_eos] total Adam steps: ', adam_t
    IF (ALLOCATED(buf_T))   DEALLOCATE(buf_T)
    IF (ALLOCATED(buf_S))   DEALLOCATE(buf_S)
    IF (ALLOCATED(buf_z))   DEALLOCATE(buf_z)
    IF (ALLOCATED(buf_rho)) DEALLOCATE(buf_rho)
    IF (ALLOCATED(buf_dT))  DEALLOCATE(buf_dT)
    IF (ALLOCATED(buf_dS))  DEALLOCATE(buf_dS)
  END SUBROUTINE nn_eos_finalize


!=======================================================================
!  PRIVATE ROUTINES
!=======================================================================

  SUBROUTINE do_train(Ntrain)
    INTEGER, INTENT(IN) :: Ntrain
    INTEGER  :: n, it
    INTEGER  :: idx(N_BATCH)
    REAL(8)  :: Xn1(N_BATCH), Xn2(N_BATCH), Xn3(N_BATCH)
    REAL(8)  :: Yn1(N_BATCH), Yn2(N_BATCH), Yn3(N_BATCH)
    REAL(8)  :: xm(3), xs(3), ym(3), ys(3)
    REAL(8)  :: h1o(NN_H), h2o(NN_H), h3o(3)
    REAL(8)  :: dh2(NN_H), dh1(NN_H)
    REAL(8)  :: dW1(NN_H,3),    db1g(NN_H)
    REAL(8)  :: dW2(NN_H,NN_H), db2g(NN_H)
    REAL(8)  :: dW3(3,NN_H),    db3g(3)
    REAL(8)  :: aW1(NN_H,3),    ab1g(NN_H)
    REAL(8)  :: aW2(NN_H,NN_H), ab2g(NN_H)
    REAL(8)  :: aW3(3,NN_H),    ab3g(3)
    REAL(8)  :: bc1, bc2, dl(3)
    INTEGER(8) :: s
    INTEGER  :: Nuse

    ! Clamp -------------------- idx 
    Nuse = MIN(N_BATCH, Ntrain)

    s = INT(adam_t+1, 8) * 998244353_8
    CALL sample_idx(idx, Ntrain, Nuse, s)

    ! Normalise inputs
    xm(1)=SUM(DBLE(buf_T(idx(1:Nuse))))/Nuse
    xm(2)=SUM(DBLE(buf_S(idx(1:Nuse))))/Nuse
    xs(1)=MAX(SQRT(SUM((DBLE(buf_T(idx(1:Nuse)))-xm(1))**2)/Nuse),1.D-4)
    xs(2)=MAX(SQRT(SUM((DBLE(buf_S(idx(1:Nuse)))-xm(2))**2)/Nuse),1.D-4)

    ! Normalise outputs
    ym(1)=SUM(DBLE(buf_rho(idx(1:Nuse))))/Nuse
    ym(2)=SUM(DBLE(buf_dT (idx(1:Nuse))))/Nuse
    ym(3)=SUM(DBLE(buf_dS (idx(1:Nuse))))/Nuse
    ys(1)=MAX(SQRT(SUM((DBLE(buf_rho(idx(1:Nuse)))-ym(1))**2)/Nuse),1.D-4)
    ys(2)=MAX(SQRT(SUM((DBLE(buf_dT (idx(1:Nuse)))-ym(2))**2)/Nuse),1.D-4)
    ys(3)=MAX(SQRT(SUM((DBLE(buf_dS (idx(1:Nuse)))-ym(3))**2)/Nuse),1.D-4)

    DO n = 1, Nuse
      Xn1(n) = (DBLE(buf_T(idx(n))) - xm(1)) / xs(1)
      Xn2(n) = (DBLE(buf_S(idx(n))) - xm(2)) / xs(2)
      Xn3(n) = 0.D0   ! unused third slot (kept for bwd1e signature)
      Yn1(n) = (DBLE(buf_rho(idx(n))) - ym(1)) / ys(1)
      Yn2(n) = (DBLE(buf_dT (idx(n))) - ym(2)) / ys(2)
      Yn3(n) = (DBLE(buf_dS (idx(n))) - ym(3)) / ys(3)
    END DO
    x_mean(1:2)=xm(1:2); x_std(1:2)=xs(1:2); y_mean=ym; y_std=ys

    aW1=0.D0; ab1g=0.D0
    aW2=0.D0; ab2g=0.D0
    aW3=0.D0; ab3g=0.D0

    DO it = 1, Nuse
      CALL fwd3(Xn1(it),Xn2(it),Xn3(it), h1o,h2o,h3o)
      ! Weighted loss: gradient accuracy more important than rho value
      dl(1) = W_RHO  * 2.D0*(h3o(1)-Yn1(it)) / DBLE(Nuse)
      dl(2) = W_GRAD * 2.D0*(h3o(2)-Yn2(it)) / DBLE(Nuse)
      dl(3) = W_GRAD * 2.D0*(h3o(3)-Yn3(it)) / DBLE(Nuse)
      CALL bwd3e(h2o, dl, dW3, db3g, dh2)
      CALL bwd2e(h1o, h2o, dh2, dW2, db2g, dh1)
      CALL bwd1e(Xn1(it),Xn2(it),Xn3(it), h1o, dh1, dW1, db1g)
      aW1=aW1+dW1; ab1g=ab1g+db1g
      aW2=aW2+dW2; ab2g=ab2g+db2g
      aW3=aW3+dW3; ab3g=ab3g+db3g
    END DO

    adam_t=adam_t+1
    bc1=1.D0-ADAM_B1**adam_t; bc2=1.D0-ADAM_B2**adam_t
    CALL adam2(W1,mW1,vW1,aW1,NN_H,3,   bc1,bc2)
    CALL adam1(b1,mb1,vb1,ab1g,NN_H,    bc1,bc2)
    CALL adam2(W2,mW2,vW2,aW2,NN_H,NN_H,bc1,bc2)
    CALL adam1(b2,mb2,vb2,ab2g,NN_H,    bc1,bc2)
    CALL adam2(W3,mW3,vW3,aW3,3,NN_H,   bc1,bc2)
    CALL adam1(b3,mb3,vb3,ab3g,3,        bc1,bc2)

    IF (MOD(adam_t,500)==0) &
      WRITE(*,'(A,I6,A)') ' [nn_eos] bg-train step=',adam_t,' (rank0)'
  END SUBROUTINE do_train


  SUBROUTINE sample_idx(idx,Ntotal,Nsample,s)
    INTEGER,    INTENT(OUT)   :: idx(Nsample)
    INTEGER,    INTENT(IN)    :: Ntotal, Nsample
    INTEGER(8), INTENT(INOUT) :: s
    INTEGER :: i,j,Nfill; REAL(8)::r
    ! Nfill = actual entries to fill (Ntotal may be < Nsample)
    Nfill = MIN(Nsample, Ntotal)
    DO i=1,Nfill; idx(i)=i; END DO
    DO i=Nfill+1,Ntotal   ! reservoir sampling
      CALL lcg(s,r); j=1+INT(r*DBLE(i))
      IF(j<=Nfill) idx(j)=i
    END DO
  END SUBROUTINE sample_idx


  SUBROUTINE init_weights(rng)
    INTEGER(8),INTENT(INOUT)::rng
    CALL init_he(W1,NN_H,2,   rng); CALL init_he(W2,NN_H,NN_H,rng)
    CALL init_he(W3,3,   NN_H,rng)
    b1=0.D0; b2=0.D0; b3=0.D0
    W1(:,3)=0.D0   ! third input slot unused; zero weights so it never fires
    x_mean(1:2)=0.D0; x_std(1:2)=1.D0; y_mean=0.D0; y_std=1.D0
    mW1=0.D0;vW1=0.D0;mb1=0.D0;vb1=0.D0
    mW2=0.D0;vW2=0.D0;mb2=0.D0;vb2=0.D0
    mW3=0.D0;vW3=0.D0;mb3=0.D0;vb3=0.D0
    adam_t=0
  END SUBROUTINE init_weights


  SUBROUTINE load_weights(fname,ios_out)
    CHARACTER(LEN=*),INTENT(IN) ::fname
    INTEGER,         INTENT(OUT)::ios_out
    INTEGER::iunit,r,c,ios; CHARACTER(LEN=512)::msg
    iunit=75
    OPEN(UNIT=iunit,FILE=TRIM(fname),STATUS='OLD', &
         ACTION='READ',IOSTAT=ios,IOMSG=msg)
    IF(ios/=0)THEN
      WRITE(*,'(2A)')'[nn_eos] load error: ',TRIM(msg)
      ios_out=ios; RETURN
    END IF
    READ(iunit,*)
    DO r=1,NN_H;  READ(iunit,*)(W1(r,c),c=1,2);    END DO
    READ(iunit,*) (b1(r),r=1,NN_H)
    DO r=1,NN_H;  READ(iunit,*)(W2(r,c),c=1,NN_H); END DO
    READ(iunit,*) (b2(r),r=1,NN_H)
    DO r=1,3;     READ(iunit,*)(W3(r,c),c=1,NN_H); END DO
    READ(iunit,*) (b3(r),r=1,3)
    READ(iunit,*) x_mean(1),x_mean(2)
    READ(iunit,*) x_std (1),x_std (2)
    READ(iunit,*) y_mean(1),y_mean(2),y_mean(3)
    READ(iunit,*) y_std (1),y_std (2),y_std (3)
    CLOSE(iunit)
    mW1=0.D0;vW1=0.D0;mb1=0.D0;vb1=0.D0
    mW2=0.D0;vW2=0.D0;mb2=0.D0;vb2=0.D0
    mW3=0.D0;vW3=0.D0;mb3=0.D0;vb3=0.D0
    adam_t=0; ios_out=0
  END SUBROUTINE load_weights


  SUBROUTINE save_weights(fname)
    CHARACTER(LEN=*),INTENT(IN)::fname
    INTEGER::iunit,r,c
    iunit=76
    OPEN(UNIT=iunit,FILE=TRIM(fname),STATUS='REPLACE',ACTION='WRITE')
    WRITE(iunit,'(I1)') 1
    DO r=1,NN_H; WRITE(iunit,'(*(ES18.10,2X))')(W1(r,c),c=1,2);    END DO
    WRITE(iunit,'(*(ES18.10,2X))')(b1(r),r=1,NN_H)
    DO r=1,NN_H; WRITE(iunit,'(*(ES18.10,2X))')(W2(r,c),c=1,NN_H); END DO
    WRITE(iunit,'(*(ES18.10,2X))')(b2(r),r=1,NN_H)
    DO r=1,3;    WRITE(iunit,'(*(ES18.10,2X))')(W3(r,c),c=1,NN_H); END DO
    WRITE(iunit,'(*(ES18.10,2X))')(b3(r),r=1,3)
    WRITE(iunit,'(*(ES18.10,2X))') x_mean(1),x_mean(2)
    WRITE(iunit,'(*(ES18.10,2X))') x_std (1),x_std (2)
    WRITE(iunit,'(*(ES18.10,2X))') y_mean(1),y_mean(2),y_mean(3)
    WRITE(iunit,'(*(ES18.10,2X))') y_std (1),y_std (2),y_std (3)
    CLOSE(iunit)
  END SUBROUTINE save_weights


  SUBROUTINE bcast_weights()
#ifdef MPI
    INTEGER :: mpi_err
    CALL MPI_Bcast(W1,    NN_H*2,   MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
    CALL MPI_Bcast(b1,    NN_H,     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
    CALL MPI_Bcast(W2,    NN_H*NN_H,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
    CALL MPI_Bcast(b2,    NN_H,     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
    CALL MPI_Bcast(W3,    3*NN_H,   MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
    CALL MPI_Bcast(b3,    3,        MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
    CALL MPI_Bcast(x_mean,2,        MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
    CALL MPI_Bcast(x_std, 2,        MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
    CALL MPI_Bcast(y_mean,3,        MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
    CALL MPI_Bcast(y_std, 3,        MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)
#endif
  END SUBROUTINE bcast_weights


  SUBROUTINE make_fname(title,fname)
    CHARACTER(LEN=*),  INTENT(IN) ::title
    CHARACTER(LEN=256),INTENT(OUT)::fname
    CHARACTER(LEN=256)::tmp; INTEGER::i,n
    tmp=ADJUSTL(TRIM(title)); n=LEN_TRIM(tmp)
    DO i=1,n; IF(tmp(i:i)==' ')tmp(i:i)='_'; END DO
    fname=TRIM(tmp)//'.txt'
  END SUBROUTINE make_fname


  SUBROUTINE fwd3(x1,x2,x3,h1o,h2o,h3o)
    REAL(8),INTENT(IN) ::x1,x2,x3
    REAL(8),INTENT(OUT)::h1o(NN_H),h2o(NN_H),h3o(3)
    REAL(8)::acc; INTEGER::r,c
    DO r=1,NN_H
      acc=b1(r)+W1(r,1)*x1+W1(r,2)*x2+W1(r,3)*x3; h1o(r)=TANH(acc)
    END DO
    DO r=1,NN_H; acc=b2(r)
      DO c=1,NN_H; acc=acc+W2(r,c)*h1o(c); END DO; h2o(r)=TANH(acc)
    END DO
    DO r=1,3; acc=b3(r)
      DO c=1,NN_H; acc=acc+W3(r,c)*h2o(c); END DO; h3o(r)=acc
    END DO
  END SUBROUTINE fwd3

  SUBROUTINE bwd3e(h2o,dout,dWo,dbo,din)
    REAL(8),INTENT(IN) ::h2o(NN_H),dout(3)
    REAL(8),INTENT(OUT)::dWo(3,NN_H),dbo(3),din(NN_H)
    INTEGER::i,j
    DO i=1,3; dbo(i)=dout(i)
      DO j=1,NN_H; dWo(i,j)=dout(i)*h2o(j); END DO; END DO
    DO j=1,NN_H; din(j)=0.D0
      DO i=1,3; din(j)=din(j)+W3(i,j)*dout(i); END DO; END DO
  END SUBROUTINE bwd3e

  SUBROUTINE bwd2e(h1o,h2o,dout,dWo,dbo,din)
    REAL(8),INTENT(IN) ::h1o(NN_H),h2o(NN_H),dout(NN_H)
    REAL(8),INTENT(OUT)::dWo(NN_H,NN_H),dbo(NN_H),din(NN_H)
    REAL(8)::delta(NN_H); INTEGER::i,j
    DO i=1,NN_H; delta(i)=dout(i)*(1.D0-h2o(i)**2); dbo(i)=delta(i)
      DO j=1,NN_H; dWo(i,j)=delta(i)*h1o(j); END DO; END DO
    DO j=1,NN_H; din(j)=0.D0
      DO i=1,NN_H; din(j)=din(j)+W2(i,j)*delta(i); END DO; END DO
  END SUBROUTINE bwd2e

  SUBROUTINE bwd1e(x1,x2,x3,h1o,dout,dWo,dbo)
    REAL(8),INTENT(IN) ::x1,x2,x3,h1o(NN_H),dout(NN_H)
    REAL(8),INTENT(OUT)::dWo(NN_H,3),dbo(NN_H)
    REAL(8)::delta(NN_H); INTEGER::i
    DO i=1,NN_H; delta(i)=dout(i)*(1.D0-h1o(i)**2); dbo(i)=delta(i)
      dWo(i,1)=delta(i)*x1; dWo(i,2)=delta(i)*x2; dWo(i,3)=delta(i)*x3
    END DO
  END SUBROUTINE bwd1e

  SUBROUTINE adam2(W,mW,vW,gW,nr,nc,bc1,bc2)
    INTEGER,INTENT(IN)::nr,nc
    REAL(8),INTENT(INOUT)::W(nr,nc),mW(nr,nc),vW(nr,nc)
    REAL(8),INTENT(IN)::gW(nr,nc),bc1,bc2; INTEGER::i,j
    DO i=1,nr; DO j=1,nc
      mW(i,j)=ADAM_B1*mW(i,j)+(1.D0-ADAM_B1)*gW(i,j)
      vW(i,j)=ADAM_B2*vW(i,j)+(1.D0-ADAM_B2)*gW(i,j)**2
      W(i,j)=W(i,j)-LR*(mW(i,j)/bc1)/(SQRT(MAX(vW(i,j)/bc2,0.D0))+ADAM_EPS)
    END DO; END DO
  END SUBROUTINE adam2

  SUBROUTINE adam1(b,mb,vb,gb,nr,bc1,bc2)
    INTEGER,INTENT(IN)::nr
    REAL(8),INTENT(INOUT)::b(nr),mb(nr),vb(nr)
    REAL(8),INTENT(IN)::gb(nr),bc1,bc2; INTEGER::i
    DO i=1,nr
      mb(i)=ADAM_B1*mb(i)+(1.D0-ADAM_B1)*gb(i)
      vb(i)=ADAM_B2*vb(i)+(1.D0-ADAM_B2)*gb(i)**2
      b(i)=b(i)-LR*(mb(i)/bc1)/(SQRT(MAX(vb(i)/bc2,0.D0))+ADAM_EPS)
    END DO
  END SUBROUTINE adam1

  SUBROUTINE init_he(W,nr,nc,s)
    INTEGER,INTENT(IN)::nr,nc; REAL(8),INTENT(OUT)::W(nr,nc)
    INTEGER(8),INTENT(INOUT)::s; REAL(8)::u1,u2,rv,sc; INTEGER::i,j
    sc=SQRT(2.D0/DBLE(nc))
    DO i=1,nr; DO j=1,nc
      CALL lcg(s,u1); u1=MAX(u1,1.D-10); CALL lcg(s,u2)
      rv=SQRT(-2.D0*LOG(u1))*COS(6.28318530718D0*u2); W(i,j)=rv*sc
    END DO; END DO
  END SUBROUTINE init_he

  SUBROUTINE lcg(s,r)
    INTEGER(8),INTENT(INOUT)::s; REAL(8),INTENT(OUT)::r
    s=MOD(6364136223846793005_8*s+1442695040888963407_8,2147483648_8)
    r=DBLE(ABS(s))/2147483648.D0
  END SUBROUTINE lcg

END MODULE nn_eos_func_mod
