!===========================================================================
! CROCO GLS Neural Network Stability Functions
! Per-dt online learning
!
!!!--------------------------------------------------------------------------
! USE: its only testing purpose of different methods (random by __sm-loop__)
!!!--------------------------------------------------------------------------
!                    ------We finalize it soon------
! ACTIVATION:  #define GLS_NN_STAB  in cppdefs.h
!============================================================================

#include "cppdefs.h"

MODULE nn_stab_func_mod

  IMPLICIT NONE
  PRIVATE

#ifdef MPI
  include 'mpif.h'
#endif

  !--------------------------------------------------------------------
  ! Tuneable parameters
  !--------------------------------------------------------------------
  INTEGER,  PARAMETER :: NN_H     = 32     ! hidden layer width
  INTEGER,  PARAMETER :: N_BATCH  = 2000   ! mini-batch per dt step
  REAL(8),  PARAMETER :: LR       = 1.D-2
  REAL(8),  PARAMETER :: ADAM_B1  = 0.9D0
  REAL(8),  PARAMETER :: ADAM_B2  = 0.999D0
  REAL(8),  PARAMETER :: ADAM_EPS = 1.D-8

  !--------------------------------------------------------------------
  ! Canuto A constants
  !--------------------------------------------------------------------
  ! Derived from Canuto A (2001)
  ! c1=5,c2=0.8,c3=1.968,c4=1.136,c6=0.4
  ! cb1=5.95,cb2=0.6,cb3=1,cb4=0,cb5=0.3333,cbb=0.72
  
  REAL(8), PARAMETER :: ca_a1=2.D0/3.D0-0.5D0*0.8D0    ! 0.26667
  REAL(8), PARAMETER :: ca_a2=1.D0-0.5D0*1.968D0        ! 0.016
  REAL(8), PARAMETER :: ca_a3=1.D0-0.5D0*1.136D0        ! 0.432
  REAL(8), PARAMETER :: ca_a5=0.5D0-0.5D0*0.4D0         ! 0.3
  REAL(8), PARAMETER :: ca_ab1=1.D0-0.6D0               ! 0.4
  REAL(8), PARAMETER :: ca_ab2=1.D0-1.D0                ! 0.0
  REAL(8), PARAMETER :: ca_ab3=2.D0*(1.D0-0.D0)         ! 2.0
  REAL(8), PARAMETER :: ca_ab5=2.D0*0.72D0*(1.D0-0.3333D0) ! 0.96005
  REAL(8), PARAMETER :: ca_nn=2.5D0,ca_nb=5.95D0
  REAL(8), PARAMETER :: sf_d0=36.D0*ca_nn**3*ca_nb**2
  REAL(8), PARAMETER :: sf_d1=84.D0*ca_a5*ca_ab3*ca_nn**2*ca_nb &
                              +36.D0*ca_ab5*ca_nn**3*ca_nb
  REAL(8), PARAMETER :: sf_d2=9.D0*(ca_ab2**2-ca_ab1**2)*ca_nn**3 &
                              -12.D0*(ca_a2**2-3.D0*ca_a3**2)*ca_nn*ca_nb**2
  REAL(8), PARAMETER :: sf_d3=12.D0*ca_a5*ca_ab3 &
                                   *(ca_a2*ca_ab1-3.D0*ca_a3*ca_ab2)*ca_nn &
                              +12.D0*ca_a5*ca_ab3*(ca_a3**2-ca_a2**2)*ca_nb &
                              +12.D0*ca_ab5*(3.D0*ca_a3**2-ca_a2**2)*ca_nn*ca_nb
  REAL(8), PARAMETER :: sf_d4=48.D0*ca_a5**2*ca_ab3**2*ca_nn &
                              +36.D0*ca_a5*ca_ab3*ca_ab5*ca_nn**2
  REAL(8), PARAMETER :: sf_d5=3.D0*(ca_a2**2-3.D0*ca_a3**2) &
                                   *(ca_ab1**2-ca_ab2**2)*ca_nn
  REAL(8), PARAMETER :: sf_n0=36.D0*ca_a1*ca_nn**2*ca_nb**2
  REAL(8), PARAMETER :: sf_n1=-12.D0*ca_a5*ca_ab3*(ca_ab1+ca_ab2)*ca_nn**2 &
                               +8.D0*ca_a5*ca_ab3 &
                                    *(6.D0*ca_a1-ca_a2-3.D0*ca_a3)*ca_nn*ca_nb &
                              +36.D0*ca_a1*ca_ab5*ca_nn**2*ca_nb
  REAL(8), PARAMETER :: sf_n2=9.D0*ca_a1*(ca_ab2**2-ca_ab1**2)*ca_nn**2
  REAL(8), PARAMETER :: sf_nb0=12.D0*ca_ab3*ca_nn**3*ca_nb
  REAL(8), PARAMETER :: sf_nb1=12.D0*ca_a5*ca_ab3**2*ca_nn**2
  REAL(8), PARAMETER :: sf_nb2=9.D0*ca_a1*ca_ab3*(ca_ab1-ca_ab2)*ca_nn**2 &
                               +(6.D0*ca_a1*(ca_a2-3.D0*ca_a3) &
                               -4.D0*(ca_a2**2-3.D0*ca_a3**2))*ca_ab3*ca_nn*ca_nb

  !--------------------------------------------------------------------
  ! Network weights (-----random init-----)
  !--------------------------------------------------------------------
  REAL(8), SAVE :: W1(NN_H,2),    b1(NN_H)
  REAL(8), SAVE :: W2(NN_H,NN_H), b2(NN_H)
  REAL(8), SAVE :: W3(2,NN_H),    b3(2)
  REAL(8), SAVE :: x_mean(2)=0.D0, x_std(2)=1.D0
  REAL(8), SAVE :: y_mean(2)=0.D0, y_std(2)=1.D0

  ! Adam 
  REAL(8), SAVE :: mW1(NN_H,2),    vW1(NN_H,2)
  REAL(8), SAVE :: mW2(NN_H,NN_H), vW2(NN_H,NN_H)
  REAL(8), SAVE :: mW3(2,NN_H),    vW3(2,NN_H)
  REAL(8), SAVE :: mb1(NN_H), vb1(NN_H)
  REAL(8), SAVE :: mb2(NN_H), vb2(NN_H)
  REAL(8), SAVE :: mb3(2),    vb3(2)
  INTEGER, SAVE :: adam_t = 0

  ! Buffer (rank 0 only, ALLOCATABLE     auto-sized from tile dims)
  REAL(4), ALLOCATABLE, SAVE :: buf_an(:), buf_am(:)
  REAL(4), ALLOCATABLE, SAVE :: buf_cmu(:), buf_cmup(:)
  INTEGER, SAVE :: buf_n    = 0
  INTEGER, SAVE :: buf_max  = 0
  INTEGER, SAVE :: last_iic = -9999

  ! State flags
  INTEGER, SAVE :: nn_ready  = 0  ! 0=use Canuto A, 1=use NN
  LOGICAL, SAVE :: buf_alloc = .FALSE.

  CHARACTER(LEN=256), SAVE :: wfile = ''

  PUBLIC :: nn_stab_init
  PUBLIC :: nn_stab_func
  PUBLIC :: nn_stab_finalize

CONTAINS

!=======================================================================
  SUBROUTINE nn_stab_init(run_title, my_rank)
!=======================================================================
    CHARACTER(LEN=*), INTENT(IN) :: run_title
    INTEGER,          INTENT(IN) :: my_rank

    INTEGER    :: ios
    INTEGER(8) :: rng
    LOGICAL    :: fexist

    CALL make_fname(run_title, wfile)

    IF (my_rank == 0) THEN
      WRITE(*,'(A,A)') ' [nn_stab] weight file: ', TRIM(wfile)
      INQUIRE(FILE=TRIM(wfile), EXIST=fexist)
      IF (fexist) THEN
        CALL load_weights(TRIM(wfile), ios)
        IF (ios == 0) THEN
          nn_ready = 1
          WRITE(*,'(A)') &
            ' [nn_stab] weights loaded -> NN active from step 1'
        ELSE
          rng = 42_8; CALL init_weights(rng)
          WRITE(*,'(A)') &
            ' [nn_stab] bad file -> Canuto A used, training in background'
        END IF
      ELSE
        rng = 42_8; CALL init_weights(rng)
        WRITE(*,'(A)') &
          ' [nn_stab] no weight file -> Canuto A this run, NN next run'
      END IF
    END IF

    ! Broadcast nn_ready flag and weights (ONCE, at startup only) for speed only
#ifdef MPI
    BLOCK
      INTEGER :: mpi_err
      CALL MPI_Bcast(nn_ready, 1, MPI_INTEGER, &
                     0, MPI_COMM_WORLD, mpi_err)
      IF (nn_ready == 1) THEN
        CALL MPI_Bcast(W1,     NN_H*2,    MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, mpi_err)
        CALL MPI_Bcast(b1,     NN_H,      MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, mpi_err)
        CALL MPI_Bcast(W2,     NN_H*NN_H, MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, mpi_err)
        CALL MPI_Bcast(b2,     NN_H,      MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, mpi_err)
        CALL MPI_Bcast(W3,     2*NN_H,    MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, mpi_err)
        CALL MPI_Bcast(b3,     2,         MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, mpi_err)
        CALL MPI_Bcast(x_mean, 2,         MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, mpi_err)
        CALL MPI_Bcast(x_std,  2,         MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, mpi_err)
        CALL MPI_Bcast(y_mean, 2,         MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, mpi_err)
        CALL MPI_Bcast(y_std,  2,         MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, mpi_err)
      END IF
    END BLOCK
#endif

  END SUBROUTINE nn_stab_init


!=======================================================================
  SUBROUTINE nn_stab_func(alpha_n, alpha_m, c_mu, c_mu_prim, &
                           my_rank, timestep,                 &
                           Istr, Iend, Jstr, Jend, Nz)
!=======================================================================
!  Call per cell inside gls_mixing_tile DO j/k/i loop.
!
!          *** NO MPI CALLS HERE. EVER. ***
!
!=======================================================================
    REAL,    INTENT(INOUT) :: alpha_n, alpha_m
    REAL,    INTENT(OUT)   :: c_mu, c_mu_prim
    INTEGER, INTENT(IN)    :: my_rank, timestep
    INTEGER, INTENT(IN)    :: Istr, Iend, Jstr, Jend, Nz

    REAL(8) :: Denom, cmu_d, cmup_d
    REAL(8) :: h1(NN_H), h2(NN_H), h3(2), acc
    INTEGER :: r, c

    !------------------------------------------------------------------
    ! RANK 0: allocate buffer on first call
    !------------------------------------------------------------------
    IF (my_rank == 0 .AND. .NOT. buf_alloc) THEN
      buf_max = (Iend-Istr+1) * (Jend-Jstr+1) * (Nz-1) + 1000
      ALLOCATE(buf_an(buf_max), buf_am(buf_max))
      ALLOCATE(buf_cmu(buf_max), buf_cmup(buf_max))
      buf_n     = 0
      buf_alloc = .TRUE.
      WRITE(*,'(A,I8)') ' [nn_stab] buffer auto-sized: ', buf_max
    END IF

    !------------------------------------------------------------------
    ! RANK 0: at timestep boundary, train on previous step buffer
    !------------------------------------------------------------------
    IF (my_rank == 0 .AND. timestep /= last_iic) THEN
      IF (buf_n >= 10) CALL do_train(buf_n)
      buf_n    = 0
      last_iic = timestep
    END IF

    !------------------------------------------------------------------
    ! Compute exact Canuto A
    !------------------------------------------------------------------
    Denom  = sf_d0 + sf_d1*DBLE(alpha_n) + sf_d2*DBLE(alpha_m) &
           + sf_d3*DBLE(alpha_n)*DBLE(alpha_m) &
           + sf_d4*DBLE(alpha_n)**2 + sf_d5*DBLE(alpha_m)**2
    cmu_d  = MAX((sf_n0+sf_n1*DBLE(alpha_n)+sf_n2*DBLE(alpha_m)) &
                 /Denom, 1.D-10)
    cmup_d = MAX((sf_nb0+sf_nb1*DBLE(alpha_n)+sf_nb2*DBLE(alpha_m)) &
                 /Denom, 1.D-10)

    !------------------------------------------------------------------
    ! RANK 0: buffer this cell
    !------------------------------------------------------------------
    IF (my_rank == 0 .AND. buf_alloc .AND. buf_n < buf_max) THEN
      buf_n           = buf_n + 1
      buf_an  (buf_n) = alpha_n
      buf_am  (buf_n) = alpha_m
      buf_cmu (buf_n) = REAL(cmu_d)
      buf_cmup(buf_n) = REAL(cmup_d)
    END IF

    !------------------------------------------------------------------
    ! Return: NN inference
    !------------------------------------------------------------------
    IF (nn_ready == 1) THEN
      ! NN forward pass
      h1(1) = (DBLE(alpha_n) - x_mean(1)) / (x_std(1) + 1.D-30)
      h1(2) = (DBLE(alpha_m) - x_mean(2)) / (x_std(2) + 1.D-30)
      DO r = 1, NN_H
        acc = b1(r) + W1(r,1)*h1(1) + W1(r,2)*h1(2)
        h2(r) = TANH(acc)
      END DO
      DO r = 1, NN_H
        acc = b2(r)
        DO c = 1, NN_H; acc = acc + W2(r,c)*h2(c); END DO
        h1(r) = TANH(acc)
      END DO
      DO r = 1, 2
        acc = b3(r)
        DO c = 1, NN_H; acc = acc + W3(r,c)*h1(c); END DO
        h3(r) = acc
      END DO
      c_mu      = REAL(EXP(MAX(MIN(h3(1)*y_std(1)+y_mean(1), &
                                    10.D0), -10.D0)))
      c_mu_prim = REAL(EXP(MAX(MIN(h3(2)*y_std(2)+y_mean(2), &
                                    10.D0), -10.D0)))
    ELSE
      ! Exact Canuto A
      c_mu      = REAL(cmu_d)
      c_mu_prim = REAL(cmup_d)
    END IF

    c_mu      = MAX(c_mu,      1.E-6)
    c_mu_prim = MAX(c_mu_prim, 1.E-6)

  END SUBROUTINE nn_stab_func


!=======================================================================
  SUBROUTINE nn_stab_finalize(my_rank)
!=======================================================================
    INTEGER, INTENT(IN) :: my_rank

    IF (my_rank /= 0) RETURN

    ! Final ***training*** on last step buffer
    IF (buf_alloc .AND. buf_n >= 10) CALL do_train(buf_n)

    ! Save weights     available for next run ???
    CALL save_weights(TRIM(wfile))
    WRITE(*,'(A,A)') ' [nn_stab] weights saved -> ', TRIM(wfile)
    WRITE(*,'(A,I8,A)') ' [nn_stab] total Adam steps: ', adam_t, &
      '  (loaded next run for instant NN)'

    IF (ALLOCATED(buf_an))   DEALLOCATE(buf_an)
    IF (ALLOCATED(buf_am))   DEALLOCATE(buf_am)
    IF (ALLOCATED(buf_cmu))  DEALLOCATE(buf_cmu)
    IF (ALLOCATED(buf_cmup)) DEALLOCATE(buf_cmup)

  END SUBROUTINE nn_stab_finalize


!=======================================================================
!  PRIVATE ROUTINES, BUT NOTHING IS PRIVATE
!=======================================================================

  SUBROUTINE do_train(Ntrain)
    INTEGER, INTENT(IN) :: Ntrain
    INTEGER  :: n, it
    INTEGER  :: idx(N_BATCH)
    REAL(8)  :: logYc(N_BATCH), logYcp(N_BATCH)
    REAL(8)  :: Xn1(N_BATCH),   Xn2(N_BATCH)
    REAL(8)  :: Yn1(N_BATCH),   Yn2(N_BATCH)
    REAL(8)  :: xm(2), xs(2), ym(2), ys(2)
    REAL(8)  :: h1o(NN_H), h2o(NN_H), h3o(2)
    REAL(8)  :: dh2(NN_H), dh1(NN_H)
    REAL(8)  :: dW1(NN_H,2),    db1g(NN_H)
    REAL(8)  :: dW2(NN_H,NN_H), db2g(NN_H)
    REAL(8)  :: dW3(2,NN_H),    db3g(2)
    REAL(8)  :: aW1(NN_H,2),    ab1g(NN_H)
    REAL(8)  :: aW2(NN_H,NN_H), ab2g(NN_H)
    REAL(8)  :: aW3(2,NN_H),    ab3g(2)
    REAL(8)  :: bc1, bc2, dl1, dl2
    INTEGER(8) :: s

    s = INT(adam_t+1, 8) * 998244353_8

    ! Reservoir sample N_BATCH from buffer
    CALL sample_idx(idx, Ntrain, N_BATCH, s)

    ! Log + normalise on mini-batch only (its mini?)
    DO n = 1, N_BATCH
      logYc(n)  = LOG(MAX(DBLE(buf_cmu (idx(n))), 1.D-10))
      logYcp(n) = LOG(MAX(DBLE(buf_cmup(idx(n))), 1.D-10))
    END DO
    xm(1)=SUM(DBLE(buf_an(idx(1:N_BATCH))))/N_BATCH
    xm(2)=SUM(DBLE(buf_am(idx(1:N_BATCH))))/N_BATCH
    xs(1)=SQRT(SUM((DBLE(buf_an(idx(1:N_BATCH)))-xm(1))**2)/N_BATCH)+1.D-30
    xs(2)=SQRT(SUM((DBLE(buf_am(idx(1:N_BATCH)))-xm(2))**2)/N_BATCH)+1.D-30
    ym(1)=SUM(logYc)/N_BATCH;  ym(2)=SUM(logYcp)/N_BATCH
    ys(1)=MAX(SQRT(SUM((logYc -ym(1))**2)/N_BATCH),1.D-4)
    ys(2)=MAX(SQRT(SUM((logYcp-ym(2))**2)/N_BATCH),1.D-4)
    DO n = 1, N_BATCH
      Xn1(n)=(DBLE(buf_an(idx(n)))-xm(1))/xs(1)
      Xn2(n)=(DBLE(buf_am(idx(n)))-xm(2))/xs(2)
      Yn1(n)=(logYc(n) -ym(1))/ys(1)
      Yn2(n)=(logYcp(n)-ym(2))/ys(2)
    END DO
    x_mean=xm; x_std=xs; y_mean=ym; y_std=ys

    ! One Adam step on mini-batch (style old 1990s...)
    aW1=0.D0; ab1g=0.D0
    aW2=0.D0; ab2g=0.D0
    aW3=0.D0; ab3g=0.D0
    DO it = 1, N_BATCH
      CALL fwd(Xn1(it),Xn2(it),h1o,h2o,h3o)
      dl1=2.D0*(h3o(1)-Yn1(it))/DBLE(N_BATCH)
      dl2=2.D0*(h3o(2)-Yn2(it))/DBLE(N_BATCH)
      CALL bwd3(h2o,[dl1,dl2],dW3,db3g,dh2)
      CALL bwd2(h1o,h2o,dh2, dW2,db2g,dh1)
      CALL bwd1(Xn1(it),Xn2(it),h1o,dh1,dW1,db1g)
      aW1=aW1+dW1; ab1g=ab1g+db1g
      aW2=aW2+dW2; ab2g=ab2g+db2g
      aW3=aW3+dW3; ab3g=ab3g+db3g
    END DO
    adam_t=adam_t+1
    bc1=1.D0-ADAM_B1**adam_t; bc2=1.D0-ADAM_B2**adam_t
    CALL adam2(W1,mW1,vW1,aW1,NN_H,2,   bc1,bc2)
    CALL adam1(b1,mb1,vb1,ab1g,NN_H,    bc1,bc2)
    CALL adam2(W2,mW2,vW2,aW2,NN_H,NN_H,bc1,bc2)
    CALL adam1(b2,mb2,vb2,ab2g,NN_H,    bc1,bc2)
    CALL adam2(W3,mW3,vW3,aW3,2,NN_H,   bc1,bc2)
    CALL adam1(b3,mb3,vb3,ab3g,2,        bc1,bc2)
    IF (MOD(adam_t,500)==0) &
      WRITE(*,'(A,I6,A)') ' [nn_stab] bg-train step=',adam_t,' (rank0)'
  END SUBROUTINE do_train


  SUBROUTINE sample_idx(idx,Ntotal,Nsample,s)
    INTEGER,    INTENT(OUT)   :: idx(Nsample)
    INTEGER,    INTENT(IN)    :: Ntotal,Nsample
    INTEGER(8), INTENT(INOUT) :: s
    INTEGER :: i,j; REAL(8)::r
    DO i=1,MIN(Nsample,Ntotal); idx(i)=i; END DO
    DO i=Nsample+1,Ntotal
      CALL lcg(s,r); j=1+INT(r*DBLE(i))
      IF(j<=Nsample) idx(j)=i
    END DO
  END SUBROUTINE sample_idx


  SUBROUTINE init_weights(rng)
    INTEGER(8),INTENT(INOUT)::rng
    CALL init_he(W1,NN_H,2,   rng); CALL init_he(W2,NN_H,NN_H,rng)
    CALL init_he(W3,2,   NN_H,rng)
    b1=0.D0; b2=0.D0; b3=0.D0
    x_mean=0.D0; x_std=1.D0; y_mean=0.D0; y_std=1.D0
    mW1=0.D0;vW1=0.D0;mb1=0.D0;vb1=0.D0
    mW2=0.D0;vW2=0.D0;mb2=0.D0;vb2=0.D0
    mW3=0.D0;vW3=0.D0;mb3=0.D0;vb3=0.D0
    adam_t=0
  END SUBROUTINE init_weights


  SUBROUTINE load_weights(fname,ios_out)
    CHARACTER(LEN=*),INTENT(IN) ::fname
    INTEGER,         INTENT(OUT)::ios_out
    INTEGER::iunit,r,c,ios; CHARACTER(LEN=512)::msg
    iunit=73
    OPEN(UNIT=iunit,FILE=TRIM(fname),STATUS='OLD', &
         ACTION='READ',IOSTAT=ios,IOMSG=msg)
    IF(ios/=0)THEN
      WRITE(*,'(2A)')'[nn_stab] load error: ',TRIM(msg)
      ios_out=ios; RETURN
    END IF
    READ(iunit,*)
    DO r=1,NN_H;  READ(iunit,*)(W1(r,c),c=1,2);    END DO
    READ(iunit,*) (b1(r),r=1,NN_H)
    DO r=1,NN_H;  READ(iunit,*)(W2(r,c),c=1,NN_H); END DO
    READ(iunit,*) (b2(r),r=1,NN_H)
    DO r=1,2;     READ(iunit,*)(W3(r,c),c=1,NN_H); END DO
    READ(iunit,*) (b3(r),r=1,2)
    READ(iunit,*) x_mean(1),x_mean(2)
    READ(iunit,*) x_std (1),x_std (2)
    READ(iunit,*) y_mean(1),y_mean(2)
    READ(iunit,*) y_std (1),y_std (2)
    CLOSE(iunit)
    mW1=0.D0;vW1=0.D0;mb1=0.D0;vb1=0.D0
    mW2=0.D0;vW2=0.D0;mb2=0.D0;vb2=0.D0
    mW3=0.D0;vW3=0.D0;mb3=0.D0;vb3=0.D0
    adam_t=0; ios_out=0
  END SUBROUTINE load_weights


  SUBROUTINE save_weights(fname)
    CHARACTER(LEN=*),INTENT(IN)::fname
    INTEGER::iunit,r,c
    iunit=74
    OPEN(UNIT=iunit,FILE=TRIM(fname),STATUS='REPLACE',ACTION='WRITE')
    WRITE(iunit,'(I1)') 1
    DO r=1,NN_H; WRITE(iunit,'(*(ES18.10,2X))')(W1(r,c),c=1,2);    END DO
    WRITE(iunit,'(*(ES18.10,2X))')(b1(r),r=1,NN_H)
    DO r=1,NN_H; WRITE(iunit,'(*(ES18.10,2X))')(W2(r,c),c=1,NN_H); END DO
    WRITE(iunit,'(*(ES18.10,2X))')(b2(r),r=1,NN_H)
    DO r=1,2;    WRITE(iunit,'(*(ES18.10,2X))')(W3(r,c),c=1,NN_H); END DO
    WRITE(iunit,'(*(ES18.10,2X))')(b3(r),r=1,2)
    WRITE(iunit,'(*(ES18.10,2X))') x_mean(1),x_mean(2)
    WRITE(iunit,'(*(ES18.10,2X))') x_std (1),x_std (2)
    WRITE(iunit,'(*(ES18.10,2X))') y_mean(1),y_mean(2)
    WRITE(iunit,'(*(ES18.10,2X))') y_std (1),y_std (2)
    CLOSE(iunit)
  END SUBROUTINE save_weights


  SUBROUTINE make_fname(title,fname)
    CHARACTER(LEN=*),  INTENT(IN) ::title
    CHARACTER(LEN=256),INTENT(OUT)::fname
    CHARACTER(LEN=256)::tmp; INTEGER::i,n
    tmp=ADJUSTL(TRIM(title)); n=LEN_TRIM(tmp)
    DO i=1,n; IF(tmp(i:i)==' ')tmp(i:i)='_'; END DO
    fname='nn_weights_'//TRIM(tmp)//'.txt'
  END SUBROUTINE make_fname


  SUBROUTINE fwd(xin1,xin2,h1o,h2o,h3o)
    REAL(8),INTENT(IN) ::xin1,xin2
    REAL(8),INTENT(OUT)::h1o(NN_H),h2o(NN_H),h3o(2)
    REAL(8)::acc; INTEGER::r,c
    DO r=1,NN_H; acc=b1(r)+W1(r,1)*xin1+W1(r,2)*xin2; h1o(r)=TANH(acc); END DO
    DO r=1,NN_H; acc=b2(r)
      DO c=1,NN_H; acc=acc+W2(r,c)*h1o(c); END DO; h2o(r)=TANH(acc); END DO
    DO r=1,2; acc=b3(r)
      DO c=1,NN_H; acc=acc+W3(r,c)*h2o(c); END DO; h3o(r)=acc; END DO
  END SUBROUTINE fwd

  SUBROUTINE bwd3(h2o,dout,dWo,dbo,din)
    REAL(8),INTENT(IN)::h2o(NN_H),dout(2)
    REAL(8),INTENT(OUT)::dWo(2,NN_H),dbo(2),din(NN_H)
    INTEGER::i,j
    DO i=1,2; dbo(i)=dout(i)
      DO j=1,NN_H; dWo(i,j)=dout(i)*h2o(j); END DO; END DO
    DO j=1,NN_H; din(j)=0.D0
      DO i=1,2; din(j)=din(j)+W3(i,j)*dout(i); END DO; END DO
  END SUBROUTINE bwd3

  SUBROUTINE bwd2(h1o,h2o,dout,dWo,dbo,din)
    REAL(8),INTENT(IN)::h1o(NN_H),h2o(NN_H),dout(NN_H)
    REAL(8),INTENT(OUT)::dWo(NN_H,NN_H),dbo(NN_H),din(NN_H)
    REAL(8)::delta(NN_H); INTEGER::i,j
    DO i=1,NN_H; delta(i)=dout(i)*(1.D0-h2o(i)**2); dbo(i)=delta(i)
      DO j=1,NN_H; dWo(i,j)=delta(i)*h1o(j); END DO; END DO
    DO j=1,NN_H; din(j)=0.D0
      DO i=1,NN_H; din(j)=din(j)+W2(i,j)*delta(i); END DO; END DO
  END SUBROUTINE bwd2

  SUBROUTINE bwd1(xin1,xin2,h1o,dout,dWo,dbo)
    REAL(8),INTENT(IN)::xin1,xin2,h1o(NN_H),dout(NN_H)
    REAL(8),INTENT(OUT)::dWo(NN_H,2),dbo(NN_H)
    REAL(8)::delta(NN_H); INTEGER::i
    DO i=1,NN_H; delta(i)=dout(i)*(1.D0-h1o(i)**2); dbo(i)=delta(i)
      dWo(i,1)=delta(i)*xin1; dWo(i,2)=delta(i)*xin2; END DO
  END SUBROUTINE bwd1

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

END MODULE nn_stab_func_mod
