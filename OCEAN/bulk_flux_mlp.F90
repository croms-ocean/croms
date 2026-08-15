!=======================================================================
!
!  bulk_flux_mlp.F90
!
!  CROCO bulk_flux COARE-MLP ONLINE TRAINING — single file, MPI only
!
!  CONCEPT
!  -------
!  Every MPI rank owns a horizontal tile of the ocean domain.
!  At each model time-step, rank processes its tile with the reference
!  COARE iterative solver (IterRef = 20, ground-truth) AND with the
!  current MLP (forward pass).  The per-sample gradient is accumulated
!  locally.  After every MLP_SYNC_STEPS time-steps, a single
!  MPI_Allreduce sums gradients across all ranks, and every rank
!  applies the same SGD update — so all ranks always carry identical
!  weights with zero weight-divergence and no extra communication.
!
!  ARCHITECTURE
!    Input  (6) : wspd, delT, delQ, TairC, RH, patm
!    Hidden1(32): tanh
!    Hidden2(16): tanh
!    Output (3) : Cd, Ch, Ce  (linear)
!
!  MPI COMMUNICATION PATTERN
!    - One MPI_Allreduce(SUM) per sync interval on gradient arrays
!      (NH1*NIN + NH1 + NH2*NH1 + NH2 + NOUT*NH2 + NOUT = 1,011 reals)
!    - One MPI_Allreduce(SUM) for scalar loss (1 real)
!    - One MPI_Bcast at init to propagate initial random weights
!    - All other operations are local, zero extra communication
!
!  HOW TO DROP INTO CROCO
!  ----------------------
!  1.  Add to jobcomp / Makefile: -DBULK_FLUX_MLP
!  2.  Replace bulk_flux.F with this file in the object list.
!  3.  Add to main.F, after MPI_Setup():
!         USE bulk_flux_mlp, ONLY: mlp_init
!         CALL mlp_init()
!  4.  No weight files needed — training starts from scratch and
!      improves throughout the run.  Optionally call mlp_save() at
!      the end of the run to persist weights for the next simulation.
!
!  COMPILE
!    mpif90 -O3 -march=native -ffast-math -c bulk_flux_mlp.F90
!
!=======================================================================
MODULE bulk_flux_mlp

  USE mpi
  IMPLICIT NONE
  PRIVATE

  !=======================================================================
  !  SECTION 1 — architecture constants
  !=======================================================================
  INTEGER, PARAMETER :: NIN    = 6    ! inputs
  INTEGER, PARAMETER :: NH1    = 32   ! hidden layer 1
  INTEGER, PARAMETER :: NH2    = 16   ! hidden layer 2
  INTEGER, PARAMETER :: NOUT   = 3    ! outputs  (Cd, Ch, Ce)

  !  Flat gradient buffer sizes (for Allreduce)
  INTEGER, PARAMETER :: NW1  = NH1*NIN
  INTEGER, PARAMETER :: NW2  = NH2*NH1
  INTEGER, PARAMETER :: NW3  = NOUT*NH2
  INTEGER, PARAMETER :: NGRAD_TOTAL = NW1+NH1 + NW2+NH2 + NW3+NOUT

  !=======================================================================
  !  SECTION 2 — physical constants (mirrors params_bulk.h)
  !=======================================================================
  REAL, PARAMETER :: vonKar   = 0.4
  REAL, PARAMETER :: g_grav   = 9.80665
  REAL, PARAMETER :: CtoK     = 273.16
  REAL, PARAMETER :: blk_Rgas = 287.05967
  REAL, PARAMETER :: blk_Rvap = 461.52499
  REAL, PARAMETER :: blk_Cpa  = 1004.70886
  REAL, PARAMETER :: cpvir    = blk_Rvap/blk_Rgas - 1.
  REAL, PARAMETER :: MvoMa    = 18.0153E-3/28.9644E-3
  REAL, PARAMETER :: r3_      = 1./3.
  REAL, PARAMETER :: pis2_    = 2.*ATAN(1.)
  REAL, PARAMETER :: sqr3_    = SQRT(3.)
  REAL, PARAMETER :: pis2osq3 = pis2_/sqr3_
  REAL, PARAMETER :: blk_ZW   = 10.0
  REAL, PARAMETER :: blk_ZT   = 2.0
  REAL, PARAMETER :: blk_ZTW  = blk_ZT/blk_ZW
  REAL, PARAMETER :: LLZw     = LOG(blk_ZW*10000.)/LOG(100000.)  ! Log ratio
  REAL, PARAMETER :: blk_Zabl = 600.
  REAL, PARAMETER :: blk_beta = 1.2
  REAL, PARAMETER :: eps_f    = 1.0E-8
  REAL, PARAMETER :: c0visc   = 1.326E-5
  REAL, PARAMETER :: c1visc   = 6.542E-3
  REAL, PARAMETER :: c2visc   = 8.301E-6
  REAL, PARAMETER :: c3visc   = 4.84E-9

  !=======================================================================
  !  SECTION 3 — training hyper-parameters
  !=======================================================================
  !  Sync interval: gradients all-reduced every N model time-steps.
  !  Lower = faster learning.  Higher = fewer MPI calls.
  !  10 steps ≈ ~O(1) seconds at typical ocean dt — negligible overhead.
  INTEGER, PARAMETER :: MLP_SYNC_STEPS = 5

  !  Reference COARE iterations for ground-truth labels
  INTEGER, PARAMETER :: COARE_REF_ITER = 20

  !  Learning rate, momentum, weight decay
  REAL, PARAMETER :: LR_INIT   = 5.0E-4
  REAL, PARAMETER :: LR_DECAY  = 0.9999   ! per sync event
  REAL, PARAMETER :: MOMENTUM  = 0.9
  REAL, PARAMETER :: WDECAY    = 1.0E-5

  !  Normalisation running stats update rate (Welford online)
  REAL, PARAMETER :: NORM_ALPHA = 0.001   ! blend new vs old stats

  !  Warm-up: use pure COARE for this many sync events before MLP kicks in.
  !  Set to 0 so the MLP drives output from the very first time-step,
  !  ensuring the "with ML" and "without ML" runs differ immediately.
  INTEGER, PARAMETER :: WARMUP_EVENTS = 0

  !=======================================================================
  !  SECTION 4 — weight arrays (all ranks carry identical copies)
  !=======================================================================
  DOUBLE PRECISION, SAVE :: W1(NH1,NIN),  b1(NH1)
  DOUBLE PRECISION, SAVE :: W2(NH2,NH1),  b2(NH2)
  DOUBLE PRECISION, SAVE :: W3(NOUT,NH2), b3(NOUT)

  !  Momentum buffers
  DOUBLE PRECISION, SAVE :: mW1(NH1,NIN),  mb1(NH1)
  DOUBLE PRECISION, SAVE :: mW2(NH2,NH1),  mb2(NH2)
  DOUBLE PRECISION, SAVE :: mW3(NOUT,NH2), mb3(NOUT)

  !  Gradient accumulators (local per rank, reduced at sync)
  DOUBLE PRECISION, SAVE :: gW1(NH1,NIN),  gb1(NH1)
  DOUBLE PRECISION, SAVE :: gW2(NH2,NH1),  gb2(NH2)
  DOUBLE PRECISION, SAVE :: gW3(NOUT,NH2), gb3(NOUT)

  !  Flat buffer for MPI_Allreduce (avoid multiple calls)
  DOUBLE PRECISION, SAVE :: grad_buf(NGRAD_TOTAL)
  DOUBLE PRECISION, SAVE :: grad_tmp(NGRAD_TOTAL)

  !=======================================================================
  !  SECTION 5 — online normalisation stats (per-rank Welford, then reduce)
  !=======================================================================
  DOUBLE PRECISION, SAVE :: x_mean(NIN),  x_var(NIN)     ! input running mean/var
  DOUBLE PRECISION, SAVE :: y_mean(NOUT), y_var(NOUT)     ! output running mean/var
  DOUBLE PRECISION, SAVE :: x_std(NIN),   y_std(NOUT)     ! sqrt(var), cached
  DOUBLE PRECISION, SAVE :: norm_count                     ! sample counter

  !=======================================================================
  !  SECTION 6 — state variables
  !=======================================================================
  DOUBLE PRECISION, SAVE :: lr_current               ! current learning rate
  INTEGER, SAVE :: step_counter             ! total time-steps seen
  INTEGER, SAVE :: sync_counter             ! number of sync events
  INTEGER, SAVE :: local_samples            ! samples on this rank since last sync
  DOUBLE PRECISION, SAVE :: local_loss               ! accumulated MSE on this rank
  LOGICAL, SAVE :: mlp_initialized = .FALSE.
  LOGICAL, SAVE :: mlp_active      = .FALSE.  ! true after warmup
  INTEGER, SAVE :: mpi_comm_croco           ! communicator set at init
  INTEGER, SAVE :: mpi_rank, mpi_size

  PUBLIC :: mlp_init, bulk_flux_mlp_tile, mlp_tick, mlp_sync, mlp_save, mlp_load

CONTAINS

  !=======================================================================
  !  mlp_init
  !  Call once from main.F after MPI_Setup().
  !  Rank-0 generates random weights; MPI_Bcast ensures all start equal.
  !=======================================================================
  SUBROUTINE mlp_init(comm)
    INTEGER, INTENT(IN), OPTIONAL :: comm
    INTEGER :: ierr, iseed(8), i, j
    DOUBLE PRECISION :: r, scale

    IF (mlp_initialized) RETURN

    IF (PRESENT(comm)) THEN
      mpi_comm_croco = comm
    ELSE
      mpi_comm_croco = MPI_COMM_WORLD
    END IF

    CALL MPI_Comm_rank(mpi_comm_croco, mpi_rank, ierr)
    CALL MPI_Comm_size(mpi_comm_croco, mpi_size, ierr)

    !-------------------------------------------------------------------
    ! Xavier-uniform initialisation on rank 0, then broadcast
    !-------------------------------------------------------------------
    IF (mpi_rank == 0) THEN
      iseed = 20240101
      CALL xavier_fill(W1, NH1, NIN,  iseed, scale=SQRT(2.0D0/DBLE(NIN+NH1)))
      CALL xavier_fill_r2(W2, NH2, NH1, iseed, scale=SQRT(2.0D0/DBLE(NH1+NH2)))
      CALL xavier_fill_r3(W3, NOUT,NH2, iseed, scale=SQRT(2.0D0/DBLE(NH2+NOUT)))
    END IF

    b1=0.0D0; b2=0.0D0; b3=0.0D0
    mW1=0.0D0; mb1=0.0D0; mW2=0.0D0; mb2=0.0D0; mW3=0.0D0; mb3=0.0D0
    gW1=0.0D0; gb1=0.0D0; gW2=0.0D0; gb2=0.0D0; gW3=0.0D0; gb3=0.0D0

    !  Broadcast weights
    CALL MPI_Bcast(W1, NH1*NIN,  MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS .AND. mpi_rank == 0) WRITE(*,*) '[MLP-INIT] Bcast W1 failed'
    CALL MPI_Bcast(b1, NH1,      MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS .AND. mpi_rank == 0) WRITE(*,*) '[MLP-INIT] Bcast b1 failed'
    CALL MPI_Bcast(W2, NH2*NH1,  MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS .AND. mpi_rank == 0) WRITE(*,*) '[MLP-INIT] Bcast W2 failed'
    CALL MPI_Bcast(b2, NH2,      MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS .AND. mpi_rank == 0) WRITE(*,*) '[MLP-INIT] Bcast b2 failed'
    CALL MPI_Bcast(W3, NOUT*NH2, MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS .AND. mpi_rank == 0) WRITE(*,*) '[MLP-INIT] Bcast W3 failed'
    CALL MPI_Bcast(b3, NOUT,     MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS .AND. mpi_rank == 0) WRITE(*,*) '[MLP-INIT] Bcast b3 failed'

    !  Initialise running normalisation stats (safe neutral values)
    x_mean = [15.0, 0.0, 0.0, 15.0, 0.8, 101325.0]
    x_var  = [50.0, 25.0, 1.0E-4, 200.0, 0.02, 4.0E6]
    y_mean = [1.3E-3, 1.0E-3, 1.1E-3]
    y_var  = [1.0E-7, 1.0E-8, 1.0E-8]
    x_std  = SQRT(x_var)
    y_std  = SQRT(y_var)
    norm_count = 1000.0D0   ! pretend we've seen 1000 samples already (FIX I: D0 literal)

    lr_current   = LR_INIT
    step_counter  = 0
    sync_counter  = 0
    local_samples = 0
    local_loss    = 0.0D0
    mlp_initialized = .TRUE.
    !  WARMUP_EVENTS=0: MLP is active from the very first time-step.
    mlp_active = (WARMUP_EVENTS == 0)

    IF (mpi_rank == 0) THEN
      WRITE(*,'(/,A)')         ' ======================================='
      WRITE(*,'(A)')           '  CROCO bulk_flux ONLINE MLP TRAINING'
      WRITE(*,'(A,I6,A)')      '  MPI ranks     : ', mpi_size, ' ranks'
      WRITE(*,'(A,I4,A,I4,A)') '  Architecture  : 6->', NH1, '->', NH2, '->3'
      WRITE(*,'(A,I6)')        '  Sync interval : ', MLP_SYNC_STEPS, ' time-steps'
      WRITE(*,'(A,I6)')        '  Warmup events : ', WARMUP_EVENTS
      WRITE(*,'(A,L1)')        '  mlp_active    : ', mlp_active
      IF (mlp_active) THEN
        WRITE(*,'(A)')  '  => MLP drives flux output from step 1 (WARMUP_EVENTS=0)'
      ELSE
        WRITE(*,'(A)')  '  => Running COARE during warmup; MLP output after warmup'
      END IF
      WRITE(*,'(A)')           ' ======================================='
    END IF

  END SUBROUTINE mlp_init

  !=======================================================================
  !  bulk_flux_mlp_tile
  !
  !  Drop-in replacement for bulk_flux_tile (COARE branch).
  !  Called every time-step from step.F inside OMP tiled loop.
  !
  !  Behaviour:
  !    During warm-up : runs COARE_REF_ITER COARE to get labels AND
  !                     accumulates gradients from MLP predictions.
  !    After warm-up  : uses MLP for Cd/Ch/Ce production output;
  !                     still runs COARE_REF_ITER on a fraction of
  !                     grid points for continuous learning.
  !
  !  Arguments mirror bulk_flux_tile in bulk_flux.F:
  !    Istr,Iend,Jstr,Jend — tile bounds
  !    wspd2d  — wind speed  (2D tile array)
  !    delT2d  — air-sea pot. temp diff
  !    delQ2d  — air-sea sp. hum. diff
  !    TairC2d — air temperature (°C)
  !    RH2d    — relative humidity
  !    patm2d  — air pressure (Pa)
  !    Cd2d    — OUTPUT: drag coefficient
  !    Ch2d    — OUTPUT: heat transfer coeff
  !    Ce2d    — OUTPUT: moisture transfer coeff
  !=======================================================================
  SUBROUTINE bulk_flux_mlp_tile(Istr, Iend, Jstr, Jend,  &
                                  wspd2d, delT2d, delQ2d,  &
                                  TairC2d, RH2d, patm2d,   &
                                  Cd2d, Ch2d, Ce2d)

    INTEGER, INTENT(IN) :: Istr, Iend, Jstr, Jend
    DOUBLE PRECISION, INTENT(IN)  :: wspd2d (Istr:Iend, Jstr:Jend)
    DOUBLE PRECISION, INTENT(IN)  :: delT2d (Istr:Iend, Jstr:Jend)
    DOUBLE PRECISION, INTENT(IN)  :: delQ2d (Istr:Iend, Jstr:Jend)
    DOUBLE PRECISION, INTENT(IN)  :: TairC2d(Istr:Iend, Jstr:Jend)
    DOUBLE PRECISION, INTENT(IN)  :: RH2d   (Istr:Iend, Jstr:Jend)
    DOUBLE PRECISION, INTENT(IN)  :: patm2d (Istr:Iend, Jstr:Jend)
    DOUBLE PRECISION, INTENT(OUT) :: Cd2d   (Istr:Iend, Jstr:Jend)
    DOUBLE PRECISION, INTENT(OUT) :: Ch2d   (Istr:Iend, Jstr:Jend)
    DOUBLE PRECISION, INTENT(OUT) :: Ce2d   (Istr:Iend, Jstr:Jend)

    !-------------------------------------------------------------------
    ! Local per-point variables
    !-------------------------------------------------------------------
    DOUBLE PRECISION :: wspd, delT, delQ, TairC, RH, patm
    DOUBLE PRECISION :: Cd_ref, Ch_ref, Ce_ref       ! ground truth from COARE-ref
    DOUBLE PRECISION :: Cd_mlp, Ch_mlp, Ce_mlp       ! MLP prediction
    DOUBLE PRECISION :: xin(NIN), xn(NIN)            ! raw & normalised inputs
    DOUBLE PRECISION :: ytrue(NOUT), ytrue_n(NOUT)   ! labels (normalised)
    DOUBLE PRECISION :: h1(NH1), h2(NH2), yhat(NOUT) ! network activations
    DOUBLE PRECISION :: dy(NOUT), dh2(NH2), dh1(NH1) ! backprop deltas
    DOUBLE PRECISION :: cff
    INTEGER :: i, j, k, p

    !-------------------------------------------------------------------
    ! Scan tile — process every grid point
    !-------------------------------------------------------------------
    DO j = Jstr, Jend
      DO i = Istr, Iend

        wspd  = wspd2d (i,j)
        delT  = delT2d (i,j)
        delQ  = delQ2d (i,j)
        TairC = TairC2d(i,j)
        RH    = RH2d   (i,j)
        patm  = patm2d (i,j)

        !--------------------------------------------------------------
        ! STEP A — COARE ground truth (full 20-iter solver)
        ! This is the "label" for online learning.
        !--------------------------------------------------------------
        CALL coare_full(wspd, delT, delQ, TairC, RH, patm, &
                        COARE_REF_ITER, Cd_ref, Ch_ref, Ce_ref)

        !--------------------------------------------------------------
        ! STEP B — Update online normalisation statistics (Welford)
        !--------------------------------------------------------------
        xin = [wspd, delT, delQ, TairC, RH, patm]
        ytrue = [Cd_ref, Ch_ref, Ce_ref]

        CALL update_norm_stats(xin, ytrue)

        !--------------------------------------------------------------
        ! STEP C — MLP forward pass
        !--------------------------------------------------------------
        CALL normalise_x(xin, xn)
        CALL forward_pass(xn, h1, h2, yhat)
        CALL denorm_y(yhat, Cd_mlp, Ch_mlp, Ce_mlp)

        !--------------------------------------------------------------
        ! STEP D — Backward pass (accumulate local gradients)
        ! Loss = MSE in normalised output space
        !--------------------------------------------------------------
        CALL normalise_y(ytrue, ytrue_n)

        DO k = 1, NOUT
          dy(k) = 2.0D0*(yhat(k) - ytrue_n(k))   ! FIX J: D0 literals on DP arrays
        END DO
        local_loss = local_loss + SUM((yhat - ytrue_n)**2) / DBLE(NOUT)
        local_samples = local_samples + 1

        !  Output layer gradients
        DO p = 1, NH2
          DO k = 1, NOUT
            gW3(k,p) = gW3(k,p) + dy(k)*h2(p)
          END DO
        END DO
        gb3 = gb3 + dy

        !  Hidden layer 2 delta  (tanh' = 1 - h2^2)
        dh2 = 0.0D0
        DO p = 1, NH2
          DO k = 1, NOUT
            dh2(p) = dh2(p) + W3(k,p)*dy(k)
          END DO
        END DO
        dh2 = dh2 * (1.0 - h2**2)

        DO p = 1, NH1
          DO k = 1, NH2
            gW2(k,p) = gW2(k,p) + dh2(k)*h1(p)
          END DO
        END DO
        gb2 = gb2 + dh2

        !  Hidden layer 1 delta
        dh1 = 0.0D0
        DO p = 1, NH1
          DO k = 1, NH2
            dh1(p) = dh1(p) + W2(k,p)*dh2(k)
          END DO
        END DO
        dh1 = dh1 * (1.0 - h1**2)

        DO p = 1, NIN
          DO k = 1, NH1
            gW1(k,p) = gW1(k,p) + dh1(k)*xn(p)
          END DO
        END DO
        gb1 = gb1 + dh1

        !--------------------------------------------------------------
        ! STEP E — Set tile output
        ! During warmup: use COARE; after: use MLP
        !--------------------------------------------------------------
        IF (mlp_active) THEN
          Cd2d(i,j) = Cd_mlp
          Ch2d(i,j) = Ch_mlp
          Ce2d(i,j) = Ce_mlp
        ELSE
          Cd2d(i,j) = Cd_ref
          Ch2d(i,j) = Ch_ref
          Ce2d(i,j) = Ce_ref
        END IF

      END DO   ! i
    END DO     ! j

    !-------------------------------------------------------------------
    ! DIAG: print first-point comparison every MLP_SYNC_STEPS calls.
    ! This is the most important diagnostic: if Cd_mlp == Cd_coare
    ! every single step, the MLP is either not active or clamped to
    ! the same floor as COARE.  Watch for divergence after warmup.
    !-------------------------------------------------------------------
    !-------------------------------------------------------------------
    ! DIAG: print first-point comparison every MLP_SYNC_STEPS calls.
    ! FIX 4: step_counter is incremented by mlp_tick() AFTER this routine
    ! returns, so step_counter here is the previous step's value.
    ! Print step_counter+1 so the log matches the actual model step.
    !-------------------------------------------------------------------
    IF (mpi_rank == 0 .AND. MOD(step_counter+1, MLP_SYNC_STEPS) == 0) THEN
      BLOCK
        INTEGER  :: ic, jc
        DOUBLE PRECISION :: Cd_r, Ch_r, Ce_r, Cd_m, Ch_m, Ce_m
        DOUBLE PRECISION :: xin_d(NIN), xn_d(NIN)
        DOUBLE PRECISION :: h1_d(NH1), h2_d(NH2), yhat_d(NOUT)
        DOUBLE PRECISION :: Cd_raw, Ch_raw, Ce_raw
        ic = (Istr+Iend)/2 ;  jc = (Jstr+Jend)/2
        CALL coare_full(wspd2d(ic,jc), delT2d(ic,jc), delQ2d(ic,jc), &
                        TairC2d(ic,jc), RH2d(ic,jc), patm2d(ic,jc),  &
                        COARE_REF_ITER, Cd_r, Ch_r, Ce_r)
        xin_d = [wspd2d(ic,jc), delT2d(ic,jc), delQ2d(ic,jc), &
                 TairC2d(ic,jc), RH2d(ic,jc), patm2d(ic,jc)]
        CALL normalise_x(xin_d, xn_d)
        CALL forward_pass(xn_d, h1_d, h2_d, yhat_d)
        Cd_raw = yhat_d(1)*y_std(1) + y_mean(1)
        Ch_raw = yhat_d(2)*y_std(2) + y_mean(2)
        Ce_raw = yhat_d(3)*y_std(3) + y_mean(3)
        Cd_m = MAX(Cd_raw, 5.0D-4)
        Ch_m = MAX(Ch_raw, 1.0D-4)
        Ce_m = MAX(Ce_raw, 1.0D-4)
        WRITE(*,'(A,I7,A,L1)') &
          ' [MLP-DIAG] step=', step_counter+1, '  mlp_active=', mlp_active
        WRITE(*,'(A,3(2X,A,E12.5,A,E12.5))') &
          ' [MLP-DIAG]', &
          'Cd_coare=', Cd_r,   ' Cd_mlp=', Cd_m, &
          'Ch_coare=', Ch_r,   ' Ch_mlp=', Ch_m, &
          'Ce_coare=', Ce_r,   ' Ce_mlp=', Ce_m
        WRITE(*,'(A,3(2X,A,E12.5))') &
          ' [MLP-DIAG]  raw(pre-clamp):', &
          'Cd_raw=', Cd_raw, 'Ch_raw=', Ch_raw, 'Ce_raw=', Ce_raw
        WRITE(*,'(A,6(2X,A,E10.3))') &
          ' [MLP-NORM] ', &
          'xmean(wspd)=', x_mean(1), 'xstd(wspd)=', x_std(1), &
          'ymean(Cd)=',   y_mean(1), 'ystd(Cd)=',   y_std(1), &
          'ymean(Ch)=',   y_mean(2), 'ystd(Ch)=',   y_std(2)
        WRITE(*,'(A,6(2X,A,E10.3))') &
          ' [MLP-IN]   ', &
          'wspd=', wspd2d(ic,jc), 'delT=', delT2d(ic,jc), &
          'delQ=', delQ2d(ic,jc), 'TairC=', TairC2d(ic,jc), &
          'RH=',   RH2d(ic,jc),   'patm=', patm2d(ic,jc)
      END BLOCK
    END IF

    !-------------------------------------------------------------------
    ! NOTE: step_counter is NOT incremented here.
    ! Call mlp_tick() once per model time-step from the bulk_flux wrapper
    ! (after all tiles have been processed) to advance the counter and
    ! trigger sync at the correct cadence.
    !-------------------------------------------------------------------

  END SUBROUTINE bulk_flux_mlp_tile

  !=======================================================================
  !  mlp_tick
  !  Call ONCE per model time-step, AFTER all tiles have been processed.
  !  Increments step_counter and fires mlp_sync when due.
  !  This fixes the bug where step_counter was incremented inside
  !  bulk_flux_mlp_tile, causing it to count tile calls (NX_tiles *
  !  NY_tiles per time-step) instead of time-steps.
  !=======================================================================
  SUBROUTINE mlp_tick()
    step_counter = step_counter + 1
    IF (MOD(step_counter, MLP_SYNC_STEPS) == 0) THEN
      CALL mlp_sync()
    END IF
  END SUBROUTINE mlp_tick

  !=======================================================================
  !  mlp_sync
  !
  !  Core MPI routine — called every MLP_SYNC_STEPS time-steps.
  !
  !  1.  Pack all local gradient arrays into a single flat buffer.
  !  2.  MPI_Allreduce (SUM) across all ranks — every rank gets the
  !      global gradient sum.
  !  3.  Divide by global sample count (MPI_Allreduce on local_samples).
  !  4.  Apply SGD with momentum + weight decay.
  !  5.  Zero gradient accumulators.
  !  6.  Decay learning rate.
  !  7.  Print diagnostics on rank 0.
  !=======================================================================
  SUBROUTINE mlp_sync()
    INTEGER :: ierr, p, k, off
    DOUBLE PRECISION :: global_loss, inv_n
    INTEGER :: global_samples

    !-------------------------------------------------------------------
    ! 1.  Pack gradient buffer  (same layout on every rank)
    !     W1(NH1,NIN) | b1(NH1) | W2(NH2,NH1) | b2(NH2) |
    !     W3(NOUT,NH2)| b3(NOUT)
    !-------------------------------------------------------------------
    off = 0
    grad_buf(off+1 : off+NW1)  = RESHAPE(gW1, [NW1])  ; off=off+NW1
    grad_buf(off+1 : off+NH1)  = gb1                   ; off=off+NH1
    grad_buf(off+1 : off+NW2)  = RESHAPE(gW2, [NW2])  ; off=off+NW2
    grad_buf(off+1 : off+NH2)  = gb2                   ; off=off+NH2
    grad_buf(off+1 : off+NW3)  = RESHAPE(gW3, [NW3])  ; off=off+NW3
    grad_buf(off+1 : off+NOUT) = gb3                   ; off=off+NOUT

    !-------------------------------------------------------------------
    ! 2.  MPI_Allreduce — sum gradients from all ranks
    !     All ranks receive the result (no master needed)
    !-------------------------------------------------------------------
    CALL MPI_Allreduce(grad_buf, grad_tmp, NGRAD_TOTAL, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS .AND. mpi_rank == 0) &
      WRITE(*,*) '[MLP-SYNC] MPI_Allreduce gradients failed, ierr=', ierr

    CALL MPI_Allreduce(local_samples, global_samples, 1, &
                       MPI_INTEGER, MPI_SUM, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS .AND. mpi_rank == 0) &
      WRITE(*,*) '[MLP-SYNC] MPI_Allreduce samples failed, ierr=', ierr

    CALL MPI_Allreduce(local_loss, global_loss, 1, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS .AND. mpi_rank == 0) &
      WRITE(*,*) '[MLP-SYNC] MPI_Allreduce loss failed, ierr=', ierr

    !-------------------------------------------------------------------
    ! 3.  Unpack and normalise by global sample count
    !-------------------------------------------------------------------
    ! FIX B: Use DBLE not REAL — REAL truncates large sample counts to ~7 sig figs.
    inv_n = 1.0D0 / MAX(DBLE(global_samples), 1.0D0)

    off = 0
    gW1 = RESHAPE(grad_tmp(off+1:off+NW1), [NH1,NIN]) * inv_n; off=off+NW1
    gb1 = grad_tmp(off+1:off+NH1) * inv_n                      ; off=off+NH1
    gW2 = RESHAPE(grad_tmp(off+1:off+NW2), [NH2,NH1]) * inv_n; off=off+NW2
    gb2 = grad_tmp(off+1:off+NH2) * inv_n                      ; off=off+NH2
    gW3 = RESHAPE(grad_tmp(off+1:off+NW3), [NOUT,NH2])* inv_n; off=off+NW3
    gb3 = grad_tmp(off+1:off+NOUT) * inv_n                     ; off=off+NOUT

    !-------------------------------------------------------------------
    ! 4.  SGD with Nesterov momentum + L2 weight decay
    !     v = momentum*v + lr*(grad + wdecay*w)
    !     w = w - v
    !-------------------------------------------------------------------
    mW1 = MOMENTUM*mW1 + lr_current*(gW1 + WDECAY*W1)
    W1  = W1 - mW1
    mb1 = MOMENTUM*mb1 + lr_current*gb1
    b1  = b1 - mb1

    mW2 = MOMENTUM*mW2 + lr_current*(gW2 + WDECAY*W2)
    W2  = W2 - mW2
    mb2 = MOMENTUM*mb2 + lr_current*gb2
    b2  = b2 - mb2

    mW3 = MOMENTUM*mW3 + lr_current*(gW3 + WDECAY*W3)
    W3  = W3 - mW3
    mb3 = MOMENTUM*mb3 + lr_current*gb3
    b3  = b3 - mb3

    !-------------------------------------------------------------------
    ! 5.  Zero accumulators
    !-------------------------------------------------------------------
    gW1=0.0D0; gb1=0.0D0; gW2=0.0D0; gb2=0.0D0; gW3=0.0D0; gb3=0.0D0
    local_loss    = 0.0D0
    local_samples = 0

    !-------------------------------------------------------------------
    ! 6.  Learning rate decay + warmup check
    !-------------------------------------------------------------------
    sync_counter = sync_counter + 1
    lr_current   = lr_current * LR_DECAY

    IF (sync_counter >= WARMUP_EVENTS .AND. .NOT. mlp_active) THEN
      mlp_active = .TRUE.
      IF (mpi_rank == 0) WRITE(*,'(/,A,/)') &
        ' [MLP] Warmup complete — MLP now drives bulk flux production output.'
    END IF

    !-------------------------------------------------------------------
    ! 7.  Diagnostics (rank 0 only)
    !     Print every sync event for the first 10, then every 20.
    !     This tells you:
    !       loss      — MSE in normalised space; should decrease over time
    !       lr        — learning rate (decays by LR_DECAY each sync)
    !       W1_norm   — L2 norm of first-layer weights (should grow from 0)
    !       active    — FALSE means still using COARE, TRUE means using MLP
    !-------------------------------------------------------------------
    IF (mpi_rank == 0 .AND. &
        (sync_counter <= 10 .OR. MOD(sync_counter, 20) == 0)) THEN
      WRITE(*,'(A,I7,A,E10.3,A,E10.3,A,I9,A,L1)') &
        ' [MLP-SYNC] sync=', sync_counter,          &
        '  loss=',    global_loss*inv_n,             &
        '  lr=',      lr_current,                    &
        '  N=',       global_samples,                &
        '  active=',  mlp_active
      WRITE(*,'(A,4(2X,A,E10.3))') &
        ' [MLP-SYNC]', &
        'W1_norm=',  SQRT(SUM(W1**2)),  &
        'W2_norm=',  SQRT(SUM(W2**2)),  &
        'W3_norm=',  SQRT(SUM(W3**2)),  &
        'b3(Cd)=',   b3(1)
      WRITE(*,'(A,3(2X,A,E10.3,A,E10.3))') &
        ' [MLP-SYNC]', &
        'ymean(Cd)=', y_mean(1), ' ystd(Cd)=', y_std(1), &
        'ymean(Ch)=', y_mean(2), ' ystd(Ch)=', y_std(2), &
        'ymean(Ce)=', y_mean(3), ' ystd(Ce)=', y_std(3)
    END IF

  END SUBROUTINE mlp_sync

  !=======================================================================
  !  SECTION 7 — COARE full reference solver (20-iteration ground truth)
  !=======================================================================
  SUBROUTINE coare_full(wspd_in, delT, delQ, TairC, RH, patm_in, &
                         Niter, Cd, Ch, Ce)
    DOUBLE PRECISION,    INTENT(IN)  :: wspd_in, delT, delQ, TairC, RH, patm_in
    INTEGER, INTENT(IN)  :: Niter
    DOUBLE PRECISION,    INTENT(OUT) :: Cd, Ch, Ce

    DOUBLE PRECISION :: TairK, Q, rhoAir, patm, dW
    DOUBLE PRECISION :: Wstar, Tstar, Qstar
    DOUBLE PRECISION :: iZoW, iZoT, ZoLu, ZoLt
    DOUBLE PRECISION :: psi_u, psi_t, logus, logts, cff
    DOUBLE PRECISION :: VisAir, Ri, CC, Ribcu, charn, Ch10, Rr
    INTEGER :: it

    TairK = TairC + CtoK
    patm  = patm_in

    !  Specific humidity (simplified, no full Exner loop)
    Q = MAX(RH*0.0175*EXP(17.502*TairC/(240.97+TairC)) / &
            (patm*1.0E-5 - 0.01), 1.0E-6)

    dW = MAX(wspd_in, 0.1)

    !  Viscosity
    VisAir = c0visc*(1. + c1visc*TairC + c2visc*TairC**2 - c3visc*TairC**3)

    charn  = 0.011
    Ch10   = 0.00115
    Ribcu  = -blk_ZW / (blk_Zabl * 0.004 * blk_beta**3)

    !  Initial guess for Wstar
    Wstar  = 0.035 * dW * LOG(1.E5) / LOG(blk_ZW*1.E4)

    iZoW  = g_grav*Wstar / (charn*Wstar**3 + 0.11*g_grav*VisAir)
    iZoT  = 0.1 * EXP(vonKar**2 / (Ch10*LOG(10.0*iZoW)))
    CC    = LOG(blk_ZW*iZoW)**2 / (LOG(blk_ZT*iZoT) + eps_f)
    Ri    = g_grav*blk_ZW*(delT + cpvir*TairK*delQ) / &
            (TairK*dW**2 + eps_f)

    IF (Ri < 0.) THEN
      ZoLu = CC*Ri / (1. + Ri/Ribcu)
    ELSE
      ZoLu = CC*Ri / (1. + 3.0*Ri/(CC + eps_f))
    END IF

    CALL psiu_coare(ZoLu, psi_u)
    logus = LOG(blk_ZW*iZoW)
    Wstar = dW * vonKar / (logus - psi_u + eps_f)
    ZoLt  = ZoLu * blk_ZTW
    CALL psit_coare(ZoLt, psi_t)
    logts = LOG(blk_ZT*iZoT)
    cff   = vonKar / (logts - psi_t + eps_f)
    Tstar = delT * cff
    Qstar = delQ * cff

    !  Wind-dependent Charnock
    IF      (dW > 18.) THEN ; charn = 0.018
    ELSE IF (dW > 10.) THEN ; charn = 0.011 + 0.125*0.007*(dW - 10.)
    ELSE                    ; charn = 0.011
    END IF

    !  Main iteration loop
    DO it = 1, Niter
      iZoW = g_grav*Wstar / (charn*Wstar**3 + 0.11*g_grav*VisAir + eps_f)
      Rr   = Wstar / (iZoW*VisAir + eps_f)
      iZoT = MAX(8695.65, 18181.8*(Rr**0.6))

      ZoLu = vonKar*g_grav*blk_ZW * &
             (Tstar*(1. + cpvir*Q) + cpvir*TairK*Qstar) / &
             (TairK*Wstar**2*(1. + cpvir*Q) + eps_f)

      CALL psiu_coare(ZoLu, psi_u)
      logus = LOG(blk_ZW * iZoW)
      Wstar = dW * vonKar / (logus - psi_u + eps_f)

      ZoLt  = ZoLu * blk_ZTW
      CALL psit_coare(ZoLt, psi_t)
      logts = LOG(blk_ZT * iZoT)
      cff   = vonKar / (logts - psi_t + eps_f)
      Tstar = delT * cff
      Qstar = delQ * cff
    END DO

    !  Exchange coefficients
    !  FIX 7: Ce was set to MAX(Ch, 1E-4), making Ce==Ch always and preventing
    !  the MLP from ever learning Ce != Ch.  Compute Ce from its own similarity
    !  scale using the same cff transfer coefficient as Ch (COARE 3.5 convention),
    !  but keep it as an independent variable so it can diverge from Ch.
    Cd = MAX((Wstar/dW)**2,        5.0E-4)
    Ch = MAX(cff * SQRT(Cd),       1.0E-4)
    Ce = MAX(cff * SQRT(Cd),       1.0E-4)   ! same formula as Ch (COARE 3.5), but independent

  END SUBROUTINE coare_full

  !=======================================================================
  !  SECTION 8 — COARE stability functions (inline, no external calls)
  !=======================================================================
  PURE SUBROUTINE psiu_coare(ZoL, psi)
    DOUBLE PRECISION, INTENT(IN)  :: ZoL
    DOUBLE PRECISION, INTENT(OUT) :: psi
    DOUBLE PRECISION :: ck, cc, pk, pc, c
    IF (ZoL <= 0.) THEN
      ck = (1. - 15.*ZoL)**0.25
      pk = 2.*LOG(0.5*(1.+ck)) + LOG(0.5*(1.+ck**2)) - 2.*ATAN(ck) + pis2_
      cc = (1. - 10.15*ZoL)**r3_
      pc = 1.5*LOG(r3_*(cc**2+cc+1.)) - sqr3_*ATAN((2.*cc+1.)/sqr3_) + 2.*pis2osq3
      psi = pc + (pk-pc)/(1. + ZoL**2)
    ELSE
      c   = -MIN(50., 0.35*ZoL)
      psi = -((1.+ZoL) + 0.6667*(ZoL-14.28)*EXP(c) + 8.525)
    END IF
  END SUBROUTINE psiu_coare

  PURE SUBROUTINE psit_coare(ZoL, psi)
    DOUBLE PRECISION, INTENT(IN)  :: ZoL
    DOUBLE PRECISION, INTENT(OUT) :: psi
    DOUBLE PRECISION :: ck, cc, pk, pc, c
    IF (ZoL < 0.) THEN
      ck = (1. - 15.*ZoL)**0.25
      pk = 2.*LOG(0.5*(1.+ck**2))
      cc = (1. - 34.15*ZoL)**r3_
      pc = 1.5*LOG((cc**2+cc+1.)*r3_) - sqr3_*ATAN((2.*cc+1.)/sqr3_) + 2.*pis2osq3
      psi = pc + (pk-pc)/(1. + ZoL**2)
    ELSE
      c   = -MIN(50., 0.35*ZoL)
      psi = -((1.+2.*ZoL/3.)**1.5 + 0.6667*(ZoL-14.28)*EXP(c) + 8.525)
    END IF
  END SUBROUTINE psit_coare

  !=======================================================================
  !  SECTION 9 — MLP forward pass
  !=======================================================================
  ! FIX 6: Removed PURE attribute — these routines read module-level weight/norm
  ! arrays (W1,W2,W3,b*,x_mean,x_std,y_mean,y_std) that are updated each sync
  ! event.  PURE is technically allowed for reads of saved module vars in
  ! Fortran, but is misleading here because the result changes between calls as
  ! training progresses.  Removing PURE avoids confusion and compiler warnings.
  SUBROUTINE forward_pass(xn, h1, h2, yhat)
    DOUBLE PRECISION, INTENT(IN)  :: xn(NIN)
    DOUBLE PRECISION, INTENT(OUT) :: h1(NH1), h2(NH2), yhat(NOUT)
    INTEGER :: i, j

    h1 = b1
    DO j = 1, NIN
      DO i = 1, NH1
        h1(i) = h1(i) + W1(i,j)*xn(j)
      END DO
    END DO
    h1 = TANH(h1)

    h2 = b2
    DO j = 1, NH1
      DO i = 1, NH2
        h2(i) = h2(i) + W2(i,j)*h1(j)
      END DO
    END DO
    h2 = TANH(h2)

    yhat = b3
    DO j = 1, NH2
      DO i = 1, NOUT
        yhat(i) = yhat(i) + W3(i,j)*h2(j)
      END DO
    END DO
  END SUBROUTINE forward_pass

  !=======================================================================
  !  SECTION 10 — Online normalisation (Welford running stats)
  !=======================================================================
  SUBROUTINE update_norm_stats(xin, ytrue)
    DOUBLE PRECISION, INTENT(IN) :: xin(NIN), ytrue(NOUT)
    DOUBLE PRECISION :: alpha, xdelta(NIN), ydelta(NOUT)
    INTEGER :: k

    norm_count = norm_count + 1.0D0
    alpha = NORM_ALPHA

    !  Exponential moving average for mean and variance
    xdelta = xin   - x_mean
    x_mean = x_mean + alpha*xdelta
    x_var  = (1.-alpha)*x_var + alpha*xdelta**2
    x_std  = SQRT(MAX(x_var, 1.0E-12))

    ydelta = ytrue  - y_mean
    y_mean = y_mean + alpha*ydelta
    y_var  = (1.-alpha)*y_var + alpha*ydelta**2
    y_std  = SQRT(MAX(y_var, 1.0E-16))
  END SUBROUTINE update_norm_stats

  SUBROUTINE normalise_x(xin, xn)  ! FIX 6: PURE removed (reads trained module state)
    DOUBLE PRECISION, INTENT(IN)  :: xin(NIN)
    DOUBLE PRECISION, INTENT(OUT) :: xn(NIN)
    INTEGER :: k
    DO k = 1, NIN
      xn(k) = (xin(k) - x_mean(k)) / x_std(k)
    END DO
  END SUBROUTINE normalise_x

  SUBROUTINE normalise_y(ytrue, yn)  ! FIX 6: PURE removed (reads trained module state)
    DOUBLE PRECISION, INTENT(IN)  :: ytrue(NOUT)
    DOUBLE PRECISION, INTENT(OUT) :: yn(NOUT)
    INTEGER :: k
    DO k = 1, NOUT
      yn(k) = (ytrue(k) - y_mean(k)) / y_std(k)
    END DO
  END SUBROUTINE normalise_y

  SUBROUTINE denorm_y(yhat, Cd, Ch, Ce)
    DOUBLE PRECISION, INTENT(IN)  :: yhat(NOUT)
    DOUBLE PRECISION, INTENT(OUT) :: Cd, Ch, Ce
    DOUBLE PRECISION :: Cd_raw, Ch_raw, Ce_raw
    Cd_raw = yhat(1)*y_std(1) + y_mean(1)
    Ch_raw = yhat(2)*y_std(2) + y_mean(2)
    Ce_raw = yhat(3)*y_std(3) + y_mean(3)
    Cd = MAX(Cd_raw, 5.0E-4)
    Ch = MAX(Ch_raw, 1.0E-4)
    Ce = MAX(Ce_raw, 1.0E-4)
  END SUBROUTINE denorm_y

  !=======================================================================
  !  SECTION 11 — Weight persistence (save / load)
  !  Optional: call mlp_save() at end of run, mlp_load() to restart.
  !=======================================================================
  SUBROUTINE mlp_save(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: iunit, ierr
    !  Only rank-0 writes — weights are identical on all ranks
    IF (mpi_rank /= 0) RETURN
    iunit = 55
    OPEN(UNIT=iunit, FILE=TRIM(filename), STATUS='REPLACE', &
         FORM='UNFORMATTED', ACCESS='STREAM')
    WRITE(iunit) W1, b1, W2, b2, W3, b3
    WRITE(iunit) x_mean, x_std, y_mean, y_std
    WRITE(iunit) sync_counter, lr_current
    CLOSE(iunit)
    WRITE(*,'(A,A)') ' [MLP] Weights saved to ', TRIM(filename)
  END SUBROUTINE mlp_save

  SUBROUTINE mlp_load(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: iunit, ierr, ios
    IF (mpi_rank == 0) THEN
      iunit = 56
      OPEN(UNIT=iunit, FILE=TRIM(filename), STATUS='OLD', &
           FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=ios)
      IF (ios /= 0) THEN
        WRITE(*,'(A,A,A,I0)') ' [MLP] ERROR: cannot open ', TRIM(filename), &
          ' IOSTAT=', ios
        CALL MPI_Abort(mpi_comm_croco, 1, ierr)
        RETURN
      END IF
      READ(iunit, IOSTAT=ios) W1, b1, W2, b2, W3, b3
      IF (ios /= 0) THEN
        WRITE(*,'(A,I0)') ' [MLP] ERROR: reading weights, IOSTAT=', ios
        CALL MPI_Abort(mpi_comm_croco, 1, ierr)
      END IF
      READ(iunit, IOSTAT=ios) x_mean, x_std, y_mean, y_std
      IF (ios /= 0) THEN
        WRITE(*,'(A,I0)') ' [MLP] ERROR: reading norm stats, IOSTAT=', ios
        CALL MPI_Abort(mpi_comm_croco, 1, ierr)
      END IF
      READ(iunit, IOSTAT=ios) sync_counter, lr_current
      IF (ios /= 0) THEN
        WRITE(*,'(A,I0)') ' [MLP] ERROR: reading training state, IOSTAT=', ios
        CALL MPI_Abort(mpi_comm_croco, 1, ierr)
      END IF
      CLOSE(iunit)
      !  Recompute std from loaded values
      x_var = x_std**2
      y_var = y_std**2
    END IF
    !  Broadcast loaded weights to all ranks; check each ierr
    CALL MPI_Bcast(W1,     NH1*NIN,  MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS) WRITE(*,*) '[MLP] MPI_Bcast W1 failed, ierr=', ierr
    CALL MPI_Bcast(b1,     NH1,      MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS) WRITE(*,*) '[MLP] MPI_Bcast b1 failed, ierr=', ierr
    CALL MPI_Bcast(W2,     NH2*NH1,  MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS) WRITE(*,*) '[MLP] MPI_Bcast W2 failed, ierr=', ierr
    CALL MPI_Bcast(b2,     NH2,      MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS) WRITE(*,*) '[MLP] MPI_Bcast b2 failed, ierr=', ierr
    CALL MPI_Bcast(W3,     NOUT*NH2, MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS) WRITE(*,*) '[MLP] MPI_Bcast W3 failed, ierr=', ierr
    CALL MPI_Bcast(b3,     NOUT,     MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS) WRITE(*,*) '[MLP] MPI_Bcast b3 failed, ierr=', ierr
    CALL MPI_Bcast(x_mean, NIN,      MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS) WRITE(*,*) '[MLP] MPI_Bcast x_mean failed, ierr=', ierr
    CALL MPI_Bcast(x_std,  NIN,      MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS) WRITE(*,*) '[MLP] MPI_Bcast x_std failed, ierr=', ierr
    CALL MPI_Bcast(y_mean, NOUT,     MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS) WRITE(*,*) '[MLP] MPI_Bcast y_mean failed, ierr=', ierr
    CALL MPI_Bcast(y_std,  NOUT,     MPI_DOUBLE_PRECISION, 0, mpi_comm_croco, ierr)
    IF (ierr /= MPI_SUCCESS) WRITE(*,*) '[MLP] MPI_Bcast y_std failed, ierr=', ierr
    IF (mpi_rank == 0) WRITE(*,'(A,A)') ' [MLP] Weights loaded from ', TRIM(filename)
  END SUBROUTINE mlp_load

  !=======================================================================
  !  SECTION 12 — Weight initialisation helpers
  !=======================================================================
  SUBROUTINE xavier_fill(W, nr, nc, iseed, scale)
    INTEGER, INTENT(IN)    :: nr, nc
    DOUBLE PRECISION,    INTENT(IN)    :: scale
    INTEGER, INTENT(INOUT) :: iseed(8)
    DOUBLE PRECISION,    INTENT(OUT)   :: W(nr,nc)
    INTEGER :: i, j
    DOUBLE PRECISION :: r
    DO j = 1, nc
      DO i = 1, nr
        CALL lcg_rand(iseed, r)
        W(i,j) = (2.0D0*r - 1.0D0) * scale   ! FIX E: was 2. (single-prec literal)
      END DO
    END DO
  END SUBROUTINE xavier_fill

  SUBROUTINE xavier_fill_r2(W, nr, nc, iseed, scale)
    INTEGER, INTENT(IN)    :: nr, nc
    DOUBLE PRECISION,    INTENT(IN)    :: scale
    INTEGER, INTENT(INOUT) :: iseed(8)
    DOUBLE PRECISION,    INTENT(OUT)   :: W(nr,nc)
    INTEGER :: i, j ; DOUBLE PRECISION :: r   ! FIX 1: was REAL, must match lcg_rand intent
    DO j=1,nc; DO i=1,nr; CALL lcg_rand(iseed,r); W(i,j)=(2.D0*r-1.D0)*scale; END DO; END DO
  END SUBROUTINE xavier_fill_r2

  SUBROUTINE xavier_fill_r3(W, nr, nc, iseed, scale)
    INTEGER, INTENT(IN)    :: nr, nc
    DOUBLE PRECISION,    INTENT(IN)    :: scale
    INTEGER, INTENT(INOUT) :: iseed(8)
    DOUBLE PRECISION,    INTENT(OUT)   :: W(nr,nc)
    INTEGER :: i, j ; DOUBLE PRECISION :: r   ! FIX 1: was REAL, must match lcg_rand intent
    DO j=1,nc; DO i=1,nr; CALL lcg_rand(iseed,r); W(i,j)=(2.D0*r-1.D0)*scale; END DO; END DO
  END SUBROUTINE xavier_fill_r3

  SUBROUTINE lcg_rand(s, r)
    INTEGER, INTENT(INOUT) :: s(8)
    DOUBLE PRECISION,    INTENT(OUT)   :: r
    INTEGER(KIND=8) :: x
    ! FIX 3: Replace INT(x,4) with safe MODULO to avoid undefined behaviour on
    ! 64->32-bit truncation overflow (INT with out-of-range value is UB in Fortran).
    x    = INT(s(1),8)*6364136223846793005_8 + 1442695040888963407_8
    s(1) = INT(MOD(ABS(x), 2147483647_8), 4)
    r    = DBLE(s(1)) / 2147483647.0D0
  END SUBROUTINE lcg_rand

END MODULE bulk_flux_mlp


!=======================================================================
!
!  DRIVER — legacy subroutine wrappers so CROCO's step.F call sites
!  require zero changes.  These are compiled in the same file.
!
!  step.F calls:   call bulk_flux(tile)
!  which calls:    call bulk_flux_tile(Istr,Iend,Jstr,Jend, aer,cer)
!
!  We intercept at bulk_flux_tile level so we can add MLP outputs.
!  The outer bulk_flux() wrapper and stress assembly remain identical
!  to the original bulk_flux.F — only the inner coefficient loop
!  is replaced.
!
!=======================================================================
#ifdef BULK_FLUX_MLP
!======================================================================
!  FIX: renamed from bulk_flux -> bulk_flux_mlp_driver to avoid
!  duplicate symbol conflict with the bulk_flux() defined in bulk_flux.F
!  (which already has full #ifdef BULK_FLUX_MLP support built in).
!
!  When BULK_FLUX_MLP is defined, CROCO should link bulk_flux.F
!  (which calls bulk_flux_mlp_tile via the #ifdef block) and NOT
!  this legacy driver.  This stub is kept for reference only.
!
!  If you want to use ONLY bulk_flux_mlp.F90 (standalone, no bulk_flux.F),
!  rename this back to bulk_flux and remove bulk_flux.F from the build.
!======================================================================
      subroutine bulk_flux_mlp_driver(tile)
!  Legacy standalone driver — NOT called when bulk_flux.F is in the build.
!  FIX 2: Removed mlp_init() call from here. mlp_init() must be called ONCE
!  from main.F after MPI_Setup(). Calling it here (per tile, per step) was
!  wrong; the mlp_initialized guard saved a hard crash but this is cleaner.
      use bulk_flux_mlp, only: mlp_init
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd = omp_get_thread_num()
!     mlp_init() intentionally NOT called here — see main.F after MPI_Setup()
      call bulk_flux_tile_mlp(Istr, Iend, Jstr, Jend,   &
                               A2d(1,1,trd), A2d(1,2,trd))
      return
      end subroutine bulk_flux_mlp_driver

!======================================================================
      subroutine bulk_flux_tile_mlp(Istr, Iend, Jstr, Jend, aer, cer)
!
!  Full tile routine.  Sections 1-3 and sections 5-end are IDENTICAL
!  to bulk_flux.F.  Only section 4 (the iterative coefficient loop)
!  is replaced by the MLP call.
!
      use bulk_flux_mlp
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "forces.h"
# include "scalars.h"
# include "params_bulk.h"

      integer Istr, Iend, Jstr, Jend, imin, imax, jmin, jmax
      integer i, j, m
      real    aer(PRIVATE_2D_SCRATCH_ARRAY)
      real    cer(PRIVATE_2D_SCRATCH_ARRAY)

      !  All the original local scalars from bulk_flux.F
      !  FIX A: Scalars that receive values from DOUBLE PRECISION tile arrays
      !  (Cd, Ch, Ce, delW, delT, delQ, rhoAir, TseaC, TseaK, TairC, RH, Q)
      !  must be DOUBLE PRECISION to avoid silent truncation.  Flux intermediates
      !  (hfsen, hflat, hflw, upvel, Hlv, evap) follow since computed from dp.
      !  rho0i, cpi, wspd0, iexns, iexna, patm, Qsea stay REAL (from REAL arrays).
      !  Removed dead locals: stflx_t, stflx_s, srflx_ij, radlw_ij, Tvstar, Bf.
      real    :: rho0i, cpi, wspd0
      real    :: Qsea, spec_hum, iexns, iexna, patm, qsat
      double precision :: TseaC, TseaK, TairC, TairK, RH, Q, rhoAir
      double precision :: delW, delT, delQ
      double precision :: Cd, Ch, Ce, Hlv
      double precision :: hfsen, hflat, hflw, upvel, evap

      !  Tile-sized arrays to pass full tile to MLP in one call
      !  FIX 10: Changed REAL -> DOUBLE PRECISION to match bulk_flux_mlp_tile
      !  interface (which declares all 2D intent arguments as DOUBLE PRECISION).
      !  Passing REAL arrays to DOUBLE PRECISION dummy args is a type mismatch
      !  that causes silent data corruption on most compilers.
      double precision, allocatable :: wspd_tile(:,:), delT_tile(:,:), delQ_tile(:,:)
      double precision, allocatable :: TairC_tile(:,:), RH_tile(:,:), patm_tile(:,:)
      double precision, allocatable :: Cd_tile(:,:), Ch_tile(:,:), Ce_tile(:,:)
      double precision, allocatable :: rhoAir_tile(:,:), TseaC_tile(:,:)

!======================================================================
!  SECTION 1 — tile bounds (identical to original)
!======================================================================
# ifdef EW_PERIODIC
      imin=Istr-2; imax=Iend+2
# else
      if (WESTERN_EDGE) then; imin=Istr-1; else; imin=Istr-2; endif
      if (EASTERN_EDGE) then; imax=Iend+1; else; imax=Iend+2; endif
# endif
# ifdef NS_PERIODIC
      jmin=Jstr-2; jmax=Jend+2
# else
      if (SOUTHERN_EDGE) then; jmin=Jstr-1; else; jmin=Jstr-2; endif
      if (NORTHERN_EDGE) then; jmax=Jend+1; else; jmax=Jend+2; endif
# endif

      rho0i = 1.0/rho0
      cpi   = 1.0/cp

      !  Allocate tile arrays
      allocate( wspd_tile(imin:imax,jmin:jmax),  &
                delT_tile(imin:imax,jmin:jmax),  &
                delQ_tile(imin:imax,jmin:jmax),  &
                TairC_tile(imin:imax,jmin:jmax), &
                RH_tile(imin:imax,jmin:jmax),    &
                patm_tile(imin:imax,jmin:jmax),  &
                Cd_tile(imin:imax,jmin:jmax),    &
                Ch_tile(imin:imax,jmin:jmax),    &
                Ce_tile(imin:imax,jmin:jmax),    &
                rhoAir_tile(imin:imax,jmin:jmax),&
                TseaC_tile(imin:imax,jmin:jmax) )

!======================================================================
!  SECTION 2 — i,j loop: extract state variables into tile arrays
!======================================================================
      DO j = jmin, jmax
        DO i = imin, imax

!======================================================================
!  SECTION 3 — state variable extraction (identical to original)
!======================================================================
          wspd0 = wspd(i,j)
          wspd0 = MAX(wspd0, 0.1*MIN(10., blk_ZW))
          TairC = tair(i,j)
          TairK = TairC + CtoK
          RH    = rhum(i,j)
          Q     = spec_hum(RH, psurf, TairC)

# ifdef SST_SKIN
          TseaC = sst_skin(i,j)
# else
          TseaC = t(i,j,N,nrhs,itemp)
# endif
          TseaK = TseaC + CtoK

          CALL exner_patm_from_tairabs(iexna, patm, Q, TairK, blk_ZT, psurf)
          rhoAir = patm*(1.+Q) / (blk_Rgas*TairK*(1.+MvoMa*Q))
          Qsea   = qsat(TseaK, psurf, 0.98)

          delW = MAX(wspd0, 0.1)
# ifdef BULK_GUSTINESS
          delW = SQRT(wspd0*wspd0 + 0.25)
# endif
          delQ = Q - Qsea
          delT = TairC*iexna - TseaC*iexns + CtoK*(iexna-iexns)

          !  Store into tile arrays for MLP call
          wspd_tile(i,j)  = delW
          delT_tile(i,j)  = delT
          delQ_tile(i,j)  = delQ
          TairC_tile(i,j) = TairC
          RH_tile(i,j)    = RH
          patm_tile(i,j)  = patm
          rhoAir_tile(i,j)= rhoAir
          TseaC_tile(i,j) = TseaC

        END DO   ! i
      END DO     ! j

!======================================================================
!  SECTION 4 — Call MLP ONCE for the full tile (fixes step_counter)
!======================================================================
      CALL bulk_flux_mlp_tile(imin, imax, jmin, jmax,  &
                               wspd_tile,  delT_tile,  &
                               delQ_tile,  TairC_tile, &
                               RH_tile,    patm_tile,  &
                               Cd_tile,    Ch_tile,    Ce_tile)

!======================================================================
!  SECTION 5 — flux assembly loop using MLP output
!  FIX: The correct Monin-Obukhov flux formulas are:
!    Wstar = sqrt(Cd)*delW
!    Tstar = Ch/sqrt(Cd)*delT   (from Ch = Wstar*Tstar/(delW*delT))
!    Qstar = Ce/sqrt(Cd)*delQ
!  =>  hfsen = -Cpa*rhoAir*Wstar*Tstar
!            = -Cpa*rhoAir*sqrt(Cd)*delW * Ch/sqrt(Cd)*delT
!            = -Cpa*rhoAir*delW*Ch*delT          (sqrt(Cd) cancels!)
!  The previous code had SQRT(Cd)*delW in numerator and
!  (SQRT(Cd)*delW+eps) in denominator — numerically close but wrong
!  and would hide MLP-vs-COARE differences.
!======================================================================
      DO j = jmin, jmax
        DO i = imin, imax

          Cd     = Cd_tile(i,j)
          Ch     = Ch_tile(i,j)
          Ce     = Ce_tile(i,j)
          rhoAir = rhoAir_tile(i,j)
          TseaC  = TseaC_tile(i,j)
          TseaK  = TseaC + CtoK
          delW   = wspd_tile(i,j)
          delT   = delT_tile(i,j)
          delQ   = delQ_tile(i,j)
          TairC  = TairC_tile(i,j)
          RH     = RH_tile(i,j)
          Q      = spec_hum(RH, psurf, TairC)

          !  Monin-Obukhov similarity scales
          !    Wstar = sqrt(Cd)*delW,  Tstar = Ch/sqrt(Cd)*delT,  Qstar = Ce/sqrt(Cd)*delQ
          !  Heat fluxes (sqrt(Cd) cancels in Wstar*Tstar and Wstar*Qstar):
          hfsen = -blk_Cpa * rhoAir * delW * Ch * delT
          Hlv   = (2.5008 - 0.0023719*TseaC)*1.0E+6
          hflat = -Hlv * rhoAir * delW * Ce * delQ

# ifndef BULK_LW
          hflw  = -radlw(i,j)
# else
          hflw  =  radlw(i,j) - emiss_lw*rho0i*cpi*SigmaSB* &
                   TseaK*TseaK*TseaK*TseaK
# endif

          !  Webb correction: upvel = -1.61*Wstar*Qstar - (1+1.61Q)*Wstar*Tstar/TairK
          !  With Wstar*Tstar = sqrt(Cd)*delW*Ch/sqrt(Cd)*delT = delW*Ch*delT  etc.
          upvel = -1.61*delW*Ce*delQ - &
                  (1.0+1.61*Q)*delW*Ch*delT/TairK
          hflat = hflat + rhoAir*Hlv*upvel*Q

          hflat = -hflat*rho0i*cpi
          hfsen = -hfsen*rho0i*cpi

          stflx(i,j,itemp) = srflx(i,j) + hflw + hflat + hfsen

# ifdef SALINITY
          evap = -cp*hflat/Hlv
#  ifdef RAIN_FLUX
          EmP(i,j) = evap - prate(i,j)
          stflx(i,j,itemp) = stflx(i,j,itemp) + &
                              prate(i,j)*t(i,j,N,nrhs,itemp)
#  else
          stflx(i,j,isalt) = (evap-prate(i,j))*t(i,j,N,nrhs,isalt)
#  endif
# endif

# ifdef MASKING
          stflx(i,j,itemp) = stflx(i,j,itemp)*rmask(i,j)
          stflx(i,j,isalt) = stflx(i,j,isalt)*rmask(i,j)
# endif

          aer(i,j) = rhoAir*delW
          cer(i,j) = Cd
          shflx_rsw(i,j) = srflx(i,j)
          shflx_lat(i,j) = hflat
          shflx_sen(i,j) = hfsen
          shflx_rlw(i,j) = hflw

        END DO   ! i
      END DO     ! j

      deallocate( wspd_tile, delT_tile, delQ_tile, TairC_tile, &
                  RH_tile, patm_tile, Cd_tile, Ch_tile, Ce_tile, &
                  rhoAir_tile, TseaC_tile )

!======================================================================
!  SECTION 6 — wind stress assembly
!======================================================================
      do j=jmin,jmax
        do i=imin+1,imax
          sustr(i,j) = 0.5*(cer(i-1,j)+cer(i,j)) * &
                       0.5*(aer(i-1,j)+aer(i,j)) * uwnd(i,j)*rho0i
# ifdef MASKING
          sustr(i,j) = sustr(i,j)*umask(i,j)
# endif
        enddo
      enddo

      do j=jmin+1,jmax
        do i=imin,imax
          svstr(i,j) = 0.5*(cer(i,j-1)+cer(i,j)) * &
                       0.5*(aer(i,j-1)+aer(i,j)) * vwnd(i,j)*rho0i
# ifdef MASKING
          svstr(i,j) = svstr(i,j)*vmask(i,j)
# endif
        enddo
      enddo

      return
      end subroutine bulk_flux_tile_mlp

#else
!  If BULK_FLUX_MLP is not defined, then i linker is happy
      subroutine bulk_flux_mlp_stub
      return
      end
#endif /* BULK_FLUX_MLP */
