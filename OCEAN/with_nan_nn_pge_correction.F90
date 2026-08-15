!=======================================================================
!  Neural-Network Pressure Gradient Error (PGE) Correction
!  for the CROCO sigma-coordinate.
!
!  Ref.
!  ----
!  Beckmann & Haidvogel (1993)      -- PGE analysis
!  Shchepetkin & McWilliams (2003)  -- Density Jacobian
!  Mellor et al. (1994)             -- sigma-coordinate PGE
!
!=======================================================================

module nn_pge_correction

  implicit none
  private

  !---------------------------------------------------------------------
  ! Change as you wish ----->>>
  !---------------------------------------------------------------------
  integer, parameter :: NIN  = 9    ! input features
  integer, parameter :: NH1  = 16   ! hidden layer 1 width
  integer, parameter :: NH2  = 8    ! hidden layer 2 width
  integer, parameter :: NOUT = 1    ! output dimension

  !---------------------------------------------------------------------
  ! Hyper-parameters
  !---------------------------------------------------------------------
  integer, parameter :: NEPOCHS    = 50       ! training epochs
  integer, parameter :: BATCH_SIZE = 256      ! too mini-batch size
  integer, parameter :: LR_STEP    = 10       ! epoch stride for LR decay
  real,    parameter :: LR_INIT    = 1.0e-3   ! Adam initial learning rate
  real,    parameter :: LR_FACTOR  = 0.5      ! LR decay multiplier
  real,    parameter :: ADAM_B1    = 0.9      ! Adam beta_1
  real,    parameter :: ADAM_B2    = 0.999    ! Adam beta_2
  real,    parameter :: ADAM_EPS   = 1.0e-8   ! Adam epsilon
  real,    parameter :: WEIGHT_DEC = 1.0e-5   ! L2 weight decay
  integer, parameter :: RAND_SEED  = 42       ! reproducible RNG seed
                                              ! if you know total number then no random

  !---------------------------------------------------------------------
  ! Weights  (nn_pge_init)
  !---------------------------------------------------------------------
  real, save :: W1(NH1, NIN),   b1(NH1)
  real, save :: W2(NH2, NH1),   b2(NH2)
  real, save :: W3(NOUT, NH2),  b3(NOUT)

  !---------------------------------------------------------------------
  ! Normalisation
  !---------------------------------------------------------------------
  real, save :: feat_mean(NIN)
  real, save :: feat_std (NIN)

  !---------------------------------------------------------------------
  ! Runtime
  !---------------------------------------------------------------------
  real,    save :: alpha_pge = 1.0    ! blending: 0=off, 1=full correction
  logical, save :: nn_loaded = .false.

  !---------------------------------------------------------------------
  ! Public interface
  !---------------------------------------------------------------------
  public :: nn_pge_init          ! load weights from nn_pge_weights.text
  public :: nn_pge_train         ! train (MPI parallel) + write .text
  public :: nn_pge_correction_u  ! apply correction to ru(i,j,k)
  public :: nn_pge_correction_v  ! apply correction to rv(i,j,k)
  public :: alpha_pge


contains

!=======================================================================
! SECTION 1  --  WEIGHT FILE I/O
!=======================================================================
  subroutine nn_pge_init(weight_file, blend_alpha)

    character(len=*), intent(in) :: weight_file
    real,             intent(in) :: blend_alpha

    integer :: ios, iunit, n
    logical :: fexist

    alpha_pge = blend_alpha
    nn_loaded = .false.

    inquire(file=trim(weight_file), exist=fexist)
    if (.not. fexist) then
      write(6,'(/,A)') &
        ' NN_PGE WARNING: nn_pge_weights.text not found.'
      write(6,'(A)') &
        ' NN_PGE WARNING: PGE correction DISABLED (safe fallback).'
      return
    end if

    iunit = 71
    open(unit=iunit, file=trim(weight_file), status='old', &
         action='read', iostat=ios)
    if (ios /= 0) then
      write(6,'(A,I0)') &
        ' NN_PGE ERROR: cannot open weight file, ios=', ios
      return
    end if

    read(iunit,*) feat_mean(1:NIN)
    read(iunit,*) feat_std (1:NIN)

    do n = 1, NH1
      read(iunit,*) W1(n, 1:NIN)
    end do
    read(iunit,*) b1(1:NH1)

    do n = 1, NH2
      read(iunit,*) W2(n, 1:NH1)
    end do
    read(iunit,*) b2(1:NH2)

    do n = 1, NOUT
      read(iunit,*) W3(n, 1:NH2)
    end do
    read(iunit,*) b3(1:NOUT)

    close(iunit)
    nn_loaded = .true.

    write(6,'(/,A)')             '  ======================================='
    write(6,'(A,A)')             '  NN_PGE  weights loaded  : ',            &
                                  trim(weight_file)
    write(6,'(A,F6.3)')          '  NN_PGE  alpha (blend)   : ', alpha_pge
    write(6,'(A,I0,A,I0,A,I0)') '  NN_PGE  arch            : ',            &
                                  NIN,' ->',NH1,' ->',NH2,' -> 1'
    write(6,'(A)')               '  ======================================='

  end subroutine nn_pge_init


!=======================================================================
!  write_weights_text
!=======================================================================
  subroutine write_weights_text(weight_file)

    character(len=*), intent(in) :: weight_file
    integer :: iunit, n

    iunit = 72
    open(unit=iunit, file=trim(weight_file), status='replace', &
         action='write')

    write(iunit,*) feat_mean(1:NIN)
    write(iunit,*) feat_std (1:NIN)

    do n = 1, NH1
      write(iunit,*) W1(n, 1:NIN)
    end do
    write(iunit,*) b1(1:NH1)

    do n = 1, NH2
      write(iunit,*) W2(n, 1:NH1)
    end do
    write(iunit,*) b2(1:NH2)

    do n = 1, NOUT
      write(iunit,*) W3(n, 1:NH2)
    end do
    write(iunit,*) b3(1:NOUT)

    close(iunit)
    write(6,'(A,A)') '  NN_PGE  weights saved   : ', trim(weight_file)

  end subroutine write_weights_text

!=======================================================================
!  INFERENCE (tanh i will change) TODO sm-loop
!=======================================================================
  pure real function forward_pass(x_raw) result(corr)

    real, intent(in) :: x_raw(NIN)

    real    :: x(NIN), h1v(NH1), h2v(NH2)
    integer :: n

    x(:) = (x_raw(:) - feat_mean(:)) / (feat_std(:) + 1.0e-8)

    do n = 1, NH1
      h1v(n) = tanh( dot_product(W1(n,:), x) + b1(n) )
    end do

    do n = 1, NH2
      h2v(n) = tanh( dot_product(W2(n,:), h1v) + b2(n) )
    end do

    corr = dot_product(W3(1,:), h2v) + b3(1)

  end function forward_pass


!=======================================================================
!  nn_pge_correction_u
!=======================================================================
  subroutine nn_pge_correction_u(Istr, Iend, Jstr, Jend, &
                                  i, j, k, N,             &
                                  rho, z_r, h,            &
                                  Hz_u, on_u,             &
                                  ru_corr)

    integer, intent(in) :: Istr, Iend, Jstr, Jend
    integer, intent(in) :: i, j, k, N
    real,    intent(in) :: rho(Istr-2:Iend+2, Jstr-2:Jend+2, N)
    real,    intent(in) :: z_r(Istr-2:Iend+2, Jstr-2:Jend+2, N)
    real,    intent(in) :: h  (Istr-2:Iend+2, Jstr-2:Jend+2)
    real,    intent(in) :: Hz_u, on_u
    real,    intent(out):: ru_corr

    real :: feats(NIN), dz_r, dz_l

    ru_corr = 0.0
    if (.not. nn_loaded) return

    if (k < N) then
      dz_r = z_r(i,   j, k+1) - z_r(i,   j, k)
      dz_l = z_r(i-1, j, k+1) - z_r(i-1, j, k)
    else
      dz_r = z_r(i,   j, k) - z_r(i,   j, k-1)
      dz_l = z_r(i-1, j, k) - z_r(i-1, j, k-1)
    end if

    feats(1) = rho(i,   j, k) - rho(i-1, j, k)
    feats(2) = rho(i,   j, k) + rho(i-1, j, k)
    feats(3) = z_r(i,   j, k) - z_r(i-1, j, k)
    feats(4) = z_r(i,   j, k) + z_r(i-1, j, k)
    feats(5) = dz_r
    feats(6) = dz_l
    feats(7) = h(i,   j)
    feats(8) = h(i-1, j)
    feats(9) = real(k) / real(N)

    ru_corr = Hz_u * on_u * alpha_pge * forward_pass(feats)

  end subroutine nn_pge_correction_u


!=======================================================================
!  nn_pge_correction_v
!=======================================================================
  subroutine nn_pge_correction_v(Istr, Iend, Jstr, Jend, &
                                  i, j, k, N,             &
                                  rho, z_r, h,            &
                                  Hz_v, om_v,             &
                                  rv_corr)

    integer, intent(in) :: Istr, Iend, Jstr, Jend
    integer, intent(in) :: i, j, k, N
    real,    intent(in) :: rho(Istr-2:Iend+2, Jstr-2:Jend+2, N)
    real,    intent(in) :: z_r(Istr-2:Iend+2, Jstr-2:Jend+2, N)
    real,    intent(in) :: h  (Istr-2:Iend+2, Jstr-2:Jend+2)
    real,    intent(in) :: Hz_v, om_v
    real,    intent(out):: rv_corr

    real :: feats(NIN), dz_n, dz_s

    rv_corr = 0.0
    if (.not. nn_loaded) return

    if (k < N) then
      dz_n = z_r(i, j,   k+1) - z_r(i, j,   k)
      dz_s = z_r(i, j-1, k+1) - z_r(i, j-1, k)
    else
      dz_n = z_r(i, j,   k) - z_r(i, j,   k-1)
      dz_s = z_r(i, j-1, k) - z_r(i, j-1, k-1)
    end if

    feats(1) = rho(i, j,   k) - rho(i, j-1, k)
    feats(2) = rho(i, j,   k) + rho(i, j-1, k)
    feats(3) = z_r(i, j,   k) - z_r(i, j-1, k)
    feats(4) = z_r(i, j,   k) + z_r(i, j-1, k)
    feats(5) = dz_n
    feats(6) = dz_s
    feats(7) = h(i, j  )
    feats(8) = h(i, j-1)
    feats(9) = real(k) / real(N)

    rv_corr = Hz_v * om_v * alpha_pge * forward_pass(feats)

  end subroutine nn_pge_correction_v

!=======================================================================
!  nn_pge_train
!=======================================================================
  subroutine nn_pge_train(rho_in, z_r_in, h_in, Nlev, Lm_in, Mm_in, &
                           mpi_comm)

    implicit none
    include 'mpif.h'

    !--- Explicit-shape using Lm_in, Mm_in, Nlev -----------------------
    !    Assumed-shape (:,:,:) breaks with CROCO's CPP-macro dimensioned
    !    arrays: the Fortran array descriptor gets wrong base/stride,
    !    making MPI receive buffers point to wrong memory (glb_nsamp_r=0).
    !    Explicit-shape with 0-based lower bound matches CROCO exactly.
    integer, intent(in) :: Nlev, Lm_in, Mm_in, mpi_comm
    real,    intent(in) :: rho_in(0:Lm_in+1, 0:Mm_in+1, Nlev)
    real,    intent(in) :: z_r_in(0:Lm_in+1, 0:Mm_in+1, Nlev)
    real,    intent(in) :: h_in  (0:Lm_in+1, 0:Mm_in+1)
    !
    !--- MPI -------------------------------------------------------------
    integer :: my_rank, n_ranks, mpi_err

    !--- local training data for this rank's tile -----------------------
    integer :: nsamp_loc
    real,    allocatable :: Xl(:,:)
    real,    allocatable :: Yl(:)
    integer, allocatable :: perm(:)

    !--- MPI datatype: MPI_DOUBLE_PRECISION matches plain 'real' under
    !    -fdefault-real-8 / -r8 (where 'real' is 8 bytes).  MPI_REAL is
    !    always 4 bytes and would cause every Allreduce to read the
    !    zero-filled low half of each element and return 0.
    !    MPI_DOUBLE_PRECISION is safe even without -r8: gfortran's
    !    MPI_DOUBLE_PRECISION is always 8 bytes regardless of flags.

    !--- Normalisation ---------------------------------------------------
    real     :: loc_sum(NIN), loc_sq(NIN)
    real     :: glb_sum(NIN), glb_sq(NIN)
    real     :: glb_nsamp_r                 ! global sample count as real

    !--- Adam updates ----------------------------------------------------
    real :: lW1(NH1,NIN),   lb1(NH1)
    real :: lW2(NH2,NH1),   lb2(NH2)
    real :: lW3(NOUT,NH2),  lb3(NOUT)

    !--- Adam moment -----------------------------------------------------
    real :: mW1(NH1,NIN),  vW1(NH1,NIN)
    real :: mb1(NH1),      vb1(NH1)
    real :: mW2(NH2,NH1),  vW2(NH2,NH1)
    real :: mb2(NH2),      vb2(NH2)
    real :: mW3(NOUT,NH2), vW3(NOUT,NH2)
    real :: mb3(NOUT),     vb3(NOUT)
!!!=====================================================================
    !--- Communication buffer  ------------------------------------------
    !    Layout (all as reals):
    !      [1 .. NH1*NIN ]             dW1 (row-major)
    !      [+NH1         ]             db1
    !      [+NH2*NH1     ]             dW2
    !      [+NH2         ]             db2
    !      [+NOUT*NH2    ]             dW3
    !      [+NOUT        ]             db3
    !      [+1           ]             loss accumulator
    !      [+1           ]             count
!!!=====================================================================
    integer, parameter :: GBUF = NH1*NIN + NH1 &
                                + NH2*NH1 + NH2 &
                                + NOUT*NH2 + NOUT + 2
    real     :: loc_buf(GBUF), glb_buf(GBUF)

    !--- per-forward / backward -----------------------------------------
    real :: xn(NIN)
    real :: h1v(NH1), h2v(NH2)
    real :: pred_s, err_s, d_out, batch_loss
    real :: d_h2(NH2), d_h1(NH1)
    real :: dW1_b(NH1,NIN), db1_b(NH1)
    real :: dW2_b(NH2,NH1), db2_b(NH2)
    real :: dW3_b(NOUT,NH2),db3_b(NOUT)

    !--- Adam scalars ---------------------------------------------------
    real    :: lr, bc1, bc2, tmp
    integer :: adam_step

    !--- loop -----------------------------------------------------------
    integer :: ep, s, si, b_start, b_end, bsz
    integer :: i, j_loc, k, kk, n1, n2, p
    integer :: i_lo, i_hi, j_lo, j_hi, k_hi, nsamp_max
    integer :: ni_iter, nj_iter, nk_iter
    real    :: dz_r, dz_l, P_r, P_l, pg_true, pg_sigma
    real    :: g_phys, rho0, HalfGRho
    real     :: loc_loss, glb_loss, glb_cnt, rmse
    integer :: nsamp_glb_int

    !--- buffer -------constants ----------------------------------------
    integer :: o_W1, o_b1, o_W2, o_b2, o_W3, o_b3, o_loss, o_cnt

    o_W1  = 1
    o_b1  = o_W1  + NH1*NIN
    o_W2  = o_b1  + NH1
    o_b2  = o_W2  + NH2*NH1
    o_W3  = o_b2  + NH2
    o_b3  = o_W3  + NOUT*NH2
    o_loss= o_b3  + NOUT
    o_cnt = o_loss + 1

    call MPI_Comm_rank(mpi_comm, my_rank, mpi_err)
    call MPI_Comm_size(mpi_comm, n_ranks, mpi_err)

    !--- Synchronise ALL ranks before any collective operation.
    !    CROCO may call nn_pge_train from inside a tile-loop that
    !    serialises ranks: fast ranks hit MPI_Allreduce before slow
    !    ranks have entered the routine, causing a collective mismatch
    !    that returns glb_nsamp_r=0 and aborts training early.
    call MPI_Barrier(mpi_comm, mpi_err)

    g_phys   = 9.81
    rho0     = 1025.0
    HalfGRho = 0.5 * g_phys / rho0

    if (my_rank == 0) then
      write(6,'(/,A)')         &
        '  ========================================================'
      write(6,'(A)')           &
        '   NN_PGE  TRAINING  (pure Fortran, MPI CPU parallel)'
      write(6,'(A,I0)')        '   MPI ranks          : ', n_ranks
      write(6,'(A,I0)')        '   Epochs             : ', NEPOCHS
      write(6,'(A,I0)')        '   Batch size         : ', BATCH_SIZE
      write(6,'(A,ES9.2)')     '   Initial LR         : ', LR_INIT
      write(6,'(A,I0,A,I0,A)') '   LR step decay      : x',             &
                                int(LR_FACTOR*10),'/10 every ',          &
                                LR_STEP,' epochs'
      write(6,'(A,I0,A,I0,A,I0)') &
                                '   Architecture       : ',              &
                                NIN,' ->',NH1,' ->',NH2,' -> 1'
      write(6,'(A)')           &
        '  ========================================================'
    end if

!!!----------------------------------------------------------------------
!!!----------------------------------------------------------------------
    !--- Array bounds from explicit-shape arguments ----------------------
    !    rho_in is declared (0:Lm_in+1, 0:Mm_in+1, Nlev) so:
    i_lo = 0
    i_hi = Lm_in + 1
    j_lo = 0
    j_hi = Mm_in + 1
    k_hi = Nlev

    !--- Loop bounds:
    !    i: from i_lo+2 to i_hi-1  (interior; need i-1 valid => start at i_lo+2=2)
    !    j: from j_lo+1 to j_hi-1  (interior; skip boundary ghost rows)
    !    k: from 1 to k_hi
    ni_iter = (i_hi - 1) - (i_lo + 2) + 1   ! = Lm_in - 1
    nj_iter = (j_hi - 1) - (j_lo + 1) + 1   ! = Mm_in - 1
    nk_iter = k_hi

    !--- Calculate maximum possible samples ------------------------------
    nsamp_max = max(1, ni_iter * nj_iter * nk_iter)

    !--- Debug output: show actual bounds on rank 0 ----------------------
    if (my_rank == 0) then
      write(6,'(A,4I6)') '   Array bounds (i_lo,i_hi,j_lo,j_hi): ', &
                          i_lo, i_hi, j_lo, j_hi
      write(6,'(A,I6)')  '   Vertical levels (k_hi from array) : ', k_hi
      write(6,'(A,I6)')  '   Nlev parameter (input)            : ', Nlev
      write(6,'(A,3I6)') '   Loop iterations (ni,nj,nk)        : ', &
                          ni_iter, nj_iter, nk_iter
      write(6,'(A,I0)')  '   Max samples per rank (nsamp_max)  : ', nsamp_max
    end if

    !--- Validate loop ranges before allocating --------------------------
    if (ni_iter < 1 .or. nj_iter < 1 .or. nk_iter < 1) then
      if (my_rank == 0) then
        write(6,'(A)')      ' NN_PGE ERROR: Invalid loop ranges!'
        write(6,'(A,3I6)')  '   ni_iter, nj_iter, nk_iter = ', ni_iter, nj_iter, nk_iter
        write(6,'(A,3I6)')  '   Lm_in, Mm_in, Nlev        = ', Lm_in, Mm_in, Nlev
        write(6,'(A)')      '   Check array dimensions passed to nn_pge_train.'
      end if
      return
    end if

    allocate( Xl(nsamp_max, NIN) )
    allocate( Yl(nsamp_max)      )
    nsamp_loc = 0

    !--- Loop over interior points only, skipping ghost/boundary ---------
    !    j_loc: from j_lo+1 to j_hi-1  (interior J, 1-based in tile)
    !    i:     from i_lo+1 to i_hi-1  (interior I; i-1 >= i_lo=1 always)
    !    k:     from 1 to k_hi
    !
    !    hydrostatic_pg is inlined to avoid passing array sections to an
    !    explicit-shape dummy -- that construct caused silent failures with
    !    several MPI+compiler combinations (temporary copy not visible to
    !    the caller; segfault or NaN before nsamp_loc increments).
    do j_loc = j_lo + 1, j_hi - 1
      do i = i_lo + 2, i_hi - 1
        do k = 1, k_hi

          if (k < k_hi) then
            dz_r = z_r_in(i,   j_loc, k+1) - z_r_in(i,   j_loc, k)
            dz_l = z_r_in(i-1, j_loc, k+1) - z_r_in(i-1, j_loc, k)
          else
            dz_r = z_r_in(i,   j_loc, k) - z_r_in(i,   j_loc, k-1)
            dz_l = z_r_in(i-1, j_loc, k) - z_r_in(i-1, j_loc, k-1)
          end if

          !--- inline: hydrostatic pressure at column (i, j_loc) ---------
          P_r = 0.0
          do kk = k_hi, k + 1, -1
            P_r = P_r + (g_phys / rho0)                                    &
                      * 0.5*(rho_in(i, j_loc, kk) + rho_in(i, j_loc, kk-1)) &
                      * abs(z_r_in(i, j_loc, kk) - z_r_in(i, j_loc, kk-1))
          end do

          !--- inline: hydrostatic pressure at column (i-1, j_loc) -------
          P_l = 0.0
          do kk = k_hi, k + 1, -1
            P_l = P_l + (g_phys / rho0)                                      &
                      * 0.5*(rho_in(i-1, j_loc, kk) + rho_in(i-1, j_loc, kk-1)) &
                      * abs(z_r_in(i-1, j_loc, kk) - z_r_in(i-1, j_loc, kk-1))
          end do

          pg_true  = P_l - P_r

          pg_sigma = -HalfGRho                                               &
                   * (rho_in(i, j_loc, k) + rho_in(i-1, j_loc, k))           &
                   * (z_r_in(i, j_loc, k) - z_r_in(i-1, j_loc, k))

          nsamp_loc = nsamp_loc + 1
          s         = nsamp_loc

          Xl(s, 1) = rho_in(i,   j_loc, k) - rho_in(i-1, j_loc, k)
          Xl(s, 2) = rho_in(i,   j_loc, k) + rho_in(i-1, j_loc, k)
          Xl(s, 3) = z_r_in(i,   j_loc, k) - z_r_in(i-1, j_loc, k)
          Xl(s, 4) = z_r_in(i,   j_loc, k) + z_r_in(i-1, j_loc, k)
          Xl(s, 5) = dz_r
          Xl(s, 6) = dz_l
          Xl(s, 7) = h_in(i,   j_loc)
          Xl(s, 8) = h_in(i-1, j_loc)
          Xl(s, 9) = real(k) / real(k_hi)
          Yl(s)    = pg_true - pg_sigma

        end do
      end do
    end do
!!!----------------------------------------------------------------------
!!!----------------------------------------------------------------------

    !--- Per-rank diagnostic: print local sample count before Allreduce --
    write(6,'(A,I0,A,I0)') '   [nn_pge] rank ', my_rank, &
                            ' local samples = ', nsamp_loc

    !--- Sum sample counts: piggyback on loc_buf(o_cnt) so the count  ---
    !    travels in the same allreduce as a zero gradient pass.          --
    !    A standalone scalar allreduce proved unreliable on several       --
    !    MPI+compiler combos (returned 0 despite MPI_SUCCESS status);    --
    !    an array-based GBUF allreduce has no aliasing or               --
    !    temporary-address issues.                                        --
    loc_buf        = 0.0
    loc_buf(o_cnt) = real(nsamp_loc)
    glb_buf        = 0.0
    call MPI_Allreduce(loc_buf, glb_buf, GBUF, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, mpi_err)
    if (mpi_err /= MPI_SUCCESS) then
      write(6,'(A,I0,A,I0)') ' NN_PGE ERROR: MPI_Allreduce(count) failed, rank=',&
                              my_rank,' mpi_err=', mpi_err
    end if
    glb_nsamp_r   = glb_buf(o_cnt)
    nsamp_glb_int = nint(glb_nsamp_r)
    loc_loss      = 0.0

    write(6,'(A,I0,A,I0)') &
      '   [nn_pge] rank ', my_rank, ' glb_nsamp_r=', nsamp_glb_int

    if (my_rank == 0) &
      write(6,'(A,I0)') '   Total training samples : ', nsamp_glb_int

    !--- Early exit if no samples collected ------------------------------
    if (nsamp_glb_int < 1) then
      if (my_rank == 0) then
        write(6,'(A)') ' NN_PGE ERROR: No training samples collected!'
        write(6,'(A)') '   Possible causes:'
        write(6,'(A)') '   - Array dimensions too small for loop bounds'
        write(6,'(A)') '   - k_hi (vertical levels) = 0'
      end if
      deallocate(Xl, Yl)
      return
    end if

    do n1 = 1, NIN
      loc_sum(n1) = sum(Xl(1:nsamp_loc, n1))
      loc_sq (n1) = sum(Xl(1:nsamp_loc, n1)**2)
    end do

    call MPI_Allreduce(loc_sum, glb_sum, NIN, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, mpi_err)
    call MPI_Allreduce(loc_sq,  glb_sq,  NIN, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, mpi_err)

    do n1 = 1, NIN
      feat_mean(n1) = glb_sum(n1) / glb_nsamp_r
      feat_std (n1) = sqrt( max( glb_sq(n1)/glb_nsamp_r        &
                               - feat_mean(n1)**2, 0.0 )        &
                           + 1.0e-8 )
    end do


    do s = 1, nsamp_loc
      do n1 = 1, NIN
        Xl(s,n1) = (Xl(s,n1) - feat_mean(n1)) / feat_std(n1)
      end do
    end do

    !-------------------------------------------------------------------
    ! Xavier weight init
    !-------------------------------------------------------------------
    if (my_rank == 0) then
      call xavier_init(lW1, NH1,  NIN,  RAND_SEED    )
      call xavier_init(lW2, NH2,  NH1,  RAND_SEED + 1)
      call xavier_init(lW3, NOUT, NH2,  RAND_SEED + 2)
      lb1 = 0.0;  lb2 = 0.0;  lb3 = 0.0
    end if

    call MPI_Bcast(lW1, NH1*NIN,  MPI_DOUBLE_PRECISION, 0, mpi_comm, mpi_err)
    call MPI_Bcast(lb1, NH1,      MPI_DOUBLE_PRECISION, 0, mpi_comm, mpi_err)
    call MPI_Bcast(lW2, NH2*NH1,  MPI_DOUBLE_PRECISION, 0, mpi_comm, mpi_err)
    call MPI_Bcast(lb2, NH2,      MPI_DOUBLE_PRECISION, 0, mpi_comm, mpi_err)
    call MPI_Bcast(lW3, NOUT*NH2, MPI_DOUBLE_PRECISION, 0, mpi_comm, mpi_err)
    call MPI_Bcast(lb3, NOUT,     MPI_DOUBLE_PRECISION, 0, mpi_comm, mpi_err)

    !--- zero Adam moments and step -----------------------------------
    mW1=0.0; vW1=0.0; mb1=0.0; vb1=0.0
    mW2=0.0; vW2=0.0; mb2=0.0; vb2=0.0
    mW3=0.0; vW3=0.0; mb3=0.0; vb3=0.0
    adam_step = 0
    lr        = LR_INIT

    !--- permutation index ---------------------------------------------
    allocate(perm(nsamp_loc))
    do s = 1, nsamp_loc
      perm(s) = s
    end do

    if (my_rank == 0) then
      write(6,'(/,A6,2X,A14,2X,A10)') 'Epoch','Global RMSE','LR'
      write(6,'(A)') '  ----------------------------------------'
    end if

    !===================================================================
    ! MAIN TRAINING
    !===================================================================
    do ep = 1, NEPOCHS

      !----- LR step decay --------------------------------------------
      if (ep > 1 .and. mod(ep-1, LR_STEP) == 0) lr = lr * LR_FACTOR

      !--------Shuffle local permutation (Fisher-Yates, LCG RNG) ------
      call fy_shuffle(perm, nsamp_loc, ep + my_rank * NEPOCHS)

      loc_loss = 0.0

      !--- 4c. Mini-batch loop -----------------------------------------
      b_start = 1
      do while (b_start <= nsamp_loc)

        b_end = min(b_start + BATCH_SIZE - 1, nsamp_loc)
        bsz   = b_end - b_start + 1
        adam_step = adam_step + 1

        !--- Gradient accumulators -------------------------
        dW1_b = 0.0;  db1_b = 0.0
        dW2_b = 0.0;  db2_b = 0.0
        dW3_b = 0.0;  db3_b = 0.0
        batch_loss = 0.0

        !--------------------------------------------------------------
        do s = b_start, b_end
          si   = perm(s)
          xn(:)= Xl(si, 1:NIN)

          !--- forward pass -------------------------------------------

          ! Layer 1
          do n1 = 1, NH1
            h1v(n1) = tanh( dot_product(lW1(n1,:), xn) + lb1(n1) )
          end do

          ! Layer 2
          do n2 = 1, NH2
            h2v(n2) = tanh( dot_product(lW2(n2,:), h1v) + lb2(n2) )
          end do

          ! Output (linear. TODO sm-loop)
          pred_s = dot_product(lW3(1,:), h2v) + lb3(1)

          err_s      = pred_s - Yl(si)
          batch_loss = batch_loss + err_s * err_s

          !--- backward pass  (MSE: dL/d_pred = err) -----------------

          ! Output
          d_out         = err_s ! dL/d_pred
          dW3_b(1,:)    = dW3_b(1,:) + d_out * h2v(:)
          db3_b(1)      = db3_b(1)   + d_out

          ! Layer 2  (tanh': d(tanh)/dx = 1 - tanh^2) --sm-loop--***
          do n2 = 1, NH2
            d_h2(n2) = d_out * lW3(1,n2) * (1.0 - h2v(n2)**2)
          end do
          do n2 = 1, NH2
            dW2_b(n2,:) = dW2_b(n2,:) + d_h2(n2) * h1v(:)
            db2_b(n2)   = db2_b(n2)   + d_h2(n2)
          end do

          ! Layer 1  (backprop through W2, then tanh')
          do n1 = 1, NH1
            d_h1(n1) = 0.0
            do n2 = 1, NH2
              d_h1(n1) = d_h1(n1) + d_h2(n2) * lW2(n2,n1)
            end do
            d_h1(n1) = d_h1(n1) * (1.0 - h1v(n1)**2)
          end do
          do n1 = 1, NH1
            dW1_b(n1,:) = dW1_b(n1,:) + d_h1(n1) * xn(:)
            db1_b(n1)   = db1_b(n1)   + d_h1(n1)
          end do

        end do

        !==============================================================
        !  Pack gradients into flat communication buffer
        !==============================================================
        loc_buf = 0.0

        p = o_W1 - 1
        do n1 = 1, NH1
          do n2 = 1, NIN
            p = p + 1;  loc_buf(p) = dW1_b(n1,n2)
          end do
        end do
        do n1 = 1, NH1
          loc_buf(o_b1 + n1 - 1) = db1_b(n1)
        end do

        p = o_W2 - 1
        do n1 = 1, NH2
          do n2 = 1, NH1
            p = p + 1;  loc_buf(p) = dW2_b(n1,n2)
          end do
        end do
        do n1 = 1, NH2
          loc_buf(o_b2 + n1 - 1) = db2_b(n1)
        end do

        p = o_W3 - 1
        do n1 = 1, NOUT
          do n2 = 1, NH2
            p = p + 1;  loc_buf(p) = dW3_b(n1,n2)
          end do
        end do
        do n1 = 1, NOUT
          loc_buf(o_b3 + n1 - 1) = db3_b(n1)
        end do

        loc_buf(o_loss) = batch_loss
        loc_buf(o_cnt)  = real(bsz)

        !==============================================================
        !  sum gradients + loss
        !==============================================================
        call MPI_Allreduce(loc_buf, glb_buf, GBUF, &
                           MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, mpi_err)

        !--- Accumulate reduced batch loss into the epoch total --------
        loc_loss = loc_loss + glb_buf(o_loss)

        !==============================================================
        !  unpack global gradients and norm
        !==============================================================
        glb_cnt  = max(glb_buf(o_cnt), 1.0)

        p = o_W1 - 1
        do n1 = 1, NH1
          do n2 = 1, NIN
            p = p + 1;  dW1_b(n1,n2) = glb_buf(p) / glb_cnt
          end do
        end do
        do n1 = 1, NH1
          db1_b(n1) = glb_buf(o_b1 + n1 - 1) / glb_cnt
        end do

        p = o_W2 - 1
        do n1 = 1, NH2
          do n2 = 1, NH1
            p = p + 1;  dW2_b(n1,n2) = glb_buf(p) / glb_cnt
          end do
        end do
        do n1 = 1, NH2
          db2_b(n1) = glb_buf(o_b2 + n1 - 1) / glb_cnt
        end do

        p = o_W3 - 1
        do n1 = 1, NOUT
          do n2 = 1, NH2
            p = p + 1;  dW3_b(n1,n2) = glb_buf(p) / glb_cnt
          end do
        end do
        do n1 = 1, NOUT
          db3_b(n1) = glb_buf(o_b3 + n1 - 1) / glb_cnt
        end do

        bc1 = 1.0 - ADAM_B1**real(adam_step)
        bc2 = 1.0 - ADAM_B2**real(adam_step)

        ! ----- W1 / b1 -----
        do n1 = 1, NH1
          do n2 = 1, NIN
            tmp           = dW1_b(n1,n2) + WEIGHT_DEC * lW1(n1,n2)
            mW1(n1,n2)    = ADAM_B1*mW1(n1,n2) + (1.0-ADAM_B1)*tmp
            vW1(n1,n2)    = ADAM_B2*vW1(n1,n2) + (1.0-ADAM_B2)*tmp**2
            lW1(n1,n2)    = lW1(n1,n2)                               &
                           - lr * (mW1(n1,n2)/bc1)                   &
                           / (sqrt(vW1(n1,n2)/bc2) + ADAM_EPS)
          end do
          tmp       = db1_b(n1)
          mb1(n1)   = ADAM_B1*mb1(n1) + (1.0-ADAM_B1)*tmp
          vb1(n1)   = ADAM_B2*vb1(n1) + (1.0-ADAM_B2)*tmp**2
          lb1(n1)   = lb1(n1) - lr*(mb1(n1)/bc1) &
                     / (sqrt(vb1(n1)/bc2) + ADAM_EPS)
        end do

        ! ----- W2 / b2 -----
        do n1 = 1, NH2
          do n2 = 1, NH1
            tmp           = dW2_b(n1,n2) + WEIGHT_DEC * lW2(n1,n2)
            mW2(n1,n2)    = ADAM_B1*mW2(n1,n2) + (1.0-ADAM_B1)*tmp
            vW2(n1,n2)    = ADAM_B2*vW2(n1,n2) + (1.0-ADAM_B2)*tmp**2
            lW2(n1,n2)    = lW2(n1,n2)                               &
                           - lr * (mW2(n1,n2)/bc1)                   &
                           / (sqrt(vW2(n1,n2)/bc2) + ADAM_EPS)
          end do
          tmp       = db2_b(n1)
          mb2(n1)   = ADAM_B1*mb2(n1) + (1.0-ADAM_B1)*tmp
          vb2(n1)   = ADAM_B2*vb2(n1) + (1.0-ADAM_B2)*tmp**2
          lb2(n1)   = lb2(n1) - lr*(mb2(n1)/bc1) &
                     / (sqrt(vb2(n1)/bc2) + ADAM_EPS)
        end do

        ! ----- W3 / b3 -----
        do n1 = 1, NOUT
          do n2 = 1, NH2
            tmp           = dW3_b(n1,n2) + WEIGHT_DEC * lW3(n1,n2)
            mW3(n1,n2)    = ADAM_B1*mW3(n1,n2) + (1.0-ADAM_B1)*tmp
            vW3(n1,n2)    = ADAM_B2*vW3(n1,n2) + (1.0-ADAM_B2)*tmp**2
            lW3(n1,n2)    = lW3(n1,n2)                               &
                           - lr * (mW3(n1,n2)/bc1)                   &
                           / (sqrt(vW3(n1,n2)/bc2) + ADAM_EPS)
          end do
          tmp       = db3_b(n1)
          mb3(n1)   = ADAM_B1*mb3(n1) + (1.0-ADAM_B1)*tmp
          vb3(n1)   = ADAM_B2*vb3(n1) + (1.0-ADAM_B2)*tmp**2
          lb3(n1)   = lb3(n1) - lr*(mb3(n1)/bc1) &
                     / (sqrt(vb3(n1)/bc2) + ADAM_EPS)
        end do

        b_start = b_start + BATCH_SIZE

      end do

      !--- Global epoch RMSE (loc_loss is already the all-rank sum from
      !    the per-batch Allreduces above; a second Allreduce would double-count)
      glb_cnt = glb_nsamp_r   ! already reduced at start of training
      rmse = sqrt(loc_loss / max(glb_cnt, 1.0))

      if (my_rank == 0) &
        write(6,'(I6,2X,ES14.6,2X,ES10.3)') ep, rmse, lr

    end do

    !-------------------------------------------------------------------
    ! cp local trained weights
    !-------------------------------------------------------------------
    W1 = lW1;  b1 = lb1
    W2 = lW2;  b2 = lb2
    W3 = lW3;  b3 = lb3
    nn_loaded = .true.

    !-------------------------------------------------------------------
    ! nn_pge_weights.text ** todo based on <title>
    !-------------------------------------------------------------------
    if (my_rank == 0) then
      call write_weights_text('nn_pge_weights.text')
      write(6,'(/,A)') '  NN_PGE  TRAINING COMPLETE'
      write(6,'(A)')   '  ========================================================'
    end if
    !--- All ranks must wait until rank 0 has flushed the weight file  ---
    !--- before any rank calls nn_pge_init (which opens the same file). -
    call MPI_Barrier(mpi_comm, mpi_err)

    deallocate(Xl, Yl, perm)

  end subroutine nn_pge_train


!=======================================================================
! HELPERS
!=======================================================================
!  Hydrostatic_pg, Hydrostatic pressure
!  P_out = (g/rho0) * integral_{0}^{-z_r(k_tgt)} rho dz'
!=======================================================================
  subroutine hydrostatic_pg(rho_col, z_col, k_tgt, Nlev, &
                             g_in, rho0_in, P_out)

    integer, intent(in) :: k_tgt, Nlev
    real,    intent(in) :: rho_col(Nlev), z_col(Nlev)
    real,    intent(in) :: g_in, rho0_in
    real,    intent(out):: P_out

    real    :: dz, rho_avg
    integer :: kk

    ! Integrate from surface level (k=Nlev) downward to k_tgt+1
    P_out = 0.0
    do kk = Nlev, k_tgt + 1, -1
      dz      = abs(z_col(kk) - z_col(kk-1))
      rho_avg = 0.5 * (rho_col(kk) + rho_col(kk-1))
      P_out   = P_out + (g_in / rho0_in) * rho_avg * dz
    end do

  end subroutine hydrostatic_pg


!=======================================================================
!  xavier_init  (TODO sm-loop)
!    w ~ Uniform[ -sqrt(6/(nrow+ncol)),  +sqrt(6/(nrow+ncol)) ]
!=======================================================================
  subroutine xavier_init(W, nrow, ncol, seed_in)

    integer, intent(in)  :: nrow, ncol, seed_in
    real,    intent(out) :: W(nrow, ncol)

    real            :: limit, u
    integer         :: i, j
    integer(kind=8) :: state

    limit = sqrt(6.0 / real(nrow + ncol))
    state = int(seed_in, kind=8)

    do i = 1, nrow
      do j = 1, ncol
        state  = mod( 1664525_8 * state + 1013904223_8, 2147483648_8 )
        u      = real(state) / 2147483648.0
        W(i,j) = limit * (2.0 * u - 1.0)
      end do
    end do

  end subroutine xavier_init


!=======================================================================
! Fisher-Yates (Knuth) with LCG
!=======================================================================
  subroutine fy_shuffle(perm, n, seed_in)

    integer, intent(inout) :: perm(n)
    integer, intent(in)    :: n, seed_in

    integer         :: i, j, tmp
    integer(kind=8) :: state

    state = int(seed_in, kind=8) + 999983_8

    do i = n, 2, -1
      state  = mod( 1664525_8 * state + 1013904223_8, 2147483648_8 )
      j      = mod( int(state), i ) + 1        ! j in [1, i]
      tmp    = perm(i)
      perm(i)= perm(j)
      perm(j)= tmp
    end do

  end subroutine fy_shuffle

end module nn_pge_correction
