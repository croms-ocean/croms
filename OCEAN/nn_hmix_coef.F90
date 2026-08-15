!=======================================================================
! CROCO – Neural-network horizontal mixing coefficient module
!=======================================================================

module nn_hmix_coef_mod
  implicit none
  private

  ! ---------------------------------------------------------------
  ! Public
  ! ---------------------------------------------------------------
  public :: nn_hmix_predict
  public :: nn_hmix_buffer
  public :: alpha_hmix

  ! ---------------------------------------------------------------
  ! Blending weight: alpha=1 is NN; alpha=0 is Smagorinsky
  !                   change if you want
  ! ---------------------------------------------------------------
  real(8), parameter :: alpha_hmix = 1.0d0

  ! ---------------------------------------------------------------
  ! Network ::
  ! ------------------------TODO-sm-loop---------------------------
  integer, parameter :: NIN  = 6
  integer, parameter :: NH1  = 32
  integer, parameter :: NH2  = 32
  integer, parameter :: NOUT = 2
  ! --------------------its just a basic net-----------------------
  ! Normalisation ::
  ! ---------------------------------------------------------------
  real(8), parameter, dimension(NIN) :: feat_mean = (/ &
       1.0d-10, 1.0d-10, 5.0d-11, 1.0d-4, 0.0d0, 20.0d0 /)
  real(8), parameter, dimension(NIN) :: feat_std  = (/ &
       1.0d-9,  1.0d-9,  5.0d-10, 1.0d-3, 1.0d0, 15.0d0 /)
  ! ---------------------------------------------------------------
  ! De-normalisation ::
  ! ---------------------------------------------------------------
  real(8), parameter :: kap_log_mean =  0.0d0
  real(8), parameter :: kap_log_std  =  2.0d0
  real(8), parameter :: nu_log_mean  =  0.0d0
  real(8), parameter :: nu_log_std   =  2.0d0
  ! ---------------------------------------------------------------
  ! Coefficients ::
  ! ---------------------------------------------------------------
  real(8), parameter :: coef_min = 0.0d0
  real(8), parameter :: coef_max = 1.0d4   ! hard upper guard
  ! ---------------------------------------------------------------
  ! Diagnostic stuff 
  ! ---------------------------------------------------------------
  logical, save :: first_call = .true.

  ! ---------------------------------------------------------------
  ! Network weights
  ! ---------------------------------------------------------------
  real(8), save :: W1(NH1, NIN)
  real(8), save :: b1(NH1)
  real(8), save :: W2(NH2, NH1)
  real(8), save :: b2(NH2)
  real(8), save :: W3(NOUT, NH2)
  real(8), save :: b3(NOUT)
  logical, save :: weights_initialised = .false.
contains
  !=================================================================
  ! Initialise network weights
  !=================================================================
  subroutine init_nn_weights()
    integer :: i, j

    W1 = 0.0d0
    do i = 1, NIN
      W1(i, i) = 0.5d0
    end do
    b1 = 0.0d0

    W2 = 0.0d0
    do i = 1, NH2
      W2(i, i) = 1.0d0
    end do
    b2 = 0.0d0

    W3       = 0.0d0
    W3(1, 1) = 1.0d0
    W3(2, 2) = 1.0d0
    b3       = 0.0d0

    weights_initialised = .true.
  end subroutine init_nn_weights

  !=================================================================
  ! Evaluate NN for one
  !======================+++TODO-sm-loop+++=========================
  subroutine nn_forward(feat, log_kap, log_nu)
    real(8), intent(in)  :: feat(NIN)
    real(8), intent(out) :: log_kap, log_nu
    real(8) :: h1(NH1), h2(NH2), out(NOUT)
    integer :: i, j
    if (.not. weights_initialised) call init_nn_weights()
    ! --- Hidden layer 1: h1 = tanh(W1 * feat + b1) ---
    do i = 1, NH1
      h1(i) = b1(i)
      do j = 1, NIN
        h1(i) = h1(i) + W1(i,j) * feat(j)
      end do
      h1(i) = tanh(h1(i))
    end do
    ! --- Hidden layer 2: h2 = tanh(W2 * h1 + b2) ---
    do i = 1, NH2
      h2(i) = b2(i)
      do j = 1, NH1
        h2(i) = h2(i) + W2(i,j) * h1(j)
      end do
      h2(i) = tanh(h2(i))
    end do
    ! --- Output layer (linear): out = W3 * h2 + b3 ---
    do i = 1, NOUT
      out(i) = b3(i)
      do j = 1, NH2
        out(i) = out(i) + W3(i,j) * h2(j)
      end do
    end do
    log_kap = out(1)
    log_nu  = out(2)
  end subroutine nn_forward

  !=================================================================
  ! nn_hmix_predict
  !=================================================================
  subroutine nn_hmix_predict(                 &
       dudx2, dvdy2,                          &
       defrate2,                              &
       bvf_ij, ri_ij, Hz_ij,                  &
       surf,                                  &
       smago_nu_ref, smago_kap_ref,           &
       kap_out, nu_out,                       &
       cfl_kap_max, my_rank, iic,             &
       Istr, Iend, Jstr, Jend, Nlev)
    real(8), intent(in)  :: dudx2, dvdy2, defrate2
    real(8), intent(in)  :: bvf_ij, ri_ij, Hz_ij
    real(8), intent(in)  :: surf
    real(8), intent(in)  :: smago_nu_ref, smago_kap_ref
    real(8), intent(out) :: kap_out, nu_out
    real(8), intent(in)  :: cfl_kap_max
    integer, intent(in)  :: my_rank, iic
    integer, intent(in)  :: Istr, Iend, Jstr, Jend, Nlev
    real(8) :: feat(NIN)
    real(8) :: log_kap, log_nu
    real(8) :: kap_nn, nu_nn
    real(8) :: ri_clip
    ! ---- Clip Richardson ----
    ri_clip = max(-5.0d0, min(5.0d0, ri_ij))
    ! ---- Build normalised vector ----
    feat(1) = (dudx2   - feat_mean(1)) / feat_std(1)
    feat(2) = (dvdy2   - feat_mean(2)) / feat_std(2)
    feat(3) = (defrate2 - feat_mean(3)) / feat_std(3)
    feat(4) = (bvf_ij  - feat_mean(4)) / feat_std(4)
    feat(5) = (ri_clip  - feat_mean(5)) / feat_std(5)
    feat(6) = (Hz_ij   - feat_mean(6)) / feat_std(6)
    ! ---- Forward pass ----
    call nn_forward(feat, log_kap, log_nu)
    ! ---- De-normalise ----
    kap_nn = exp(kap_log_mean + kap_log_std * log_kap)
    nu_nn  = exp(nu_log_mean  + nu_log_std  * log_nu )

    !**************************************************************
    ! ---- Blend with Smagorinsky ref.----
    kap_nn = alpha_hmix * kap_nn + (1.0d0 - alpha_hmix) * smago_kap_ref
    nu_nn  = alpha_hmix * nu_nn  + (1.0d0 - alpha_hmix) * smago_nu_ref
    ! ---- floor and CFL cap ----
    kap_out = max(coef_min, min(kap_nn, min(coef_max, cfl_kap_max)))
    nu_out  = max(coef_min, min(nu_nn,  min(coef_max, cfl_kap_max)))
  end subroutine nn_hmix_predict

  !=================================================================
  ! nn_hmix_buffer
  !=================================================================
  !*****************************************************************
  !     You may try first_call = .true. , default is .false.
  !*****************************************************************
  subroutine nn_hmix_buffer(my_rank, iic)
    integer, intent(in) :: my_rank, iic

    if (first_call .and. my_rank == 0) then
      write(6,'(A,I8)') &
           '[NN_HMIX] nn_hmix_buffer: first activation at iic=', iic
      write(6,'(A,F6.3)') &
           '[NN_HMIX] blending weight alpha_hmix = ', alpha_hmix
      first_call = .false.
    end if
  end subroutine nn_hmix_buffer

end module nn_hmix_coef_mod
