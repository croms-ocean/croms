!=======================================================================
!  CROCO  --  TGRU_VMIX
!  Transformer + GRU Vertical Mixing Scheme
!-----------------------------------------------------------------------
!  Implemented by :: S. Maishal, 2026
!  Indian Institute of Technology Kharagpur, India
!=======================================================================

module tgru_vmix_mod
  implicit none
  private
!**********************************************************************
  !=====  Architecture Hyper-Parameters ===============================
  integer, parameter :: d_feat  = 6    ! Features per interface
  integer, parameter :: d_model = 8    ! Embedding dim
  integer, parameter :: n_head  = 2    ! Attention in my heads
  integer, parameter :: d_head  = 4    ! = d_model / n_head
  integer, parameter :: d_gru   = 8    ! GRU hidden size
  integer, parameter :: d_out   = 2    ! Akv_scale, Akt_scale

  !=====  Physical constants ============================================
  real(8), parameter :: grav     = 9.81d0
  real(8), parameter :: rho0_ref = 1025.d0
  real(8), parameter :: eps_dbl  = 1.d-14

  !=====  Mixing bounds =================================================
  real(8), parameter :: Akv_background = 1.d-5
  real(8), parameter :: Akt_background = 1.d-6
  real(8), parameter :: Akv_maximum    = 5.d-2
  real(8), parameter :: Akt_maximum    = 5.d-2
  real(8), parameter :: Akv_convective = 1.d-1
  real(8), parameter :: Akt_convective = 1.d-2
  real(8), parameter :: Ri_crit        = 0.7d0
  real(8), parameter :: nu0m           = 5.d-2

  !The internal-wave mixing (very stable conditions)
  real(8), parameter :: wave_bkg_Akv   = 1.d-4
  real(8), parameter :: wave_bkg_Akt   = 1.d-5

  ! Galperin mixing-length cap
  real(8), parameter :: galp           = 0.53d0
  real(8), parameter :: N2_scale       = 1.d-4

  ! Blinded GRU to shear
  real(8), parameter :: Sh2_scale      = 1.d-6

  !             !!!  YOU MAY CHANGE  !!!
  ! Spin-up ramp: 6-hour exponential ramp to prevent cold-start shock
  real(8), parameter :: ramp_tau       = 3600.d0 * 6.d0

  ! Initial mixing boost parameters
  real(8), parameter :: boost_tau       = 3600.d0 * 12.d0 ! boost decays over 12 hours
  real(8), parameter :: Akv_boost_init  = 2.d-2           ! initial boost
  real(8), parameter :: Akt_boost_init  = 5.d-3           ! initial boost for tracers OK

  ! Time-smoothing coefficient (EMA)
  real(8), parameter :: smooth_alpha   = 0.5d0

  ! Persistent Akv/Akt for time smoothing
  real(8), allocatable, save :: Akv_old(:,:,:)
  real(8), allocatable, save :: Akt_old(:,:,:)
  logical, save :: smooth_arrays_allocated = .false.

  ! Tighter convective blending
  real(8), parameter :: N2_conv_width  = 1.d-7

  ! Shift sigmoid midpoint
  real(8), parameter :: N2_conv_offset = 2.944d-7

  ! Normalisation ceiling
  real(8), parameter :: log_N2_max     = 7.d0

  ! ------------------------TODO sm-loop --------------------------------
  real(8), parameter :: dz_min         = 1.d-3
  real(8), save :: W_enc(d_model, d_feat)
  real(8), save :: b_enc(d_model)
  real(8), save :: W_Q(d_head, d_model, n_head)
  real(8), save :: W_K(d_head, d_model, n_head)
  real(8), save :: W_V(d_head, d_model, n_head)
  real(8), save :: W_O(d_model, d_model)
  real(8), save :: Wz(d_gru, d_model), Uz(d_gru, d_gru), bz(d_gru)
  real(8), save :: Wr(d_gru, d_model), Ur(d_gru, d_gru), br(d_gru)
  real(8), save :: Wn(d_gru, d_model), Un(d_gru, d_gru), bn(d_gru)
  real(8), save :: W_dec(d_out, d_gru)
  real(8), save :: b_dec(d_out)

  logical, save :: tgru_initialised = .false.

  !--- Public API -------------------------------------------------------
  public :: tgru_vmix_init_impl
  public :: tgru_vmix_column
  public :: tgru_vmix_tile_exec

contains

!======================================================================
!  Physics-informed weight initialisation
!======================================================================
  subroutine tgru_vmix_init_impl()
    implicit none
    integer :: i, j, h
    real(8) :: sc

    if (tgru_initialised) return

    ! --- Encoder -------------------------------------------------------
    W_enc = 0.d0 ;  b_enc = 0.d0
    W_enc(1,1) =  2.0d0 ;  W_enc(2,1) = -1.0d0
    W_enc(3,2) =  2.0d0 ;  W_enc(4,2) = -1.0d0
    W_enc(5,3) =  1.0d0 ;  W_enc(6,3) = -1.0d0
    W_enc(7,4) = -2.0d0 ! Richardson number
    W_enc(8,5) =  1.0d0
    W_enc(1,6) =  0.2d0

    ! --- Attention matrices-----------------------------------------------
    sc = 1.d0 / sqrt(dble(d_head))
    do h = 1, n_head
      W_Q(:,:,h) = 0.d0
      W_K(:,:,h) = 0.d0
      W_V(:,:,h) = 0.d0
      do i = 1, d_head
        do j = 1, d_model
          if (i == mod(j-1, d_head)+1) then
            W_Q(i,j,h) = sc
            W_K(i,j,h) = sc
            W_V(i,j,h) = sc
          end if
        end do
      end do
    end do
    W_O = 0.d0
    do i = 1, d_model
      W_O(i,i) = 1.d0
    end do

    ! --- GRU weights ---------------------------------------------------
    Wz = 0.d0 ;  bz = 0.5d0
    Wr = 0.d0 ;  br = 1.0d0
    Wn = 0.d0 ;  bn = 0.d0

    Uz = 0.d0 ;  Ur = 0.d0 ;  Un = 0.d0
    do i = 1, d_gru
      Uz(i,i) = 0.1d0    ! [FIX-4] was 0
      Ur(i,i) = 0.1d0    ! [FIX-4] was 0
      Un(i,i) = 0.1d0    ! [FIX-4] was 0
    end do

    do i = 1, d_gru
      Wn(i, mod(i-1, d_model)+1) = 1.d0
    end do
    Wz(:,7) = -0.5d0

    ! --- Decoder -------------------------------------------------------
    W_dec    = 0.d0
    b_dec(1) = log(1.d-4)
    b_dec(2) = log(1.d-2)
    W_dec(1,7) = -0.75d0
    W_dec(2,7) = -0.60d0
    W_dec(1,3) =  0.8d0
    W_dec(2,3) =  0.8d0

    tgru_initialised = .true.
  end subroutine tgru_vmix_init_impl

!======================================================================
! Water column
!======================================================================
  subroutine tgru_vmix_column(Nz, model_time, Hz_col, u_col, v_col, rho_col,   &
                               z_r_col, z_w_col,                     &
                               Akv_col, Akt_col, bvf_col)
    implicit none
    integer,  intent(in)  :: Nz
    real(8),  intent(in)  :: model_time
    real(8),  intent(in)  :: Hz_col(1:Nz)
    real(8),  intent(in)  :: u_col(1:Nz), v_col(1:Nz)
    real(8),  intent(in)  :: rho_col(1:Nz)
    real(8),  intent(in)  :: z_r_col(1:Nz), z_w_col(0:Nz)
    real(8),  intent(out) :: Akv_col(0:Nz), Akt_col(0:Nz,2)
    real(8),  intent(out) :: bvf_col(0:Nz)

    real(8) :: N2(0:Nz), Sh2(0:Nz), Ri_k(0:Nz)
    real(8) :: feat(d_feat, 0:Nz)
    real(8) :: enc(d_model, 0:Nz), attn(d_model, 0:Nz)
    real(8) :: h_gru(d_gru), h_new(d_gru)
    real(8) :: z_gate(d_gru), r_gate(d_gru), n_gate(d_gru)
    real(8) :: dec_out(d_out)
    real(8) :: H_total, dz, drho, du, dv, shear2, n2k, dz_k
    real(8) :: Akv_ml, Akt_ml, w_conv, Akv_ri, Akt_ri
    real(8) :: spin_ramp
    real(8) :: Akv_nn, Akt_nn
    real(8) :: mix_boost
    real(8) :: Akv_with_boost, Akt_with_boost
    real(8) :: nu_sx, ratio_ri
    real(8) :: Akv_shear, Akt_shear
    real(8) :: Akv_gal
    integer :: k, f

    H_total    = max(sum(Hz_col), eps_dbl)
    mix_boost = Akv_boost_init * exp(-max(model_time, 0.d0) / boost_tau)
    spin_ramp  = 1.d0 - exp(-max(model_time, 0.d0) / ramp_tau)
    N2(0)      = 0.d0 ;  Sh2(0) = 0.d0
    Ri_k(0)    = 1.d10 ;  bvf_col(0) = 0.d0

    do k = 1, Nz-1
      dz = z_r_col(k+1) - z_r_col(k)
      if (abs(dz) < dz_min) dz = sign(dz_min, dz)
      drho       = rho_col(k+1) - rho_col(k)
      n2k        = -grav / rho0_ref * drho / dz
      N2(k)      = n2k
      bvf_col(k) = n2k
      du         = u_col(k+1) - u_col(k)
      dv         = v_col(k+1) - v_col(k)
      shear2     = (du*du + dv*dv) / (dz*dz)
      Sh2(k)     = shear2
      if (shear2 > eps_dbl) then
        Ri_k(k) = max(n2k, 0.d0) / shear2
      else
        Ri_k(k) = 1.d10
      end if
    end do

    if (Nz >= 2) then
      dz = 0.5d0*(z_w_col(Nz) - z_r_col(Nz)) + &
           0.5d0*(z_r_col(Nz) - z_r_col(Nz-1))
      if (abs(dz) < dz_min) dz = sign(dz_min, dz)
      drho       = rho_col(Nz) - rho_col(Nz-1)
      n2k        = -grav / rho0_ref * drho / dz
      N2(Nz)     = n2k
      bvf_col(Nz)= n2k
      du         = u_col(Nz) - u_col(Nz-1)
      dv         = v_col(Nz) - v_col(Nz-1)
      shear2     = (du*du + dv*dv) / (dz*dz)
      Sh2(Nz)    = shear2
      if (shear2 > eps_dbl) then
        Ri_k(Nz) = max(n2k, 0.d0) / shear2
      else
        Ri_k(Nz) = 1.d10
      end if
    else
      N2(Nz)      = N2(Nz-1)
      Sh2(Nz)     = Sh2(Nz-1)
      Ri_k(Nz)    = Ri_k(Nz-1)
      bvf_col(Nz) = N2(Nz)
    end if

    if (Nz >= 1) then
      dz = z_r_col(1) - z_w_col(0)
      if (abs(dz) < dz_min) dz = sign(dz_min, dz)
      N2(0)       = 0.d0
      bvf_col(0)  = 0.d0
      Sh2(0)      = 0.d0
      Ri_k(0)     = 1.d10
    end if

    ! --- Feature construction ------------------------------------------
    do k = 0, Nz
      if (k == 0) then
        dz_k = max(Hz_col(1),   eps_dbl)
      else if (k == Nz) then
        dz_k = max(Hz_col(Nz),  eps_dbl)
      else
        dz_k = max(0.5d0*(Hz_col(k) + Hz_col(k+1)), eps_dbl)
      end if

      feat(1,k) = tanh(sign(1.d0, N2(k)+eps_dbl) * &
                       log(1.d0 + abs(N2(k)) / N2_scale) / log_N2_max)

      feat(2,k) = tanh(Sh2(k) / Sh2_scale)
      feat(3,k) = z_w_col(k) / H_total
      feat(4,k) = -(min(Ri_k(k), Ri_crit) / Ri_crit)
      feat(5,k) = sign(1.d0, N2(k)+eps_dbl) * &
                  log(1.d0 + abs(N2(k)) / N2_scale)
      feat(6,k) = log(dz_k)
    end do

    ! --- Encoder: TODO sm-loop----------------------------
    do k = 0, Nz
      enc(:,k) = b_enc
      do f = 1, d_feat
        enc(:,k) = enc(:,k) + W_enc(:,f) * feat(f,k)
      end do
      enc(:,k) = tanh(enc(:,k))
    end do

    call mhsa(Nz, enc, attn)
    attn = attn + enc

    ! --- GRU sweep-----------------------------------------
    h_gru = 0.d0
    do k = Nz, 0, -1
      z_gate  = sigmoid_v(d_gru, matmul(Wz,attn(:,k)) + matmul(Uz,h_gru) + bz)
      r_gate  = sigmoid_v(d_gru, matmul(Wr,attn(:,k)) + matmul(Ur,h_gru) + br)
      n_gate  = tanh(matmul(Wn,attn(:,k)) + matmul(Un,r_gate*h_gru) + bn)
      h_new   = (1.d0 - z_gate) * h_gru + z_gate * n_gate
      h_gru   = h_new

      dec_out = b_dec + matmul(W_dec, h_gru)
      dec_out(1) = max(-13.8d0, min(-2.3d0, dec_out(1)))
      dec_out(2) = max(-13.8d0, min(-0.7d0, dec_out(2)))

      if (Ri_k(k) <= 0.d0) then
        nu_sx = 1.d0
      else if (Ri_k(k) < Ri_crit) then
        ratio_ri = Ri_k(k) / Ri_crit
        nu_sx    = (1.d0 - ratio_ri*ratio_ri)**3
      else
        nu_sx = 0.d0
      end if

    ! wave
      Akv_shear = exp(dec_out(1)) * nu_sx * nu0m
      Akt_shear = exp(dec_out(2)) * nu_sx * nu0m
      Akv_ri = Akv_background + wave_bkg_Akv + Akv_shear
      Akt_ri = Akt_background + wave_bkg_Akt + Akt_shear
      
      if (Ri_k(k) > eps_dbl .and. N2(k) > eps_dbl) then
        Akv_gal = galp * galp * Akv_ri / max(Ri_k(k), eps_dbl)
        Akv_ri  = min(Akv_ri, max(Akv_gal, Akv_background + wave_bkg_Akv))
        Akt_ri  = min(Akt_ri, max(Akv_gal, Akt_background + wave_bkg_Akt))
      end if

      w_conv = 1.d0 / (1.d0 + exp(min((N2(k)+N2_conv_offset) / N2_conv_width, 30.d0)))

      Akv_nn = (1.d0 - w_conv) * Akv_ri + w_conv * Akv_convective
      Akt_nn = (1.d0 - w_conv) * Akt_ri + w_conv * Akt_convective
      Akv_ml = spin_ramp * Akv_nn + (1.d0 - spin_ramp) * Akv_background
      Akt_ml = spin_ramp * Akt_nn + (1.d0 - spin_ramp) * Akt_background
      Akv_with_boost = Akv_ml + mix_boost
      Akt_with_boost = Akt_ml + mix_boost * (Akt_boost_init / Akv_boost_init)
      Akv_col(k)   = max(Akv_background, min(Akv_maximum, Akv_with_boost))
      Akt_col(k,1) = max(Akt_background, min(Akt_maximum, Akt_with_boost))
      Akt_col(k,2) = Akt_col(k,1)
    end do

    Akv_col(0)   = Akv_background
    Akt_col(0,1) = Akt_background
    Akt_col(0,2) = Akt_background

    if (Nz >= 2) then
      Akv_col(Nz)   = max(1.5d0*Akv_col(Nz-1) - 0.5d0*Akv_col(Nz-2), &
                          Akv_background)
      Akt_col(Nz,1) = max(1.5d0*Akt_col(Nz-1,1) - 0.5d0*Akt_col(Nz-2,1), &
                          Akt_background)
      Akt_col(Nz,2) = Akt_col(Nz,1)
    else
      Akv_col(Nz)   = Akv_background
      Akt_col(Nz,1) = Akt_background
      Akt_col(Nz,2) = Akt_background
    end if

  end subroutine tgru_vmix_column

!======================================================================
!  ---Multi-Head Scaled ***Dot-Product*** Self-Attention---
!======================================================================
  subroutine mhsa(Nz, X, out)
    implicit none
    integer, intent(in)  :: Nz
    real(8), intent(in)  :: X(d_model, 0:Nz)
    real(8), intent(out) :: out(d_model, 0:Nz)
    integer :: h, q, k
    real(8) :: Q_h(d_head,0:Nz), K_h(d_head,0:Nz), V_h(d_head,0:Nz)
    real(8) :: score(0:Nz), attn_w(0:Nz), cat(d_model,0:Nz)
    real(8) :: scale, rmax, rsum

    scale = 1.d0 / sqrt(dble(d_head))
    cat   = 0.d0 ;  out = 0.d0
    do h = 1, n_head
      do q = 0, Nz
        Q_h(:,q) = matmul(W_Q(:,:,h), X(:,q))
        K_h(:,q) = matmul(W_K(:,:,h), X(:,q))
        V_h(:,q) = matmul(W_V(:,:,h), X(:,q))
      end do
      do q = 0, Nz
        do k = 0, Nz
          score(k) = dot_product(Q_h(:,q), K_h(:,k)) * scale
        end do
        rmax = maxval(score(0:Nz)) ;  rsum = 0.d0
        do k = 0, Nz
          attn_w(k) = exp(score(k) - rmax)
          rsum      = rsum + attn_w(k)
        end do
        attn_w(0:Nz) = attn_w(0:Nz) / max(rsum, eps_dbl)
        do k = 0, Nz
          cat((h-1)*d_head+1:h*d_head, q) =                          &
            cat((h-1)*d_head+1:h*d_head, q) + attn_w(k)*V_h(:,k)
        end do
      end do
    end do
    do q = 0, Nz
      out(:,q) = matmul(W_O, cat(:,q))
    end do
  end subroutine mhsa

!======================================================================
!  Legacy smooth Richardson-number-based amplitude
!======================================================================
!!!--------------------------------------------------------------------
  pure function ri_mixing(Ri_val, log_scale) result(Ak)
    implicit none
    real(8), intent(in) :: Ri_val, log_scale
    real(8) :: Ak, alpha
    if (Ri_val <= 0.d0) then
      alpha = 1.d0
    else if (Ri_val < Ri_crit) then
      alpha = (1.d0 - Ri_val/Ri_crit)**3
    else
      alpha = 0.d0
    end if
    Ak = alpha * exp(max(min(log_scale, 0.d0), -15.d0))
  end function ri_mixing

!======================================================================
!  Element-Wise
!======================================================================
  pure function sigmoid_v(m, x) result(y)
    implicit none
    integer, intent(in) :: m
    real(8), intent(in) :: x(m)
    real(8) :: y(m)
    y = 1.d0 / (1.d0 + exp(-min(max(x, -30.d0), 30.d0)))
  end function sigmoid_v
!!!--------------------------------------------------------------------
!======================================================================
!  tgru_vmix_tile_exec
!======================================================================
  subroutine tgru_vmix_tile_exec(Istr, Iend, Jstr, Jend,            &
                                  Nz, nstep, i_temp, i_salt,          &
                                  model_time,                          &
                                  Lmx, Mmx,                           &
                                  Hz_a,                               &
                                  u_a, v_a,                           &
                                  rho_a, z_r_a, z_w_a,               &
                                  Akv_a, Akt_a)
    implicit none
    integer, intent(in) :: Istr, Iend, Jstr, Jend
    integer, intent(in) :: Nz, nstep, i_temp, i_salt, Lmx, Mmx
    real(8), intent(in) :: model_time
    real, intent(in)    :: Hz_a (0:Lmx+1, 0:Mmx+1, Nz)
    real, intent(in)    :: u_a  (0:Lmx+1, 0:Mmx+1, Nz, 3)
    real, intent(in)    :: v_a  (0:Lmx+1, 0:Mmx+1, Nz, 3)
    real, intent(in)    :: rho_a(0:Lmx+1, 0:Mmx+1, Nz)
    real, intent(in)    :: z_r_a(0:Lmx+1, 0:Mmx+1, Nz)
    real, intent(in)    :: z_w_a(0:Lmx+1, 0:Mmx+1, 0:Nz)
    real, intent(inout) :: Akv_a(0:Lmx+1, 0:Mmx+1, 0:Nz)
    real, intent(inout) :: Akt_a(0:Lmx+1, 0:Mmx+1, 0:Nz, 2)

    integer :: i, j, k
    real(8), allocatable :: Hz_col(:), u_col(:), v_col(:)
    real(8), allocatable :: rho_col(:), z_r_col(:), z_w_col(:)
    real(8), allocatable :: Akv_col(:), Akt_col(:,:), bvf_col(:)

    if (.not. tgru_initialised) call tgru_vmix_init_impl()

    if (.not. smooth_arrays_allocated) then
      allocate(Akv_old(0:Lmx+1, 0:Mmx+1, 0:Nz))
      allocate(Akt_old(0:Lmx+1, 0:Mmx+1, 0:Nz))
      Akv_old = Akv_background
      Akt_old = Akt_background
      smooth_arrays_allocated = .true.
    end if

    allocate(Hz_col(Nz),  u_col(Nz),   v_col(Nz))
    allocate(rho_col(Nz), z_r_col(Nz), z_w_col(0:Nz))
    allocate(Akv_col(0:Nz), Akt_col(0:Nz,2), bvf_col(0:Nz))

    do j = Jstr, Jend
      do i = Istr, Iend

        do k = 1, Nz
          Hz_col(k)  = dble(Hz_a(i, j, k))
          u_col(k)   = 0.5d0*(dble(u_a(i,  j,k,nstep)) +            &
                               dble(u_a(i+1,j,k,nstep)))
          v_col(k)   = 0.5d0*(dble(v_a(i,j,  k,nstep)) +            &
                               dble(v_a(i,j+1,k,nstep)))

! *******************---5-point cross-stencil---**************************
           rho_col(k) = (dble(rho_a(i-1,j,  k)) +                    &
                        dble(rho_a(i+1,j,  k)) +                     &
                        dble(rho_a(i,  j-1,k)) +                     &
                        dble(rho_a(i,  j+1,k)) +                     &
                        dble(rho_a(i,  j,  k))) / 5.d0
! *******************---5-point cross-stencil---**************************

          z_r_col(k) = dble(z_r_a(i, j, k))
        end do
        do k = 0, Nz
          z_w_col(k) = dble(z_w_a(i, j, k))
        end do

        call tgru_vmix_column(Nz, model_time, Hz_col, u_col, v_col, rho_col,    &
                               z_r_col, z_w_col,                      &
                               Akv_col, Akt_col, bvf_col)

        do k = 0, Nz
          Akv_col(k)   = smooth_alpha * Akv_col(k)   + (1.d0-smooth_alpha) * Akv_old(i,j,k)
          Akt_col(k,1) = smooth_alpha * Akt_col(k,1) + (1.d0-smooth_alpha) * Akt_old(i,j,k)
          Akt_col(k,2) = Akt_col(k,1)
          Akv_old(i,j,k) = Akv_col(k)
          Akt_old(i,j,k) = Akt_col(k,1)

          Akv_a(i,j,k)        = real(Akv_col(k))
          Akt_a(i,j,k,i_temp) = real(Akt_col(k,1))
          Akt_a(i,j,k,i_salt) = real(Akt_col(k,2))
        end do

      end do
    end do

    deallocate(Hz_col, u_col, v_col, rho_col, z_r_col, z_w_col)
    deallocate(Akv_col, Akt_col, bvf_col)
  end subroutine tgru_vmix_tile_exec

end module tgru_vmix_mod
