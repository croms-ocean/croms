!========================================================================*!
!*************************************************************************!
!*   CR*** - OBC_M3ORLANSKI_LSTM_PSC                                     *!
!*-----------------------------------------------------------------------*!
!*   Implemented by :: S. Maishal, 2026                                  *!
!*   Indian Institute of Technology Kharagpur, India                     *!
!======================================================================= *!
!*   DESCRIPTION:                                                        *!
!*    ML-based Phase Speed Correction (PSC) + LSTM Temporal Boundary     *!
!*    Memory for 3-D momentum open boundary conditions.                  *!
!=========================================================================!
!*    This module replaces / augments the classic Orlanski radiation     *!
!*    formula for the M3 (3-D baroclinic) velocity OBC by:               *!
!*                                                                       *!
!*    (1) ML-PSC  - A compact, single-hidden-layer feed-forward network  *!
!*        that corrects the instantaneous Orlanski phase speed (cx, cy)  *!
!*        using the local velocity shear, density gradient, and Froude   *!
!*        number as predictors.                                          *!
!*                                                                       *!
!*    (2) LSTM-TBM - A per-boundary-point Long Short-Term Memory cell    *!
!*        (single layer, hidden size = LSTM_HS) that accumulates a       *!
!*        temporal memory of the boundary state over time.  Hidden and   *!
!*        cell states are MPI-aware.                                     *!
!*                                                                       *!
!*    (3) Blending:                                                      *!
!*          u_obc = alpha * u_lstm  +  (1-alpha) * u_orlanski_psc        *!
!*        where alpha ramps over LSTM_WARMUP_STEPS.                      *!
!*                                                                       *!
!*          *** RUNTIME ONLINE and OFFLINE WEIGHT LOADING ***            *!
!========================================================================*!
!*  DESIGN PRINCIPLES:                                                   *!
!*    - Everything lives here.                                           *!
!*    - F90 module with SAVE: LSTM state persists across calls.          *!
!*    - MPI collective: MPI_Allreduce after each baroclinic step.        *!
!*    - Zero external dependencies: no BLAS, no NetCDF, no Python, etc.  *!
!*    - Thread-safe: OpenMP CRITICAL guards on LSTM state updates.       *!
!*  LICENSE: https://www.gnu.org/licenses/gpl-3.0.en.html                *!
!*************************************************************************!
!=========================================================================!

module obc_orlanski_lstm_psc_m3

  implicit none
  save

  include 'mpif.h'

  !--------------------------------------------------------------------
  ! Tunable hyper-parameters ( !!!  -YOU MAY CHANGE HERE-  !!! )
  !--------------------------------------------------------------------
  integer, parameter :: LSTM_HS           = 8   ! LSTM hidden size
  integer, parameter :: LSTM_SEQ          = 5   ! temporal window (INFO-4U)
  integer, parameter :: LSTM_WARMUP_STEPS = 10  ! steps before full blending
  integer, parameter :: ML_HIDDEN         = 6   ! PSC net hidden size
  integer, parameter :: ML_IN             = 3   ! PSC net input size
  integer, parameter :: LSTM_IN           = 2   ! LSTM input size (for u and v)
  real,    parameter :: LSTM_BLEND_MIN    = 0.0
  real,    parameter :: LSTM_BLEND_MAX    = 1.0
  real,    parameter :: EPS_LSTM          = 1.E-20

  !--------------------------------------------------------------------
  ! Step counter
  !--------------------------------------------------------------------
  integer :: lstm_step_count = 0

  !--------------------------------------------------------------------
  ! LSTM state (wall, boundary_pt, level, hidden_unit)
  ! Wall: 1=West 2=East 3=South 4=North
  ! Increase MAX_BDRY_PTS if your grid exceeds 4096 points on one edge.
  ! *******************************************************************
  !--------------------------------------------------------------------
  integer, parameter :: MAX_BDRY_PTS = 4096 !Your Pts (sm-loop change it)
  integer, parameter :: MAX_N        = 100

  real :: lstm_h(4, MAX_BDRY_PTS, MAX_N, LSTM_HS)
  real :: lstm_c(4, MAX_BDRY_PTS, MAX_N, LSTM_HS)

  !--------------------------------------------------------------------
  ! Flag: Have the weights been loaded this run ???
  !--------------------------------------------------------------------
  logical :: lstm_weights_loaded = .false.

  !--------------------------------------------------------------------
  ! PSC feed-forward weights (-*-mini-*-)
  !
  !   Network:  x(ML_IN) -> Dense(ML_HIDDEN, tanh) -> Dense(2, tanh)
  !   Outputs:  [delta_cx, delta_cy] – additive corrections to Orlanski
  !             phase speeds, normalised to [-1,1].
  !--------------------------------------------------------------------
  real :: W1(ML_HIDDEN, ML_IN)
  real :: b1(ML_HIDDEN)
  real :: W2(2, ML_HIDDEN)
  real :: b2(2)

  !--------------------------------------------------------------------
  !   LSTM weight matrices
  !--------------------------------------------------------------------
  real :: LSTM_WI_X(LSTM_HS, LSTM_IN)
  real :: LSTM_WI_H(LSTM_HS, LSTM_HS)
  real :: LSTM_BI(LSTM_HS)

  real :: LSTM_WF_X(LSTM_HS, LSTM_IN)
  real :: LSTM_WF_H(LSTM_HS, LSTM_HS)
  real :: LSTM_BF(LSTM_HS)

  real :: LSTM_WG_X(LSTM_HS, LSTM_IN)
  real :: LSTM_WG_H(LSTM_HS, LSTM_HS)
  real :: LSTM_BG(LSTM_HS)

  real :: LSTM_WO_X(LSTM_HS, LSTM_IN)
  real :: LSTM_WO_H(LSTM_HS, LSTM_HS)
  real :: LSTM_BO(LSTM_HS)

  real :: LSTM_WOUT(LSTM_HS)
  real :: LSTM_BOUT

contains

  !=======================================================================
  !  Activation ::
  !=======================================================================

  pure elemental real function sigmoid(x)
    real, intent(in) :: x
    sigmoid = 1.0 / (1.0 + exp(-x))
  end function sigmoid

  pure elemental real function tanh_act(x)
    real, intent(in) :: x
    real :: ex2
    ex2 = exp(-2.0*x)
    tanh_act = (1.0 - ex2) / (1.0 + ex2 + EPS_LSTM)
  end function tanh_act

  !=======================================================================
  !  lstm_weights_init_defaults
  !  (Xavier-uniform)
  !=======================================================================

  subroutine lstm_weights_init_defaults()

    real :: scale_psc_W1, scale_psc_W2
    real :: scale_xi, scale_xh, scale_out
    integer :: i, j

    ! Xavier uniform scale factors
    !-----------------------------
    scale_psc_W1 = sqrt(6.0 / real(ML_IN    + ML_HIDDEN))
    scale_psc_W2 = sqrt(6.0 / real(ML_HIDDEN + 2))
    scale_xi     = sqrt(6.0 / real(LSTM_IN  + LSTM_HS))
    scale_xh     = sqrt(6.0 / real(LSTM_HS  + LSTM_HS))
    scale_out    = sqrt(6.0 / real(LSTM_HS  + 1))

    ! PSC network  --  Xavier bounds
    !-------------------------------
    do i = 1, ML_HIDDEN
      do j = 1, ML_IN
        W1(i,j) = scale_psc_W1 * (2.0*real(mod(i*ML_IN+j,7))/6.0 - 1.0)
      enddo
      b1(i) = 0.0
    enddo
    do i = 1, 2
      do j = 1, ML_HIDDEN
        W2(i,j) = scale_psc_W2 * (2.0*real(mod(i*ML_HIDDEN+j,5))/4.0 - 1.0)
      enddo
      b2(i) = 0.0
    enddo

    ! LSTM gates --  Xavier pattern
    !------------------------------
    do i = 1, LSTM_HS
      do j = 1, LSTM_IN
        LSTM_WI_X(i,j) = scale_xi*(2.0*real(mod(i*LSTM_IN +j, 7))/6.0-1.0)
        LSTM_WF_X(i,j) = scale_xi*(2.0*real(mod(i*LSTM_IN +j+1,7))/6.0-1.0)
        LSTM_WG_X(i,j) = scale_xi*(2.0*real(mod(i*LSTM_IN +j+2,7))/6.0-1.0)
        LSTM_WO_X(i,j) = scale_xi*(2.0*real(mod(i*LSTM_IN +j+3,7))/6.0-1.0)
      enddo
      do j = 1, LSTM_HS
        LSTM_WI_H(i,j) = scale_xh*(2.0*real(mod(i*LSTM_HS+j,   11))/10.0-1.0)
        LSTM_WF_H(i,j) = scale_xh*(2.0*real(mod(i*LSTM_HS+j+1, 11))/10.0-1.0)
        LSTM_WG_H(i,j) = scale_xh*(2.0*real(mod(i*LSTM_HS+j+2, 11))/10.0-1.0)
        LSTM_WO_H(i,j) = scale_xh*(2.0*real(mod(i*LSTM_HS+j+3, 11))/10.0-1.0)
      enddo
      LSTM_BI(i)   = 0.0
      LSTM_BF(i)   = 1.0
      LSTM_BG(i)   = 0.0
      LSTM_BO(i)   = 0.0
      LSTM_WOUT(i) = scale_out*(2.0*real(mod(i,5))/4.0 - 1.0)
    enddo
    LSTM_BOUT = 0.0

  end subroutine lstm_weights_init_defaults

  !=======================================================================
  !  lstm_weights_load
  !  TODO sm-loop for CUDA--
  !  Read all PSC and LSTM weights from a plain-text ASCII file at runtime.
  !=======================================================================

  subroutine lstm_weights_load(comm, mynode)

    integer, intent(in) :: comm, mynode

    character(len=256) :: wfile
    character(len=64)  :: tag
    integer :: iunit, ios, ierr
    logical :: file_exists

    ! ------------------------------------------------------------------
    ! Initialise defaults
    ! ------------------------------------------------------------------
    call lstm_weights_init_defaults()

    ! ------------------------------------------------------------------
    !         Rank 0 determines the file path.
    !         We use get_environment_variable (F2003, gfortran, ifort/ifx,
    !                                                   pgf90, nvfortran).
    ! ------------------------------------------------------------------
    if (mynode .eq. 0) then

      call get_environment_variable('LSTM_WEIGHT_FILE', wfile, status=ios)
      if (ios .ne. 0 .or. len_trim(wfile) .eq. 0) then
        wfile = 'lstm_weights.txt'
      endif

      inquire(file=trim(wfile), exist=file_exists)

      if (.not. file_exists) then
        write(6,'(A)') '===================================================='
        write(6,'(A)') ' WARNING: LSTM weight file not found:'
        write(6,'(2A)')' WARNING:   ', trim(wfile)
        write(6,'(A)') ' WARNING: Using case-neutral Xavier defaults.'
        write(6,'(A)') ' WARNING: Train on your grid and set LSTM_WEIGHT_FILE'
        write(6,'(A)') ' WARNING: or place lstm_weights.txt in run directory.'
        write(6,'(A)') '===================================================='
        ! Broadcast flag and return
        ios = 0
        call MPI_Bcast(ios, 1, MPI_INTEGER, 0, comm, ierr)
        lstm_weights_loaded = .false.
        return
      endif

      write(6,'(2A)') ' LSTM: loading weights from ', trim(wfile)

      iunit = 99
      open(unit=iunit, file=trim(wfile), status='old', action='read', &
           iostat=ios)
      if (ios .ne. 0) then
        write(6,'(A,I4)') ' LSTM: ERROR opening weight file, iostat=', ios
        ios = 0
        call MPI_Bcast(ios, 1, MPI_INTEGER, 0, comm, ierr)
        lstm_weights_loaded = .false.
        return
      endif

      do
        read(iunit, *, iostat=ios) tag
        if (ios .ne. 0) exit
        if (tag(1:1) .eq. '#' .or. len_trim(tag) .eq. 0) cycle
        select case (trim(tag))

          case ('PSC_W1')
            read(iunit, *, iostat=ios) W1
          case ('PSC_b1')
            read(iunit, *, iostat=ios) b1
          case ('PSC_W2')
            read(iunit, *, iostat=ios) W2
          case ('PSC_b2')
            read(iunit, *, iostat=ios) b2

          case ('LSTM_WI_X')
            read(iunit, *, iostat=ios) LSTM_WI_X
          case ('LSTM_WI_H')
            read(iunit, *, iostat=ios) LSTM_WI_H
          case ('LSTM_BI')
            read(iunit, *, iostat=ios) LSTM_BI

          case ('LSTM_WF_X')
            read(iunit, *, iostat=ios) LSTM_WF_X
          case ('LSTM_WF_H')
            read(iunit, *, iostat=ios) LSTM_WF_H
          case ('LSTM_BF')
            read(iunit, *, iostat=ios) LSTM_BF

          case ('LSTM_WG_X')
            read(iunit, *, iostat=ios) LSTM_WG_X
          case ('LSTM_WG_H')
            read(iunit, *, iostat=ios) LSTM_WG_H
          case ('LSTM_BG')
            read(iunit, *, iostat=ios) LSTM_BG

          case ('LSTM_WO_X')
            read(iunit, *, iostat=ios) LSTM_WO_X
          case ('LSTM_WO_H')
            read(iunit, *, iostat=ios) LSTM_WO_H
          case ('LSTM_BO')
            read(iunit, *, iostat=ios) LSTM_BO

          case ('LSTM_WOUT')
            read(iunit, *, iostat=ios) LSTM_WOUT
          case ('LSTM_BOUT')
            read(iunit, *, iostat=ios) LSTM_BOUT

          case default
            ! I dont know who are. 
        end select

        if (ios .ne. 0) then
          write(6,'(2A,I4)') ' LSTM: read error on tag ', trim(tag), ios
          exit
        endif

      enddo

      close(iunit)
      write(6,'(A)') ' LSTM: weight file loaded successfully.'
      ios = 1

    endif

    ! ------------------------------------------------------------------
    ! Broadcast load-success //
    ! ------------------------------------------------------------------
    call MPI_Bcast(ios, 1, MPI_INTEGER, 0, comm, ierr)

    ! ------------------------------------------------------------------
    ! Broadcast ALL weight //
    ! ------------------------------------------------------------------
    call MPI_Bcast(W1,         ML_HIDDEN*ML_IN,         MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(b1,         ML_HIDDEN,               MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(W2,         2*ML_HIDDEN,             MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(b2,         2,                       MPI_DOUBLE_PRECISION, 0, comm, ierr)

    call MPI_Bcast(LSTM_WI_X,  LSTM_HS*LSTM_IN,         MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(LSTM_WI_H,  LSTM_HS*LSTM_HS,         MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(LSTM_BI,    LSTM_HS,                 MPI_DOUBLE_PRECISION, 0, comm, ierr)

    call MPI_Bcast(LSTM_WF_X,  LSTM_HS*LSTM_IN,         MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(LSTM_WF_H,  LSTM_HS*LSTM_HS,         MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(LSTM_BF,    LSTM_HS,                 MPI_DOUBLE_PRECISION, 0, comm, ierr)

    call MPI_Bcast(LSTM_WG_X,  LSTM_HS*LSTM_IN,         MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(LSTM_WG_H,  LSTM_HS*LSTM_HS,         MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(LSTM_BG,    LSTM_HS,                 MPI_DOUBLE_PRECISION, 0, comm, ierr)

    call MPI_Bcast(LSTM_WO_X,  LSTM_HS*LSTM_IN,         MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(LSTM_WO_H,  LSTM_HS*LSTM_HS,         MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(LSTM_BO,    LSTM_HS,                 MPI_DOUBLE_PRECISION, 0, comm, ierr)

    call MPI_Bcast(LSTM_WOUT,  LSTM_HS,                 MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(LSTM_BOUT,  1,                       MPI_DOUBLE_PRECISION, 0, comm, ierr)

    lstm_weights_loaded = (ios .eq. 1)
    lstm_h          = 0.0
    lstm_c          = 0.0
    lstm_step_count = 0

  end subroutine lstm_weights_load

  !=======================================================================
  !  Phase Speed Correction ::
  !=======================================================================
  subroutine psc_correct(dft, dfx, dfy, cff, cx_orl, cy_orl, &
                         cx_ml, cy_ml)

    real, intent(in)  :: dft, dfx, dfy, cff
    real, intent(in)  :: cx_orl, cy_orl
    real, intent(out) :: cx_ml, cy_ml

    real    :: x_in(ML_IN), h1(ML_HIDDEN), out(2)
    real    :: cff_norm
    real,    parameter :: delta_scale = 0.20
    integer :: i, j

    cff_norm = max(sqrt(cff), EPS_LSTM)
    x_in(1)  = max(-3.0, min(3.0, dft / cff_norm))
    x_in(2)  = max(-3.0, min(3.0, dfx / cff_norm))
    x_in(3)  = max(-3.0, min(3.0, dfy / cff_norm))

    do i = 1, ML_HIDDEN
      h1(i) = b1(i)
      do j = 1, ML_IN
        h1(i) = h1(i) + W1(i,j) * x_in(j)
      enddo
      h1(i) = tanh_act(h1(i))
    enddo

    do i = 1, 2
      out(i) = b2(i)
      do j = 1, ML_HIDDEN
        out(i) = out(i) + W2(i,j) * h1(j)
      enddo
      out(i) = tanh_act(out(i))
    enddo

    cx_ml = cx_orl + delta_scale * out(1)
    cy_ml = cy_orl + delta_scale * out(2)

  end subroutine psc_correct

  !=======================================================================
  !  lstm_update
  !  LSTM cell
  !  Returns scalar velocity correction u_corr [m/s].
  !=======================================================================

  subroutine lstm_update(iwall, j_idx, k_idx, u_val, v_val, u_corr)

    integer, intent(in)  :: iwall, j_idx, k_idx
    real,    intent(in)  :: u_val, v_val
    real,    intent(out) :: u_corr

    real :: x_in(LSTM_IN)
    real :: h_prev(LSTM_HS), c_prev(LSTM_HS)
    real :: gate_i(LSTM_HS), gate_f(LSTM_HS)
    real :: gate_g(LSTM_HS), gate_o(LSTM_HS)
    real :: c_new(LSTM_HS),  h_new(LSTM_HS)
    integer :: i, p

    h_prev(:) = lstm_h(iwall, j_idx, k_idx, :)
    c_prev(:) = lstm_c(iwall, j_idx, k_idx, :)

    x_in(1) = max(-5.0, min(5.0, u_val))
    x_in(2) = max(-5.0, min(5.0, v_val))

    do i = 1, LSTM_HS
      gate_i(i) = LSTM_BI(i)
      gate_f(i) = LSTM_BF(i)
      gate_g(i) = LSTM_BG(i)
      gate_o(i) = LSTM_BO(i)
      do p = 1, LSTM_IN
        gate_i(i) = gate_i(i) + LSTM_WI_X(i,p) * x_in(p)
        gate_f(i) = gate_f(i) + LSTM_WF_X(i,p) * x_in(p)
        gate_g(i) = gate_g(i) + LSTM_WG_X(i,p) * x_in(p)
        gate_o(i) = gate_o(i) + LSTM_WO_X(i,p) * x_in(p)
      enddo
      do p = 1, LSTM_HS
        gate_i(i) = gate_i(i) + LSTM_WI_H(i,p) * h_prev(p)
        gate_f(i) = gate_f(i) + LSTM_WF_H(i,p) * h_prev(p)
        gate_g(i) = gate_g(i) + LSTM_WG_H(i,p) * h_prev(p)
        gate_o(i) = gate_o(i) + LSTM_WO_H(i,p) * h_prev(p)
      enddo
      gate_i(i) = sigmoid(gate_i(i))
      gate_f(i) = sigmoid(gate_f(i))
      gate_g(i) = tanh_act(gate_g(i))
      gate_o(i) = sigmoid(gate_o(i))
    enddo

    do i = 1, LSTM_HS
      c_new(i) = gate_f(i)*c_prev(i) + gate_i(i)*gate_g(i)
      h_new(i) = gate_o(i) * tanh_act(c_new(i))
    enddo

    lstm_h(iwall, j_idx, k_idx, :) = h_new(:)
    lstm_c(iwall, j_idx, k_idx, :) = c_new(:)

    u_corr = LSTM_BOUT
    do i = 1, LSTM_HS
      u_corr = u_corr + LSTM_WOUT(i) * h_new(i)
    enddo
    u_corr = max(-0.5, min(0.5, u_corr))

  end subroutine lstm_update

  !=======================================================================
  !  lstm_mpi_sync
  !  MPI_Allreduce to synchronise LSTM
  !=======================================================================

  subroutine lstm_mpi_sync(comm)

    integer, intent(in) :: comm

    real, save :: buf_h(4, MAX_BDRY_PTS, MAX_N, LSTM_HS)
    real, save :: buf_c(4, MAX_BDRY_PTS, MAX_N, LSTM_HS)
    integer :: ierr, ntot

    ntot = 4 * MAX_BDRY_PTS * MAX_N * LSTM_HS

    call MPI_Allreduce(lstm_h, buf_h, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    call MPI_Allreduce(lstm_c, buf_c, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

    lstm_h = buf_h
    lstm_c = buf_c

    lstm_step_count = lstm_step_count + 1

  end subroutine lstm_mpi_sync

  !=======================================================================
  !  lstm_reset
  !  Fresh weights on restart.
  !=======================================================================
  subroutine lstm_reset()
    lstm_h(:,:,:,:) = 0.0
    lstm_c(:,:,:,:) = 0.0
    lstm_step_count = 0
  end subroutine lstm_reset

  !=======================================================================
  !  obc_orlanski_lstm_psc_m3_u
  !  Apply PSC + LSTM to XI-velocity u at all 4 open boundaries
  !  Called from u3dbc_tile when OBC_M3ORLANSKI_LSTM_PSC is defined.
  !=======================================================================

  subroutine obc_orlanski_lstm_psc_m3_u(Istr, Iend, Jstr, Jend, &
                             nstp, nnew, N_lev,        &
                             eps_in, tau_in, tau_out,   &
                             u_arr, v_arr,              &
                             grad,                      &
                             umask_arr, pmask_arr,       &
                             uclm_arr)

    implicit none

    integer, intent(in) :: Istr, Iend, Jstr, Jend
    integer, intent(in) :: nstp, nnew, N_lev
    real,    intent(in) :: eps_in, tau_in, tau_out

    real, intent(inout) :: u_arr(0:,0:,:,:)
    real, intent(in)    :: v_arr(0:,0:,:,:)
    real, intent(inout) :: grad(0:,0:)
    real, intent(in)    :: umask_arr(0:,0:)
    real, intent(in)    :: pmask_arr(0:,0:)
    real, intent(in)    :: uclm_arr(0:,0:,:)
!-----------------------
#ifdef M3_FRC_BRY
# include "boundary.h"
#endif
!-----------------------
    integer :: i, j, k, j_idx
    real    :: dft, dfx, dfy, cff
    real    :: cx_orl, cy_orl, cx_ml, cy_ml
    real    :: tau, u_orl, u_lstm_corr, alpha, u_final
    if (lstm_step_count >= LSTM_WARMUP_STEPS) then
      alpha = LSTM_BLEND_MAX
    else
      alpha = LSTM_BLEND_MIN + &
              (LSTM_BLEND_MAX - LSTM_BLEND_MIN) * &
              real(lstm_step_count) / real(max(1, LSTM_WARMUP_STEPS))
    endif

    !==================================================================
    ! WESTERN BOUNDARY
    !==================================================================
    if (.true.) then

      do k = 1, N_lev
        do j = Jstr, Jend+1
          grad(Istr  ,j) = (u_arr(Istr  ,j,k,nstp) - &
                            u_arr(Istr  ,j-1,k,nstp)) * pmask_arr(Istr,j)
          grad(Istr+1,j) = (u_arr(Istr+1,j,k,nstp) - &
                            u_arr(Istr+1,j-1,k,nstp)) * pmask_arr(Istr+1,j)
        enddo

        do j = Jstr, Jend
          j_idx = j - Jstr + 1

          dft = u_arr(Istr+1,j,k,nstp) - u_arr(Istr+1,j,k,nnew)
          dfx = u_arr(Istr+1,j,k,nnew) - u_arr(Istr+2,j,k,nnew)

          tau = tau_out
          if (dfx*dft .lt. 0.0) then ; dft = 0.0 ; tau = tau_in ; endif

          if (dft*(grad(Istr+1,j)+grad(Istr+1,j+1)) .gt. 0.0) then
            dfy = grad(Istr+1,j)
          else
            dfy = grad(Istr+1,j+1)
          endif

          cff    = max(dfx*dfx + dfy*dfy, eps_in)
          cx_orl = dft*dfx
          cy_orl = min(cff, max(dft*dfy, -cff))

          call psc_correct(dft, dfx, dfy, cff, cx_orl, cy_orl, cx_ml, cy_ml)

          u_orl = (  cff  * u_arr(Istr  ,j,k,nstp)     &
                   + cx_ml* u_arr(Istr+1,j,k,nnew)     &
                   - max(cy_ml,0.0)*grad(Istr,j  )      &
                   - min(cy_ml,0.0)*grad(Istr,j+1)      &
                  ) / (cff + cx_ml + eps_in)

          !$OMP CRITICAL (lstm_west_u)
          call lstm_update(1, j_idx, k, &
                           u_arr(Istr,j,k,nstp), v_arr(Istr,j,k,nstp), &
                           u_lstm_corr)
          !$OMP END CRITICAL (lstm_west_u)

          u_final = (1.0 - alpha)*u_orl + alpha*(u_orl + u_lstm_corr)
!--------------------------------------------------------------------
#ifdef M3_FRC_BRY
          u_final = (1.0-tau)*u_final + tau*ubry_west(j,k)
#elif defined M3CLIMATOLOGY
          u_final = (1.0-tau)*u_final + tau*uclm_arr(Istr,j,k)
#endif
!--------------------------------------------------------------------
          u_arr(Istr,j,k,nnew) = u_final * umask_arr(Istr,j)

        enddo
      enddo

    endif

    !==================================================================
    ! EASTERN BOUNDARY
    !==================================================================
    if (.true.) then

      do k = 1, N_lev
        do j = Jstr, Jend+1
          grad(Iend  ,j) = (u_arr(Iend  ,j,k,nstp) - &
                            u_arr(Iend  ,j-1,k,nstp)) * pmask_arr(Iend,j)
          grad(Iend+1,j) = (u_arr(Iend+1,j,k,nstp) - &
                            u_arr(Iend+1,j-1,k,nstp)) * pmask_arr(Iend+1,j)
        enddo

        do j = Jstr, Jend
          j_idx = j - Jstr + 1

          dft = u_arr(Iend,j,k,nstp) - u_arr(Iend  ,j,k,nnew)
          dfx = u_arr(Iend,j,k,nnew) - u_arr(Iend-1,j,k,nnew)

          tau = tau_out
          if (dfx*dft .lt. 0.0) then ; dft = 0.0 ; tau = tau_in ; endif

          if (dft*(grad(Iend,j)+grad(Iend,j+1)) .gt. 0.0) then
            dfy = grad(Iend,j)
          else
            dfy = grad(Iend,j+1)
          endif

          cff    = max(dfx*dfx + dfy*dfy, eps_in)
          cx_orl = dft*dfx
          cy_orl = min(cff, max(dft*dfy, -cff))

          call psc_correct(dft, dfx, dfy, cff, cx_orl, cy_orl, cx_ml, cy_ml)

          u_orl = (  cff  * u_arr(Iend+1,j,k,nstp)     &
                   + cx_ml* u_arr(Iend  ,j,k,nnew)      &
                   - max(cy_ml,0.0)*grad(Iend+1,j  )    &
                   - min(cy_ml,0.0)*grad(Iend+1,j+1)    &
                  ) / (cff + cx_ml + eps_in)

          !$OMP CRITICAL (lstm_east_u)
          call lstm_update(2, j_idx, k, &
                           u_arr(Iend+1,j,k,nstp), v_arr(Iend+1,j,k,nstp), &
                           u_lstm_corr)
          !$OMP END CRITICAL (lstm_east_u)

          u_final = (1.0 - alpha)*u_orl + alpha*(u_orl + u_lstm_corr)
!--------------------------------------------------------------------
#ifdef M3_FRC_BRY
          u_final = (1.0-tau)*u_final + tau*ubry_east(j,k)
#elif defined M3CLIMATOLOGY
          u_final = (1.0-tau)*u_final + tau*uclm_arr(Iend+1,j,k)
#endif
!--------------------------------------------------------------------
          u_arr(Iend+1,j,k,nnew) = u_final * umask_arr(Iend+1,j)

        enddo
      enddo

    endif

    !==================================================================
    ! SOUTHERN BOUNDARY
    !==================================================================
    if (.true.) then

      do k = 1, N_lev
        do i = Istr-1, Iend
          grad(i,Jstr-1) = u_arr(i+1,Jstr-1,k,nstp) - u_arr(i,Jstr-1,k,nstp)
          grad(i,Jstr  ) = u_arr(i+1,Jstr  ,k,nstp) - u_arr(i,Jstr  ,k,nstp)
        enddo

        do i = Istr, Iend
          j_idx = i - Istr + 1

          dft = u_arr(i,Jstr,k,nstp) - u_arr(i,Jstr  ,k,nnew)
          dfx = u_arr(i,Jstr,k,nnew) - u_arr(i,Jstr+1,k,nnew)

          tau = tau_out
          if (dfx*dft .lt. 0.0) then ; dft = 0.0 ; tau = tau_in ; endif

          if (dft*(grad(i-1,Jstr)+grad(i,Jstr)) .gt. 0.0) then
            dfy = grad(i-1,Jstr)
          else
            dfy = grad(i  ,Jstr)
          endif

          cff    = max(dfx*dfx + dfy*dfy, eps_in)
          cx_orl = dft*dfx
          cy_orl = min(cff, max(dft*dfy, -cff))

          call psc_correct(dft, dfx, dfy, cff, cx_orl, cy_orl, cx_ml, cy_ml)

          u_orl = (  cff  * u_arr(i,Jstr-1,k,nstp)     &
                   + cx_ml* u_arr(i,Jstr  ,k,nnew)      &
                   - max(cy_ml,0.0)*grad(i-1,Jstr-1)    &
                   - min(cy_ml,0.0)*grad(i  ,Jstr-1)    &
                  ) / (cff + cx_ml + eps_in)

          !$OMP CRITICAL (lstm_south_u)
          call lstm_update(3, j_idx, k, &
                           u_arr(i,Jstr-1,k,nstp), v_arr(i,Jstr-1,k,nstp), &
                           u_lstm_corr)
          !$OMP END CRITICAL (lstm_south_u)

          u_final = (1.0 - alpha)*u_orl + alpha*(u_orl + u_lstm_corr)
!--------------------------------------------------------------------
#ifdef M3_FRC_BRY
          u_final = (1.0-tau)*u_final + tau*ubry_south(i,k)
#elif defined M3CLIMATOLOGY
          u_final = (1.0-tau)*u_final + tau*uclm_arr(i,Jstr-1,k)
#endif
!--------------------------------------------------------------------
          u_arr(i,Jstr-1,k,nnew) = u_final * umask_arr(i,Jstr-1)

        enddo
      enddo

    endif  ! lsouth

    !==================================================================
    ! NORTHERN BOUNDARY
    !==================================================================
    if (.true.) then

      do k = 1, N_lev
        do i = Istr-1, Iend
          grad(i,Jend  ) = u_arr(i+1,Jend  ,k,nstp) - u_arr(i,Jend  ,k,nstp)
          grad(i,Jend+1) = u_arr(i+1,Jend+1,k,nstp) - u_arr(i,Jend+1,k,nstp)
        enddo

        do i = Istr, Iend
          j_idx = i - Istr + 1

          dft = u_arr(i,Jend,k,nstp) - u_arr(i,Jend  ,k,nnew)
          dfx = u_arr(i,Jend,k,nnew) - u_arr(i,Jend-1,k,nnew)

          tau = tau_out
          if (dfx*dft .lt. 0.0) then ; dft = 0.0 ; tau = tau_in ; endif

          if (dft*(grad(i-1,Jend)+grad(i,Jend)) .gt. 0.0) then
            dfy = grad(i-1,Jend)
          else
            dfy = grad(i  ,Jend)
          endif

          cff    = max(dfx*dfx + dfy*dfy, eps_in)
          cx_orl = dft*dfx
          cy_orl = min(cff, max(dft*dfy, -cff))

          call psc_correct(dft, dfx, dfy, cff, cx_orl, cy_orl, cx_ml, cy_ml)

          u_orl = (  cff  * u_arr(i,Jend+1,k,nstp)     &
                   + cx_ml* u_arr(i,Jend  ,k,nnew)      &
                   - max(cy_ml,0.0)*grad(i-1,Jend+1)    &
                   - min(cy_ml,0.0)*grad(i  ,Jend+1)    &
                  ) / (cff + cx_ml + eps_in)

          !$OMP CRITICAL (lstm_north_u)
          call lstm_update(4, j_idx, k, &
                           u_arr(i,Jend+1,k,nstp), v_arr(i,Jend+1,k,nstp), &
                           u_lstm_corr)
          !$OMP END CRITICAL (lstm_north_u)

          u_final = (1.0 - alpha)*u_orl + alpha*(u_orl + u_lstm_corr)
!--------------------------------------------------------------------
#ifdef M3_FRC_BRY
          u_final = (1.0-tau)*u_final + tau*ubry_north(i,k)
#elif defined M3CLIMATOLOGY
          u_final = (1.0-tau)*u_final + tau*uclm_arr(i,Jend+1,k)
#endif
!--------------------------------------------------------------------
          u_arr(i,Jend+1,k,nnew) = u_final * umask_arr(i,Jend+1)

        enddo
      enddo

    endif  ! lnorth

  end subroutine obc_orlanski_lstm_psc_m3_u

  !=======================================================================
  !  obc_orlanski_lstm_psc_m3_v
  !  PSC + LSTM to ETA-velocity v at all 4 open boundaries.
  !  Mirrors obc_orlanski_lstm_psc_m3_u with **transposed** stencil dir.
  !=======================================================================

  subroutine obc_orlanski_lstm_psc_m3_v(Istr, Iend, Jstr, Jend,    &
                             nstp, nnew, N_lev,                    &
                             eps_in, tau_in, tau_out,              &
                             u_arr, v_arr,                         &
                             grad,                                 & 
                             vmask_arr, pmask_arr,                 &
                             vclm_arr)
    implicit none
    integer, intent(in) :: Istr, Iend, Jstr, Jend
    integer, intent(in) :: nstp, nnew, N_lev
    real,    intent(in) :: eps_in, tau_in, tau_out
    real, intent(in)    :: u_arr(0:,0:,:,:)
    real, intent(inout) :: v_arr(0:,0:,:,:)
    real, intent(inout) :: grad(0:,0:)
    real, intent(in)    :: vmask_arr(0:,0:)
    real, intent(in)    :: pmask_arr(0:,0:)
    real, intent(in)    :: vclm_arr(0:,0:,:)
!--------------------------------------------------------------------
#ifdef M3_FRC_BRY
# include "boundary.h"
#endif
!--------------------------------------------------------------------
    integer :: i, j, k, j_idx
    real    :: dft, dfx, dfy, cff
    real    :: cx_orl, cy_orl, cx_ml, cy_ml
    real    :: tau, v_orl, u_lstm_corr, alpha, v_final

    if (lstm_step_count >= LSTM_WARMUP_STEPS) then
      alpha = LSTM_BLEND_MAX
    else
      alpha = LSTM_BLEND_MIN + &
              (LSTM_BLEND_MAX - LSTM_BLEND_MIN) * &
              real(lstm_step_count) / real(max(1, LSTM_WARMUP_STEPS))
    endif

    !==================================================================
    ! WESTERN boundary
    !==================================================================
    if (.true.) then

      do k = 1, N_lev
        do j = Jstr, Jend+1
          grad(Istr-1,j) = (v_arr(Istr-1,j,k,nstp) - v_arr(Istr-1,j-1,k,nstp)) &
                           * pmask_arr(Istr-1,j)
          grad(Istr  ,j) = (v_arr(Istr  ,j,k,nstp) - v_arr(Istr  ,j-1,k,nstp)) &
                           * pmask_arr(Istr  ,j)
        enddo

        do j = Jstr, Jend
          j_idx = j - Jstr + 1
          dft = v_arr(Istr  ,j,k,nstp) - v_arr(Istr  ,j,k,nnew)
          dfx = v_arr(Istr  ,j,k,nnew) - v_arr(Istr+1,j,k,nnew)

          tau = tau_out
          if (dfx*dft .lt. 0.0) then ; dft = 0.0 ; tau = tau_in ; endif

          if (dft*(grad(Istr,j)+grad(Istr,j+1)) .gt. 0.0) then
            dfy = grad(Istr,j)
          else
            dfy = grad(Istr,j+1)
          endif

          cff    = max(dfx*dfx+dfy*dfy, eps_in)
          cx_orl = dft*dfx
          cy_orl = min(cff, max(dft*dfy,-cff))

          call psc_correct(dft,dfx,dfy,cff,cx_orl,cy_orl,cx_ml,cy_ml)

          v_orl = ( cff  *v_arr(Istr-1,j,k,nstp)       &
                   +cx_ml*v_arr(Istr  ,j,k,nnew)       &
                   -max(cy_ml,0.0)*grad(Istr-1,j  )    &
                   -min(cy_ml,0.0)*grad(Istr-1,j+1)    &
                  ) / (cff+cx_ml+eps_in)

          !$OMP CRITICAL (lstm_west_v)
          call lstm_update(1, j_idx, k, &
                           u_arr(Istr-1,j,k,nstp), v_arr(Istr-1,j,k,nstp), &
                           u_lstm_corr)
          !$OMP END CRITICAL (lstm_west_v)

          v_final = (1.0-alpha)*v_orl + alpha*(v_orl+u_lstm_corr)
!--------------------------------------------------------------------
#ifdef M3_FRC_BRY
          v_final = (1.0-tau)*v_final + tau*vbry_west(j,k)
#elif defined M3CLIMATOLOGY
          v_final = (1.0-tau)*v_final + tau*vclm_arr(Istr-1,j,k)
#endif
!--------------------------------------------------------------------
          v_arr(Istr-1,j,k,nnew) = v_final * vmask_arr(Istr-1,j)
        enddo
      enddo

    endif

    !==================================================================
    ! EASTERN boundary
    !==================================================================
    if (.true.) then

      do k = 1, N_lev
        do j = Jstr, Jend+1
          grad(Iend  ,j) = (v_arr(Iend  ,j,k,nstp)-v_arr(Iend  ,j-1,k,nstp)) &
                           * pmask_arr(Iend,j)
          grad(Iend+1,j) = (v_arr(Iend+1,j,k,nstp)-v_arr(Iend+1,j-1,k,nstp)) &
                           * pmask_arr(Iend+1,j)
        enddo

        do j = Jstr, Jend
          j_idx = j - Jstr + 1
          dft = v_arr(Iend  ,j,k,nstp) - v_arr(Iend  ,j,k,nnew)
          dfx = v_arr(Iend  ,j,k,nnew) - v_arr(Iend-1,j,k,nnew)

          tau = tau_out
          if (dfx*dft .lt. 0.0) then ; dft = 0.0 ; tau = tau_in ; endif

          if (dft*(grad(Iend,j)+grad(Iend,j+1)) .gt. 0.0) then
            dfy = grad(Iend,j)
          else
            dfy = grad(Iend,j+1)
          endif

          cff    = max(dfx*dfx+dfy*dfy, eps_in)
          cx_orl = dft*dfx
          cy_orl = min(cff, max(dft*dfy,-cff))

          call psc_correct(dft,dfx,dfy,cff,cx_orl,cy_orl,cx_ml,cy_ml)

          v_orl = ( cff  *v_arr(Iend+1,j,k,nstp)      &
                   +cx_ml*v_arr(Iend  ,j,k,nnew)       &
                   -max(cy_ml,0.0)*grad(Iend+1,j  )    &
                   -min(cy_ml,0.0)*grad(Iend+1,j+1)    &
                  ) / (cff+cx_ml+eps_in)

          !$OMP CRITICAL (lstm_east_v)
          call lstm_update(2, j_idx, k, &
                           u_arr(Iend+1,j,k,nstp), v_arr(Iend+1,j,k,nstp), &
                           u_lstm_corr)
          !$OMP END CRITICAL (lstm_east_v)

          v_final = (1.0-alpha)*v_orl + alpha*(v_orl+u_lstm_corr)
!--------------------------------------------------------------------
#ifdef M3_FRC_BRY
          v_final = (1.0-tau)*v_final + tau*vbry_east(j,k)
#elif defined M3CLIMATOLOGY
          v_final = (1.0-tau)*v_final + tau*vclm_arr(Iend+1,j,k)
#endif
!--------------------------------------------------------------------
          v_arr(Iend+1,j,k,nnew) = v_final * vmask_arr(Iend+1,j)
        enddo
      enddo

    endif

    !==================================================================
    ! SOUTHERN boundary
    !==================================================================
    if (.true.) then

      do k = 1, N_lev
        do i = Istr, Iend+1
          grad(i,Jstr-1) = (v_arr(i,Jstr-1,k,nstp)-v_arr(i-1,Jstr-1,k,nstp)) &
                           * pmask_arr(i,Jstr-1)
          grad(i,Jstr  ) = (v_arr(i,Jstr  ,k,nstp)-v_arr(i-1,Jstr  ,k,nstp)) &
                           * pmask_arr(i,Jstr  )
        enddo

        do i = Istr, Iend
          j_idx = i - Istr + 1
          dft = v_arr(i,Jstr  ,k,nstp) - v_arr(i,Jstr  ,k,nnew)
          dfx = v_arr(i,Jstr  ,k,nnew) - v_arr(i,Jstr+1,k,nnew)

          tau = tau_out
          if (dfx*dft .lt. 0.0) then ; dft = 0.0 ; tau = tau_in ; endif

          if (dft*(grad(i,Jstr)+grad(i+1,Jstr)) .gt. 0.0) then
            dfy = grad(i  ,Jstr)
          else
            dfy = grad(i+1,Jstr)
          endif

          cff    = max(dfx*dfx+dfy*dfy, eps_in)
          cx_orl = dft*dfx
          cy_orl = min(cff, max(dft*dfy,-cff))

          call psc_correct(dft,dfx,dfy,cff,cx_orl,cy_orl,cx_ml,cy_ml)

          v_orl = ( cff  *v_arr(i,Jstr-1,k,nstp)       &
                   +cx_ml*v_arr(i,Jstr  ,k,nnew)       &
                   -max(cy_ml,0.0)*grad(i  ,Jstr-1)    &
                   -min(cy_ml,0.0)*grad(i+1,Jstr-1)    &
                  ) / (cff+cx_ml+eps_in)

          !$OMP CRITICAL (lstm_south_v)
          call lstm_update(3, j_idx, k, &
                           u_arr(i,Jstr-1,k,nstp), v_arr(i,Jstr-1,k,nstp), &
                           u_lstm_corr)
          !$OMP END CRITICAL (lstm_south_v)

          v_final = (1.0-alpha)*v_orl + alpha*(v_orl+u_lstm_corr)
!--------------------------------------------------------------------
#ifdef M3_FRC_BRY
          v_final = (1.0-tau)*v_final + tau*vbry_south(i,k)
#elif defined M3CLIMATOLOGY
          v_final = (1.0-tau)*v_final + tau*vclm_arr(i,Jstr-1,k)
#endif
!--------------------------------------------------------------------
          v_arr(i,Jstr-1,k,nnew) = v_final * vmask_arr(i,Jstr-1)
        enddo
      enddo

    endif

    !==================================================================
    ! NORTHERN boundary
    !==================================================================
    if (.true.) then

      do k = 1, N_lev
        do i = Istr, Iend+1
          grad(i,Jend  ) = (v_arr(i,Jend  ,k,nstp)-v_arr(i-1,Jend  ,k,nstp)) &
                           * pmask_arr(i,Jend  )
          grad(i,Jend+1) = (v_arr(i,Jend+1,k,nstp)-v_arr(i-1,Jend+1,k,nstp)) &
                           * pmask_arr(i,Jend+1)
        enddo

        do i = Istr, Iend
          j_idx = i - Istr + 1
          dft = v_arr(i,Jend  ,k,nstp) - v_arr(i,Jend  ,k,nnew)
          dfx = v_arr(i,Jend  ,k,nnew) - v_arr(i,Jend-1,k,nnew)

          tau = tau_out
          if (dfx*dft .lt. 0.0) then ; dft = 0.0 ; tau = tau_in ; endif

          if (dft*(grad(i,Jend)+grad(i+1,Jend)) .gt. 0.0) then
            dfy = grad(i  ,Jend)
          else
            dfy = grad(i+1,Jend)
          endif

          cff    = max(dfx*dfx+dfy*dfy, eps_in)
          cx_orl = dft*dfx
          cy_orl = min(cff, max(dft*dfy,-cff))

          call psc_correct(dft,dfx,dfy,cff,cx_orl,cy_orl,cx_ml,cy_ml)

          v_orl = ( cff  *v_arr(i,Jend+1,k,nstp)       &
                   +cx_ml*v_arr(i,Jend  ,k,nnew)       &
                   -max(cy_ml,0.0)*grad(i  ,Jend+1)    &
                   -min(cy_ml,0.0)*grad(i+1,Jend+1)    &
                  ) / (cff+cx_ml+eps_in)

          !$OMP CRITICAL (lstm_north_v)
          call lstm_update(4, j_idx, k, &
                           u_arr(i,Jend+1,k,nstp), v_arr(i,Jend+1,k,nstp), &
                           u_lstm_corr)
          !$OMP END CRITICAL (lstm_north_v)

          v_final = (1.0-alpha)*v_orl + alpha*(v_orl+u_lstm_corr)
!--------------------------------------------------------------------
#ifdef M3_FRC_BRY
          v_final = (1.0-tau)*v_final + tau*vbry_north(i,k)
#elif defined M3CLIMATOLOGY
          v_final = (1.0-tau)*v_final + tau*vclm_arr(i,Jend+1,k)
#endif
!--------------------------------------------------------------------
          v_arr(i,Jend+1,k,nnew) = v_final * vmask_arr(i,Jend+1)
        enddo
      enddo

    endif

  end subroutine obc_orlanski_lstm_psc_m3_v

  !=======================================================================
  !  ONLINE TRAINING 
  !=======================================================================

#ifdef OBC_LSTM_ONLINE_TRAIN

  !--------------------------------------------------------------------
  ! Training hyper-parameters
  !--------------------------------------------------------------------
  real,    parameter :: TRAIN_LR        = 1.0e-3   ! Adam learning rate
  real,    parameter :: TRAIN_BETA1     = 0.9      ! Adam 1st moment decay
  real,    parameter :: TRAIN_BETA2     = 0.999    ! Adam 2nd moment decay
  real,    parameter :: TRAIN_EPS_ADAM  = 1.0e-8   ! Adam denominator floor
  real,    parameter :: TRAIN_GRAD_CLIP = 1.0      ! gradient clipping norm
  integer, parameter :: TRAIN_SEQ_LEN   = 5         ! circular buffer depth
  integer, parameter :: TRAIN_FREQ      = 1         ! update every N steps

  !--------------------------------------------------------------------
  ! Adam step 
  !--------------------------------------------------------------------
  integer :: adam_t = 0

  !--------------------------------------------------------------------
  ! Circular observation buffer  (wall, bdry_pt, level, time, uv)
  ! Kept small: only TRAIN_SEQ_LEN recent steps are needed.
  !--------------------------------------------------------------------
  real    :: train_buf(4, MAX_BDRY_PTS, MAX_N, TRAIN_SEQ_LEN, LSTM_IN)
  integer :: train_buf_ptr = 0           ! current write pointer (1-based)
  integer :: train_buf_count = 0         ! how many steps have been buffered

  !--------------------------------------------------------------------
  ! PSC Adam moment buffers
  !--------------------------------------------------------------------
  real :: adam_m_W1(ML_HIDDEN, ML_IN),   adam_v_W1(ML_HIDDEN, ML_IN)
  real :: adam_m_b1(ML_HIDDEN),          adam_v_b1(ML_HIDDEN)
  real :: adam_m_W2(2, ML_HIDDEN),       adam_v_W2(2, ML_HIDDEN)
  real :: adam_m_b2(2),                  adam_v_b2(2)

  !--------------------------------------------------------------------
  ! LSTM Adam moment buffers  (input gate)
  !--------------------------------------------------------------------
  real :: adam_m_WI_X(LSTM_HS, LSTM_IN), adam_v_WI_X(LSTM_HS, LSTM_IN)
  real :: adam_m_WI_H(LSTM_HS, LSTM_HS), adam_v_WI_H(LSTM_HS, LSTM_HS)
  real :: adam_m_BI(LSTM_HS),            adam_v_BI(LSTM_HS)

  !--------------------------------------------------------------------
  ! LSTM Adam moment buffers  (forget gate)
  !--------------------------------------------------------------------
  real :: adam_m_WF_X(LSTM_HS, LSTM_IN), adam_v_WF_X(LSTM_HS, LSTM_IN)
  real :: adam_m_WF_H(LSTM_HS, LSTM_HS), adam_v_WF_H(LSTM_HS, LSTM_HS)
  real :: adam_m_BF(LSTM_HS),            adam_v_BF(LSTM_HS)

  !--------------------------------------------------------------------
  ! LSTM Adam moment buffers  (cell gate)
  !--------------------------------------------------------------------
  real :: adam_m_WG_X(LSTM_HS, LSTM_IN), adam_v_WG_X(LSTM_HS, LSTM_IN)
  real :: adam_m_WG_H(LSTM_HS, LSTM_HS), adam_v_WG_H(LSTM_HS, LSTM_HS)
  real :: adam_m_BG(LSTM_HS),            adam_v_BG(LSTM_HS)

  !--------------------------------------------------------------------
  ! LSTM Adam moment buffers  (output gate)
  !--------------------------------------------------------------------
  real :: adam_m_WO_X(LSTM_HS, LSTM_IN), adam_v_WO_X(LSTM_HS, LSTM_IN)
  real :: adam_m_WO_H(LSTM_HS, LSTM_HS), adam_v_WO_H(LSTM_HS, LSTM_HS)
  real :: adam_m_BO(LSTM_HS),            adam_v_BO(LSTM_HS)

  !--------------------------------------------------------------------
  ! LSTM output projection Adam moment buffers
  !--------------------------------------------------------------------
  real :: adam_m_WOUT(LSTM_HS), adam_v_WOUT(LSTM_HS)
  real :: adam_m_BOUT,          adam_v_BOUT

  !=======================================================================
  !  lstm_train_init
  !=======================================================================

  subroutine lstm_train_init()

    adam_t          = 0
    train_buf_ptr   = 0
    train_buf_count = 0
    train_buf       = 0.0

    ! PSC Adam moments
    adam_m_W1 = 0.0 ; adam_v_W1 = 0.0
    adam_m_b1 = 0.0 ; adam_v_b1 = 0.0
    adam_m_W2 = 0.0 ; adam_v_W2 = 0.0
    adam_m_b2 = 0.0 ; adam_v_b2 = 0.0

    ! LSTM Adam moments – input gate
    adam_m_WI_X = 0.0 ; adam_v_WI_X = 0.0
    adam_m_WI_H = 0.0 ; adam_v_WI_H = 0.0
    adam_m_BI   = 0.0 ; adam_v_BI   = 0.0

    ! LSTM Adam moments – forget gate
    adam_m_WF_X = 0.0 ; adam_v_WF_X = 0.0
    adam_m_WF_H = 0.0 ; adam_v_WF_H = 0.0
    adam_m_BF   = 0.0 ; adam_v_BF   = 0.0

    ! LSTM Adam moments – cell gate
    adam_m_WG_X = 0.0 ; adam_v_WG_X = 0.0
    adam_m_WG_H = 0.0 ; adam_v_WG_H = 0.0
    adam_m_BG   = 0.0 ; adam_v_BG   = 0.0

    ! LSTM Adam moments – output gate
    adam_m_WO_X = 0.0 ; adam_v_WO_X = 0.0
    adam_m_WO_H = 0.0 ; adam_v_WO_H = 0.0
    adam_m_BO   = 0.0 ; adam_v_BO   = 0.0

    ! Output projection
    adam_m_WOUT = 0.0 ; adam_v_WOUT = 0.0
    adam_m_BOUT = 0.0 ; adam_v_BOUT = 0.0

    write(6,'(A)') ' LSTM TRAIN: online training initialised (Adam optimiser).'

  end subroutine lstm_train_init

  !=======================================================================
  !  adam_update_scalar
  !=======================================================================

  subroutine adam_update_scalar(param, grad_val, m, v, t_adam)

    real,    intent(inout) :: param, m, v
    real,    intent(in)    :: grad_val
    integer, intent(in)    :: t_adam

    real :: m_hat, v_hat, g

    g = grad_val
    ! Gradient clipping (element-wise)
    g = max(-TRAIN_GRAD_CLIP, min(TRAIN_GRAD_CLIP, g))

    m     = TRAIN_BETA1 * m + (1.0 - TRAIN_BETA1) * g
    v     = TRAIN_BETA2 * v + (1.0 - TRAIN_BETA2) * g * g
    m_hat = m / (1.0 - TRAIN_BETA1**t_adam)
    v_hat = v / (1.0 - TRAIN_BETA2**t_adam)

    param = param - TRAIN_LR * m_hat / (sqrt(v_hat) + TRAIN_EPS_ADAM)

  end subroutine adam_update_scalar

  !=======================================================================
  !  psc_forward_and_grad
  !=======================================================================

  subroutine psc_forward_and_grad(x_in,            &
                                  dL_dW1, dL_db1,  &
                                  dL_dW2, dL_db2)

    real, intent(in)  :: x_in(ML_IN)
    real, intent(out) :: dL_dW1(ML_HIDDEN, ML_IN)
    real, intent(out) :: dL_db1(ML_HIDDEN)
    real, intent(out) :: dL_dW2(2, ML_HIDDEN)
    real, intent(out) :: dL_db2(2)

    real :: h1(ML_HIDDEN), pre1(ML_HIDDEN)
    real :: out(2),         pre2(2)
    real :: dL_dpre2(2),   dL_dpre1(ML_HIDDEN)
    integer :: i, j

    ! --- Forward pass ---------------------------------------------------
    ! Hidden layer: pre1 = W1*x + b1,  h1 = tanh(pre1)
    do i = 1, ML_HIDDEN
      pre1(i) = b1(i)
      do j = 1, ML_IN
        pre1(i) = pre1(i) + W1(i,j) * x_in(j)
      enddo
      h1(i) = tanh_act(pre1(i))
    enddo

    ! Output layer --------------------------------------
    do i = 1, 2
      pre2(i) = b2(i)
      do j = 1, ML_HIDDEN
        pre2(i) = pre2(i) + W2(i,j) * h1(j)
      enddo
      out(i) = tanh_act(pre2(i))
    enddo

    ! --- Backward pass -------------
    do i = 1, 2
      dL_dpre2(i) = out(i) * (1.0 - out(i)*out(i))
    enddo

    do i = 1, 2
      dL_db2(i) = dL_dpre2(i)
      do j = 1, ML_HIDDEN
        dL_dW2(i,j) = dL_dpre2(i) * h1(j)
      enddo
    enddo

    do j = 1, ML_HIDDEN
      dL_dpre1(j) = 0.0
      do i = 1, 2
        dL_dpre1(j) = dL_dpre1(j) + dL_dpre2(i) * W2(i,j)
      enddo
      dL_dpre1(j) = dL_dpre1(j) * (1.0 - h1(j)*h1(j))
    enddo

    do i = 1, ML_HIDDEN
      dL_db1(i) = dL_dpre1(i)
      do j = 1, ML_IN
        dL_dW1(i,j) = dL_dpre1(i) * x_in(j)
      enddo
    enddo

  end subroutine psc_forward_and_grad

  !=======================================================================
  !  lstm_forward_and_grad
  !=======================================================================

  subroutine lstm_forward_and_grad(x_in, h_prev, c_prev, y_target, &
                                   dL_dWI_X, dL_dWI_H, dL_dBI,    &
                                   dL_dWF_X, dL_dWF_H, dL_dBF,    &
                                   dL_dWG_X, dL_dWG_H, dL_dBG,    &
                                   dL_dWO_X, dL_dWO_H, dL_dBO,    &
                                   dL_dWOUT, dL_dBOUT)

    real, intent(in)  :: x_in(LSTM_IN)
    real, intent(in)  :: h_prev(LSTM_HS), c_prev(LSTM_HS)
    real, intent(in)  :: y_target

    real, intent(out) :: dL_dWI_X(LSTM_HS, LSTM_IN)
    real, intent(out) :: dL_dWI_H(LSTM_HS, LSTM_HS)
    real, intent(out) :: dL_dBI(LSTM_HS)
    real, intent(out) :: dL_dWF_X(LSTM_HS, LSTM_IN)
    real, intent(out) :: dL_dWF_H(LSTM_HS, LSTM_HS)
    real, intent(out) :: dL_dBF(LSTM_HS)
    real, intent(out) :: dL_dWG_X(LSTM_HS, LSTM_IN)
    real, intent(out) :: dL_dWG_H(LSTM_HS, LSTM_HS)
    real, intent(out) :: dL_dBG(LSTM_HS)
    real, intent(out) :: dL_dWO_X(LSTM_HS, LSTM_IN)
    real, intent(out) :: dL_dWO_H(LSTM_HS, LSTM_HS)
    real, intent(out) :: dL_dBO(LSTM_HS)
    real, intent(out) :: dL_dWOUT(LSTM_HS)
    real, intent(out) :: dL_dBOUT

    ! --- Cached gate pre-activations and activations --------------------
    real :: gi(LSTM_HS), gf(LSTM_HS), gg(LSTM_HS), go(LSTM_HS)
    real :: ai(LSTM_HS), af(LSTM_HS), ag(LSTM_HS), ao(LSTM_HS)
    real :: c_new(LSTM_HS), h_new(LSTM_HS), tanh_c(LSTM_HS)

    ! --- Backward intermediates -----------------------------------------
    real :: y_pred, dL_dy, dL_dh(LSTM_HS)
    real :: dL_dao(LSTM_HS), dL_dc(LSTM_HS)
    real :: dL_dai(LSTM_HS), dL_daf(LSTM_HS), dL_dag(LSTM_HS)
    real :: dL_dgo(LSTM_HS), dL_dgi(LSTM_HS), dL_dgf(LSTM_HS), dL_dgg(LSTM_HS)
    integer :: ii, p

    ! ---------------------------------------------------------------
    ! Forward pass
    ! ---------------------------------------------------------------
    do ii = 1, LSTM_HS
      gi(ii) = LSTM_BI(ii)
      gf(ii) = LSTM_BF(ii)
      gg(ii) = LSTM_BG(ii)
      go(ii) = LSTM_BO(ii)
      do p = 1, LSTM_IN
        gi(ii) = gi(ii) + LSTM_WI_X(ii,p) * x_in(p)
        gf(ii) = gf(ii) + LSTM_WF_X(ii,p) * x_in(p)
        gg(ii) = gg(ii) + LSTM_WG_X(ii,p) * x_in(p)
        go(ii) = go(ii) + LSTM_WO_X(ii,p) * x_in(p)
      enddo
      do p = 1, LSTM_HS
        gi(ii) = gi(ii) + LSTM_WI_H(ii,p) * h_prev(p)
        gf(ii) = gf(ii) + LSTM_WF_H(ii,p) * h_prev(p)
        gg(ii) = gg(ii) + LSTM_WG_H(ii,p) * h_prev(p)
        go(ii) = go(ii) + LSTM_WO_H(ii,p) * h_prev(p)
      enddo
      ai(ii) = sigmoid(gi(ii))
      af(ii) = sigmoid(gf(ii))
      ag(ii) = tanh_act(gg(ii))
      ao(ii) = sigmoid(go(ii))
    enddo

    do ii = 1, LSTM_HS
      c_new(ii)  = af(ii)*c_prev(ii) + ai(ii)*ag(ii)
      tanh_c(ii) = tanh_act(c_new(ii))
      h_new(ii)  = ao(ii) * tanh_c(ii)
    enddo

    ! Output
    y_pred = LSTM_BOUT
    do ii = 1, LSTM_HS
      y_pred = y_pred + LSTM_WOUT(ii) * h_new(ii)
    enddo
    y_pred = max(-0.5, min(0.5, y_pred))

    ! ---------------------------------------------------------------
    ! Backward pass
    ! ---------------------------------------------------------------
    dL_dy = y_pred - y_target
    dL_dBOUT = dL_dy
    do ii = 1, LSTM_HS
      dL_dWOUT(ii) = dL_dy * h_new(ii)
    enddo

    do ii = 1, LSTM_HS
      dL_dh(ii) = dL_dy * LSTM_WOUT(ii)
    enddo

    do ii = 1, LSTM_HS
      dL_dao(ii) = dL_dh(ii) * tanh_c(ii)
      dL_dgo(ii) = dL_dao(ii) * ao(ii) * (1.0 - ao(ii))
    enddo

    do ii = 1, LSTM_HS
      dL_dc(ii) = dL_dh(ii) * ao(ii) * (1.0 - tanh_c(ii)*tanh_c(ii))
    enddo

    do ii = 1, LSTM_HS
      dL_dai(ii) = dL_dc(ii) * ag(ii)
      dL_dgi(ii) = dL_dai(ii) * ai(ii) * (1.0 - ai(ii))

      dL_daf(ii) = dL_dc(ii) * c_prev(ii)
      dL_dgf(ii) = dL_daf(ii) * af(ii) * (1.0 - af(ii))

      dL_dag(ii) = dL_dc(ii) * ai(ii)
      dL_dgg(ii) = dL_dag(ii) * (1.0 - ag(ii)*ag(ii))
    enddo

    ! --- Accumulate weight gradients ------------------------------------
    do ii = 1, LSTM_HS
      dL_dBI(ii) = dL_dgi(ii)
      dL_dBF(ii) = dL_dgf(ii)
      dL_dBG(ii) = dL_dgg(ii)
      dL_dBO(ii) = dL_dgo(ii)
      do p = 1, LSTM_IN
        dL_dWI_X(ii,p) = dL_dgi(ii) * x_in(p)
        dL_dWF_X(ii,p) = dL_dgf(ii) * x_in(p)
        dL_dWG_X(ii,p) = dL_dgg(ii) * x_in(p)
        dL_dWO_X(ii,p) = dL_dgo(ii) * x_in(p)
      enddo
      do p = 1, LSTM_HS
        dL_dWI_H(ii,p) = dL_dgi(ii) * h_prev(p)
        dL_dWF_H(ii,p) = dL_dgf(ii) * h_prev(p)
        dL_dWG_H(ii,p) = dL_dgg(ii) * h_prev(p)
        dL_dWO_H(ii,p) = dL_dgo(ii) * h_prev(p)
      enddo
    enddo

  end subroutine lstm_forward_and_grad

  !=======================================================================
  !  lstm_train_step
  !=======================================================================

  subroutine lstm_train_step(iwall, j_idx, k_idx, u_obs, v_obs)

    integer, intent(in) :: iwall, j_idx, k_idx
    real,    intent(in) :: u_obs, v_obs

    ! PSC gradient arrays
    real :: dL_dW1(ML_HIDDEN, ML_IN), dL_db1(ML_HIDDEN)
    real :: dL_dW2(2, ML_HIDDEN),     dL_db2(2)

    ! LSTM gradient arrays
    real :: dL_dWI_X(LSTM_HS, LSTM_IN), dL_dWI_H(LSTM_HS, LSTM_HS)
    real :: dL_dBI(LSTM_HS)
    real :: dL_dWF_X(LSTM_HS, LSTM_IN), dL_dWF_H(LSTM_HS, LSTM_HS)
    real :: dL_dBF(LSTM_HS)
    real :: dL_dWG_X(LSTM_HS, LSTM_IN), dL_dWG_H(LSTM_HS, LSTM_HS)
    real :: dL_dBG(LSTM_HS)
    real :: dL_dWO_X(LSTM_HS, LSTM_IN), dL_dWO_H(LSTM_HS, LSTM_HS)
    real :: dL_dBO(LSTM_HS)
    real :: dL_dWOUT(LSTM_HS), dL_dBOUT

    real :: x_psc(ML_IN), x_lstm(LSTM_IN)
    real :: h_prev(LSTM_HS), c_prev(LSTM_HS)
    real :: u_prev, dft_l, dfx_l, dfy_l, cff_l
    real :: y_target
    integer :: ptr_new, ptr_prev, ii, p

    ! ------------------------------------------------------------------+
    !*******************************************************************+
    ! Advance the global circular buffer pointer.                       +
    ! This is a module-level pointer shared by all points.              +
    ! ******************************************************************+
    ! ------------------------------------------------------------------+

    ptr_new  = train_buf_ptr
    ptr_prev = mod(ptr_new - 2 + TRAIN_SEQ_LEN, TRAIN_SEQ_LEN) + 1

    ! Early exit if buffer not yet full--
    if (train_buf_count < TRAIN_SEQ_LEN) return

    ! Early exit if not yet at training frequency--
    if (mod(lstm_step_count, TRAIN_FREQ) .ne. 0) return

    ! ------------------------------------------------------------------
    ! Build PSC feature vector from circular buffer.
    ! ------------------------------------------------------------------
    u_prev = train_buf(iwall, j_idx, k_idx, ptr_prev, 1)
    dft_l  = u_obs - u_prev
    dfx_l  = u_obs - u_prev
    dfy_l  = v_obs - train_buf(iwall, j_idx, k_idx, ptr_prev, 2)

    cff_l  = max(sqrt(dfx_l*dfx_l + dfy_l*dfy_l), EPS_LSTM)
    x_psc(1) = max(-3.0, min(3.0, dft_l / cff_l))
    x_psc(2) = max(-3.0, min(3.0, dfx_l / cff_l))
    x_psc(3) = max(-3.0, min(3.0, dfy_l / cff_l))

    ! ------------------------------------------------------------------
    ! PSC weight update
    ! ------------------------------------------------------------------
    call psc_forward_and_grad(x_psc, dL_dW1, dL_db1, dL_dW2, dL_db2)

    adam_t = adam_t + 1

    do ii = 1, ML_HIDDEN
      do p = 1, ML_IN
        call adam_update_scalar(W1(ii,p), dL_dW1(ii,p), &
                                adam_m_W1(ii,p), adam_v_W1(ii,p), adam_t)
      enddo
      call adam_update_scalar(b1(ii), dL_db1(ii), &
                              adam_m_b1(ii), adam_v_b1(ii), adam_t)
    enddo
    do ii = 1, 2
      do p = 1, ML_HIDDEN
        call adam_update_scalar(W2(ii,p), dL_dW2(ii,p), &
                                adam_m_W2(ii,p), adam_v_W2(ii,p), adam_t)
      enddo
      call adam_update_scalar(b2(ii), dL_db2(ii), &
                              adam_m_b2(ii), adam_v_b2(ii), adam_t)
    enddo

    ! ------------------------------------------------------------------
    ! LSTM self-supervised target and weight update
    ! ------------------------------------------------------------------
    y_target = max(-0.5, min(0.5, dft_l))

    x_lstm(1) = max(-5.0, min(5.0, u_obs))
    x_lstm(2) = max(-5.0, min(5.0, v_obs))

    h_prev(:) = lstm_h(iwall, j_idx, k_idx, :)
    c_prev(:) = lstm_c(iwall, j_idx, k_idx, :)

    call lstm_forward_and_grad(x_lstm, h_prev, c_prev, y_target, &
                               dL_dWI_X, dL_dWI_H, dL_dBI,       &
                               dL_dWF_X, dL_dWF_H, dL_dBF,       &
                               dL_dWG_X, dL_dWG_H, dL_dBG,       &
                               dL_dWO_X, dL_dWO_H, dL_dBO,       &
                               dL_dWOUT, dL_dBOUT)

    ! Adam update – LSTM input gate
    do ii = 1, LSTM_HS
      do p = 1, LSTM_IN
        call adam_update_scalar(LSTM_WI_X(ii,p), dL_dWI_X(ii,p), &
                                adam_m_WI_X(ii,p), adam_v_WI_X(ii,p), adam_t)
      enddo
      do p = 1, LSTM_HS
        call adam_update_scalar(LSTM_WI_H(ii,p), dL_dWI_H(ii,p), &
                                adam_m_WI_H(ii,p), adam_v_WI_H(ii,p), adam_t)
      enddo
      call adam_update_scalar(LSTM_BI(ii), dL_dBI(ii), &
                              adam_m_BI(ii), adam_v_BI(ii), adam_t)
    enddo

    ! Adam update – LSTM forget gate
    do ii = 1, LSTM_HS
      do p = 1, LSTM_IN
        call adam_update_scalar(LSTM_WF_X(ii,p), dL_dWF_X(ii,p), &
                                adam_m_WF_X(ii,p), adam_v_WF_X(ii,p), adam_t)
      enddo
      do p = 1, LSTM_HS
        call adam_update_scalar(LSTM_WF_H(ii,p), dL_dWF_H(ii,p), &
                                adam_m_WF_H(ii,p), adam_v_WF_H(ii,p), adam_t)
      enddo
      call adam_update_scalar(LSTM_BF(ii), dL_dBF(ii), &
                              adam_m_BF(ii), adam_v_BF(ii), adam_t)
    enddo

    ! Adam update – LSTM cell gate
    do ii = 1, LSTM_HS
      do p = 1, LSTM_IN
        call adam_update_scalar(LSTM_WG_X(ii,p), dL_dWG_X(ii,p), &
                                adam_m_WG_X(ii,p), adam_v_WG_X(ii,p), adam_t)
      enddo
      do p = 1, LSTM_HS
        call adam_update_scalar(LSTM_WG_H(ii,p), dL_dWG_H(ii,p), &
                                adam_m_WG_H(ii,p), adam_v_WG_H(ii,p), adam_t)
      enddo
      call adam_update_scalar(LSTM_BG(ii), dL_dBG(ii), &
                              adam_m_BG(ii), adam_v_BG(ii), adam_t)
    enddo

    ! Adam update – LSTM output gate
    do ii = 1, LSTM_HS
      do p = 1, LSTM_IN
        call adam_update_scalar(LSTM_WO_X(ii,p), dL_dWO_X(ii,p), &
                                adam_m_WO_X(ii,p), adam_v_WO_X(ii,p), adam_t)
      enddo
      do p = 1, LSTM_HS
        call adam_update_scalar(LSTM_WO_H(ii,p), dL_dWO_H(ii,p), &
                                adam_m_WO_H(ii,p), adam_v_WO_H(ii,p), adam_t)
      enddo
      call adam_update_scalar(LSTM_BO(ii), dL_dBO(ii), &
                              adam_m_BO(ii), adam_v_BO(ii), adam_t)
    enddo

    ! Adam update – output projection
    do ii = 1, LSTM_HS
      call adam_update_scalar(LSTM_WOUT(ii), dL_dWOUT(ii), &
                              adam_m_WOUT(ii), adam_v_WOUT(ii), adam_t)
    enddo
    call adam_update_scalar(LSTM_BOUT, dL_dBOUT, &
                            adam_m_BOUT, adam_v_BOUT, adam_t)

  end subroutine lstm_train_step

  !=======================================================================
  !  lstm_train_buffer_push
  !=======================================================================

  subroutine lstm_train_buffer_push(iwall, j_idx, k_idx, u_obs, v_obs)

    integer, intent(in) :: iwall, j_idx, k_idx
    real,    intent(in) :: u_obs, v_obs

    if (iwall .eq. 1 .and. j_idx .eq. 1 .and. k_idx .eq. 1) then
      train_buf_ptr   = mod(train_buf_ptr, TRAIN_SEQ_LEN) + 1
      train_buf_count = min(train_buf_count + 1, TRAIN_SEQ_LEN)
    endif

    train_buf(iwall, j_idx, k_idx, train_buf_ptr, 1) = u_obs
    train_buf(iwall, j_idx, k_idx, train_buf_ptr, 2) = v_obs

  end subroutine lstm_train_buffer_push

  !=======================================================================
  !  lstm_train_save_weight
  !=======================================================================

  subroutine lstm_train_save_weights(outfile, mynode)

    character(len=*), intent(in) :: outfile
    integer,          intent(in) :: mynode

    integer :: iunit, ii, p, ios

    if (mynode .ne. 0) return

    iunit = 98
    open(unit=iunit, file=trim(outfile), status='replace', action='write', &
         iostat=ios)
    if (ios .ne. 0) then
      write(6,'(2A,I4)') ' LSTM TRAIN: ERROR opening output file ', &
                           trim(outfile), ios
      return
    endif

    write(iunit,'(A)') '# CROMS LSTM OBC weight file'
    write(iunit,'(A)') '# Generated by obc_orlanski_lstm_psc_m3 online training'
    write(iunit,'(A,I0,A,I0,A,I0,A,I0)') &
          '# LSTM_HS=', LSTM_HS, '  ML_HIDDEN=', ML_HIDDEN, &
          '  ML_IN=', ML_IN, '  LSTM_IN=', LSTM_IN
    write(iunit,'(A)') '#'

    ! --- PSC weights ---------------------------------------------------
    write(iunit,'(A)') '# PSC feed-forward network'
    write(iunit,'(A)') 'PSC_W1'
    write(iunit,'(*(ES16.8E2,1X))') ((W1(ii,p), p=1,ML_IN), ii=1,ML_HIDDEN)
    write(iunit,'(A)') 'PSC_b1'
    write(iunit,'(*(ES16.8E2,1X))') (b1(ii), ii=1,ML_HIDDEN)
    write(iunit,'(A)') 'PSC_W2'
    write(iunit,'(*(ES16.8E2,1X))') ((W2(ii,p), p=1,ML_HIDDEN), ii=1,2)
    write(iunit,'(A)') 'PSC_b2'
    write(iunit,'(*(ES16.8E2,1X))') (b2(ii), ii=1,2)
    write(iunit,'(A)') '#'

    ! --- LSTM weights --------------------------------------------------
    write(iunit,'(A)') '# LSTM weights'
    write(iunit,'(A)') 'LSTM_WI_X'
    write(iunit,'(*(ES16.8E2,1X))') &
          ((LSTM_WI_X(ii,p), p=1,LSTM_IN), ii=1,LSTM_HS)
    write(iunit,'(A)') 'LSTM_WI_H'
    write(iunit,'(*(ES16.8E2,1X))') &
          ((LSTM_WI_H(ii,p), p=1,LSTM_HS), ii=1,LSTM_HS)
    write(iunit,'(A)') 'LSTM_BI'
    write(iunit,'(*(ES16.8E2,1X))') (LSTM_BI(ii), ii=1,LSTM_HS)

    write(iunit,'(A)') 'LSTM_WF_X'
    write(iunit,'(*(ES16.8E2,1X))') &
          ((LSTM_WF_X(ii,p), p=1,LSTM_IN), ii=1,LSTM_HS)
    write(iunit,'(A)') 'LSTM_WF_H'
    write(iunit,'(*(ES16.8E2,1X))') &
          ((LSTM_WF_H(ii,p), p=1,LSTM_HS), ii=1,LSTM_HS)
    write(iunit,'(A)') 'LSTM_BF'
    write(iunit,'(*(ES16.8E2,1X))') (LSTM_BF(ii), ii=1,LSTM_HS)

    write(iunit,'(A)') 'LSTM_WG_X'
    write(iunit,'(*(ES16.8E2,1X))') &
          ((LSTM_WG_X(ii,p), p=1,LSTM_IN), ii=1,LSTM_HS)
    write(iunit,'(A)') 'LSTM_WG_H'
    write(iunit,'(*(ES16.8E2,1X))') &
          ((LSTM_WG_H(ii,p), p=1,LSTM_HS), ii=1,LSTM_HS)
    write(iunit,'(A)') 'LSTM_BG'
    write(iunit,'(*(ES16.8E2,1X))') (LSTM_BG(ii), ii=1,LSTM_HS)

    write(iunit,'(A)') 'LSTM_WO_X'
    write(iunit,'(*(ES16.8E2,1X))') &
          ((LSTM_WO_X(ii,p), p=1,LSTM_IN), ii=1,LSTM_HS)
    write(iunit,'(A)') 'LSTM_WO_H'
    write(iunit,'(*(ES16.8E2,1X))') &
          ((LSTM_WO_H(ii,p), p=1,LSTM_HS), ii=1,LSTM_HS)
    write(iunit,'(A)') 'LSTM_BO'
    write(iunit,'(*(ES16.8E2,1X))') (LSTM_BO(ii), ii=1,LSTM_HS)

    write(iunit,'(A)') 'LSTM_WOUT'
    write(iunit,'(*(ES16.8E2,1X))') (LSTM_WOUT(ii), ii=1,LSTM_HS)
    write(iunit,'(A)') 'LSTM_BOUT'
    write(iunit,'(ES16.8E2)') LSTM_BOUT

    close(iunit)
    write(6,'(2A)') ' LSTM TRAIN: weights saved to ', trim(outfile)

  end subroutine lstm_train_save_weights

#endif /* OBC_LSTM_ONLINE_TRAIN */

end module obc_orlanski_lstm_psc_m3
