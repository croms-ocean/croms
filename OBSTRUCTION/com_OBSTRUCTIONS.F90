#include "cppdefs.h"

MODULE com_OBSTRUCTIONS

#ifdef OBSTRUCTION
   !!==========================================================================
   !!                   ***  MODULE com_OBSTRUCTIONS  ***
   !!
   !! ** Purpose : Declaration of variables used by obstructions module
   !!
   !!==========================================================================

   IMPLICIT NONE

   ! default
   PUBLIC

   ! -------------------------------------------------------------------------
   ! Definition of rsh, rlg, riosh, riolg, lchain
   ! -------------------------------------------------------------------------
   INTEGER, PARAMETER :: riosh = 8, riolg = 8, rlg = 8, rsh = 8  ! working precision
   INTEGER, PARAMETER :: lchain = 200
   REAL(KIND=rsh), PARAMETER :: pi = 3.14159265358979323846
   REAL(KIND=rsh), PARAMETER :: obst_c2turb = 1.92 ! obst_c2turb = beta2 in gls kepsilon

   INTEGER :: obst_kmax ! number of layers in water
   INTEGER :: ierrorlog, iwarnlog, iscreenlog ! logging IO

   ! Variables nomenclature :
   ! obst_*    : variables for obstructions
   ! obst_i_*  : variables used for initialization
   ! obst_c_*  : constant parameters
   ! obst_fn_* : filename variables
   ! obst_l_* or l_obst* : logical variables

   !--------------------------------------------------------------------------
   ! * VARIABLES WITHOUT DEPENDANCE ON THE NUMBER OF OBSTRUCTIONS VARIABLES  :
   !--------------------------------------------------------------------------

   ! * Variables for inputs purpose
   CHARACTER(LEN=lchain) :: obst_fn_position                         ! Name of the input file for the obstruction position within the domain

   ! * Variables for outputs purpose
   LOGICAL :: l_obstout_pos                                          ! Write obstruction position (iv,i,j)
   LOGICAL :: l_obstout_height_f                                     ! Write 2D obstruction height (forcing) (iv,i,j)
   LOGICAL :: l_obstout_height_e                                     ! Write 2D obstruction height (effective) (iv,i,j)
   LOGICAL :: l_obstout_dens_f                                       ! Write 2D obstruction density (forcing) (iv,i,j)
   LOGICAL :: l_obstout_dens_e                                       ! Write 3D obstruction density (3D effective) (iv,k,i,j)
   LOGICAL :: l_obstout_dens_2d                                      ! Write 2D obstruction density (2D effective) (i,j)
   LOGICAL :: l_obstout_width_f                                      ! Write 2D obstruction width (forcing) (iv,i,j)
   LOGICAL :: l_obstout_width_e                                      ! Write 3D obstruction width (3D effective) (iv,k,i,j)
   LOGICAL :: l_obstout_width_2d                                     ! Write 2D obstruction width (2D effective) (i,j)
   LOGICAL :: l_obstout_thick_f                                      ! Write 2D obstruction thick (forcing )(iv,i,j)
   LOGICAL :: l_obstout_thick_e                                      ! Write 3D obstruction thick (3D effective) (iv,k,i,j)
   LOGICAL :: l_obstout_theta                                        ! Write 3D obstruction bending angle, (iv,k,i,j)
   LOGICAL :: l_obstout_frac_xy                                      ! Write 2D obstruction coverage, (iv,i,j)
   LOGICAL :: l_obstout_frac_z                                       ! Write 3D obstruction fraction of sigma layer, (iv,k,i,j)
   LOGICAL :: l_obstout_fuzvz                                        ! Write 3D obstruction friction force (k,i,j)
   LOGICAL :: l_obstout_a2d                                          ! Write 2D obstruction horizontal area (iv+3,i,j)
   LOGICAL :: l_obstout_a3d                                          ! Write 3D obstruction horizontal area (iv+3,k,i,j)
   LOGICAL :: l_obstout_s2d                                          ! Write 2D obstruction vertical area (iv+3,i,j)
   LOGICAL :: l_obstout_s3d                                          ! Write 3D obstruction vertical area (iv+3,k,i,j)
   LOGICAL :: l_obstout_drag                                         ! Write 3D obstruction drag coefficient (iv,k,i,j)
   LOGICAL :: l_obstout_tau                                          ! Write 3D obstruction turbulence dissipation scale (k,i,j)

   INTEGER, DIMENSION(:), ALLOCATABLE :: hisObst                     ! Output identifier
   INTEGER, DIMENSION(:), ALLOCATABLE :: avgObst                     ! Output identifier
   LOGICAL, DIMENSION(:), ALLOCATABLE :: outObst                     ! To choose which variable is outputed
   LOGICAL, DIMENSION(:), ALLOCATABLE :: out2DObst                   ! To indicate if the output variable is 2D or 3D
   CHARACTER(LEN=75), DIMENSION(:, :), ALLOCATABLE :: vname_obst     ! Vector of characteristics of each outputed variables

   CHARACTER(LEN=lchain) :: obst_nout_pos                             ! Name obstruction position (iv,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_height_f                        ! Name 2D obstruction height (forcing) (iv,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_height_e                        ! Name 2D obstruction height (effective) (iv,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_dens_f                          ! Name 2D obstruction density (forcing) (iv,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_dens_e                          ! Name 3D obstruction density (3D effective) (iv,k,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_width_f                         ! Name 2D obstruction width (forcing) (iv,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_width_e                         ! Name 3D obstruction width (3D effective) (iv,k,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_thick_f                         ! Name 2D obstruction thick (forcing) (iv,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_thick_e                         ! Name 3D obstruction thick (3D effective) (iv,k,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_theta                           ! Name 3D obstruction bending angle (iv,k,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_frac_xy                         ! Name 2D obstruction coverage, (iv,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_frac_z                          ! Name 3D obstruction fraction of sigma layer, (iv,k,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_fuzvz                           ! Name 3D obstruction friction force (k,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_a2d                             ! Name 2D obstruction horizontal area (iv+3,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_a3d                             ! Name 3D obstruction horizontal area (iv+3,k,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_s2d                             ! Name 2D obstruction vertical area (iv+3,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_s3d                             ! Name 3D obstruction vertical area (iv+3,k,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_drag                            ! Name 3D obstruction drag coefficient (k,i,j)
   CHARACTER(LEN=lchain) :: obst_nout_tau                             ! Name 3D obstruction turbulence dissipation scale (k,i,j)

   ! Other variables/parameters
   INTEGER  :: obst_nbvar                       ! The total number of obstruction variables
   INTEGER  :: obst_nv_up                       ! Number of variable for upward obstructions
   INTEGER  :: obst_nv_do                       ! Number of variable for downward obstruction
   INTEGER  :: obst_nv_3d                       ! Number of variable for full 3d obstructions
   INTEGER  :: obst_nv_turb                     ! Number of variable for obstructions using full turbulence procedure
   INTEGER  :: obst_nv_noturb                   ! Number of variable for obstructions using simplified (roughness length) procedure
   INTEGER  :: obst_nv_rigid_up                 ! Number of variable for rigid and upward (from sea-bed) obstructions
   INTEGER  :: obst_nv_rigid_do                 ! Number of variable for rigid and downward (from sea-surface) obstructions
   INTEGER  :: obst_nv_flexi_up                 ! Number of variable for flexible and upward (from sea-bed) obstructions
   INTEGER  :: obst_nv_flexi_do                 ! Number of variable for flexible and downward (from sea-surface) obstructions

   INTEGER  :: obst_nb_max_hnorm                ! Maximum number of vertical steps from the distribution file

   LOGICAL  :: obst_l_z0bstress_tot             ! IF ONLY ONE obstruction VARIABLE USED Z0SED

   REAL(KIND=rsh)  :: obst_c_paramhuv           ! The coefficient of obstruction height for computation of velocity

   !--------------------------------------------------------------------------
   ! * VARIABLES DEPENDING ONLY ON THE NUMBER OF OBSTRUCTIONS VARIABLES (iv) :
   !--------------------------------------------------------------------------

   INTEGER, DIMENSION(:), ALLOCATABLE :: obst_fracxy_type              ! The type of correction for horizontal coverage (grid cell not completely fill with obstructions) (iv)
   INTEGER, DIMENSION(:), ALLOCATABLE :: obst_nbhnorm                  ! Number of vertical steps from the distribution file (iv)
   INTEGER, DIMENSION(:), ALLOCATABLE :: obst_c_abdel_nmax             ! Number of segments for Abdlerhman method (bending) (iv)

   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_filetimeserie          ! For reading time-series of obstructions characteristics from a file (iv)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_filedistri             ! For reading the vertical distribution of obstructions from a file (iv)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_init_spatial           ! For reading spatial file for density, height, width and thickness (iv)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_flexible               ! For obstructions flexibility (iv)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_cylinder               ! For obstruction shape (cylinder, ellispe / parallelepipeds), (iv)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_noturb                 ! For use of simplified formulation (roughness length) instead of full turbulent formulation (iv)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_abdelrough_cste        ! For use constant Abdelrhman 2003 coefficient (iv)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_fracxy                 ! For horizontal coverage correction (grid cell not completely fill with obstructions) (iv)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_abdelposture           ! For computation of obstructions posture following abdelrhman 2007 (iv)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_param_height           ! For use a parameterization (from velocity) for obstruction height (iv)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_drag_cste              ! For constant/variable drag coefficient (or corrected from bending angle) (iv)
   LOGICAL, DIMENSION(:), ALLOCATABLE :: obst_l_z0bstress              ! For using a roughness length for obstructions (iv)

   INTEGER, DIMENSION(:), ALLOCATABLE :: obst_z0bstress_option         ! For using various parameterizations of bottom roughness

   CHARACTER(LEN=lchain), DIMENSION(:), ALLOCATABLE :: obst_varname        ! Name of obstructions variables (iv)
   CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE :: obst_type                ! Type of obstruction, one of 'UP', 'DO' or '3D' (iv)
   CHARACTER(LEN=lchain), DIMENSION(:), ALLOCATABLE :: obst_fn_var         ! Name of the input file for obstruction variables
   CHARACTER(LEN=lchain), DIMENSION(:), ALLOCATABLE :: obst_fn_timeserie   ! Name of the file for times-series of obstructions characteristics (iv)
   CHARACTER(LEN=lchain), DIMENSION(:), ALLOCATABLE :: obst_fn_initspatial ! Name of the file for spatial obstructions characteristics (iv)
   CHARACTER(LEN=lchain), DIMENSION(:), ALLOCATABLE :: obst_fn_distrib     ! Name of the file for the vertical distribution of obstructions (iv)

   ! Timeseries variables
   INTEGER, DIMENSION(:), ALLOCATABLE :: obst_ts_tmax                  ! number maximum of value in timeserie
   INTEGER, DIMENSION(:), ALLOCATABLE :: obst_ts_tbefore               ! time index before current time (for interpolation)
   INTEGER, DIMENSION(:), ALLOCATABLE :: obst_ts_tafter                ! time index after current time (for interpolation)
   REAL(KIND=rlg), DIMENSION(:, :), ALLOCATABLE :: obst_ts_time        ! time axis of timeserie
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_ts_height      ! height of timeserie
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_ts_dens        ! height of timeserie
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_ts_width       ! height of timeserie
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_ts_thick       ! height of timeserie

   ! Initialization variables
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_i_height          ! Initial height of obstructions (iv), [m]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_i_width           ! Initial width of obstructions element (iv), [m]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_i_thick           ! Initial thickness of obstructions element (iv), [m]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_i_dens            ! Initial density of obstructions (iv), [m-2]

   ! Constant parameters
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_rho             ! Initial density of obstructions (iv), [kg.m-3]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_drag            ! Drag coefficient for obstructions (iv), [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_lift            ! Lift coefficient for obstructions (iv), [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_fracxy_k0       ! First parameter for horizontal coverage correction (iv), [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_fracxy_k1       ! Second parameter for horizontal coverage correction (iv), [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_fracxy_l        ! Third parameter for horizontal coverage correction (iv), [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_lz              ! obstructions spacing coefficient (iv), [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_shelter         ! sheltering coefficient (iv) [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_crough_x0       ! First coefficient (or cst value) for Abdelrhman (2003) c value [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_crough_x1       ! Second coefficient for Abdelrhman (2003) c value [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_height_x0       ! First parameter for obstruction or Zostera noltii height (iv), [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_height_x1       ! Second parameter for obstruction or Zostera noltii height (iv), [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_z0bstress       ! roughness length for obstructions (iv), [m]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_z0bstress_x0    ! First parameter for roughness length parameterization (iv), [-]
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: obst_c_z0bstress_x1    ! Second parameter for roughness length parameterization (iv), [-]

   !--------------------------------------------------------------------------
   ! * VARIABLES DEPENDING ONLY ON GRID (i,j) or (k,i,j)                     :
   !--------------------------------------------------------------------------

   ! Variables on (i,j)
   !----------------------
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_z0bed           ! z0 for bed without obstructions from source code (used where no obstructions are present) (i,j), [m]
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_z0bstress       ! z0Sed from sedimento (i,j) but modified by obstructions, [m]

   ! Variables on (i,j,k)
   !---------------------
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_fuz           ! Sink term in the momentum equation for 3D formulation (x direction) (k,i,j), [N.m-2]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_fvz           ! Sink term in the momentum equation for 3D formulation (y direction) (k,i,j), [N.m-2]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_t             ! Source term T(z) in y equation (k,i,j), [N.m-2]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_tau           ! Source term tau_veg in y equation (k,i,j), [N.m-2]

   !--------------------------------------------------------------------------
   ! * VARIABLES DEPENDING ON BOTH (iv,i,j) or (iv,k,i,j) or other           :
   !--------------------------------------------------------------------------
   ! Variables on (iv,i,j)
   !----------------------
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_dens_inst       ! Instantaneous obstruction real density (iv,i,j), [m-2]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_width_inst      ! Instantaneous obstruction real width (iv,i,j), [m]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_thick_inst      ! Instantaneous obstruction real thickness (iv,i,j), [m]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_height_inst     ! Instantaneous obstruction real unbend height (iv,i,j), [m]

   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_position      ! The table where the position of the differents variables is defined through occupation rate (iv,i,j), [-]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_height        ! Obstruction height within the domain (iv,i,j), [m]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_fracxy        ! Obstruction correction factor for horizontal coverage (iv,i,j), [-]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_a2d           ! Horizontal area occupied by obstructions per unit area (iv+3,i,j), [-]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_s2d           ! Vertical area occupied by obstructions per unit area (iv+3,i,j), [-]
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: obst_z0obst        ! z0 for each obstructions + for no turb (iv+3,i,j), [m]

   ! Variables on (iv,k,i,j)
   !------------------------
   REAL(KIND=rsh), DIMENSION(:, :, :, :), ALLOCATABLE :: obst_dens3d      ! Obstructions density for use with 3D formulations (iv,k,i,j), [m-2]
   REAL(KIND=rsh), DIMENSION(:, :, :, :), ALLOCATABLE :: obst_width3d     ! Obstructions width (in flow direction) for use with 3D formulations (iv,k,i,j), [m]
   REAL(KIND=rsh), DIMENSION(:, :, :, :), ALLOCATABLE :: obst_thick3d     ! Obstructions thickness (perpendicular to flow direction) for use with 3D formulations (iv,k,i,j), [m]
   REAL(KIND=rsh), DIMENSION(:, :, :, :), ALLOCATABLE :: obst_theta3d     ! Obstruction bending angle (from vertical) (iv,k,i,j), [rad]
   REAL(KIND=rsh), DIMENSION(:, :, :, :), ALLOCATABLE :: obst_fracz3d     ! Fraction of sigma layer occupied by obstruction (iv,k,i,j) [0:1]
   REAL(KIND=rsh), DIMENSION(:, :, :, :), ALLOCATABLE :: obst_drag3d      ! Obstructions drag coefficient used by 3D formulations (iv,k,i,j), [-]

   REAL(KIND=rsh), DIMENSION(:, :, :, :), ALLOCATABLE :: obst_a3d         ! Horizontal area occupied by obstructions at height k per unit area (iv+3,k,i,j), [-]
   REAL(KIND=rsh), DIMENSION(:, :, :, :), ALLOCATABLE :: obst_s3d         ! Vertical area occupied by obstructions at height k per unit area (iv+3,k,i,j), [-]

   ! Variables on (iv,kk)
   !---------------------
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_dens_norm       ! Normalized density from the distribution file (iv,kk) [-]
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_height_norm     ! Normalized heigth from the distribution file (iv,kk) [-]

   ! Variables on (iv,t)
   !--------------------
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_dens_t          ! Instantaneous obstruction density (iv,t), allocated within obst_readfile_char subroutine, [-]
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_width_t         ! Instantaneous obstruction width (iv,t), allocated within obst_readfile_char subroutine, [m]
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_thick_t         ! Instantaneous obstruction thickness (iv,t), allocated within obst_readfile_char subroutine, [m]
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: obst_height_t        ! Instantaneous obstruction heights (iv,t), allocated within obst_readfile_char subroutine, [m]

CONTAINS
!!=============================================================================

#endif

END MODULE com_OBSTRUCTIONS
