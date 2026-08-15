MODULE OBSTRUCTIONS1DV
   !
   ! This module manage the obstruction computation (in one cell "i,j").
   !
   ! No use, no commons, no include this module is independant.
   !
   ! Authors : F.Ganthy
   ! adaptation from MARS version to an independant module S.Le Gac, 2023
   !
   ! PUBLIC SUBROUTINES, VARIABLES and TYPES names of this module are precede by "o1dv_"
   ! All variables are private except o1dv_obst_param
   !
   ! Vertical axis is from 1 at bottom to o1dv_kmax at sea surface
   !

   IMPLICIT NONE
   PRIVATE
   ! VARIABLE
   PUBLIC :: o1dv_obst_param
   ! SUBROUTINE
   PUBLIC :: o1dv_init
   PUBLIC :: o1dv_alloc
   PUBLIC :: o1dv_main
   ! FUNCTION
   PUBLIC :: o1dv_comp_z0sedim
   ! TYPE
   PUBLIC :: o1dv_param_type
   PUBLIC :: o1dv_out_type ! use of structure output to add flexibility for using it in the hydrodynamic code part

   ! Declaration

   ! parameters
   INTEGER, PARAMETER :: lchain = 200 ! max length for string
   INTEGER, PARAMETER :: rsh = 8 ! working precision
   REAL(KIND=rsh), PARAMETER :: pi = 3.14159265358979323846 ! pi value in the module
   REAL(KIND=rsh), PARAMETER :: o1dv_c2turb = 1.92 ! obst_c2turb = beta2 in kepsilon
   REAL(KIND=rsh), PARAMETER :: o1dv_p_hmin = 0.05_rsh ! Minimum coefficient value for obstruction height under bending

   TYPE o1dv_param_type
      ! characteristic of an obstruction
      CHARACTER(LEN=lchain) :: name        ! obstruction name (arbitrary, user defined)
      CHARACTER(LEN=2)  :: type            ! type of obstruction, one of 'UP', 'DO' or '3D'
      LOGICAL       :: l_flexible          ! True if obstruction is flexible
      LOGICAL       :: l_abdelposture      ! True for computation of obstructions posture following abdelrhman 2007
      LOGICAL       :: l_param_height      ! True for use a parameterization (from velocity) for obstruction height
      LOGICAL       :: l_cylinder          ! True for cylinder/ellipse obstruction shape , False for parallelepiped
      LOGICAL       :: l_noturb            ! True for use of simplified formulation (roughness length) instead of full turbulent formulation
      LOGICAL       :: l_drag_cste         ! True for constant drag coefficient, False, for correction from bending angle
      LOGICAL       :: l_abdelrough_cste   ! True for use constant Abdelrhman 2003 coefficient
      LOGICAL       :: l_fracxy            ! True for horizontal coverage correction (grid cell not completely fill with obstructions)
      LOGICAL       :: l_z0bstress         ! True for using a roughness length for obstructions (for sedim module)
      INTEGER       :: fracxy_type         ! The type of correction for horizontal coverage (grid cell not completely fill with obstructions)
      INTEGER       :: c_abdel_nmax        ! Number of segments for Abdlerhman method (bending)
      INTEGER       :: z0bstress_option    ! For using various parameterizations of bottom roughness
      REAL(KIND=rsh) :: c_rho              ! Initial density of obstructions [kg.m-3]
      REAL(KIND=rsh) :: c_height_x0        ! First parameter for obstruction height [-]
      REAL(KIND=rsh) :: c_height_x1        ! Second parameter for obstruction height [-]
      REAL(KIND=rsh) :: c_shelter          ! sheltering coefficient [-]
      REAL(KIND=rsh) :: c_lift             ! Lift coefficient for obstructions [-]
      REAL(KIND=rsh) :: c_drag             ! Drag coefficient for obstructions [-]
      REAL(KIND=rsh) :: c_lz               ! obstructions spacing coefficient [-]
      REAL(KIND=rsh) :: c_crough_x0        ! First coefficient (or cst value) for Abdelrhman (2003) c value [-]
      REAL(KIND=rsh) :: c_crough_x1        ! Second coefficient for Abdelrhman (2003) c value [-]
      REAL(KIND=rsh) :: c_fracxy_k0        ! First parameter for horizontal coverage correction [-]
      REAL(KIND=rsh) :: c_fracxy_k1        ! Second parameter for horizontal coverage correction [-]
      REAL(KIND=rsh) :: c_fracxy_l         ! Third parameter for horizontal coverage correction [-]
      REAL(KIND=rsh) :: c_z0bstress        ! roughness length for obstructions [m]
      REAL(KIND=rsh) :: c_z0bstress_x0     ! First parameter for roughness length parameterization [-]
      REAL(KIND=rsh) :: c_z0bstress_x1     ! Second parameter for roughness length parameterization [-]
      LOGICAL        :: l_filedistri       ! For using a vertical distribution of obstructions
      INTEGER        :: nbhnorm            ! Number of vertical steps for the vertical distribution
      REAL(KIND=rsh), ALLOCATABLE :: height_norm(:) ! Normalized heigth from the vertical distribution [-]
      REAL(KIND=rsh), ALLOCATABLE :: dens_norm(:)   ! Normalized density from the vertical distribution [-]
   END TYPE

   TYPE o1dv_out_type
      ! output of obstruction module
      REAL(KIND=rsh), ALLOCATABLE :: height(:)     ! iv
      REAL(KIND=rsh), ALLOCATABLE :: dens3d(:, :)  ! iv,k
      REAL(KIND=rsh), ALLOCATABLE :: width3d(:, :) ! iv,k
      REAL(KIND=rsh), ALLOCATABLE :: thick3d(:, :) ! iv,k
      REAL(KIND=rsh), ALLOCATABLE :: theta3d(:, :) ! iv,k
      REAL(KIND=rsh), ALLOCATABLE :: fracz3d(:, :) ! iv,k
      REAL(KIND=rsh), ALLOCATABLE :: drag3d(:, :)  ! iv,k
      REAL(KIND=rsh), ALLOCATABLE :: fracxy(:)     ! iv
      REAL(KIND=rsh), ALLOCATABLE :: fuz(:)        ! k
      REAL(KIND=rsh), ALLOCATABLE :: fvz(:)        ! k
      REAL(KIND=rsh), ALLOCATABLE :: tau(:)        ! k
      REAL(KIND=rsh), ALLOCATABLE :: t(:)          ! k
      REAL(KIND=rsh), ALLOCATABLE :: a3d(:, :)     ! iv+3,k
      REAL(KIND=rsh), ALLOCATABLE :: s3d(:, :)     ! iv+3,k
      REAL(KIND=rsh), ALLOCATABLE :: a2d(:)        ! iv+3
      REAL(KIND=rsh), ALLOCATABLE :: s2d(:)        ! iv+3
      REAL(KIND=rsh), ALLOCATABLE :: z0obst(:)     ! iv+3
      REAL(KIND=rsh) :: height_mean
      REAL(KIND=rsh) :: dens_mean
      REAL(KIND=rsh) :: width_mean
   END TYPE

   ! variables of the module

   INTEGER :: o1dv_nbvar    ! Total number of obstruction variables
   INTEGER :: o1dv_kmax     ! Number of vertical layer in water (from hydro model)
   INTEGER :: o1dv_nb_max_hnorm ! Maximum number of vertical steps from the vertical distribution
   INTEGER :: o1dv_abdel_nmax ! Number of segments for Abdlerhman method (bending)
   INTEGER :: o1dv_iwarnlog ! Logging io
   REAL(KIND=rsh) :: o1dv_h0fond ! Threshold, no computation for water depth below this range
   REAL(KIND=rsh) :: o1dv_c_paramhuv ! The coefficient of obstruction height for computation of velocity
   TYPE(o1dv_param_type), DIMENSION(:), ALLOCATABLE :: o1dv_obst_param ! Obstruction parameters
   TYPE(o1dv_out_type) :: o1dv_output_reset
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: o1dv_zz
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: o1dv_uvz
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: o1dv_abdel_fx         ! Horizontal force acting on leaves seagment (Abdelrhman method for bending) (iv,o1dv_c_abdel_nmax)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: o1dv_abdel_fz         ! Vertical force acting on leaves seagment (Abdelrhman method for bending) (iv,o1dv_c_abdel_nmax)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: o1dv_abdel_zcent      ! Z position of segment centres (Abdelrhman method for bending) (iv,o1dv_c_abdel_nmax)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: o1dv_abdel_t0cent     ! First bending angles of segment centres (Abdelrhman method for bending) (iv,o1dv_c_abdel_nmax)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: o1dv_abdel_t1cent     ! Second bending angles of segment centres (Abdelrhman method for bending) (iv,o1dv_c_abdel_nmax)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: o1dv_abdel_dtheta     ! Difference in bending angles of segment centres between two iterations (Abdelrhman method for bending) (iv,o1dv_c_abdel_nmax)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: o1dv_abdel_uvcent     ! Local velocity at segment centres (Abdelrhman method for bending) (iv,o1dv_c_abdel_nmax)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: o1dv_abdel_zn         ! Height of end of each segment (Abdelrhman method for bending) (iv,0:o1dv_abdel_c_nmax+1)
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: o1dv_abdel_tn

CONTAINS

   !==============================================================================
   SUBROUTINE o1dv_alloc(nbvar, kmax, nb_max_hnorm)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_alloc  ***
      !!
      !! ** Purpose : Initialize the dimensions
      !!---------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nbvar  ! number of obstructions
      INTEGER, INTENT(IN) :: kmax   ! number of layers in water column
      INTEGER, INTENT(IN) :: nb_max_hnorm   ! maximum number of vertical steps from the vertical distribution

      INTEGER :: iv

      o1dv_nbvar = nbvar
      o1dv_kmax = kmax
      o1dv_nb_max_hnorm = nb_max_hnorm
      ALLOCATE (o1dv_obst_param(nbvar))
      DO iv = 1, nbvar
         ALLOCATE (o1dv_obst_param(iv)%height_norm(nb_max_hnorm))
         ALLOCATE (o1dv_obst_param(iv)%dens_norm(nb_max_hnorm))
      END DO

      RETURN
   END SUBROUTINE o1dv_alloc

   !==============================================================================
   SUBROUTINE o1dv_init(h0fond, paramhuv, &
                        o1dv_varname, o1dv_type, o1dv_l_flexible, o1dv_l_abdelposture, &
                        o1dv_l_param_height, o1dv_l_cylinder, o1dv_l_noturb, o1dv_l_drag_cste, &
                        o1dv_l_abdelrough_cste, o1dv_l_fracxy, o1dv_l_z0bstress, &
                        o1dv_fracxy_type, o1dv_c_abdel_nmax, o1dv_z0bstress_option, &
                        o1dv_c_rho, o1dv_c_height_x0, o1dv_c_height_x1, o1dv_c_shelter, &
                        o1dv_c_lift, o1dv_c_drag, o1dv_c_lz, o1dv_c_crough_x0, o1dv_c_crough_x1, &
                        o1dv_c_fracxy_k0, o1dv_c_fracxy_k1, o1dv_c_fracxy_l, o1dv_c_z0bstress, &
                        o1dv_c_z0bstress_x0, o1dv_c_z0bstress_x1, &
                        o1dv_l_filedistri, o1dv_nbhnorm, o1dv_height_norm, o1dv_dens_norm, &
                        o1dv_stdout)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_init  ***
      !!
      !! ** Purpose : Allocate the global arrays of the module.
      !!---------------------------------------------------------------------
      IMPLICIT NONE

      REAL, INTENT(IN) :: h0fond
      REAL, INTENT(IN) :: paramhuv
      CHARACTER(LEN=lchain), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_varname
      CHARACTER(LEN=2), DIMENSION(o1dv_nbvar), INTENT(IN)  :: o1dv_type
      LOGICAL, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_l_flexible
      LOGICAL, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_l_abdelposture
      LOGICAL, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_l_param_height
      LOGICAL, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_l_cylinder
      LOGICAL, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_l_noturb
      LOGICAL, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_l_drag_cste
      LOGICAL, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_l_abdelrough_cste
      LOGICAL, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_l_fracxy
      LOGICAL, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_l_z0bstress
      INTEGER, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_fracxy_type
      INTEGER, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_c_abdel_nmax
      INTEGER, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_z0bstress_option
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_rho
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_height_x0
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_height_x1
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_shelter
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_lift
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_drag
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_lz
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_crough_x0
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_crough_x1
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_fracxy_k0
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_fracxy_k1
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_fracxy_l
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_z0bstress
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_z0bstress_x0
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_c_z0bstress_x1
      LOGICAL, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_l_filedistri
      INTEGER, DIMENSION(o1dv_nbvar), INTENT(IN)       :: o1dv_nbhnorm
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar, o1dv_nb_max_hnorm), INTENT(IN) :: o1dv_height_norm
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar, o1dv_nb_max_hnorm), INTENT(IN) :: o1dv_dens_norm
      INTEGER, INTENT(IN) :: o1dv_stdout

      INTEGER :: iv

      o1dv_h0fond = h0fond
      o1dv_iwarnlog = o1dv_stdout
      o1dv_c_paramhuv = paramhuv
      o1dv_abdel_nmax = MAXVAL(o1dv_c_abdel_nmax)

      DO iv = 1, o1dv_nbvar
         o1dv_obst_param(iv)%name = o1dv_varname(iv)
         o1dv_obst_param(iv)%type = o1dv_type(iv)
         o1dv_obst_param(iv)%l_flexible = o1dv_l_flexible(iv)
         o1dv_obst_param(iv)%l_abdelposture = o1dv_l_abdelposture(iv)
         o1dv_obst_param(iv)%l_param_height = o1dv_l_param_height(iv)
         o1dv_obst_param(iv)%l_cylinder = o1dv_l_cylinder(iv)
         o1dv_obst_param(iv)%l_noturb = o1dv_l_noturb(iv)
         o1dv_obst_param(iv)%l_drag_cste = o1dv_l_drag_cste(iv)
         o1dv_obst_param(iv)%l_abdelrough_cste = o1dv_l_abdelrough_cste(iv)
         o1dv_obst_param(iv)%l_fracxy = o1dv_l_fracxy(iv)
         o1dv_obst_param(iv)%l_z0bstress = o1dv_l_z0bstress(iv)
         o1dv_obst_param(iv)%fracxy_type = o1dv_fracxy_type(iv)
         o1dv_obst_param(iv)%c_abdel_nmax = o1dv_c_abdel_nmax(iv)
         o1dv_obst_param(iv)%z0bstress_option = o1dv_z0bstress_option(iv)
         o1dv_obst_param(iv)%c_rho = o1dv_c_rho(iv)
         o1dv_obst_param(iv)%c_height_x0 = o1dv_c_height_x0(iv)
         o1dv_obst_param(iv)%c_height_x1 = o1dv_c_height_x1(iv)
         o1dv_obst_param(iv)%c_shelter = o1dv_c_shelter(iv)
         o1dv_obst_param(iv)%c_lift = o1dv_c_lift(iv)
         o1dv_obst_param(iv)%c_drag = o1dv_c_drag(iv)
         o1dv_obst_param(iv)%c_lz = o1dv_c_lz(iv)
         o1dv_obst_param(iv)%c_crough_x0 = o1dv_c_crough_x0(iv)
         o1dv_obst_param(iv)%c_crough_x1 = o1dv_c_crough_x1(iv)
         o1dv_obst_param(iv)%c_fracxy_k0 = o1dv_c_fracxy_k0(iv)
         o1dv_obst_param(iv)%c_fracxy_k1 = o1dv_c_fracxy_k1(iv)
         o1dv_obst_param(iv)%c_fracxy_l = o1dv_c_fracxy_l(iv)
         o1dv_obst_param(iv)%c_z0bstress = o1dv_c_z0bstress(iv)
         o1dv_obst_param(iv)%c_z0bstress_x0 = o1dv_c_z0bstress_x0(iv)
         o1dv_obst_param(iv)%c_z0bstress_x1 = o1dv_c_z0bstress_x1(iv)
         o1dv_obst_param(iv)%l_filedistri = o1dv_l_filedistri(iv)
         o1dv_obst_param(iv)%nbhnorm = o1dv_nbhnorm(iv)
         o1dv_obst_param(iv)%height_norm(:) = o1dv_height_norm(iv, :)
         o1dv_obst_param(iv)%dens_norm(:) = o1dv_dens_norm(iv, :)
      END DO

      o1dv_output_reset = o1dv_init_output()
      IF (ANY(o1dv_l_abdelposture)) THEN
         ALLOCATE (o1dv_zz(0:o1dv_kmax + 1))
         ALLOCATE (o1dv_uvz(0:o1dv_kmax + 1))
         ALLOCATE (o1dv_abdel_fx(1:o1dv_abdel_nmax))
         ALLOCATE (o1dv_abdel_fz(1:o1dv_abdel_nmax))
         ALLOCATE (o1dv_abdel_zcent(1:o1dv_abdel_nmax))
         ALLOCATE (o1dv_abdel_t0cent(1:o1dv_abdel_nmax))
         ALLOCATE (o1dv_abdel_t1cent(1:o1dv_abdel_nmax))
         ALLOCATE (o1dv_abdel_dtheta(1:o1dv_abdel_nmax))
         ALLOCATE (o1dv_abdel_uvcent(1:o1dv_abdel_nmax))
         ALLOCATE (o1dv_abdel_zn(0:o1dv_abdel_nmax + 1))
         ALLOCATE (o1dv_abdel_tn(0:o1dv_abdel_nmax + 1))
      END IF

      RETURN
   END SUBROUTINE o1dv_init

   !==============================================================================
   FUNCTION o1dv_init_output() result(o1dv_output)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_init_output  ***
      !!
      !! ** Purpose : Initialize output by allocating arrays and fill with zeros
      !!
      !! ** Description : Allocate with o1dv_kmax and o1dv_nbvar
      !!---------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      TYPE(o1dv_out_type) :: o1dv_output

      ALLOCATE (o1dv_output%height(o1dv_nbvar))
      o1dv_output%height(:) = 0.
      ALLOCATE (o1dv_output%dens3d(o1dv_nbvar, o1dv_kmax))
      o1dv_output%dens3d(:, :) = 0.
      ALLOCATE (o1dv_output%width3d(o1dv_nbvar, o1dv_kmax))
      o1dv_output%width3d(:, :) = 0.
      ALLOCATE (o1dv_output%thick3d(o1dv_nbvar, o1dv_kmax))
      o1dv_output%thick3d(:, :) = 0.
      ALLOCATE (o1dv_output%theta3d(o1dv_nbvar, o1dv_kmax))
      o1dv_output%theta3d(:, :) = 0.
      ALLOCATE (o1dv_output%fracz3d(o1dv_nbvar, o1dv_kmax))
      o1dv_output%fracz3d(:, :) = 0.
      ALLOCATE (o1dv_output%drag3d(o1dv_nbvar, o1dv_kmax))
      o1dv_output%drag3d(:, :) = 0.
      ALLOCATE (o1dv_output%fracxy(o1dv_nbvar))
      o1dv_output%fracxy(:) = 0.
      ALLOCATE (o1dv_output%fuz(o1dv_kmax))
      o1dv_output%fuz(:) = 0.
      ALLOCATE (o1dv_output%fvz(o1dv_kmax))
      o1dv_output%fvz(:) = 0.
      ALLOCATE (o1dv_output%tau(o1dv_kmax))
      o1dv_output%tau(:) = 0.
      ALLOCATE (o1dv_output%t(o1dv_kmax))
      o1dv_output%t(:) = 0.
      ALLOCATE (o1dv_output%a3d(o1dv_nbvar + 3, o1dv_kmax))
      o1dv_output%a3d(:, :) = 0.
      ALLOCATE (o1dv_output%s3d(o1dv_nbvar + 3, o1dv_kmax))
      o1dv_output%s3d(:, :) = 0.
      ALLOCATE (o1dv_output%a2d(o1dv_nbvar + 3))
      o1dv_output%a2d(:) = 0.
      ALLOCATE (o1dv_output%s2d(o1dv_nbvar + 3))
      o1dv_output%s2d(:) = 0.
      ALLOCATE (o1dv_output%z0obst(o1dv_nbvar + 3))
      o1dv_output%z0obst(:) = 0.
      o1dv_output%height_mean = 0.
      o1dv_output%dens_mean = 0.
      o1dv_output%width_mean = 0.

   END FUNCTION o1dv_init_output

   !==============================================================================
   SUBROUTINE o1dv_main(o1dv_hwat, o1dv_cmu, o1dv_z0bed, o1dv_uz, o1dv_vz, o1dv_dz, o1dv_zc, &
                        o1dv_height_f, o1dv_dens_f, o1dv_thick_f, o1dv_width_f, o1dv_position, o1dv_height_p, &
                        o1dv_output)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_main  ***
      !!
      !! ** Purpose : From hydrodynamic variables and OBSTRUCTIONS parameters,
      !!              compute OBSTRUCTIONS output variable
      !!
      !! ** Description : Main subroutine of this module, entry point for 
      !!                  hydrodynamic model
      !!---------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      REAL(KIND=rsh), INTENT(IN) :: o1dv_hwat
      REAL(KIND=rsh), INTENT(IN) :: o1dv_cmu
      REAL(KIND=rsh), INTENT(IN) :: o1dv_z0bed
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_uz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_vz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_dz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_zc
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_height_f
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_dens_f
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_thick_f
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_width_f
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_position
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_height_p ! previous value of height for o1dv_comp_height subroutine
      TYPE(o1dv_out_type), INTENT(INOUT) :: o1dv_output

      TYPE(o1dv_param_type) :: obst
      INTEGER :: iv

      ! initialize output with zero
      o1dv_output = o1dv_output_reset

      IF (o1dv_hwat .gt. o1dv_h0fond) THEN

         DO iv = 1, o1dv_nbvar
            obst = o1dv_obst_param(iv)
            IF (o1dv_position(iv) .gt. 0.0_rsh) THEN
               IF (obst%l_flexible) THEN
                  IF (obst%l_abdelposture) THEN
                     CALL o1dv_abdelposture(obst, o1dv_hwat, &
                                            o1dv_uz, o1dv_vz, o1dv_dz, o1dv_zc, &
                                            o1dv_height_f(iv), o1dv_dens_f(iv), o1dv_thick_f(iv), o1dv_width_f(iv), &
                                            o1dv_output%height(iv), o1dv_output%theta3d(iv, :))
                  ELSE
                     o1dv_output%height(iv) = o1dv_height_p(iv)

                     CALL o1dv_comp_height(obst, o1dv_hwat, &
                                           o1dv_uz, o1dv_vz, o1dv_dz, o1dv_zc, &
                                           o1dv_height_f(iv), o1dv_height_p(iv), o1dv_output%height(iv))

                     CALL o1dv_comp_theta(obst, o1dv_hwat, o1dv_dz, o1dv_zc, &
                                          o1dv_height_f(iv), o1dv_output%height(iv), o1dv_output%theta3d(iv, :))

                  END IF
               ELSE ! Obstruction is not flexible
                  o1dv_output%height(iv) = o1dv_height_f(iv)
                  o1dv_output%theta3d(iv, :) = 0.0_rsh
               END IF

               CALL o1dv_comp_distrib(obst, o1dv_hwat, o1dv_dz, o1dv_zc, &
                                      o1dv_dens_f(iv), o1dv_output%height(iv), o1dv_output%dens3d(iv, :))

               CALL o1dv_comp_fracz(obst, o1dv_hwat, o1dv_dz, o1dv_zc, o1dv_output%height(iv), &
                                    o1dv_output%dens3d(iv, :), o1dv_output%fracz3d(iv, :))

               CALL o1dv_comp_diam(obst, o1dv_dz, &
                                   o1dv_height_f(iv), o1dv_thick_f(iv), o1dv_width_f(iv), &
                                   o1dv_output%dens3d(iv, :), o1dv_output%theta3d(iv, :), o1dv_output%fracz3d(iv, :), &
                                   o1dv_output%width3d(iv, :), o1dv_output%thick3d(iv, :))

               CALL o1dv_comp_fracxy(obst, o1dv_position(iv), o1dv_output%fracxy(iv))

            ELSE ! No obstruction within this cell
               o1dv_output%height(iv) = 0.0_rsh
               o1dv_output%theta3d(iv, :) = 0.0_rsh
            END IF ! END test on o1dv_position (presence of obstructions)

         END DO ! nbvar

         CALL o1dv_comp_projarea(o1dv_position, o1dv_dz, o1dv_output)

         CALL o1dv_comp_obstroughness(o1dv_position, o1dv_z0bed, &
                                      o1dv_height_f, o1dv_dens_f, o1dv_width_f, o1dv_output)

         CALL o1dv_comp_hydroparam(o1dv_position, o1dv_cmu, o1dv_uz, o1dv_vz, o1dv_output)

      ELSE ! Not enough water
         DO iv = 1, o1dv_nbvar
            IF (o1dv_position(iv) .gt. 0.0_rsh) THEN
               o1dv_output%height(iv) = 0.99_rsh*o1dv_height_f(iv)
            END IF
         END DO
         o1dv_output%z0obst(:) = o1dv_z0bed
      END IF

   END SUBROUTINE o1dv_main

   !==============================================================================
   SUBROUTINE o1dv_abdelposture(obst, o1dv_hwat, o1dv_uz, o1dv_vz, o1dv_dz, o1dv_zc, &
                                o1dv_height_inst, o1dv_dens_inst, o1dv_thick_inst, o1dv_width_inst, &
                                o1dv_height, o1dv_theta3d)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_abdelposture  ***
      !!
      !! ** Purpose : Computes obstruction posture (height, diameters and bending angle)
      !!
      !! ** Description : Based on the balance between forces of drag, lift, friction,
      !!                  weight and buoyancy of single obstruction element discretized
      !!                  along the vertical, the Abdelrhman's model defined the bending
      !!                  of obstruction elements
      !!
      !! ** Reference : Abdelrhman M.A. (2007) Modeling coupling between eelgrass
      !!                Zostera marina and water flow. Mar. Ecol. Prog. Ser. 338:81-96
      !!
      !!---------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      TYPE(o1dv_param_type), INTENT(IN) :: obst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_hwat
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_uz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_vz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_dz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_zc
      REAL(KIND=rsh), INTENT(IN) :: o1dv_height_inst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_dens_inst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_thick_inst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_width_inst
      REAL(KIND=rsh), INTENT(INOUT) :: o1dv_height
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(INOUT) :: o1dv_theta3d

      !! * Local declaration
      LOGICAL                     :: l_theta
      INTEGER                     :: k, niter, niter_eff, s
      INTEGER                     :: k0, k1
      REAL(KIND=rsh)              :: z0, z1, zc, zt0, zt1, dz, udz, uv, htmp
      REAL(KIND=rsh)              :: lseg, cd, cl, cf, dtheti, phi, re, cshelt
      REAL(KIND=rsh)              :: th0, th0p, th0m, th1, th1p, th1m, sm0, sm0p, sm0m, sm1, sm1p, sm1m
      REAL(KIND=rsh)              :: drag_xz_0, drag_xz_1m, drag_xz_1p
      REAL(KIND=rsh)              :: lift_xz_0, lift_xz_1m, lift_xz_1p
      REAL(KIND=rsh)              :: buoy_xz_0, buoy_xz_1m, buoy_xz_1p
      REAL(KIND=rsh)              :: drag_x, fric_x, lift_z, buoy_z, fric_z

      REAL(KIND=rsh), PARAMETER   :: g = 9.81_rsh
      REAL(KIND=rsh), PARAMETER   :: rw = 1025.0_rsh
      REAL(KIND=rsh), PARAMETER   :: nu = 1.4E-6_rsh
      REAL(KIND=rsh), PARAMETER   :: gamma = 1.0e-6_rsh
      INTEGER, PARAMETER          :: niter_max = 25                    ! Maximum number of iterations
      REAL(KIND=rsh), PARAMETER   :: dtheta_max = 0.25_rsh              ! Angle variation (degree) to reach stability
      REAL(KIND=rsh), PARAMETER   :: dtheta_iter = 10.0_rsh*pi/180.0_rsh ! Searching step
      REAL(KIND=rsh), PARAMETER   :: phi_max = 10.0_rsh*pi/180.0_rsh ! Max angle for lift
      !!----------------------------------------------------------------------
      !! * Executable part

      IF (o1dv_hwat <= o1dv_height_inst) THEN
         o1dv_height = o1dv_hwat
         o1dv_theta3d(:) = ACOS(o1dv_height/o1dv_height_inst)
      ELSE
         !--------------------------
         ! **** Initializations ****
         !--------------------------
         lseg = o1dv_height_inst/REAL(obst%c_abdel_nmax, rsh)
         o1dv_abdel_zcent(:) = 0.0_rsh
         o1dv_abdel_t0cent(:) = 0.0_rsh
         o1dv_abdel_t1cent(:) = 0.0_rsh
         o1dv_abdel_dtheta(:) = 0.0_rsh
         o1dv_abdel_uvcent(:) = 0.0_rsh
         o1dv_abdel_zn(:) = 0.0_rsh
         o1dv_abdel_tn(:) = 0.0_rsh
         !---------------------------------------
         ! **** Computes velocity k=0:kmax+1 ****
         !---------------------------------------
         ! For k=0
         k = 0
         o1dv_zz(k) = 0.0_rsh
         o1dv_uvz(k) = 0.0_rsh
         ! For k=1,o1dv_kmax
         DO k = 1, o1dv_kmax
            o1dv_uvz(:) = SQRT(o1dv_uz(k)**2.0_rsh + o1dv_vz(k)**2.0_rsh)
            o1dv_zz(k) = o1dv_zc(k) + o1dv_dz(k)*0.5_rsh
         END DO
         ! For k=kmax+1
         k = o1dv_kmax + 1
         o1dv_zz(k) = o1dv_hwat
         o1dv_uvz(k) = o1dv_uvz(k - 1)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! **** STARTING ITERATIVE PROCEDURE OVER ALL SEGMENTS **** !
         ! ****      TO REACH A STABLE OBSTRUCTION STATE       **** !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         niter = 0
         DO WHILE (niter .LT. niter_max)
            !-------------------------
            ! Some re-initializations
            !------------------------
            o1dv_abdel_fx(:) = 0.0_rsh
            o1dv_abdel_fz(:) = 0.0_rsh
            !---------------------------------------------------------------------------!
            ! *** COMPUTE HEIGHT AT EACH HALF-SEGMENT THROUGOUT AN UPWARD PROCEDURE *** !
            !---------------------------------------------------------------------------!
            DO s = 1, obst%c_abdel_nmax
               ! Height at half segment
               IF (s .EQ. 1) THEN ! First segment
                  o1dv_abdel_zcent(s) = (lseg/2.0_rsh)*COS(o1dv_abdel_t0cent(s))
               ELSE             ! Other segments
                  o1dv_abdel_zcent(s) = (o1dv_abdel_zcent(s - 1) + (lseg/2.0)*COS(o1dv_abdel_t0cent(s - 1))) + &
                                        (lseg/2.0_rsh)*COS(o1dv_abdel_t0cent(s))
               END IF
               ! Velocity at height segment
               k = 0
               DO WHILE ((k .LT. o1dv_kmax) .AND. (o1dv_zz(k) .LE. o1dv_abdel_zcent(s)))
                  k = k + 1
               END DO
               k0 = MAX(MIN(k - 1, o1dv_kmax - 1), 0)
               k1 = MIN(MAX(k, 1), o1dv_kmax)
               o1dv_abdel_uvcent(s) = o1dv_uvz(k0) + (o1dv_abdel_zcent(s) - o1dv_zz(k0))* &
                                      (o1dv_uvz(k1) - o1dv_uvz(k0))/(o1dv_zz(k1) - o1dv_zz(k0))
            END DO
            !----------------------------------------------------------!
            ! *** GENERAL DOWNWARD PROCEDURE FOR EACH LEAF SEGMENT *** !
            !----------------------------------------------------------!
            DO s = obst%c_abdel_nmax, 1, -1 ! Downward loop on each leaf segment
               !-------------------------------------------------------------!
               ! * Starting iterative procedure to search theta by solving * !
               ! *  iteratively Eq. 7 or Eq. 10 of Abdelrhman 2007 with    * !
               ! * with values starting at theta=0 and ending at pi/2 with * !
               ! * its symetry at -pi/2, respectively for th0m,th0p,th1m,  * !
               ! *              th1p and sm0p,sm0m,sm1p,sm1m               * !
               !-------------------------------------------------------------!
               th0 = 0.0_rsh ! Starting angle = 0 (vertical segment)
               l_theta = .TRUE.  ! Condition to stop iterative research of theta value
               !------------------------------------------------------------!
               ! Computing xz contributions and summing all moments for th0 !
               !------------------------------------------------------------!
               phi = ABS((pi/2.0_rsh) - th0)
               cshelt = MIN(MAX(1.0_rsh, obst%c_shelter*o1dv_dens_inst* &
                                o1dv_width_inst*o1dv_thick_inst* &
                                lseg/(lseg*COS(th0)) - 1.0_rsh), 4.0_rsh)
               ! For drag
               cd = (phi*obst%c_drag/(pi/2.0_rsh))/cshelt
               drag_xz_0 = 0.5_rsh*cd*rw*(o1dv_abdel_uvcent(s)**2.0_rsh)* &
                           o1dv_width_inst*((lseg**2.0_rsh)/2.0_rsh)*(COS(th0)**2.0_rsh)
               ! For lift
               IF (phi .LE. phi_max) THEN
                  cl = (obst%c_lift*(phi/phi_max))/cshelt
               ELSE
                  cl = 0.0_rsh
               END IF
               lift_xz_0 = 0.5_rsh*cl*rw*(o1dv_abdel_uvcent(s)**2.0_rsh)* &
                           o1dv_width_inst*((lseg**2.0_rsh)/2.0_rsh)*(SIN(th0)**2.0_rsh)
               ! For buoyancy
               buoy_xz_0 = (rw - obst%c_rho)*g*o1dv_width_inst* &
                           o1dv_thick_inst*((lseg**2.0_rsh)/2.0_rsh)*SIN(th0)
               ! Now summing all moments
               IF (s .EQ. obst%c_abdel_nmax) THEN
                  sm0p = drag_xz_0 + lift_xz_0 - buoy_xz_0
               ELSE
                  sm0p = drag_xz_0 + lift_xz_0 - buoy_xz_0 + &
                         o1dv_abdel_fx(s + 1)*lseg*COS(th0) - &
                         o1dv_abdel_fz(s + 1)*lseg*SIN(th0)
               END IF
               sm0m = sm0p
               th0m = th0
               th0p = th0
               IF (o1dv_abdel_uvcent(s) .LE. 0.001) THEN
                  !------------------------------!
                  ! Not enough velocity, theta=0 !
                  !------------------------------!
                  o1dv_abdel_t1cent(s) = 0.0_rsh
                  o1dv_abdel_fx(s) = 0.0_rsh
                  o1dv_abdel_fz(s) = 0.0_rsh
               ELSE
                  !----------------------------------------!
                  ! Effective start of iterative procedure !
                  !----------------------------------------!
                  DO WHILE (l_theta)
                     ! Incrementation of new theta values
                     th1m = th0m - dtheta_iter
                     th1p = th0p + dtheta_iter
                     !-------------------------------------------------------------!
                     ! Computing xz contributions and summing all moments for th1m !
                     !-------------------------------------------------------------!
                     phi = ABS((pi/2.0_rsh) - th1m)
                     cshelt = MIN(MAX(1.0_rsh, obst%c_shelter*o1dv_dens_inst* &
                                      o1dv_width_inst*o1dv_thick_inst* &
                                      lseg/(lseg*COS(th1m)) - 1.0_rsh), 4.0_rsh)
                     ! For drag
                     cd = (phi*obst%c_drag/(pi/2.0_rsh))/cshelt
                     drag_xz_1m = 0.5_rsh*cd*rw*(o1dv_abdel_uvcent(s)**2.0_rsh)* &
                                  o1dv_width_inst*((lseg**2.0_rsh)/2.0_rsh)*(COS(th1m)**2.0_rsh)
                     ! For lift
                     IF (phi .LE. phi_max) THEN
                        cl = (obst%c_lift*(phi/phi_max))/cshelt
                     ELSE
                        cl = 0.0_rsh
                     END IF
                     lift_xz_1m = 0.5_rsh*cl*rw*(o1dv_abdel_uvcent(s)**2.0_rsh)* &
                                  o1dv_width_inst*((lseg**2.0_rsh)/2.0_rsh)*(SIN(th1m)**2.0_rsh)
                     ! For buoyancy
                     buoy_xz_1m = (rw - obst%c_rho)*g*o1dv_width_inst* &
                                  o1dv_thick_inst*((lseg**2.0_rsh)/2.0_rsh)*SIN(th1m)
                     ! Now summing all moments
                     IF (s .EQ. obst%c_abdel_nmax) THEN
                        sm1m = drag_xz_1m + lift_xz_1m - buoy_xz_1m
                     ELSE
                        sm1m = drag_xz_1m + lift_xz_1m - buoy_xz_1m + &
                               o1dv_abdel_fx(s + 1)*lseg*COS(th1m) - &
                               o1dv_abdel_fz(s + 1)*lseg*SIN(th1m)
                     END IF
                     !-------------------------------------------------------------!
                     ! Computing xz contributions and summing all moments for th1p !
                     !-------------------------------------------------------------!
                     phi = ABS((pi/2.0_rsh) - th1p)
                     cshelt = MIN(MAX(1.0_rsh, obst%c_shelter*o1dv_dens_inst* &
                                      o1dv_width_inst*o1dv_thick_inst* &
                                      lseg/(lseg*COS(th1p)) - 1.0_rsh), 4.0_rsh)
                     ! For drag
                     cd = (phi*obst%c_drag/(pi/2.0_rsh))/cshelt
                     drag_xz_1p = 0.5_rsh*cd*rw*(o1dv_abdel_uvcent(s)**2.0_rsh)* &
                                  o1dv_width_inst*((lseg**2.0_rsh)/2.0_rsh)*(COS(th1p)**2.0_rsh)
                     ! For lift
                     IF (phi .LE. phi_max) THEN
                        cl = (obst%c_lift*(phi/phi_max))/cshelt
                     ELSE
                        cl = 0.0_rsh
                     END IF
                     lift_xz_1p = 0.5_rsh*cl*rw*(o1dv_abdel_uvcent(s)**2.0_rsh)* &
                                  o1dv_width_inst*((lseg**2.0_rsh)/2.0_rsh)*(SIN(th1p)**2.0_rsh)
                     ! For buoyancy
                     buoy_xz_1p = (rw - obst%c_rho)*g*o1dv_width_inst* &
                                  o1dv_thick_inst*((lseg**2.0_rsh)/2.0_rsh)*SIN(th1p)
                     ! Now summing all moments
                     IF (s .EQ. obst%c_abdel_nmax) THEN
                        sm1p = drag_xz_1p + lift_xz_1p - buoy_xz_1p
                     ELSE
                        sm1p = drag_xz_1p + lift_xz_1p - buoy_xz_1p + &
                               o1dv_abdel_fx(s + 1)*lseg*COS(th1p) - &
                               o1dv_abdel_fz(s + 1)*lseg*SIN(th1p)
                     END IF
                     !------------------------------------------!
                     ! Now testing if there sign change between !
                     ! sm0p and sm1p or between sm0m and sm1m   !
                     !------------------------------------------!
                     IF (((sm0m .GE. 0.0_rsh) .AND. (sm1m .LE. 0.0_rsh)) .OR. ((sm0m .LE. 0.0_rsh) .AND. (sm1m .GE. 0.0_rsh))) THEN
                        l_theta = .FALSE.
                        th0 = ABS(th0m)
                        th1 = ABS(th1m)
                        sm0 = sm0m
                        sm1 = sm1m
                     ELSEIF (((sm0p .GE. 0.0_rsh) .AND. (sm1p .LE. 0.0_rsh)) .OR. &
                             ((sm0p .LE. 0.0_rsh) .AND. (sm1p .GE. 0.0_rsh))) THEN
                        l_theta = .FALSE.
                        th0 = ABS(th0p)
                        th1 = ABS(th1p)
                        sm0 = sm0p
                        sm1 = sm1p
                     ELSEIF (th1p .GT. (pi/2.0_rsh) + dtheta_iter) THEN
                        IF (o1dv_abdel_uvcent(s) .LE. 0.1_rsh) THEN
                           th0 = -(pi/4.0_rsh)
                           th1 = (pi/4.0_rsh)
                           sm0 = -1.0_rsh
                           sm1 = 1.0_rsh
                        ELSE
                           th0 = (pi/2.0_rsh) - (pi/4.0_rsh)
                           th1 = (pi/2.0_rsh) + (pi/4.0_rsh)
                           sm0 = -1.0_rsh
                           sm1 = 1.0_rsh
                        END IF
                        l_theta = .FALSE.
                        WRITE (o1dv_iwarnlog, *) ' '
                        WRITE (o1dv_iwarnlog, *) ' '
                        WRITE (o1dv_iwarnlog, *) '**************************************************************************'
                        WRITE (o1dv_iwarnlog, *) '***** module OBSTRUCTIONS, subroutine o1dv_COMP_ABDELPOSTURE *****'
                        WRITE (o1dv_iwarnlog, *) ' WARNING : no solution was found for the sum of moment'
                        WRITE (o1dv_iwarnlog, *) '           At obst', obst%name, 's', s
                        WRITE (o1dv_iwarnlog, *) '           o1dv_hwat', o1dv_hwat, 'o1dv_height_inst', o1dv_height_inst
                        WRITE (o1dv_iwarnlog, *) '           uzvz', o1dv_abdel_uvcent(s)
                        WRITE (o1dv_iwarnlog, *) ' --> Th=0 or Th=pi/2 applied depending on local velocity !!! '
                        WRITE (o1dv_iwarnlog, *) '**************************************************************************'
                     ELSE
                        th0m = th1m
                        th0p = th1p
                        sm0m = sm1m
                        sm0p = sm1p
                     END IF
                  END DO ! End of iterative loop for theta value
                  !-------------------------------------------!
                  ! Linear interpolation of final theta value !
                  !-------------------------------------------!
                  o1dv_abdel_t1cent(s) = th0 - sm0*(th1 - th0)/(sm1 - sm0)
                  !-----------------------------------------------------!
                  ! Apply this good theta value to compute forces which !
                  !     will acts to the next obstruction segment       !
                  !-----------------------------------------------------!
                  IF (obst%c_abdel_nmax .GT. 1) THEN
                     phi = ABS((pi/2.0_rsh) - o1dv_abdel_t1cent(s))
                     cshelt = MIN(MAX(1.0_rsh, obst%c_shelter*o1dv_dens_inst* &
                                      o1dv_width_inst*o1dv_thick_inst* &
                                      lseg/(lseg*COS(o1dv_abdel_t1cent(s))) - 1.0_rsh), 4.0_rsh)
                     ! Computing drag
                     cd = (phi*obst%c_drag/(pi/2.0_rsh))/cshelt
                     drag_x = 0.5_rsh*cd*rw*(o1dv_abdel_uvcent(s)**2.0_rsh)* &
                              o1dv_width_inst*lseg*COS(o1dv_abdel_t1cent(s))
                     ! Computing friction
                     re = (o1dv_abdel_uvcent(s)*o1dv_width_inst)/nu ! reynolds number
                     cf = (0.074_rsh*(re**(-1.0_rsh/5.0_rsh)))/cshelt
                     fric_x = cf*rw*(o1dv_abdel_uvcent(s)**2.0_rsh)* &
                              o1dv_width_inst*lseg* &
                              (SIN(o1dv_abdel_t1cent(s))**3.0_rsh)
                     fric_z = cf*rw*(o1dv_abdel_uvcent(s)**2.0_rsh)* &
                              o1dv_width_inst*lseg* &
                              (SIN(o1dv_abdel_t1cent(s))**2.0_rsh)*COS(o1dv_abdel_t1cent(s))
                     ! Computing lift
                     IF (phi .LE. phi_max) THEN
                        cl = (obst%c_lift*(phi/phi_max))/cshelt
                     ELSE
                        cl = 0.0_rsh
                     END IF
                     lift_z = 0.5*cl*rw*(o1dv_abdel_uvcent(s)**2.0_rsh)* &
                              o1dv_width_inst*lseg*SIN(o1dv_abdel_t1cent(s))
                     ! Computing buoyancy
                     buoy_z = (rw - obst%c_rho)*g*o1dv_width_inst* &
                              o1dv_thick_inst*lseg
                     ! Now summing all forces in x and z directions
                     IF (s .EQ. obst%c_abdel_nmax) THEN
                        o1dv_abdel_fx(s) = -(drag_x + fric_x)
                        o1dv_abdel_fz(s) = -(-lift_z + fric_z + buoy_z)
                     ELSE
                        o1dv_abdel_fx(s) = -(drag_x + fric_x + o1dv_abdel_fx(s + 1))
                        o1dv_abdel_fz(s) = -(-lift_z + fric_z + buoy_z + o1dv_abdel_fz(s + 1))
                     END IF
                  END IF
               END IF ! End test on enough velocity
            END DO ! End of downward loop on each leaf segment
            !---------------------------------------------------!
            ! *** TEST FOR CONVERGENCE OF OBSTRUCTION STATE *** !
            !---------------------------------------------------!
            ! Computing difference in bending angle between current iteration and previous one
            DO s = 1, obst%c_abdel_nmax
               o1dv_abdel_dtheta(s) = ABS(o1dv_abdel_t0cent(s) - o1dv_abdel_t1cent(s))*180.0_rsh/pi
            END DO
            ! Stop iterative procedure if maximum bendig difference is less than dtheta_max
            IF (MAXVAL(o1dv_abdel_dtheta(:)) .LE. dtheta_max) THEN
               niter_eff = niter
               niter = niter_max
            ELSE
               niter_eff = niter
               niter = niter + 1
            END IF
            !-------------------------------------------------------------------------!
            ! *** APPLY NEW THETA VALUE TO OLD ONE FOR NEXT ITERATION OR AVERAGED *** !
            ! ***     THETA BETWEEN LAST ITERATION AND PREVIOUS ONE IF STABLE     *** !
            ! ***   IF STABLE OBSTRUCTION STATE WAS NOT ACHIEVED (NO CONVERGENCE) *** !
            !-------------------------------------------------------------------------!
            IF (niter_eff .EQ. niter_max - 1) THEN
               o1dv_abdel_t0cent(:) = 0.5_rsh*(o1dv_abdel_t0cent(:) + o1dv_abdel_t1cent(:))
            ELSE
               o1dv_abdel_t0cent(:) = o1dv_abdel_t1cent(:)
            END IF
         END DO ! End of iterative loop for osbtruction state
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! ****                   FINALIZATION                 **** !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !------------------------------------------------------------------!
         ! *** WRITING INTO WARNING LOG IF CONVERGENCE WAS NOT ACHIEVED *** !
         !------------------------------------------------------------------!
         IF (niter_eff .EQ. niter_max - 1) THEN
            WRITE (o1dv_iwarnlog, *) ' '
            WRITE (o1dv_iwarnlog, *) ' '
            WRITE (o1dv_iwarnlog, *) '************************************************************************'
            WRITE (o1dv_iwarnlog, *) '**** module OBSTRUCTIONS, subroutine o1dv_COMP_ABDELPOSTURE ****'
            WRITE (o1dv_iwarnlog, *) ' WARNING : convergence for flow-obstruction coupling was not reached'
            WRITE (o1dv_iwarnlog, *) '           At obst', obst%name
            WRITE (o1dv_iwarnlog, *) '           Niter', niter_eff, 'Niter_max', niter_max
            WRITE (o1dv_iwarnlog, *) '           Dtheta', MAXVAL(o1dv_abdel_dtheta(:)), 'Dtheta_max', dtheta_max
            WRITE (o1dv_iwarnlog, *) '           This may significantly alter hydrodynamic results and'
            WRITE (o1dv_iwarnlog, *) '           produce canopy oscillations...'
            WRITE (o1dv_iwarnlog, *) '           To prevent this, bending angles obtain at the last'
            WRITE (o1dv_iwarnlog, *) '           iteration and the preivous one are averaged'
            WRITE (o1dv_iwarnlog, *) '************************************************************************'
         END IF
         !--------------------------------------!
         ! *** Updating obstructions height *** !
         !--------------------------------------!
         o1dv_height = 0.0_rsh
         DO s = 1, obst%c_abdel_nmax
            o1dv_height = o1dv_height + lseg*COS(o1dv_abdel_t0cent(s))
         END DO
         o1dv_height = MIN(MAX(o1dv_height, o1dv_p_hmin*o1dv_height_inst), &
                           0.99_rsh*o1dv_height_inst)
         !--------------------------------------------------------------!
         ! **** INTERPOLATION OF REAL BENDING ANGLES ON SIGMA GRID **** !
         !--------------------------------------------------------------!
         o1dv_theta3d(:) = 0.0_rsh
         s = 0
         o1dv_abdel_zn(s) = 0.0_rsh
         o1dv_abdel_tn(s) = 0.0_rsh
         DO s = 1, obst%c_abdel_nmax
            o1dv_abdel_zn(s) = o1dv_abdel_zcent(s)
            o1dv_abdel_tn(s) = o1dv_abdel_t0cent(s)
         END DO
         s = obst%c_abdel_nmax + 1
         o1dv_abdel_zn(s) = o1dv_height
         o1dv_abdel_tn(s) = o1dv_abdel_t0cent(s - 1)
         DO k = 1, o1dv_kmax
            z0 = o1dv_zc(k) - o1dv_dz(k)/2.0_rsh
            z1 = o1dv_zc(k) + o1dv_dz(k)/2.0_rsh
            IF (z1 .LE. o1dv_height) THEN
            DO s = 0, obst%c_abdel_nmax - 1
               zt0 = o1dv_abdel_zn(s)
               zt1 = o1dv_abdel_zn(s + 1)
               th0 = o1dv_abdel_tn(s)
               th1 = o1dv_abdel_tn(s + 1)
               IF ((o1dv_zc(k) .GE. zt0) .AND. (o1dv_zc(k) .LE. zt1)) THEN
                  o1dv_theta3d(k) = th0 + (o1dv_zc(k) - zt0)*(th1 - th0)/(zt1 - zt0)
               END IF
            END DO ! * END LOO ON S
            ELSEIF ((z0 .LE. o1dv_height) .AND. (z1 .GT. o1dv_height)) THEN
            o1dv_theta3d(k) = o1dv_abdel_t0cent(obst%c_abdel_nmax)
            END IF ! * END TEST ON Z0 AND Z1
         END DO ! * END LOOP ON K
         ! Additional check
         DO k = 2, o1dv_kmax - 1
            IF (o1dv_theta3d(k) .EQ. 0.0_rsh) THEN
            IF ((o1dv_theta3d(k - 1) .GT. 0.0_rsh) .AND. (o1dv_theta3d(k + 1) .GT. 0.0_rsh)) THEN
               o1dv_theta3d(k) = 0.5*(o1dv_theta3d(k - 1) + o1dv_theta3d(k + 1))
            END IF
            END IF
         END DO
      END IF ! * END TEST ON HWAT AND o1dv_HEIGHT

   END SUBROUTINE o1dv_abdelposture

   !!==========================================================================================================

   SUBROUTINE o1dv_comp_height(obst, o1dv_hwat, o1dv_uz, o1dv_vz, o1dv_dz, o1dv_zc, &
                               o1dv_height_inst, o1dv_height_p, o1dv_height)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_comp_height  ***
      !!
      !! ** Purpose : Computes obstructions height depending on choosen
      !!              parameterization
      !!
      !! ** Description : For rigid obstructions, o1dv_height is not changed
      !!                  For flexible obstructions, partial depth averaged velocity
      !!                  (from bottom or surface to o1dv_c_paramhuv*previous o1dv_height)
      !!                  corresponding to unconfined canopy is computed and used
      !!                  to compute o1dv_height. 
      !!
      !!---------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      TYPE(o1dv_param_type), INTENT(IN) :: obst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_hwat
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_uz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_vz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_dz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_zc
      REAL(KIND=rsh), INTENT(IN) :: o1dv_height_inst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_height_p
      REAL(KIND=rsh), INTENT(INOUT) :: o1dv_height

      !! * Local declaration
      INTEGER            :: k
      REAL(KIND=rsh)   :: htmp, dz, udz
      REAL(KIND=rsh)   :: uv

      !! * Executable part

      IF (obst%l_param_height) THEN ! * USE EXPONENTIAL PARAMETERIZATION

         !-------------------------------------------------
         ! * COMPUTATION OF PARTIAL DEPTH-AVERAGED VELOCITY
         !-------------------------------------------------
         IF (o1dv_height_p .EQ. 0.0_rsh) THEN
            htmp = o1dv_c_paramhuv*o1dv_height_inst ! at init
         ELSE
            htmp = o1dv_c_paramhuv*o1dv_height_p
         END IF
         dz = 0.0_rsh
         udz = 0.0_rsh
         IF (obst%type == "DO") THEN ! * DOWNWARD
            k = o1dv_kmax
            DO WHILE ((k .GE. 1) .AND. (o1dv_hwat - (o1dv_zc(k) - o1dv_dz(k)/2.0_rsh) .LT. htmp))
               dz = dz + o1dv_dz(k)
               udz = udz + SQRT(o1dv_uz(k)**2.0_rsh + o1dv_vz(k)**2.0_rsh)*o1dv_dz(k)
               k = k - 1
            END DO
         ELSE ! * UPWARD
            k = 1
            DO WHILE ((k .LE. o1dv_kmax) .AND. (o1dv_zc(k) + o1dv_dz(k)/2.0_rsh .LT. htmp))
               dz = dz + o1dv_dz(k)
               udz = udz + SQRT(o1dv_uz(k)**2.0_rsh + o1dv_vz(k)**2.0_rsh)*o1dv_dz(k)
               k = k + 1
            END DO
         END IF ! END TEST ON UPWARD/DOWNWARD
         uv = udz/dz
         !-----------------------------------------
         ! * COMPUTATION OF BENT OBSTRUCTION HEIGHT
         !-----------------------------------------
         htmp = obst%c_height_x0*o1dv_height_inst*EXP(obst%c_height_x1*uv)
         o1dv_height = MIN(MAX(htmp, o1dv_p_hmin*o1dv_height_inst), &
                           0.99_rsh*o1dv_height_inst)

      ELSE ! * USE ONLY FIRST COEFFICIENT
         !-------------------------------------------
         ! * COMPUTATION OF UNBENT OBSTRUCTION HEIGHT
         !-------------------------------------------
         o1dv_height = o1dv_height_inst*obst%c_height_x0

      END IF ! * END TEST ON PARAMETERIZATION

   END SUBROUTINE o1dv_comp_height

   !!==========================================================================================================

   SUBROUTINE o1dv_comp_theta(obst, o1dv_hwat, o1dv_dz, o1dv_zc, &
                              o1dv_height_inst, o1dv_height, o1dv_theta3d)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_comp_theta  ***
      !!
      !! ** Purpose : Computes obstructions bending angla depending on choosen
      !!              parameterization
      !!
      !! ** Description : For rigid obstructions, o1dv_theta is zero 
      !!                  For flexible obstructions, o1dv_height is used
      !!
      !!---------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      TYPE(o1dv_param_type), INTENT(IN) :: obst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_hwat
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_dz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_zc
      REAL(KIND=rsh), INTENT(IN) :: o1dv_height_inst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_height
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(INOUT) :: o1dv_theta3d

      !! * Local declaration
      INTEGER            :: k

      !! * Executable part

      !*********************************
      ! *** COMPUTATION OF BENDING ANGLE
      !*********************************
      IF (obst%type == "DO") THEN ! * DOWNWARD
         k = o1dv_kmax
         DO WHILE ((k .GE. 1) .AND. (o1dv_hwat - (o1dv_zc(k) + o1dv_dz(k)/2.0_rsh) .LT. o1dv_height))
            o1dv_theta3d(k) = ACOS(o1dv_height/o1dv_height_inst)
            k = k - 1
         END DO
      ELSE ! * UPWARD
         k = 1
         DO WHILE ((k .LE. o1dv_kmax) .AND. (o1dv_zc(k) - o1dv_dz(k)/2.0_rsh .GE. o1dv_height))
            o1dv_theta3d(k) = ACOS(o1dv_height/o1dv_height_inst)
            k = k + 1
         END DO
      END IF ! * END TEST ON UPWARD/DOWNWARD

   END SUBROUTINE o1dv_comp_theta

   !==========================================================================================================

   SUBROUTINE o1dv_comp_distrib(obst, o1dv_hwat, o1dv_dz, o1dv_zc, &
                                o1dv_dens_inst, o1dv_height, o1dv_dens3d)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_comp_distrib  ***
      !!
      !! ** Purpose : Computes vertical distribution of obstructions densities
      !!
      !!---------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      TYPE(o1dv_param_type), INTENT(IN) :: obst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_hwat
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_dz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_zc
      REAL(KIND=rsh), INTENT(IN) :: o1dv_dens_inst
      REAL(KIND=rsh), INTENT(INOUT) :: o1dv_height
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(INOUT) :: o1dv_dens3d

      !! * Local declaration
      INTEGER           :: k, kk
      REAL(KIND=rsh)  :: z0, z1, zc, zn0, zn1, dn0, dn1

      !! * Executable part
      IF (obst%l_filedistri) THEN ! * VARIABLE DISTRIBUTION
         !----------------------------------------
         ! *** COMPUTATION OF DENSITY DISTRIBUTION
         !----------------------------------------
         IF (obst%type == "DO") THEN ! * DOWNWARD
            DO k = o1dv_kmax, 1, -1
               z0 = MAX(o1dv_hwat - (o1dv_zc(k) - o1dv_dz(k)/2.0_rsh), 0.0_rsh)
               z1 = MIN(o1dv_hwat - (o1dv_zc(k) + o1dv_dz(k)/2.0_rsh), o1dv_hwat)
               IF ((z1 .LE. o1dv_height) .OR. ((z0 .LE. o1dv_height) .AND. (z1 .GT. o1dv_height))) THEN
                  IF (z1 .LE. o1dv_height) THEN
                     zc = o1dv_hwat - o1dv_zc(k)
                  ELSEIF ((z0 .LE. o1dv_height) .AND. (z1 .GT. o1dv_height)) THEN
                     zc = z0 + 0.5_rsh*(o1dv_height - z0)
                  END IF
                  DO kk = 1, obst%nbhnorm - 1
                     !--------------------------------------------------
                     ! *** CONVERSION OF NORMALIZED HEIGHT AND DENSITY
                     !--------------------------------------------------
                     zn0 = MIN(MAX(obst%height_norm(kk)*o1dv_height/100.0_rsh, 0.0_rsh), o1dv_height)
                     zn1 = MIN(MAX(obst%height_norm(kk + 1)*o1dv_height/100.0_rsh, 0.0_rsh), o1dv_height)
                     dn0 = MIN(MAX(obst%dens_norm(kk)*o1dv_dens_inst/100.0_rsh, 0.0_rsh), o1dv_dens_inst)
                     dn1 = MIN(MAX(obst%dens_norm(kk + 1)*o1dv_dens_inst/100.0_rsh, 0.0_rsh), o1dv_dens_inst)
                     IF ((zc .GE. zn0) .AND. (zc .LE. zn1)) THEN
                        o1dv_dens3d(k) = dn0 + (zc - zn0)*(dn1 - dn0)/(zn1 - zn0)
                     END IF
                  END DO
               END IF
            END DO
         ELSE ! * UPWARD
            DO k = 1, o1dv_kmax
               z0 = MAX(o1dv_zc(k) - o1dv_dz(k)/2.0_rsh, 0.0_rsh)
               z1 = MIN(o1dv_zc(k) + o1dv_dz(k)/2.0_rsh, o1dv_hwat)
               IF ((z1 .LE. o1dv_height) .OR. ((z0 .LE. o1dv_height) .AND. (z1 .GT. o1dv_height))) THEN
                  IF (z1 .LE. o1dv_height) THEN
                     zc = o1dv_zc(k)
                  ELSEIF ((z0 .LE. o1dv_height) .AND. (z1 .GT. o1dv_height)) THEN
                     zc = z0 + 0.5_rsh*(o1dv_height - z0)
                  END IF
                  DO kk = 1, obst%nbhnorm - 1
                     zn0 = MIN(MAX(obst%height_norm(kk)*o1dv_height/100.0_rsh, 0.0_rsh), o1dv_height)
                     zn1 = MIN(MAX(obst%height_norm(kk + 1)*o1dv_height/100.0_rsh, 0.0_rsh), o1dv_height)
                     dn0 = MIN(MAX(obst%dens_norm(kk)*o1dv_dens_inst/100.0_rsh, 0.0_rsh), o1dv_dens_inst)
                     dn1 = MIN(MAX(obst%dens_norm(kk + 1)*o1dv_dens_inst/100.0_rsh, 0.0_rsh), o1dv_dens_inst)
                     IF ((zc .GE. zn0) .AND. (zc .LE. zn1)) THEN
                        o1dv_dens3d(k) = dn0 + (zc - zn0)*(dn1 - dn0)/(zn1 - zn0)
                     END IF
                  END DO
               END IF
            END DO
         END IF ! * END TEST ON UPWARD/DOWNWARD
      ELSE ! * CONSTANT DISTRIBUTION
         IF (obst%type == "DO") THEN ! * DOWNWARD
            DO k = o1dv_kmax, 1, -1
               z0 = MAX(o1dv_hwat - (o1dv_zc(k) - o1dv_dz(k)/2.0_rsh), 0.0_rsh)
               z1 = MIN(o1dv_hwat - (o1dv_zc(k) + o1dv_dz(k)/2.0_rsh), o1dv_hwat)
               IF ((z1 .LE. o1dv_height) .OR. ((z0 .LE. o1dv_height) .AND. (z1 .GT. o1dv_height))) THEN
                  o1dv_dens3d(k) = o1dv_dens_inst
               END IF
            END DO
         ELSE ! * UPWARD
            DO k = 1, o1dv_kmax
               z0 = MAX(o1dv_zc(k) - o1dv_dz(k)/2.0_rsh, 0.0_rsh)
               z1 = MIN(o1dv_zc(k) + o1dv_dz(k)/2.0_rsh, o1dv_hwat)
               IF ((z1 .LE. o1dv_height) .OR. ((z0 .LE. o1dv_height) .AND. (z1 .GT. o1dv_height))) THEN
                  o1dv_dens3d(k) = o1dv_dens_inst
               END IF
            END DO
         END IF ! * END TEST ON UPWARD/DOWNWARD
      END IF ! * END TEST ON CONSTANT DISTRIBUTION

   END SUBROUTINE o1dv_comp_distrib

   !==========================================================================================================

   SUBROUTINE o1dv_comp_fracz(obst, o1dv_hwat, o1dv_dz, o1dv_zc, o1dv_height, &
                              o1dv_dens3d, o1dv_fracz3d)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_comp_fracz ***
      !!
      !! ** Purpose : Computes the fraction of sigma layer occupied by obstruction
      !!
      !!---------------------------------------------------------------------
      IMPLICIT NONE
      !! * Arguments
      TYPE(o1dv_param_type), INTENT(IN) :: obst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_hwat
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_dz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_zc
      REAL(KIND=rsh), INTENT(IN) :: o1dv_height
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_dens3d
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(INOUT) :: o1dv_fracz3d

      !! * Local declaration
      INTEGER           :: k
      REAL(KIND=rsh)  :: z0, z1

      !! * Executable part
      select case (obst%type)

      case ('DO') ! downward
         DO k = o1dv_kmax, 1, -1
            z0 = o1dv_hwat - (o1dv_zc(k) - o1dv_dz(k)/2.0_rsh)
            z1 = o1dv_hwat - (o1dv_zc(k) + o1dv_dz(k)/2.0_rsh)
            IF (z1 .LE. o1dv_height) THEN
               o1dv_fracz3d(k) = 1.0_rsh
            ELSEIF ((z0 .LE. o1dv_height) .AND. (z1 .GT. o1dv_height)) THEN
               o1dv_fracz3d(k) = (o1dv_height - z0)/(z1 - z0)
            END IF
         END DO

      case ('UP') ! upward
         DO k = 1, o1dv_kmax
            z0 = o1dv_zc(k) - o1dv_dz(k)/2.0_rsh
            z1 = o1dv_zc(k) + o1dv_dz(k)/2.0_rsh
            IF (z1 .LE. o1dv_height) THEN
               o1dv_fracz3d(k) = 1.0_rsh
            ELSEIF ((z0 .LE. o1dv_height) .AND. (z1 .GT. o1dv_height)) THEN
               o1dv_fracz3d(k) = (o1dv_height - z0)/(z1 - z0)
            END IF
         END DO

      case ('3D') ! 3D obst
         DO k = 1, o1dv_kmax
            z0 = o1dv_zc(k) - o1dv_dz(k)/2.0_rsh
            z1 = o1dv_zc(k) + o1dv_dz(k)/2.0_rsh
            IF (z1 .LE. o1dv_height) THEN
               IF (o1dv_dens3d(k) .EQ. 0.0_rsh) THEN
                  o1dv_fracz3d(k) = 0.0_rsh
               ELSE
                  o1dv_fracz3d(k) = 1.0_rsh
               END IF
            ELSEIF ((z0 .LE. o1dv_height) .AND. (z1 .GT. o1dv_height)) THEN
               o1dv_fracz3d(k) = (o1dv_height - z0)/(z1 - z0)
            END IF
         END DO

      case default
         WRITE (o1dv_iwarnlog, *) "Invalid obstruction type"

      end select

   END SUBROUTINE o1dv_comp_fracz

   !!==========================================================================================================

   SUBROUTINE o1dv_comp_diam(obst, o1dv_dz, &
                             o1dv_height_inst, o1dv_thick_inst, o1dv_width_inst, &
                             o1dv_dens3d, o1dv_theta3d, o1dv_fracz3d, o1dv_width3d, o1dv_thick3d)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_comp_diam  ***
      !!
      !! ** Purpose : Computes obstruction diameters (width and thickness)
      !!
      !! ** Description : Computes obstruction width and thickness according to
      !!                  bending angle (for flexible obstructions)
      !!
      !!---------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      TYPE(o1dv_param_type), INTENT(IN) :: obst
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_dz
      REAL(KIND=rsh), INTENT(IN) :: o1dv_height_inst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_thick_inst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_width_inst
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_dens3d
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_theta3d
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_fracz3d
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(INOUT) :: o1dv_width3d
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(INOUT) :: o1dv_thick3d

      !! * Local declaration
      INTEGER          :: k
      REAL(KIND=rsh) :: lsegl, ltotl

      !! * Executable part
      IF (obst%type == "DO") THEN ! * DOWNWARD
         ltotl = 0.0_rsh
         DO k = o1dv_kmax, 1, -1
            IF (o1dv_dens3d(k) .GT. 0.0_rsh) THEN
               o1dv_width3d(k) = o1dv_width_inst
               IF (obst%l_cylinder) THEN ! * CYLINDER
                  IF (o1dv_theta3d(k) .EQ. 0.0_rsh) THEN
                     o1dv_thick3d(k) = o1dv_width_inst
                  ELSE
                     lsegl = MIN(MAX(o1dv_fracz3d(k)*o1dv_dz(k)/COS(o1dv_theta3d(k)), &
                                     o1dv_fracz3d(k)*o1dv_dz(k)), o1dv_height_inst - ltotl)
                     o1dv_thick3d(k) = MIN(lsegl, MAX(o1dv_width_inst, lsegl*SIN(o1dv_theta3d(k))))
                     ltotl = ltotl + lsegl
                  END IF
               ELSE ! * PARALELEPIPED
                  IF (o1dv_theta3d(k) .EQ. 0.0_rsh) THEN
                     o1dv_thick3d(k) = o1dv_thick_inst
                  ELSE
                     lsegl = MIN(MAX(o1dv_fracz3d(k)*o1dv_dz(k)/COS(o1dv_theta3d(k)), &
                                     o1dv_fracz3d(k)*o1dv_dz(k)), o1dv_height_inst - ltotl)
                     o1dv_thick3d(k) = MIN(lsegl, MAX(o1dv_thick_inst, lsegl*SIN(o1dv_theta3d(k))))
                     ltotl = ltotl + lsegl
                  END IF
               END IF ! * END TEST CYLINDER
            END IF
         END DO
      ELSE ! * UPWARD
         ltotl = 0.0_rsh
         DO k = 1, o1dv_kmax
            IF (o1dv_dens3d(k) .GT. 0.0_rsh) THEN
               o1dv_width3d(k) = o1dv_width_inst
               IF (obst%l_cylinder) THEN ! * CYLINDER
                  IF (o1dv_theta3d(k) .EQ. 0.0_rsh) THEN
                     o1dv_thick3d(k) = o1dv_width_inst
                  ELSE
                     lsegl = MIN(MAX(o1dv_fracz3d(k)*o1dv_dz(k)/COS(o1dv_theta3d(k)), &
                                     o1dv_fracz3d(k)*o1dv_dz(k)), o1dv_height_inst - ltotl)
                     o1dv_thick3d(k) = MIN(lsegl, MAX(o1dv_width_inst, lsegl*SIN(o1dv_theta3d(k))))
                     ltotl = ltotl + lsegl
                  END IF
               ELSE ! * PARALELEPIPED
                  IF (o1dv_theta3d(k) .EQ. 0.0_rsh) THEN
                     o1dv_thick3d(k) = o1dv_thick_inst
                  ELSE
                     lsegl = MIN(MAX(o1dv_fracz3d(k)*o1dv_dz(k)/COS(o1dv_theta3d(k)), &
                                     o1dv_fracz3d(k)*o1dv_dz(k)), o1dv_height_inst - ltotl)
                     o1dv_thick3d(k) = MIN(lsegl, MAX(o1dv_thick_inst, lsegl*SIN(o1dv_theta3d(k))))
                     ltotl = ltotl + lsegl
                  END IF
               END IF ! * END TEST CYLINDER
            END IF
         END DO
      END IF ! * END TEST ON UPWARD/DOWNWARD

   END SUBROUTINE o1dv_comp_diam

   !!==========================================================================================================

   SUBROUTINE o1dv_comp_fracxy(obst, o1dv_position, o1dv_fracxy)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_comp_fracxy  ***
      !!
      !! ** Purpose : Computes obstructions correction term for coverage (fragmentation)
      !!              within one single grid cell
      !!---------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      TYPE(o1dv_param_type), INTENT(IN) :: obst
      REAL(KIND=rsh), INTENT(IN) :: o1dv_position
      REAL(KIND=rsh), INTENT(INOUT) :: o1dv_fracxy

      !! * Local declaration
      REAL(KIND=rsh)     :: kv

      !! * Executable part
      IF (obst%l_fracxy) THEN
         IF (obst%fracxy_type .EQ. 0) THEN
            !------------------
            ! Linear correction
            !------------------
            o1dv_fracxy = o1dv_position
         ELSEIF (obst%fracxy_type .EQ. 1) THEN
            !------------------------------
            ! Simple exponential correction
            !------------------------------
            o1dv_fracxy = (EXP(o1dv_position*obst%c_fracxy_k0) - 1.0_rsh)/ &
                          (EXP(obst%c_fracxy_k0) - 1.0_rsh)
         ELSEIF (obst%fracxy_type .EQ. 2) THEN
            !-------------------------------
            ! Complex exponential correction
            !-------------------------------
            kv = obst%c_fracxy_k0 + obst%c_fracxy_k1*EXP(-(1.0_rsh - o1dv_position)*obst%c_fracxy_l)
            o1dv_fracxy = (EXP(o1dv_position*kv) - 1.0_rsh)/ &
                          (EXP(kv) - 1.0_rsh)
         ELSEIF (obst%fracxy_type .EQ. 3) THEN
            !---------------------------------------------------------
            ! Constant correction factor (for parameterization purpose)
            !---------------------------------------------------------
            o1dv_fracxy = obst%c_fracxy_k0*o1dv_position
         END IF ! * END TEST ON CORRECTION TYPE
      ELSE
         !--------------
         ! No correction
         !--------------
         o1dv_fracxy = 1.0_rsh
      END IF ! * END TEST ON CORRECTION TYPE

   END SUBROUTINE o1dv_comp_fracxy

   !!==========================================================================================================

   SUBROUTINE o1dv_comp_projarea(o1dv_position, o1dv_dz, o1dv_output)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_comp_projarea  ***
      !!
      !! ** Purpose : Computes obstructions projected horizontal and vertical area
      !!---------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_position
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN)  :: o1dv_dz
      TYPE(o1dv_out_type), INTENT(INOUT) :: o1dv_output

      !! * Local declaration
      INTEGER                    :: iv, k
      REAL(KIND=rsh)           :: dzas, spos, wtmp, dtmp, sfz
      REAL(KIND=rsh), PARAMETER :: gamma = 1.0_rsh - 1.0e-6_rsh

      !! * Executable part
      !*************************

      !----------------------------------
      ! 3D VARIABLES !
      !----------------------------------
      DO k = 1, o1dv_kmax
         !----------------------
         ! * LOOP ON o1dv_NB_VAR
         !----------------------
         DO iv = 1, o1dv_nbvar
            IF ((o1dv_position(iv) .GT. 0.0_rsh) .AND. (o1dv_output%dens3d(iv, k) .GT. 0.0_rsh)) THEN
               !------------------------------------------
               ! * Horizontal surface area of obstructions
               !------------------------------------------
               IF (o1dv_obst_param(iv)%l_cylinder) THEN ! Cylindric/Ellipse obstruction
                  o1dv_output%a3d(iv, k) = MIN(pi*(0.5_rsh*o1dv_output%width3d(iv, k)*0.5_rsh*o1dv_output%thick3d(iv, k))* &
                                               o1dv_output%dens3d(iv, k)*o1dv_output%fracxy(iv), gamma)
               ELSE
                  o1dv_output%a3d(iv, k) = MIN(o1dv_output%width3d(iv, k)*o1dv_output%thick3d(iv, k)* &
                                               o1dv_output%dens3d(iv, k)*o1dv_output%fracxy(iv), gamma)
               END IF
               !----------------------------------------
               ! * Vertical surface area of obstructions
               !----------------------------------------
               o1dv_output%s3d(iv, k) = MIN(o1dv_output%fracz3d(iv, k)*o1dv_output%width3d(iv, k)* &
                                            SQRT(o1dv_output%dens3d(iv, k))*o1dv_output%fracxy(iv), gamma)
            END IF ! * END TEST ON POSITION AND DENSITY
         END DO ! * END LOOP ON o1dv_NBVAR
         !------------------------------------------------------------------------
         ! * Computes o1dv_output%a3d and o1dv_output%s3d for NoTurb/Turb and Total obstructions
         !------------------------------------------------------------------------
         DO iv = 1, o1dv_nbvar
            IF (o1dv_obst_param(iv)%l_noturb) THEN
               !-------------------------------------------
               ! First for NoTurb (simplified obstructions)
               !-------------------------------------------
               o1dv_output%a3d(o1dv_nbvar + 1, k) = MIN(o1dv_output%a3d(o1dv_nbvar + 1, k) + o1dv_output%a3d(iv, k), gamma)
               o1dv_output%s3d(o1dv_nbvar + 1, k) = MIN(o1dv_output%s3d(o1dv_nbvar + 1, k) + o1dv_output%s3d(iv, k), gamma)
            ELSE
               !---------------------------------------
               ! Then for Turb (turbulent obstructions)
               !---------------------------------------
               o1dv_output%a3d(o1dv_nbvar + 2, k) = MIN(o1dv_output%a3d(o1dv_nbvar + 2, k) + o1dv_output%a3d(iv, k), gamma)
               o1dv_output%s3d(o1dv_nbvar + 2, k) = MIN(o1dv_output%s3d(o1dv_nbvar + 2, k) + o1dv_output%s3d(iv, k), gamma)
            END IF
            !--------------------------------
            ! Then for Tot (all obstructions)
            !--------------------------------
            o1dv_output%a3d(o1dv_nbvar + 3, k) = MIN(o1dv_output%a3d(o1dv_nbvar + 3, k) + o1dv_output%a3d(iv, k), gamma)
            o1dv_output%s3d(o1dv_nbvar + 3, k) = MIN(o1dv_output%s3d(o1dv_nbvar + 3, k) + o1dv_output%s3d(iv, k), gamma)
         END DO ! * END LOOP ON o1dv_NBVAR
      END DO ! * END LOOP ON K

      !----------------------------------
      ! Now computes depth-averaged value
      !----------------------------------
      DO iv = 1, o1dv_nbvar
         IF (o1dv_position(iv) .GT. 0.0_rsh) THEN
            dzas = 0.0_rsh
            DO k = 1, o1dv_kmax
               !-----------
               ! For o1dv_a
               !-----------
               o1dv_output%a2d(iv) = o1dv_output%a2d(iv) + (o1dv_output%a3d(iv, k)*o1dv_dz(k))
               !-----------
               ! For o1dv_s
               !-----------
               o1dv_output%s2d(iv) = o1dv_output%s2d(iv) + (o1dv_output%s3d(iv, k)*o1dv_dz(k))
               !------------
               dzas = dzas + o1dv_dz(k)
            END DO ! * END LOOP ON K
            IF (dzas .GT. 0.0_rsh) THEN
               o1dv_output%a2d(iv) = MIN(o1dv_output%a2d(iv)/dzas, gamma)
               o1dv_output%s2d(iv) = MIN(o1dv_output%s2d(iv)/dzas, gamma)
            END IF
         END IF ! * END TEST ON POSITION
      END DO ! * END LOOP ON IV
      !------------------------------
      ! Then For NoTurb, Turb and Tot
      !------------------------------
      DO iv = 1, o1dv_nbvar
         IF (o1dv_position(iv) .GT. 0.0_rsh) THEN
            IF (o1dv_obst_param(iv)%l_noturb) THEN
               !-----------
               ! For NoTurb
               !-----------
               o1dv_output%a2d(o1dv_nbvar + 1) = o1dv_output%a2d(o1dv_nbvar + 1) + o1dv_output%a2d(iv)
               o1dv_output%s2d(o1dv_nbvar + 1) = o1dv_output%s2d(o1dv_nbvar + 1) + o1dv_output%s2d(iv)
            ELSE
               !---------
               ! For Turb
               !---------
               o1dv_output%a2d(o1dv_nbvar + 2) = o1dv_output%a2d(o1dv_nbvar + 2) + o1dv_output%a2d(iv)
               o1dv_output%s2d(o1dv_nbvar + 2) = o1dv_output%s2d(o1dv_nbvar + 2) + o1dv_output%s2d(iv)
            END IF ! * END TEST ON NOTURB
            !----------
            ! For Total
            !----------
            o1dv_output%a2d(o1dv_nbvar + 3) = o1dv_output%a2d(o1dv_nbvar + 3) + o1dv_output%a2d(iv)
            o1dv_output%s2d(o1dv_nbvar + 3) = o1dv_output%s2d(o1dv_nbvar + 3) + o1dv_output%s2d(iv)
         END IF ! * END TEST ON POSITION
      END DO ! * END LOOP ON IV

      o1dv_output%a2d(o1dv_nbvar + 1) = MIN(o1dv_output%a2d(o1dv_nbvar + 1), gamma)
      o1dv_output%a2d(o1dv_nbvar + 2) = MIN(o1dv_output%a2d(o1dv_nbvar + 2), gamma)
      o1dv_output%a2d(o1dv_nbvar + 3) = MIN(o1dv_output%a2d(o1dv_nbvar + 3), gamma)
      o1dv_output%s2d(o1dv_nbvar + 1) = MIN(o1dv_output%s2d(o1dv_nbvar + 1), gamma)
      o1dv_output%s2d(o1dv_nbvar + 2) = MIN(o1dv_output%s2d(o1dv_nbvar + 2), gamma)
      o1dv_output%s2d(o1dv_nbvar + 3) = MIN(o1dv_output%s2d(o1dv_nbvar + 3), gamma)

      !---------------------------
      ! * Averaged characteristics
      !---------------------------
      spos = 0.0_rsh
      DO iv = 1, o1dv_nbvar
         o1dv_output%height_mean = o1dv_output%height_mean + o1dv_position(iv)*o1dv_output%height(iv)
         wtmp = 0.0_rsh
         dtmp = 0.0_rsh
         sfz = 0.0_rsh
         DO k = 1, o1dv_kmax
            wtmp = wtmp + o1dv_output%fracz3d(iv, k)*o1dv_output%width3d(iv, k)
            dtmp = dtmp + o1dv_output%fracz3d(iv, k)*o1dv_output%dens3d(iv, k)
            sfz = sfz + o1dv_output%fracz3d(iv, k)
         END DO
         IF (sfz .GT. 0.0_rsh) THEN
            wtmp = wtmp/sfz
            dtmp = dtmp/sfz
         END IF
         o1dv_output%dens_mean = o1dv_output%dens_mean + o1dv_position(iv)*dtmp
         o1dv_output%width_mean = o1dv_output%width_mean + o1dv_position(iv)*wtmp
         spos = spos + o1dv_position(iv)
      END DO

      IF (spos .GT. 0.0_rsh) THEN
         o1dv_output%dens_mean = o1dv_output%dens_mean/spos
         o1dv_output%width_mean = o1dv_output%width_mean/spos
         o1dv_output%height_mean = o1dv_output%height_mean/spos
      END IF

      !-------------------------------------
   END SUBROUTINE o1dv_comp_projarea

   !!==========================================================================================================

   SUBROUTINE o1dv_comp_obstroughness(o1dv_position, o1dv_z0bed, &
                                      o1dv_height_f, o1dv_dens_f, o1dv_width_f, o1dv_output)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE o1dv_comp_obstroughness  ***
      !!
      !! ** Purpose : Computes obstructions bottom roughness
      !!---------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_position
      REAL(KIND=rsh), INTENT(IN) ::  o1dv_z0bed
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_height_f
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_dens_f
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_width_f
      TYPE(o1dv_out_type), INTENT(INOUT) :: o1dv_output

      !! * Local declaration
      INTEGER         :: iv
      REAL(KIND=rsh)  :: a, d, z0a, coef, z0tot, ctot
      !!----------------------------------------------------------------------
      !! * Executable part
      !-----------------------------------------

      o1dv_output%z0obst(:) = o1dv_z0bed

      ! ******************************************************* !
      ! *** Computation for all upward and not 3d variables *** !
      ! ******************************************************* !
      DO iv = 1, o1dv_nbvar
         !--------------------------------
         ! *** Abdelhrman parameterization
         !--------------------------------
         IF (o1dv_obst_param(iv)%type == "UP") THEN ! * UPWARD (not downward and not 3D)
            IF (o1dv_position(iv) .GT. 0.0_rsh) THEN
               IF (o1dv_obst_param(iv)%l_abdelrough_cste) THEN
                  coef = o1dv_obst_param(iv)%c_crough_x0
               ELSE
                  coef = o1dv_obst_param(iv)%c_crough_x1 + o1dv_obst_param(iv)%c_crough_x0*(o1dv_height_f(iv)**2.0_rsh)/ &
                         ((o1dv_output%height(iv)**2.0_rsh)*o1dv_width_f(iv)*o1dv_dens_f(iv))
               END IF
               a = o1dv_width_f(iv)*o1dv_height_f(iv)*(o1dv_output%height(iv)/o1dv_height_f(iv))
               d = (coef*o1dv_output%height(iv)*o1dv_output%height(iv)*o1dv_width_f(iv))/ &
                   (a + coef*o1dv_width_f(iv)*o1dv_output%height(iv))
               z0a = (0.5_rsh*o1dv_width_f(iv)*o1dv_output%height(iv)*o1dv_output%height(iv)*a)/ &
                     (a + coef*o1dv_width_f(iv)*o1dv_output%height(iv))**2.0_rsh
               z0a = MAX(z0a + d, o1dv_z0bed)
               ! Now apply fraction
               o1dv_output%z0obst(iv) = (z0a*o1dv_position(iv)) + ((1.0_rsh - o1dv_position(iv))*o1dv_z0bed)
            END IF
         END IF
      END DO ! * END LOOP o1dv_NBVAR
      ! ************************************************* !
      ! *** Total computation for NoTurb Formulations *** !
      ! ************************************************* !
      ! As a weighted-average
      !----------------------
      z0tot = 0.0_rsh
      ctot = 0.0_rsh
      DO iv = 1, o1dv_nbvar
         IF (o1dv_obst_param(iv)%l_noturb) THEN
            z0tot = z0tot + o1dv_output%z0obst(iv)
            ctot = ctot + 1.0_rsh
         END IF
      END DO
      IF (z0tot .EQ. 0.0_rsh) THEN
         o1dv_output%z0obst(o1dv_nbvar + 1) = o1dv_z0bed
      ELSE
         o1dv_output%z0obst(o1dv_nbvar + 1) = MAX(z0tot/ctot, o1dv_z0bed)
      END IF
      ! *********************************************** !
      ! *** Total computation for Turb Formulations *** !
      ! *********************************************** !
      ! As a weighted-average
      !----------------------
      z0tot = 0.0_rsh
      ctot = 0.0_rsh
      DO iv = 1, o1dv_nbvar
         IF (.NOT. o1dv_obst_param(iv)%l_noturb) THEN
            z0tot = z0tot + o1dv_output%z0obst(iv)
            ctot = ctot + 1.0_rsh
         END IF
      END DO
      IF (z0tot .EQ. 0.0_rsh) THEN
         o1dv_output%z0obst(o1dv_nbvar + 2) = o1dv_z0bed
      ELSE
         o1dv_output%z0obst(o1dv_nbvar + 2) = MAX(z0tot/ctot, o1dv_z0bed)
      END IF
      ! ************************* !
      ! *** Total computation *** !
      ! ************************* !
      ! As a weighted-average
      !----------------------
      z0tot = 0.0_rsh
      ctot = 0.0_rsh
      DO iv = 1, o1dv_nbvar
         z0tot = z0tot + o1dv_output%z0obst(iv)
         ctot = ctot + 1.0_rsh
      END DO
      IF (z0tot .EQ. 0.0_rsh) THEN
         o1dv_output%z0obst(o1dv_nbvar + 3) = o1dv_z0bed
      ELSE
         o1dv_output%z0obst(o1dv_nbvar + 3) = MAX(z0tot/ctot, o1dv_z0bed)
      END IF

      !-----------------------------------------
   END SUBROUTINE o1dv_comp_obstroughness

   !!==========================================================================================================

   SUBROUTINE o1dv_comp_hydroparam(o1dv_position, o1dv_cmu, o1dv_uz, o1dv_vz, o1dv_output)
   !!---------------------------------------------------------------------
   !!                 *** SUBROUTINE o1dv_comp_hydroparam  ***
   !!
   !! ** Purpose : Computes obstructions parameters used for hydrodynamics
   !!---------------------------------------------------------------------

      IMPLICIT NONE

   !! * Arguments
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: o1dv_position
      REAL(KIND=rsh), INTENT(IN) ::  o1dv_cmu
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_uz
      REAL(KIND=rsh), DIMENSION(o1dv_kmax), INTENT(IN) :: o1dv_vz
      TYPE(o1dv_out_type), INTENT(INOUT) :: o1dv_output

   !! * Local declaration
      INTEGER                     :: iv, k
      REAL(KIND=rsh)            :: phi, ntot, clz, nclz, lz
      REAL(KIND=rsh)            :: uzvz, fuv, tuz, tvz
      REAL(KIND=rsh), PARAMETER  :: gamma = 1.0E-5
   !!----------------------------------------------------------------------
   !! * Executable part

      ! *************************************************** !
      ! *************** TURBULENT VARIABLES *************** !
      ! *************************************************** !
      DO k = 1, o1dv_kmax
         uzvz = SQRT(o1dv_uz(k)**2.0_rsh + o1dv_vz(k)**2.0_rsh)
         IF (uzvz .GT. 0.0001_rsh) THEN
            ! *************************************************** !
            ! *********** TURBULENT VARIABLES 1ST PART ********** !
            ! *************************************************** !
            ntot = 0.0_rsh
            tuz = 0.0_rsh
            tvz = 0.0_rsh
            clz = 0.0_rsh
            nclz = 0.0_rsh
            DO iv = 1, o1dv_nbvar
            IF ((o1dv_position(iv) .GT. 0.0_rsh) .AND. (.NOT. o1dv_obst_param(iv)%l_noturb)) THEN
               IF (o1dv_output%dens3d(iv, k) .GT. 0.0_rsh) THEN
                  !----------------------------
                  ! * Total obstruction density
                  !----------------------------
                  ntot = ntot + o1dv_output%dens3d(iv, k)*o1dv_output%fracxy(iv)
                  !-------------------
                  ! * Drag coefficient
                  !-------------------
                  IF (o1dv_obst_param(iv)%l_drag_cste) THEN
                     o1dv_output%drag3d(iv, k) = o1dv_obst_param(iv)%c_drag
                  ELSE
                     phi = (pi/2.0_rsh) - o1dv_output%theta3d(iv, k)
                     o1dv_output%drag3d(iv, k) = gamma + ABS(phi)*(o1dv_obst_param(iv)%c_drag - gamma)/(pi/2.0_rsh)
                  END IF
                  !-------------------
                  ! * Resistance force
                  !-------------------
                  ! Here, fuv has unit:
                  ! [fuv] = - * - * m * m.s-1 * m-2 * - * -
                  ! [fuv] = s-1
                  ! Here o1dv_fuz_e (o1dv_fvz_e) has unit:
                  ! [o1dv_fuz_e]     = [fuv] * [o1dv_uz]
                  ! [o1dv_fuz_e]     = s-1   * m.s-1
                  ! [o1dv_fuz_e]     = m.s-2
                  ! [o1dv_fuz_e*rho] = kg.m-3 * m.s-2
                  ! [o1dv_fuz_e*rho] = N.m-3
                  !-------------------
                  fuv = 0.5_rsh*o1dv_output%drag3d(iv, k)*o1dv_output%width3d(iv, k)*uzvz* &
                        o1dv_output%dens3d(iv, k)*o1dv_output%fracxy(iv)*o1dv_output%fracz3d(iv, k)
                  o1dv_output%fuz(k) = o1dv_output%fuz(k) + (fuv*o1dv_uz(k))
                  o1dv_output%fvz(k) = o1dv_output%fvz(k) + (fuv*o1dv_vz(k))
                  !--------------------------
                  ! * Work spent by the fluid
                  !--------------------------
                  ! Here, tuz has unit:
                  ! [tuz] = s-1 * m.s-1 * m.s-1
                  ! [tuz] = m2.s-3
                  tuz = tuz + fuv*(o1dv_uz(k)**2.0_rsh)
                  tvz = tvz + fuv*(o1dv_vz(k)**2.0_rsh)
                  !-------------------------------
                  ! * Averaging coefficient for lz
                  !-------------------------------
                  clz = clz + o1dv_obst_param(iv)%c_lz
                  nclz = nclz + 1.0_rsh
               END IF ! * END TEST ON OBSTRUCTION DENSITY
            END IF ! * END TEST ON POSITION AND TURB VARIABLE
            END DO ! * END LOOP ON NBVAR
            ! *************************************************** !
            ! *********** TURBULENT VARIABLES 2ND PART ********** !
            ! *************************************************** !
            IF (ntot .GT. 0.0_rsh) THEN
               !--------------------------
               ! * Work spent by the fluid
               !--------------------------
               ! Here, o1dv_t has unit:
               ! [o1dv_t] = m2.s-3
               o1dv_output%t(k) = SQRT(tuz**2.0_rsh + tvz**2.0_rsh)
               !-------------------------------------
               ! * Smallest distance between elements
               !-------------------------------------
               clz = clz/nclz
               lz = clz*SQRT((1.0_rsh - o1dv_output%a3d(o1dv_nbvar + 2, k))/ntot)
               !-----------------------------------------------------
               ! * Dissipation timescale of eddies between structures
               ! ----------------------------------------------------
               IF ((o1dv_output%t(k) /= 0.0_rsh) .AND. (o1dv_c2turb /= 0.0_rsh) .AND. (o1dv_cmu /= 0.0_rsh)) THEN
                  o1dv_output%tau(k) = 1.0_rsh/(1.0_rsh/(o1dv_c2turb*sqrt(o1dv_cmu))* &
                                                ((lz*lz)/(o1dv_output%t(k)))**(1.0_rsh/3.0_rsh))
               END IF
            END IF ! * END TEST ON OBSTRUCTION DENSITY
         END IF ! * END TEST ON ENOUGH VELOCITY
      END DO ! * END LOOP ON k

      !-------------------------------------
   END SUBROUTINE o1dv_comp_hydroparam

   !!==========================================================================================================
   FUNCTION o1dv_comp_z0sedim(position, height, width, dens, z0sedim) result(o1dv_z0sedim)
   !!---------------------------------------------------------------------
   !!                 *** FUNCTION o1dv_comp_z0sedim ***
   !!
   !! ** Purpose : Computes z0sedim in presence of obstructions
   !!
   !!---------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: position
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: height
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: width
      REAL(KIND=rsh), DIMENSION(o1dv_nbvar), INTENT(IN) :: dens
      REAL(KIND=rsh), INTENT(IN) :: z0sedim
      REAL(KIND=rsh) :: o1dv_z0sedim
      !! * Local declaration
      INTEGER :: iv
      REAL(KIND=rsh) :: z0tmp, stmp, z00, oah
      REAL(KIND=rsh), PARAMETER :: epsi = 1E-6
      !!----------------------------------------------------------------------
      !! * Executable part
      !!------------------
      z0tmp = 0.0_rsh
      stmp = 0.0_rsh
      DO iv = 1, o1dv_nbvar
         IF ((o1dv_obst_param(iv)%type == "UP") .AND. (o1dv_obst_param(iv)%l_z0bstress)) THEN
            IF (position(iv) .GT. 0.0_rsh) THEN
               IF (o1dv_obst_param(iv)%z0bstress_option == 0) THEN
                  !-----------------------
                  ! Constant value is used
                  !-----------------------
                  z00 = o1dv_obst_param(iv)%c_z0bstress
               ELSE
                  !-------------------------
                  ! Parameterization is used
                  !-------------------------
                  oah = dens(iv)*width(iv)*height(iv)
                  z00 = o1dv_obst_param(iv)%c_z0bstress_x0*oah**o1dv_obst_param(iv)%c_z0bstress_x1
                  z00 = MAX(MIN(z00, 0.01_rsh), epsi)
               END IF ! END test on parameterization
               z0tmp = z0tmp + (position(iv)*z00) + ((1.0_rsh - position(iv))*z0sedim)
               stmp = stmp + 1.0_rsh
            ELSE
               z0tmp = z0tmp + z0sedim
               stmp = stmp + 1.0_rsh
            END IF
         ELSE
            z0tmp = z0tmp + z0sedim
            stmp = stmp + 1.0_rsh
         END IF
      END DO ! END LOOP obst_nbvar
      !----------
      ! Averaging
      !----------
      o1dv_z0sedim = MAX(epsi, z0tmp/stmp)

   !!**********************************************************************
   END FUNCTION o1dv_comp_z0sedim

   !==================================================================================================
END MODULE OBSTRUCTIONS1DV
