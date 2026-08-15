#include "cppdefs.h"

MODULE oce_trc

#if defined key_pisces

   use scalars
   use ncscrum
   USE par_pisces
   USE in_out_manager
   USE lib_mpp

   IMPLICIT NONE
   PUBLIC

#include "grid.h"
#include "ocean3d.h"
#include "forces.h"
#include "mixing.h"
#include "diagnostics.h"

   !! * Substitutions
#  include "ocean2pisces.h90"
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"


   PUBLIC   trc_oce_rgb        ! routine called by traqsr.F90
   PUBLIC   trc_oce_rgb_read   ! routine called by traqsr.F90
   PUBLIC   trc_oce_ext_lev    ! function called by traqsr.F90 at least
   PUBLIC   trc_oce_alloc      ! function called by nemogcm.F90
   PUBLIC   ocean_2_pisces 

   PUBLIC   glob_sum
   PUBLIC   tracer_stat
# ifdef NEMO
   PUBLIC   fld_read, fld_fill
# endif
   INTERFACE glob_sum
      MODULE PROCEDURE  glob_sum_2d, glob_sum_3d
   END INTERFACE   

   LOGICAL , PUBLIC ::   l_co2cpl  = .false.   !: atmospheric pco2 recieved from oasis
   LOGICAL , PUBLIC ::   l_offline = .false.   !: offline passive tracers flag
   REAL(wp), PUBLIC ::   r_si2                 !: largest depth of extinction (blue & 0.01 mg.m-3)  (RGB)
   LOGICAL , PUBLIC ::   ln_trcdc2dm           !: Diurnal cycle for TOP
   INTEGER, PUBLIC  ::   nksr    !: =nkV, i.e. maximum level of light extinction
   !
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) ::   tra   !: traceur concentration for next time step
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:)   ::   tmask
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:)   ::   etot3     !: light absortion coefficient
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)     ::   oce_co2   !: ocean carbon flux
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)     ::   atm_co2   !: ocean carbon flux
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)     ::   fr_i      !: ocean carbon flux
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)     ::   e1e2t
   REAL(wp) , PUBLIC, DIMENSION(3,61)   ::   rkrgb    ! tabulated attenuation coefficients for RGB absorption

# if defined NEMO
   TYPE, PUBLIC ::   FLD_N      !: Namelist field informations
      CHARACTER(len = 256) ::   clname      ! generic name of the NetCDF flux file
      REAL(wp)             ::   freqh       ! frequency of each flux file
      CHARACTER(len = 34)  ::   clvar       ! generic name of the variable in the NetCDF flux file
      LOGICAL              ::   ln_tint     ! time interpolation or not (T/F)
      LOGICAL              ::   ln_clim     ! climatology or not (T/F)
      CHARACTER(len = 8)   ::   clftyp      ! type of data file 'daily', 'monthly' or yearly'
      CHARACTER(len = 256) ::   wname       ! generic name of a NetCDF weights file to be used, blank if not
      CHARACTER(len = 34)  ::   vcomp       ! symbolic component name if a vector that needs rotation
      !                                     ! a string starting with "U" or "V" for each component
      !                                     ! chars 2 onwards identify which components go together
      CHARACTER(len = 34)  ::   lname       ! generic name of a NetCDF land/sea mask file to be used, blank if not
      !                                     ! 0=sea 1=land
   END TYPE FLD_N

 TYPE, PUBLIC ::   FLD        !: Input field related variables
      CHARACTER(len = 256)            ::   clrootname   ! generic name of the NetCDF file
      CHARACTER(len = 256)            ::   clname       ! current name of the NetCDF file
      REAL(wp)                        ::   freqh        ! frequency of each flux file
      CHARACTER(len = 34)             ::   clvar        ! generic name of the variable in the NetCDF flux file
      LOGICAL                         ::   ln_tint      ! time interpolation or not (T/F)
      REAL(wp)                        ::   rec_shft     ! record shift (from -0.5 to +0.5) when ln_tint=T
      LOGICAL                         ::   ln_clim      ! climatology or not (T/F)
      CHARACTER(len = 8)              ::   clftyp       ! type of data file 'daily', 'monthly' or yearly'
      CHARACTER(len = 1)              ::   cltype       ! nature of grid-points: T, U, V...
      REAL(wp)                        ::   zsgn         ! -1. the sign change across the north fold, =  1. otherwise
      INTEGER                         ::   num          ! iom id of the jpfld files to be read
      INTEGER , DIMENSION(2,2)        ::   nrec         ! before/after record (1: index, 2: second since Jan. 1st 00h of yr nit000)
      INTEGER                         ::   nbb          ! index of before values
      INTEGER                         ::   naa          ! index of after  values
      INTEGER , ALLOCATABLE, DIMENSION(:) ::   nrecsec  ! records time in seconds
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::   fnow   ! input fields interpolated to now time step
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   fdta   ! 2 consecutive record of input fields
      CHARACTER(len = 256)            ::   wgtname      ! current name of the NetCDF weight file acting as a key
      !                                                 ! into the WGTLIST structure
      CHARACTER(len = 34)             ::   vcomp        ! symbolic name for a vector component that needs rotation
      LOGICAL, DIMENSION(2)           ::   rotn         ! flag to indicate whether before/after field has been rotated
      INTEGER                         ::   nreclast     ! last record to be read in the current file
      CHARACTER(len = 256)            ::   lsmname      ! current name of the NetCDF mask file acting as a key
      !                                                 !
      !                                                 ! Variables related to BDY
      INTEGER                         ::   igrd         !   grid type for bdy data
      INTEGER                         ::   ibdy         !   bdy set id number
      INTEGER, POINTER, DIMENSION(:)  ::   imap         !   Array of integer pointers to 1D arrays
      LOGICAL                         ::   ltotvel      !   total velocity or not (T/F)
      LOGICAL                         ::   lzint        !   T if it requires a vertical interpolation
   END TYPE FLD

   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trsfac    ! multiplicative factor for SBC tracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trcsbc    ! structure of data input SBC (file informations, fields read)
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trcfac    ! multiplicative factor for CBC tracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trccbc    ! structure of data input CBC (file informations, fields read)
# endif

  REAL, DIMENSION(GLOBAL_2D_ARRAY,N+1,3) :: gdepw   ! W-depht
  REAL, DIMENSION(GLOBAL_2D_ARRAY,N+1,3) :: e3w     ! W-vertical scale factor

  INTEGER, PUBLIC :: Istrp,Iendp,Jstrp,Jendp
!$OMP threadprivate(Istrp,Iendp)
!$OMP threadprivate(Jstrp,Jendp)

  LOGICAL :: ln_ctl     = .false.
  LOGICAL :: ln_qsr_bio = .false.
  REAL(wp)    :: rn_si0     = 0.35
!  INTEGER :: numout
  INTEGER :: jpdom_data = 1
  INTEGER :: jpdom_global = 1
  CHARACTER(len = 8)  ::   cn_cfg = "BENGUELA"      ! current name of the NetCDF mask file acting as a key
   

CONTAINS

   SUBROUTINE ocean_2_pisces( kt, Istr, Iend, Jstr, Jend )
   
      INTEGER, intent(in) :: kt, Istr, Iend, Jstr, Jend
      INTEGER :: ji, jj, jk, jl

      IF( kt == nit000 ) THEN
         Istrp = Istr
         Iendp = Iend
         Jstrp = Jstr
         Jendp = Jend
      ENDIF

      DO jl = 1, 3
         DO jk = 1, N+1
            DO jj =  Jstrp, Jendp 
               DO ji =  Istrp, Iendp 
                  gdepw(ji,jj,N+2-jk,jl) = -(z_w(ji,jj,jk-1)-z_w(ji,jj,N))
               END DO
            END DO
         END DO

         DO jk = 2, N
            DO jj =  Jstrp, Jendp 
               DO ji =  Istrp, Iendp 
                  e3w(ji,jj,jk,jl) = -z_r(ji,jj,N+1-jk) + z_r(ji,jj,N+2-jk)
              END DO
            END DO
         END DO

         DO jj =  Jstrp, Jendp 
            DO ji =  Istrp, Iendp 
               e3w(ji,jj,1  ,jl) = -2 * z_r(ji,jj,N)
               e3w(ji,jj,N+1,jl) = 2 * ( -z_w(ji,jj,0) + z_r(ji,jj,1) )
            END DO
         END DO
         !
      ENDDO
      !
   END SUBROUTINE ocean_2_pisces

   INTEGER FUNCTION trc_oce_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  trc_oce_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( tra(A2D(0),jpk,jptra),tmask(A2D(0),jpk), e1e2t(A2D(0)), &
         &      etot3(A2D(0),jpk), oce_co2(A2D(0)), &
         &     atm_co2(A2D(0)), fr_i(A2D(0)),STAT=trc_oce_alloc )

      IF( trc_oce_alloc /= 0 )   CALL ctl_warn('trc_oce_alloc: failed to allocate etot3 array')
      !
   END FUNCTION trc_oce_alloc


   SUBROUTINE trc_oce_rgb( prgb )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of of the optical scheme
      !!
      !! ** Method  :   Set a look up table for the optical coefficients
      !!                i.e. the attenuation coefficient for R-G-B light
      !!                tabulated in Chlorophyll class (from JM Andre)
      !!
      !! ** Action  :   prgb(3,61) tabulated R-G-B attenuation coef.
      !!
      !! Reference  : Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(3,61), INTENT(out) ::   prgb   ! tabulated attenuation coefficient
      !
      INTEGER  ::   jc     ! dummy loop indice
      INTEGER  ::   irgb   ! temporary integer
      REAL(wp) ::   zchl   ! temporary scalar
      REAL(wp), DIMENSION(4,61) ::   zrgb   ! tabulated attenuation coefficient (formerly read in 'kRGB61.txt')
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   trc_oce_rgb : Initialisation of the optical look-up table'
         WRITE(numout,*) '   ~~~~~~~~~~~ '
      ENDIF
      !
      !  Chlorophyll        !     Blue attenuation     !     Green attenuation    !     Red attenuation      !
      zrgb(1, 1) =  0.010   ;   zrgb(2, 1) = 0.01618   ;   zrgb(3, 1) = 0.07464   ;   zrgb(4, 1) = 0.37807
      zrgb(1, 2) =  0.011   ;   zrgb(2, 2) = 0.01654   ;   zrgb(3, 2) = 0.07480   ;   zrgb(4, 2) = 0.37823
      zrgb(1, 3) =  0.013   ;   zrgb(2, 3) = 0.01693   ;   zrgb(3, 3) = 0.07499   ;   zrgb(4, 3) = 0.37840
      zrgb(1, 4) =  0.014   ;   zrgb(2, 4) = 0.01736   ;   zrgb(3, 4) = 0.07518   ;   zrgb(4, 4) = 0.37859
      zrgb(1, 5) =  0.016   ;   zrgb(2, 5) = 0.01782   ;   zrgb(3, 5) = 0.07539   ;   zrgb(4, 5) = 0.37879
      zrgb(1, 6) =  0.018   ;   zrgb(2, 6) = 0.01831   ;   zrgb(3, 6) = 0.07562   ;   zrgb(4, 6) = 0.37900
      zrgb(1, 7) =  0.020   ;   zrgb(2, 7) = 0.01885   ;   zrgb(3, 7) = 0.07586   ;   zrgb(4, 7) = 0.37923
      zrgb(1, 8) =  0.022   ;   zrgb(2, 8) = 0.01943   ;   zrgb(3, 8) = 0.07613   ;   zrgb(4, 8) = 0.37948
      zrgb(1, 9) =  0.025   ;   zrgb(2, 9) = 0.02005   ;   zrgb(3, 9) = 0.07641   ;   zrgb(4, 9) = 0.37976
      zrgb(1,10) =  0.028   ;   zrgb(2,10) = 0.02073   ;   zrgb(3,10) = 0.07672   ;   zrgb(4,10) = 0.38005
      zrgb(1,11) =  0.032   ;   zrgb(2,11) = 0.02146   ;   zrgb(3,11) = 0.07705   ;   zrgb(4,11) = 0.38036
      zrgb(1,12) =  0.035   ;   zrgb(2,12) = 0.02224   ;   zrgb(3,12) = 0.07741   ;   zrgb(4,12) = 0.38070
      zrgb(1,13) =  0.040   ;   zrgb(2,13) = 0.02310   ;   zrgb(3,13) = 0.07780   ;   zrgb(4,13) = 0.38107
      zrgb(1,14) =  0.045   ;   zrgb(2,14) = 0.02402   ;   zrgb(3,14) = 0.07821   ;   zrgb(4,14) = 0.38146
      zrgb(1,15) =  0.050   ;   zrgb(2,15) = 0.02501   ;   zrgb(3,15) = 0.07866   ;   zrgb(4,15) = 0.38189
      zrgb(1,16) =  0.056   ;   zrgb(2,16) = 0.02608   ;   zrgb(3,16) = 0.07914   ;   zrgb(4,16) = 0.38235
      zrgb(1,17) =  0.063   ;   zrgb(2,17) = 0.02724   ;   zrgb(3,17) = 0.07967   ;   zrgb(4,17) = 0.38285
      zrgb(1,18) =  0.071   ;   zrgb(2,18) = 0.02849   ;   zrgb(3,18) = 0.08023   ;   zrgb(4,18) = 0.38338
      zrgb(1,19) =  0.079   ;   zrgb(2,19) = 0.02984   ;   zrgb(3,19) = 0.08083   ;   zrgb(4,19) = 0.38396
      zrgb(1,20) =  0.089   ;   zrgb(2,20) = 0.03131   ;   zrgb(3,20) = 0.08149   ;   zrgb(4,20) = 0.38458
      zrgb(1,21) =  0.100   ;   zrgb(2,21) = 0.03288   ;   zrgb(3,21) = 0.08219   ;   zrgb(4,21) = 0.38526
      zrgb(1,22) =  0.112   ;   zrgb(2,22) = 0.03459   ;   zrgb(3,22) = 0.08295   ;   zrgb(4,22) = 0.38598
      zrgb(1,23) =  0.126   ;   zrgb(2,23) = 0.03643   ;   zrgb(3,23) = 0.08377   ;   zrgb(4,23) = 0.38676
      zrgb(1,24) =  0.141   ;   zrgb(2,24) = 0.03842   ;   zrgb(3,24) = 0.08466   ;   zrgb(4,24) = 0.38761
      zrgb(1,25) =  0.158   ;   zrgb(2,25) = 0.04057   ;   zrgb(3,25) = 0.08561   ;   zrgb(4,25) = 0.38852
      zrgb(1,26) =  0.178   ;   zrgb(2,26) = 0.04289   ;   zrgb(3,26) = 0.08664   ;   zrgb(4,26) = 0.38950
      zrgb(1,27) =  0.200   ;   zrgb(2,27) = 0.04540   ;   zrgb(3,27) = 0.08775   ;   zrgb(4,27) = 0.39056
      zrgb(1,28) =  0.224   ;   zrgb(2,28) = 0.04811   ;   zrgb(3,28) = 0.08894   ;   zrgb(4,28) = 0.39171
      zrgb(1,29) =  0.251   ;   zrgb(2,29) = 0.05103   ;   zrgb(3,29) = 0.09023   ;   zrgb(4,29) = 0.39294
      zrgb(1,30) =  0.282   ;   zrgb(2,30) = 0.05420   ;   zrgb(3,30) = 0.09162   ;   zrgb(4,30) = 0.39428
      zrgb(1,31) =  0.316   ;   zrgb(2,31) = 0.05761   ;   zrgb(3,31) = 0.09312   ;   zrgb(4,31) = 0.39572
      zrgb(1,32) =  0.355   ;   zrgb(2,32) = 0.06130   ;   zrgb(3,32) = 0.09474   ;   zrgb(4,32) = 0.39727
      zrgb(1,33) =  0.398   ;   zrgb(2,33) = 0.06529   ;   zrgb(3,33) = 0.09649   ;   zrgb(4,33) = 0.39894
      zrgb(1,34) =  0.447   ;   zrgb(2,34) = 0.06959   ;   zrgb(3,34) = 0.09837   ;   zrgb(4,34) = 0.40075
      zrgb(1,35) =  0.501   ;   zrgb(2,35) = 0.07424   ;   zrgb(3,35) = 0.10040   ;   zrgb(4,35) = 0.40270
      zrgb(1,36) =  0.562   ;   zrgb(2,36) = 0.07927   ;   zrgb(3,36) = 0.10259   ;   zrgb(4,36) = 0.40480
      zrgb(1,37) =  0.631   ;   zrgb(2,37) = 0.08470   ;   zrgb(3,37) = 0.10495   ;   zrgb(4,37) = 0.40707
      zrgb(1,38) =  0.708   ;   zrgb(2,38) = 0.09056   ;   zrgb(3,38) = 0.10749   ;   zrgb(4,38) = 0.40952
      zrgb(1,39) =  0.794   ;   zrgb(2,39) = 0.09690   ;   zrgb(3,39) = 0.11024   ;   zrgb(4,39) = 0.41216
      zrgb(1,40) =  0.891   ;   zrgb(2,40) = 0.10374   ;   zrgb(3,40) = 0.11320   ;   zrgb(4,40) = 0.41502
      zrgb(1,41) =  1.000   ;   zrgb(2,41) = 0.11114   ;   zrgb(3,41) = 0.11639   ;   zrgb(4,41) = 0.41809
      zrgb(1,42) =  1.122   ;   zrgb(2,42) = 0.11912   ;   zrgb(3,42) = 0.11984   ;   zrgb(4,42) = 0.42142
      zrgb(1,43) =  1.259   ;   zrgb(2,43) = 0.12775   ;   zrgb(3,43) = 0.12356   ;   zrgb(4,43) = 0.42500
      zrgb(1,44) =  1.413   ;   zrgb(2,44) = 0.13707   ;   zrgb(3,44) = 0.12757   ;   zrgb(4,44) = 0.42887
      zrgb(1,45) =  1.585   ;   zrgb(2,45) = 0.14715   ;   zrgb(3,45) = 0.13189   ;   zrgb(4,45) = 0.43304
      zrgb(1,46) =  1.778   ;   zrgb(2,46) = 0.15803   ;   zrgb(3,46) = 0.13655   ;   zrgb(4,46) = 0.43754
      zrgb(1,47) =  1.995   ;   zrgb(2,47) = 0.16978   ;   zrgb(3,47) = 0.14158   ;   zrgb(4,47) = 0.44240
      zrgb(1,48) =  2.239   ;   zrgb(2,48) = 0.18248   ;   zrgb(3,48) = 0.14701   ;   zrgb(4,48) = 0.44765
      zrgb(1,49) =  2.512   ;   zrgb(2,49) = 0.19620   ;   zrgb(3,49) = 0.15286   ;   zrgb(4,49) = 0.45331
      zrgb(1,50) =  2.818   ;   zrgb(2,50) = 0.21102   ;   zrgb(3,50) = 0.15918   ;   zrgb(4,50) = 0.45942
      zrgb(1,51) =  3.162   ;   zrgb(2,51) = 0.22703   ;   zrgb(3,51) = 0.16599   ;   zrgb(4,51) = 0.46601
      zrgb(1,52) =  3.548   ;   zrgb(2,52) = 0.24433   ;   zrgb(3,52) = 0.17334   ;   zrgb(4,52) = 0.47313
      zrgb(1,53) =  3.981   ;   zrgb(2,53) = 0.26301   ;   zrgb(3,53) = 0.18126   ;   zrgb(4,53) = 0.48080
      zrgb(1,54) =  4.467   ;   zrgb(2,54) = 0.28320   ;   zrgb(3,54) = 0.18981   ;   zrgb(4,54) = 0.48909
      zrgb(1,55) =  5.012   ;   zrgb(2,55) = 0.30502   ;   zrgb(3,55) = 0.19903   ;   zrgb(4,55) = 0.49803
      zrgb(1,56) =  5.623   ;   zrgb(2,56) = 0.32858   ;   zrgb(3,56) = 0.20898   ;   zrgb(4,56) = 0.50768
      zrgb(1,57) =  6.310   ;   zrgb(2,57) = 0.35404   ;   zrgb(3,57) = 0.21971   ;   zrgb(4,57) = 0.51810
      zrgb(1,58) =  7.079   ;   zrgb(2,58) = 0.38154   ;   zrgb(3,58) = 0.23129   ;   zrgb(4,58) = 0.52934
      zrgb(1,59) =  7.943   ;   zrgb(2,59) = 0.41125   ;   zrgb(3,59) = 0.24378   ;   zrgb(4,59) = 0.54147
      zrgb(1,60) =  8.912   ;   zrgb(2,60) = 0.44336   ;   zrgb(3,60) = 0.25725   ;   zrgb(4,60) = 0.55457
      zrgb(1,61) = 10.000   ;   zrgb(2,61) = 0.47804   ;   zrgb(3,61) = 0.27178   ;   zrgb(4,61) = 0.56870
      !
      prgb(:,:) = zrgb(2:4,:)
      !
      r_si2 = 1.e0 / zrgb(2, 1)        ! blue with the smallest chlorophyll concentration)
      IF(lwp) WRITE(numout,*) '      RGB longest depth of extinction    r_si2 = ', r_si2
      !
      DO jc = 1, 61                         ! check
         zchl = zrgb(1,jc)
         irgb = NINT( 41 + 20.* LOG10( zchl ) + 1.e-15 )
         IF( irgb /= jc ) THEN
            IF(lwp) WRITE(numout,*) '    jc =', jc, '  Chl = ', zchl, '  Chl class = ', irgb
            CALL ctl_stop( 'STOP', 'trc_oce_rgb : inconsistency in Chl tabulated attenuation coeff.' )
         ENDIF
      END DO
      !
   END SUBROUTINE trc_oce_rgb


   SUBROUTINE trc_oce_rgb_read( prgb )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of of the optical scheme
      !!
      !! ** Method  :   read the look up table for the optical coefficients
      !!
      !! ** input   :   xkrgb(61) precomputed array corresponding to the
      !!                          attenuation coefficient (from JM Andre)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(3,61), INTENT(out) ::   prgb   ! tabulated attenuation coefficient
      !
      INTEGER  ::   jc, jb ! dummy loop indice
      INTEGER  ::   irgb   ! temporary integer
      REAL(wp) ::   zchl   ! temporary scalar
      INTEGER  ::   numlight
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' trc_oce_rgb_read : optical look-up table read in kRGB61.txt file'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~'
         WRITE(numout,*)
      ENDIF
      !
      CALL ctl_opn( numlight, 'kRGB61.txt', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
      DO jc = 1, 61
         READ(numlight,*) zchl, ( prgb(jb,jc), jb = 1, 3 )
         irgb = NINT( 41 + 20.* LOG10( zchl ) + 1.e-15 )
         IF(lwp) WRITE(numout,*) '    jc =', jc, '  Chl = ', zchl, '  irgb = ', irgb
         IF( irgb /= jc ) THEN
            IF(lwp) WRITE(numout,*) '    jc =', jc, '  Chl = ', zchl, '  Chl class = ', irgb
            CALL ctl_stop( 'STOP','trc_oce_rgb_read : inconsistency in Chl tabulated attenuation coeff.' )
         ENDIF
      END DO
      CLOSE( numlight )
      !
      r_si2 = 1.e0 / prgb(1, 1)      ! blue with the smallest chlorophyll concentration)
      IF(lwp) WRITE(numout,*) '      RGB longest depth of extinction    r_si2 = ', r_si2
      !
   END SUBROUTINE trc_oce_rgb_read


   FUNCTION trc_oce_ext_lev( prldex, pqsr_frc ) RESULT( pjl )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE trc_oce_ext_lev  ***
      !!
      !! ** Purpose :   compute max. level for light penetration
      !!
      !! ** Method  :   the function provides the level at which irradiance
      !!                becomes negligible (i.e. = 1.e-15 W/m2) for 3 or 2 bands light
      !!                penetration: I(z) = pqsr_frc * EXP(hext/prldex) = 1.e-15 W/m2
      !!                # prldex is the longest depth of extinction:
      !!                   - prldex = 23 m (2 bands case)
      !!                   - prldex = 62 m (3 bands case: blue waveband & 0.01 mg/m2 for the chlorophyll)
      !!                # pqsr_frc is the fraction of solar radiation which penetrates,
      !!                considering Qsr=240 W/m2 and rn_abs = 0.58:
      !!                   - pqsr_frc = Qsr * (1-rn_abs)   = 1.00e2 W/m2 (2 bands case)
      !!                   - pqsr_frc = Qsr * (1-rn_abs)/3 = 0.33e2 W/m2 (3 bands case & equi-partition)
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::   prldex    ! longest depth of extinction
      REAL(wp), INTENT(in) ::   pqsr_frc  ! frac. solar radiation which penetrates
      !
      INTEGER  ::   jk, pjl            ! levels
      REAL(wp) ::   zhext              ! deepest level till which light penetrates
      REAL(wp) ::   zprec = 15._wp     ! precision to reach -LOG10(1.e-15)
      REAL(wp) ::   zem                ! temporary scalar
      !!----------------------------------------------------------------------
      !
      ! It is not necessary to compute anything below the following depth
      zhext = prldex * ( LOG(10._wp) * zprec + LOG(pqsr_frc) )
      !
      ! Level of light extinction
      pjl = jpkm1
      DO jk = jpkm1, 1, -1
         IF(SUM(tmask(A2D(0),jk)) > 0 ) THEN
!            zem = MAXVAL( gdepw_1d(jk+1) * tmask(:,:,jk) )
            zem = MAXVAL( gdepw(A2D(0),jk+1,1) * tmask(A2D(0),jk) )
            IF( zem >= zhext )   pjl = jk                       ! last T-level reached by Qsr
         ELSE
            pjl = jk                                            ! or regional sea-bed depth
         ENDIF
      END DO
      !
   END FUNCTION trc_oce_ext_lev

      FUNCTION glob_sum_2d( cdname, ptab )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_sum_2d ***
      !!
      !! ** Purpose : perform a sum in calling DDPDD routine
      !!----------------------------------------------------------------------
      CHARACTER(len=*),  INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(wp), INTENT(in), DIMENSION(A2D(0)) ::   ptab
      REAL(wp)                             ::   glob_sum_2d   ! global masked sum
      !!
      REAL(wp)   ::   ztmp
      INTEGER    ::   ji, jj   ! dummy loop indices
      !!-----------------------------------------------------------------------
      !
      ztmp = 0.e0
      DO_2D( 0, 0, 0, 0 )
         ztmp =  ptab(ji,jj) * rmask(ji,jj)
      END_2D
      IF( lk_mpp )   CALL mpp_sum( ztmp )   ! sum over the global domain
      glob_sum_2d = REAL(ztmp,wp)
      !
   END FUNCTION glob_sum_2d

      FUNCTION glob_sum_3d( cdname, ptab )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_sum_3d ***
      !!
      !! ** Purpose : perform a sum on a 3D array in calling DDPDD routine
      !!----------------------------------------------------------------------
      CHARACTER(len=*),  INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(wp), INTENT(in), DIMENSION(A2D(0),jpk) ::   ptab
      REAL(wp)                               ::   glob_sum_3d   ! global masked sum
      !!
      REAL(wp)   ::   ztmp
      INTEGER    ::   ji, jj, jk   ! dummy loop indices
      !!-----------------------------------------------------------------------
      !
      !
      ztmp = 0.e0
      DO_3D( 0, 0, 0, 0, 1, jpk)
         ztmp =  ptab(ji,jj,jk) * rmask(ji,jj)
      END_3D
      IF( lk_mpp )   CALL mpp_sum( ztmp )   ! sum over the global domain
      glob_sum_3d = REAL(ztmp,wp)
      !
   END FUNCTION glob_sum_3d

   SUBROUTINE tracer_stat( kt, clname )
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_stat  ***
      !!
      !! ** purpose  :   Compute tracers statistics
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)  :: kt
      CHARACTER(len=20), DIMENSION(jptra), INTENT(in) ::   clname
      !
      INTEGER  :: ji, jj, jk, jn
      REAL(wp) :: ztra, zmin, zmax, zmean, areatot, zcoef
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY,jptra)  :: ptra
      REAL(wp), DIMENSION(PRIV_3D_BIOARRAY)        :: zmask, zvol
      !!----------------------------------------------------------------------

      IF( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' TRACER STAT at time-step kt = ', kt
         WRITE(numout,*)
      ENDIF
      !
! to have coherent units when calling tracer_stat
      IF( kt .eq. nit000 ) THEN
        zcoef = 1.e-6
      ELSE
        zcoef = 1.
      ENDIF

      DO jn = 1, jptra
         DO_3D( 0, 0, 0, 0, 1, jpk)
            ptra(ji,jj,jk,jn) = tr(ji,jj,jk,jn,nnew) * zcoef
         END_3D
      ENDDO
      areatot = 0.                                                           ! total volume of the ocean
      DO_3D( 0, 0, 0, 0, 1, jpk)
          zvol(ji,jj,jk)  = e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk,Kmm) * tmask(ji,jj,jk)
          zmask(ji,jj,jk) = tmask(ji,jj,jk) * tmask_i(ji,jj)
          areatot         = areatot + zvol(ji,jj,jk)
     END_3D
     IF( lk_mpp )   CALL mpp_sum( areatot )     ! sum over the global domain

     DO jn = 1, jptra
         ztra = 0.
         DO_3D( 0, 0, 0, 0, 1, jpk)
            ztra  = ztra + ptra(ji,jj,jk,jn) * zvol(ji,jj,jk)
         END_3D
         zmin  = MINVAL( ptra(:,:,:,jn), mask= ( zmask(:,:,:) /= 0. ) )
         zmax  = MAXVAL( ptra(:,:,:,jn), mask= ( zmask(:,:,:) /= 0. ) )
         IF( lk_mpp ) THEN
            CALL mpp_sum( ztra )      ! min over the global domain
            CALL mpp_min( zmin )      ! min over the global domain
            CALL mpp_max( zmax )      ! max over the global domain
         END IF
         zmean  = ztra / areatot
         IF(lwp) WRITE(numout,9000) jn, TRIM( clname(jn) ), zmean, zmin, zmax
      END DO
      WRITE(numout,*)
9000  FORMAT(' tracer nb :',i2,'    name :',a10,'    mean :',e18.10,'    min :',e18.10, '    max :',e18.10 )
      !
   END SUBROUTINE tracer_stat

# ifdef NEMO
   SUBROUTINE fld_read( kt, kn_fsbc, sd, kit, pt_offset, Kmm )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_read  ***
      !!                   
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!                (momentum, heat, freshwater and runoff) 
      !!
      !! ** Method  :   READ each input fields in NetCDF files using IOM
      !!      and intepolate it to the model time-step.
      !!         Several assumptions are made on the input file:
      !!      blahblahblah....
      !!----------------------------------------------------------------------
      INTEGER  , INTENT(in   )               ::   kt        ! ocean time step
      INTEGER  , INTENT(in   )               ::   kn_fsbc   ! sbc computation period (in time step) 
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd        ! input field related variables
      INTEGER  , INTENT(in   ), OPTIONAL     ::   kit       ! subcycle timestep for timesplitting option
      REAL(wp) , INTENT(in   ), OPTIONAL     ::   pt_offset ! provide fields at time other than "now"
      INTEGER  , INTENT(in   ), OPTIONAL     ::   Kmm       ! ocean time level index
      !!
      INTEGER  ::   imf          ! size of the structure sd
      INTEGER  ::   jf           ! dummy indices
      INTEGER  ::   isecsbc      ! number of seconds between Jan. 1st 00h of nit000 year and the middle of sbc time step
      INTEGER  ::   ibb, iaa     ! shorter name for sd(jf)%nbb and sd(jf)%naa
      LOGICAL  ::   ll_firstcall ! true if this is the first call to fld_read for this set of fields
      REAL(wp) ::   zt_offset    ! local time offset variable
      REAL(wp) ::   ztinta       ! ratio applied to after  records when doing time interpolation
      REAL(wp) ::   ztintb       ! ratio applied to before records when doing time interpolation
      CHARACTER(LEN=1000) ::   clfmt  ! write format
      !!---------------------------------------------------------------------
   END SUBROUTINE fld_read
   

   SUBROUTINE fld_fill( sdf, sdf_n, cdir, cdcaller, cdtitle, cdnam, knoprint )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_fill  ***
      !!
      !! ** Purpose :   fill the data structure (sdf) with the associated information
      !!              read in namelist (sdf_n) and control print
      !!----------------------------------------------------------------------
      TYPE(FLD)  , DIMENSION(:)          , INTENT(inout) ::   sdf        ! structure of input fields (file informations, fields read)
      TYPE(FLD_N), DIMENSION(:)          , INTENT(in   ) ::   sdf_n      ! array of namelist information structures
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdir       ! Root directory for location of flx files
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdcaller   ! name of the calling routine
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdtitle    ! description of the calling routine
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdnam      ! name of the namelist from which sdf_n comes
      INTEGER                  , OPTIONAL, INTENT(in   ) ::   knoprint   ! no calling routine information printed
      !
   END SUBROUTINE fld_fill
# endif
#endif
END MODULE oce_trc

