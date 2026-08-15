#include "cppdefs.h"

MODULE p4zbc
   !!======================================================================
   !!                         ***  MODULE p4zbc  ***
   !! TOP :   PISCES surface boundary conditions of external inputs of nutrients
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, C. Ethe) Original code
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_bc        :  Read and interpolate time-varying nutrients fluxes
   !!   p4z_bc_init   :  Initialization of p4z_bc
   !!----------------------------------------------------------------------
   USE oce_trc
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE iom             !  I/O manager
#ifdef AGRIF
      USE param, ONLY : Lmmpi,Mmmpi
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_bc
   PUBLIC   p4z_bc_init   

   !! * Module variables
   LOGICAL , PUBLIC ::   ln_dust      !: boolean for dust input from the atmosphere
   LOGICAL , PUBLIC ::   ln_ndepo     !: boolean for atmospheric deposition of N
   LOGICAL , PUBLIC ::   ln_ironsed   !: boolean for Fe input from sediments
   REAL(wp), PUBLIC ::   sedfeinput   !: Coastal release of Iron
   REAL(wp), PUBLIC ::   mfrac        !: Mineral Content of the dust
   REAL(wp), PUBLIC ::   wdust        !: Sinking speed of the dust 
   REAL(wp), PUBLIC ::   lgw_rath     !: Weak ligand ratio from sed hydro sources
   LOGICAL , PUBLIC ::   ll_bc
   LOGICAL , PUBLIC ::   ll_dust      !: boolean for dust input from the atmosphere

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   dust             !: dust fields
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ironsed          !: Coastal supply of iron
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   dustmo, ferdepmo !: 2 consecutive set of dust fields 
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   po4depmo, sildepmo !: 2 consecutive set of dust fields 
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   no3depmo, nh4depmo !: 2 consecutive set of dust fields 

   REAL(wp), PUBLIC :: sedsilfrac, sedcalfrac
   REAL(wp), PUBLIC :: year2daydta

   LOGICAL  :: l_dia_iron, l_dia_dust, l_dia_ndep

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsbc.F90 10868 2019-04-15 10:32:56Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_bc( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_bc  ***
      !!
      !! ** purpose :   read and interpolate the external sources of nutrients
      !!
      !! ** method  :   read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt              ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level index     
      !
      INTEGER :: ji, jj, jk
      INTEGER, PARAMETER :: jpmois = 12
      INTEGER  :: irec1, irec2, i15
      REAL(wp) :: zpdtan, zpdtmo, zdemi, zt
      REAL(wp) :: zxy, zjulian, zsec
      !
      REAL(wp)   :: zdust, zwdust, ztra
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zno3dep, znh4dep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zpo4dep, zsildep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zirondep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw2d
      !!---------------------------------------------------------------------
      !
      ! Compute dust at nit000 or only if there is more than 1 time record in dust file
      IF( kt == nit000 ) THEN
        l_dia_iron   = iom_use( "Ironsed" ) 
        l_dia_dust   = iom_use( "Irondep" ) .OR. iom_use( "pdust" ) &
           &            .OR. iom_use( "Sildep" ) .OR. iom_use( "Po4dep" )
        l_dia_ndep   = iom_use( "No3dep" ) .OR. iom_use( "Nh4dep" )
        IF( lwp ) THEN
           WRITE(numout,*) ' '
           WRITE(numout,*) ' Number of days per year in file year2daydta = ', year2daydta
           WRITE(numout,*) ' '
        ENDIF
      ENDIF
      !
      zpdtan = ( year2daydta * day2sec ) / rdt
      zpdtmo = zpdtan / float( jpmois )
      zdemi  = zpdtmo / 2.
      zt     = ( float( kt ) + zdemi) / zpdtmo


      !  recherche de l'indice des enregistrements
      !  du modele dynamique encadrant le pas de temps kt.
      !  --------------------------------------------------
      irec1 = int( zt )
      irec2 = irec1 + 1
      irec1 = MOD( irec1, jpmois )
      IF ( irec1 == 0 ) irec1 = jpmois
      irec2 = MOD( irec2, jpmois )
      IF ( irec2 == 0 ) irec2 = jpmois

      zxy = zt - float(int ( zt ) )
      !
      IF( ln_dust ) THEN
         ALLOCATE( zirondep(A2D(0), jpk) )
         !
         DO_2D( 0, 0, 0, 0 )
            dust(ji,jj) = ( 1. - zxy ) * dustmo(ji,jj,irec1) + zxy * dustmo(ji,jj,irec2)
         END_2D
         !
         DO_2D( 0, 0, 0, 0 )
            zirondep(ji,jj,1) =  ( ( 1. - zxy ) * ferdepmo(ji,jj,irec1) &
                    &                   + zxy   * ferdepmo(ji,jj,irec2) ) &
               &                   * rfact2 / e3t(ji,jj,1,Kmm) 
         END_2D
         !                                              ! Iron solubilization of particles in the water column
         !                                              ! dust in kg/m2/s ---> 1/55.85 to put in mol/Fe ;  wdust in m/j
         zwdust = 0.03  / ( wdust / rday ) / ( 250. * rday )
         DO_3D( 0, 0, 0, 0, 2, jpk)
            zirondep(ji,jj,jk) = ( dust(ji,jj) * mfrac / mMass_Fe ) * zwdust &
                   &        * rfact * EXP( -gdept(ji,jj,jk,Kmm) /( 250. * wdust ) )
         END_3D
         !                                              ! Iron solubilization of particles in the water column
         DO_3D( 0, 0, 0, 0, 1, jpk)
            tr(ji,jj,jk,jpfer,Krhs) = tr(ji,jj,jk,jpfer,Krhs) + zirondep(ji,jj,jk)
         END_3D
         !
         IF( lk_iomput .AND. l_dia_dust ) THEN
            ALLOCATE( zw3d(GLOBAL_2D_ARRAY,jpk) )   ;   zw3d(:,:,:) = 0._wp
            ALLOCATE( zw2d(GLOBAL_2D_ARRAY) )       ;   zw2d(:,:) = 0._wp
            DO_3D( 0, 0, 0, 0, 1, jpk )
               zw3d(ji,jj,jkR) = zirondep(ji,jj,jk) &
                    &         * 1.e+3 * rfact2r * e3t(ji,jj,jk,Kmm) * tmask(ji,jj,jk)
            END_3D
            CALL iom_put( "Irondep", zw3d )  ! surface downward dust depo of iron
            !
            DO_2D( 0, 0, 0, 0 )
               zw2d(ji,jj) = dust(ji,jj) / ( wdust * rday ) * tmask(ji,jj,1) ! dust concentration at surface
            END_2D
            CALL iom_put( "pdust", zw2d ) ! dust concentration at surface
            DEALLOCATE( zw2d, zw3d )
            !
         ENDIF
#if defined key_trc_diaadd
         DO_3D( 0, 0, 0, 0, 1, jpk )
            trc3d(ji,jj,jkR,jp_irondep)  = zirondep(ji,jj,jk) * 1.e+3 * tmask(ji,jj,jk)
         END_3D
# endif
         DEALLOCATE( zirondep )
         !                                              
#if ! defined key_pisces_npzd
         ! Atmospheric input of PO4 and Si dissolves in the water column
         ALLOCATE( zpo4dep(A2D(0)), zsildep(A2D(0)) )
         DO_2D( 0, 0, 0, 0 )
            zpo4dep(ji,jj) = ( 1. - zxy ) * po4depmo(ji,jj,irec1) &
               &                  + zxy   * po4depmo(ji,jj,irec2)
            ztra   = zpo4dep(ji,jj) * rfact2 / e3t(ji,jj,1,Kmm)
            tr(ji,jj,1,jppo4,Krhs) = tr(ji,jj,1,jppo4,Krhs) + ztra
         END_2D
         !
         DO_2D( 0, 0, 0, 0 )
            zsildep(ji,jj) = ( 1. - zxy ) * sildepmo(ji,jj,irec1) &
               &                  + zxy   * sildepmo(ji,jj,irec2)
            ztra   = zsildep(ji,jj) * rfact2 / e3t(ji,jj,1,Kmm)
            tr(ji,jj,1,jpsil,Krhs) = tr(ji,jj,1,jpsil,Krhs) + ztra
         END_2D
#if defined key_trc_diaadd
         DO_2D( 0, 0, 0, 0 )
            trc2d(ji,jj,jp_sildep) = zsildep(ji,jj) * 1.e+3 * tmask(ji,jj,1)         ! Si surface deposition
            trc2d(ji,jj,jp_po4dep) = zpo4dep(ji,jj) * po4r * 1.e+3 * tmask(ji,jj,1)  ! PO4 surface deposition
         END_2D
# endif

         DO_3D( 0, 0, 0, 0, 2, jpk )
            zdust = dust(ji,jj) * zwdust * rfact2 &
               &     * EXP( -gdept(ji,jj,jk,Kmm) /( 250. * wdust ) )
            tr(ji,jj,jk,jppo4,Krhs) = tr(ji,jj,jk,jppo4,Krhs) &
                    &                 + zdust * 1.e-3 / mMass_P
            tr(ji,jj,jk,jpsil,Krhs) = tr(ji,jj,jk,jpsil,Krhs) &
                    &                 + zdust * 0.269 / mMass_Si
         END_3D
         !
         IF( lk_iomput .AND. l_dia_dust ) THEN
            !
            ALLOCATE( zw2d(GLOBAL_2D_ARRAY) )   ;   zw2d(:,:) = 0.
            DO_2D( 0, 0, 0, 0 )
               zw2d(ji,jj) = zpo4dep(ji,jj) * 1.e+3 * po4r * tmask(ji,jj,1)        ! Si surface deposition
            END_2D
            CALL iom_put( "Po4dep", zw2d ) ! PO4 concentration at surface

            DO_2D( 0, 0, 0, 0 )
               zw2d(ji,jj) = zsildep(ji,jj) * 1.e+3 * tmask(ji,jj,1)        ! Si surface deposition
            END_2D
            CALL iom_put( "Sildep", zw2d ) ! Si concentration at surface
            DEALLOCATE( zw2d )
         ENDIF
         !
         DEALLOCATE( zpo4dep, zsildep )
# endif
         !
      ENDIF
      !
      IF( ln_ndepo ) THEN
         ALLOCATE( zno3dep(A2D(0)) )
         !
         DO_2D( 0, 0, 0, 0 )
            zno3dep(ji,jj) = ( 1. - zxy ) * no3depmo(ji,jj,irec1) &
               &                  + zxy   * no3depmo(ji,jj,irec2)
            ztra  = zno3dep(ji,jj) * rfact2 / e3t(ji,jj,1,Kmm)
            tr(ji,jj,1,jpno3,Krhs) = tr(ji,jj,1,jpno3,Krhs) + ztra
            tr(ji,jj,1,jptal,Krhs) = tr(ji,jj,1,jptal,Krhs) - rno3 * ztra
         END_2D
         !
         IF( lk_iomput .AND. l_dia_ndep ) THEN
             ALLOCATE( zw2d(GLOBAL_2D_ARRAY) )   ;   zw2d(:,:) = 0.
             DO_2D( 0, 0, 0, 0 )
                zw2d(ji,jj) = zno3dep(ji,jj) * 1.e+3 * rno3 * tmask(ji,jj,1)
             END_2D
             CALL iom_put( "No3dep", zw2d )
             DEALLOCATE( zw2d )
         ENDIF
#if defined key_trc_diaadd
         DO_2D( 0, 0, 0, 0 )
            trc2d(ji,jj,jp_no3dep )  = zno3dep(ji,jj) * 1.e+3 * rno3 * tmask(ji,jj,1)
         END_2D
# endif
         DEALLOCATE( zno3dep )

#if ! defined key_pisces_npzd
         ALLOCATE( znh4dep(A2D(0)) )
         !
         DO_2D( 0, 0, 0, 0 )
            znh4dep(ji,jj) = ( 1. - zxy ) * nh4depmo(ji,jj,irec1) &
                &                 + zxy   * nh4depmo(ji,jj,irec2)
            ztra  = znh4dep(ji,jj) * rfact2 / e3t(ji,jj,1,Kmm)
            tr(ji,jj,1,jpnh4,Krhs) = tr(ji,jj,1,jpnh4,Krhs) + ztra 
            tr(ji,jj,1,jptal,Krhs) = tr(ji,jj,1,jptal,Krhs) + rno3 * ztra
         END_2D
         !
         IF( lk_iomput .AND. l_dia_ndep ) THEN
             ALLOCATE( zw2d(GLOBAL_2D_ARRAY) )   ;   zw2d(:,:) = 0.
             DO_2D( 0, 0, 0, 0 )
                zw2d(ji,jj) = znh4dep(ji,jj) * 1.e+3 * rno3 * tmask(ji,jj,1)
             END_2D
             CALL iom_put( "Nh4dep", zw2d )
             DEALLOCATE( zw2d )
         ENDIF
#if defined key_trc_diaadd
         DO_2D( 0, 0, 0, 0 )
            trc2d(ji,jj,jp_nh4dep ) = znh4dep(ji,jj) * 1.e+3 * rno3 * tmask(ji,jj,1)
         END_2D
# endif
         DEALLOCATE( znh4dep )
# endif
      ENDIF
      ! Add the external input of iron from sediment mobilization
      ! ------------------------------------------------------
      IF( ln_ironsed ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpk )
           tr(ji,jj,jk,jpfer,Krhs) = tr(ji,jj,jk,jpfer,Krhs) &
               &                   + ironsed(ji,jj,jk) * rfact2
         END_3D
         !
         IF( lk_iomput .AND. l_dia_iron ) THEN
            ALLOCATE( zw3d(GLOBAL_2D_ARRAY,jpk) )   ;   zw3d(:,:,:) = 0.
            DO_3D( 0, 0, 0, 0, 1, jpk )
               zw3d(ji,jj,jk) = ironsed(ji,jj,jk) * 1.e+3 * tmask(ji,jj,jk)
            END_3D
            CALL iom_put( "Ironsed", zw3d )  ! iron inputs from sediments
            DEALLOCATE( zw3d )
         ENDIF
#if defined key_trc_diaadd
         DO_3D( 0, 0, 0, 0, 1, jpk )
            trc3d(ji,jj,jkR,jp_ironsed ) = ironsed(ji,jj,jk) * 1e+3 * tmask(ji,jj,jk)  ! iron from  sediment
         END_3D
#endif
      ENDIF

   END SUBROUTINE p4z_bc


   SUBROUTINE p4z_bc_init( Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_bc_init  ***
      !!
      !! ** purpose :   initialization of the external sources of nutrients
      !!
      !! ** method  :   read the files and compute the budget
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)  :: Kmm
#if defined key_pisces      
# include "netcdf.inc"
      INTEGER  :: ji, jj, jk, irec
      INTEGER :: ncid, varid, dimid, ierr, lstr, lenstr, nf_fread, nrec_dust, nrec_ndep
      INTEGER :: vartype, nvatts, latt, nvdims
      INTEGER :: vdims(5)
      INTEGER  :: ios                 ! Local integer output status for namelist read
      CHARACTER(len=16) :: varname, dimname, attname
      REAL(wp) :: zexpide, zdenitide, zmaskt

      REAL     ::  cycle_length
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  dustmp
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  zcmask
#ifdef MPI
#define LOCALLM Lmmpi
#define LOCALMM Mmmpi
#else
#define LOCALLM Lm
#define LOCALMM Mm
#endif
      !
      NAMELIST/nampisbc/ln_dust, ln_ndepo,  ln_ironsed, &
        &                sedfeinput, wdust, mfrac, lgw_rath
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_bc_init : initialization of the external sources of nutrients '
         WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF
      !                            !* set file information
      READ_NML_REF(numnatp,nampisbc)
      READ_NML_CFG(numnatp,nampisbc)
      IF(lwm) WRITE ( numonp, nampisbc )      
 
      IF(lwp) THEN
         WRITE(numout,*) '   Namelist : nampissbc '
         WRITE(numout,*) '      dust input from the atmosphere           ln_dust     = ', ln_dust
         WRITE(numout,*) '      atmospheric deposition of n              ln_ndepo    = ', ln_ndepo
         WRITE(numout,*) '      Fe input from sediments                  ln_ironsed  = ', ln_ironsed
         IF( ln_ironsed ) THEN
            WRITE(numout,*) '      coastal release of iron                  sedfeinput  = ', sedfeinput
         ENDIF
         IF( ln_ligand ) THEN
            WRITE(numout,*) '      Weak ligand ratio from sed hydro sources  lgw_rath   = ', lgw_rath
         ENDIF
         IF( ln_dust ) THEN
            WRITE(numout,*) '      Mineral Fe content of the dust           mfrac       = ', mfrac
            WRITE(numout,*) '      sinking speed of the dust                wdust       = ', wdust
         ENDIF
      ENDIF

      IF( ln_dust .OR. ln_ndepo .OR. ln_ironsed ) THEN   ;   ll_bc = .TRUE.
      ELSE                                               ;   ll_bc = .FALSE.
      ENDIF

      ll_dust = ln_dust .OR. ln_sediment

      ! coastal and island masks
      ! ------------------------
      IF( ln_ironsed ) THEN     
         ALLOCATE( zcmask(A2D(0),jpk), ironsed(A2D(0),jpk) )
         zcmask(:,:,:) = 0.0
         DO_2D( 0, 0, 0, 0 )
            zcmask(ji,jj,jpk) = 1
         END_2D

         DO_3D( 0, 0, 0, 0, 1, jpk )
            IF( tmask_i(ji,jj) /= 0. ) THEN
               zmaskt = tmask_i(ji+1,jj) * tmask_i(ji-1,jj  ) * tmask_i(ji,jj+1)    &
                 &                      * tmask_i(ji  ,jj-1) * tmask_i(ji,jj  )
               IF( zmaskt == 0. )   zcmask(ji,jj,jk ) = 0.1
            ENDIF
         END_3D
         !
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zexpide   = MIN( 8.,( gdept(ji,jj,jk,Kmm) / 500. )**(-1.5) )
            zdenitide = -0.9543 + 0.7662 * LOG( zexpide ) - 0.235 * LOG( zexpide )**2
            zcmask(ji,jj,jk) = zcmask(ji,jj,jk) * MIN( 1., EXP( zdenitide ) / 0.5 )
         END_3D
         ! Coastal supply of iron
         ! -------------------------
         ironsed(:,:,jpk) = 0.
         DO_3D( 0, 0, 0, 0, 1, jpk )
            ironsed(ji,jj,jk) = sedfeinput * zcmask(ji,jj,jk) / ( e3t(ji,jj,jk,Kmm) * rday )
         END_3D
         DEALLOCATE( zcmask)
      ENDIF
      !
      !
      !    READ DUST INPUT FROM ATMOSPHERE
      !    -------------------------------------
      IF( ln_dust .OR. ln_ndepo ) THEN
         lstr = lenstr(bioname)
         ierr = nf_open (bioname(1:lstr), nf_nowrite, ncid)
         IF (ierr .eq. nf_noerr ) THEN
            IF( lwp) WRITE(numout,*) ' Read atmospheric deposition in file = ', TRIM(bioname)
         ELSE
            IF( lwp) CALL ctl_stop( 'STOP', 'p4z_bc_ini : Needed input file for atmospheric deposition &
                & or put ln_dust and ln_ndepo to false' )
         ENDIF
         ierr = nf_inq_varid(ncid,"dust_time",varid)
! bug if compilation with gfortran
!         ierr =nf_inq_var (ncid, varid, varname, vartype, nvdims,  vdims,  nvatts) 
         ierr =nf_inq_varnatts (ncid, varid, nvatts) 
         year2daydta = year2day
         DO ji = 1, nvatts
            ierr = nf_inq_attname (ncid, varid, ji, attname)
            IF (ierr == nf_noerr) THEN
               latt = lenstr(attname)
               IF (attname(1:latt) == 'cycle_length') THEN
                  ierr = nf_get_att_FTYPE (ncid, varid, attname(1:latt), cycle_length)
                  IF (ierr == nf_noerr) THEN
                     year2daydta = cycle_length
                  ELSE
                     IF (lwp) write(numout,'(/1x,4A/)') 'SET_CYCLE ERROR while ', &
                     &        'reading attribute ''', attname(1:latt), '''.'
                  ENDIF
               ENDIF
            ENDIF
         END DO
      ENDIF

      IF ( ln_dust ) THEN
         lstr = lenstr(bioname)
         ierr = nf_open (bioname(1:lstr), nf_nowrite, ncid)
         IF ( ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,4) bioname
         ENDIF
         ierr = nf_inq_varid (ncid,"dust",varid)
         IF (ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,5) "dust", bioname
         ENDIF
         ierr = nf_inq_dimid(ncid,"dust_time",dimid)
         ierr = nf_inq_dimlen(ncid,dimid,nrec_dust)
         !
         IF (lwp) WRITE(numout,*) ' nrec_dust = ', nrec_dust
         ALLOCATE( dustmp(GLOBAL_2D_ARRAY,nrec_dust), dustmo(GLOBAL_2D_ARRAY,nrec_dust) )
         ALLOCATE( ferdepmo(GLOBAL_2D_ARRAY,nrec_dust), dust(PRIV_2D_BIOARRAY) )
         !
         DO irec = 1, nrec_dust
            ierr = nf_fread(dustmp(START_2D_ARRAY,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. lwp ) THEN
               WRITE(numout,6) "dust", irec
            ENDIF
         END DO
         !
         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,'(6x,A,1x,I4)') &
#ifdef MPI
         &                   'TRCINI_PISCES -- Read dust deposition ', mynode
#else
         &                   'TRCINI_PISCES -- Read dust deposition '
#endif
         DO irec = 1, nrec_dust
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  dustmo(ji,jj,irec) = MAX( 0., dustmp(ji,jj,irec) )
               ENDDO
            ENDDO
         ENDDO
         !
         ierr = nf_inq_varid (ncid,"dustfer",varid)
         IF (ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,5) "dustfer", bioname
         ENDIF
         DO irec = 1, nrec_dust
            ierr = nf_fread(dustmp(START_2D_ARRAY,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. lwp ) THEN
               WRITE(numout,6) "dustfer", irec
            ENDIF
         END DO
         !
         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,'(6x,A,1x,I4)') &
#ifdef MPI
         &                   'TRCINI_PISCES -- Read iron deposition ', mynode
#else
         &                   'TRCINI_PISCES -- Read iron deposition '
#endif
         DO irec = 1, nrec_dust
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  ferdepmo(ji,jj,irec) = MAX( 0., dustmp(ji,jj,irec) ) * mfrac / mMass_Fe
               ENDDO
            ENDDO
         ENDDO
         !
#if ! defined key_pisces_npzd
         ALLOCATE( sildepmo(GLOBAL_2D_ARRAY,nrec_dust), po4depmo(GLOBAL_2D_ARRAY,nrec_dust) )
         !
         ierr = nf_inq_varid (ncid,"dustpo4",varid)
         IF (ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,5) "dustpo4", bioname
         ENDIF
         DO irec = 1, nrec_dust
            ierr = nf_fread(dustmp(START_2D_ARRAY,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. lwp ) THEN
               WRITE(numout,6) "dustpo4", irec
            ENDIF
         END DO
         !
         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,'(6x,A,1x,I4)') &
#ifdef MPI
         &                   'TRCINI_PISCES -- Read phosphorus deposition ', mynode
#else
         &                   'TRCINI_PISCES -- Read phosphorus deposition '
#endif
         DO irec = 1, nrec_dust
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  po4depmo(ji,jj,irec) = MAX( 0., dustmp(ji,jj,irec) ) * 1e-3 / mMass_P / po4r
               ENDDO
            ENDDO
         ENDDO
         !
         ierr = nf_inq_varid (ncid,"dustsi",varid)
         IF (ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,5) "dustsi", bioname
         ENDIF
         DO irec = 1, nrec_dust
            ierr = nf_fread(dustmp(START_2D_ARRAY,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. lwp ) THEN
               WRITE(numout,6) "dustsi", irec
            ENDIF
         END DO
         !
         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,'(6x,A,1x,I4)') &
#ifdef MPI
         &                   'TRCINI_PISCES -- Read silicate deposition ', mynode
#else
         &                   'TRCINI_PISCES -- Read silicate deposition '
#endif
         DO irec = 1, nrec_dust
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  sildepmo(ji,jj,irec) = MAX( 0., dustmp(ji,jj,irec) ) * 0.269 / mMass_Si
               ENDDO
            ENDDO
         ENDDO
         !
#endif
         ierr = nf_close(ncid)
         DEALLOCATE(dustmp)
         !
      ELSE
         ALLOCATE( dust(PRIV_2D_BIOARRAY) )
         dust(:,:) = 0.
      ENDIF
!
!    READ N DEPOSITION FROM ATMOSPHERE (use dust_time for time)
!    -------------------------------------
!
      IF (ln_ndepo) THEN
         lstr = lenstr(bioname)
         ierr = nf_open (bioname(1:lstr), nf_nowrite, ncid)
         IF (ierr .NE. nf_noerr .AND. lwp) THEN
            WRITE(numout,4) bioname
         ENDIF
         ierr = nf_inq_dimid(ncid,"ndepo_time",dimid)
         ierr = nf_inq_dimlen(ncid,dimid,nrec_ndep)
         ALLOCATE( dustmp(GLOBAL_2D_ARRAY,nrec_ndep) )
#if ! defined key_pisces_npzd
         ALLOCATE( no3depmo(GLOBAL_2D_ARRAY,nrec_ndep) )
         ierr = nf_inq_varid (ncid,"noyndepo",varid)
         IF (ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,5) "noyndepo", bioname
         ENDIF
         !
         DO irec = 1, nrec_ndep
            ierr = nf_fread(dustmp(START_2D_ARRAY,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. lwp ) THEN
               WRITE(numout,6) "noyndepo", irec
            ENDIF
         END DO
         !
         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,'(6x,A,1x,I4)') &
#ifdef MPI
         &                   'TRCINI_PISCES -- Read Nitrate deposition ', mynode
#else
         &                   'TRCINI_PISCES -- Read Nitrate deposition '
#endif

         DO irec = 1, nrec_ndep
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  no3depmo(ji,jj,irec) = MAX( 0., dustmp(ji,jj,irec) ) / rno3 / mMass_N 
               END DO
            END DO
         END DO
         !
         ALLOCATE( nh4depmo(GLOBAL_2D_ARRAY,nrec_ndep) )
         ierr = nf_inq_varid (ncid,"nhxndepo",varid)
         IF (ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,5) "nhxndepo", bioname
         ENDIF
         DO irec = 1, nrec_ndep
            ierr = nf_fread(dustmp(START_2D_ARRAY,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. lwp ) THEN
               WRITE(numout,6) "nhxndepo", irec
            ENDIF
         END DO
         !
         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,'(6x,A,1x,I4)') &
#ifdef MPI
         &                   'TRCINI_PISCES -- Read Ammoniun deposition ', mynode
#else
         &                   'TRCINI_PISCES -- Read Ammoniun deposition '
#endif

         DO irec = 1, nrec_ndep
            DO jj= 1, LOCALMM
               DO ji =1, LOCALLM
                  nh4depmo(ji,jj,irec) = MAX( 0., dustmp(ji,jj,irec) ) / rno3 / mMass_N 
               END DO
            END DO
         END DO
         !
#else         
         ALLOCATE( no3depmo(GLOBAL_2D_ARRAY,nrec_ndep) )
         ierr = nf_inq_varid (ncid,"ndep",varid)
         IF (ierr .NE. nf_noerr .AND. lwp ) THEN
            WRITE(numout,5) "ndep", bioname
         ENDIF
         !
         DO irec = 1, nrec_ndep
            ierr = nf_fread(dustmp(START_2D_ARRAY,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. lwp ) THEN
               WRITE(numout,6) "ndep", irec
            ENDIF
         END DO
         !
         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,'(6x,A,1x,I4)') &
#ifdef MPI
         &                   'TRCINI_PISCES -- Read Nitrate deposition ', mynode
#else
         &                   'TRCINI_PISCES -- Read Nitrate deposition '
#endif

         DO irec = 1, nrec_ndep
            DO jj = 1, LOCALMM
               DO ji = 1, LOCALLM
                  no3depmo(ji,jj,irec) = MAX( 0., dustmp(ji,jj,irec) ) / rno3 / mMass_N 
               END DO
            END DO
         END DO
         !
#endif
         ierr = nf_close(ncid)
         DEALLOCATE( dustmp )

      ENDIF

  4   FORMAT(/,' TRCINI_PISCES - unable to open forcing netCDF ',1x,A)
  5   FORMAT(/,' TRCINI_PISCES - unable to find forcing variable: ',A, &
     &                               /,14x,'in forcing netCDF  ',A)
  6   FORMAT(/,' TRCINI_PISCES - error while reading variable: ',A,2x, &
     &                                           ' at TIME index = ',i4)
#endif

   END SUBROUTINE p4z_bc_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sbc                        ! Empty routine
   END SUBROUTINE p4z_sbc
#endif


   !!======================================================================
END MODULE p4zbc
