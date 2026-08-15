#include "cppdefs.h"

MODULE p4zsms
   !!======================================================================
   !!                         ***  MODULE p4zsms  ***
   !! TOP :   PISCES Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_sms        : Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zbio          ! Biological model
   USE p4zche          ! Chemical model
   USE p4zlys          ! Calcite saturation
   USE p4zflx          ! Gas exchange
   USE p4zbc           ! External source of nutrients
   USE p4zdiaz         !  Diazotrophy
   USE p4zsed          ! Sedimentation
   USE p4zint          ! time interpolation
   USE p4zrem          ! remineralisation
   USE iom             ! I/O manager
!   USE trd_oce         ! Ocean trends variables
!   USE trdtrc          ! TOP trends variables
   USE sedmodel        ! Sediment model
!   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sms_init   ! called in trcini_pisces.F90
   PUBLIC   p4z_sms        ! called in trcsms_pisces.F90

   INTEGER ::    numco2, numnut, numnit      ! logical unit for co2 budget
   REAL(wp) ::   xfact, xfact3

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xnegtr     ! Array used to indicate negative tracer values
   LOGICAL :: l_budget

   INTEGER  ::  rstph, rstfe, rstszn, rstszd, rstszp
   INTEGER  ::  rstthet, rstxksi, rstxksim, rstpisstep, rstpistime, rstpistime2
   INTEGER  ::  rstpoc, rstgoc

   !! * Substitutions
#  include "ocean2pisces.h90"   
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsms.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sms( kt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sms  ***
      !!
      !! ** Purpose :   Managment of the call to Biological sources and sinks 
      !!              routines of PISCES bio-model
      !!
      !! ** Method  : - calls the various SMS subroutines
      !!              - calls the sediment module (if ln_sediment)  
      !!              - several calls of bio and sed (possible time-splitting)
      !!              - handles the potential negative concentrations (xnegtr)
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt              ! ocean time-step index      
      INTEGER, INTENT( in ) ::   Kbb, Kmm, Krhs  ! time level index
      !!
      INTEGER ::   ji, jj, jk, jnt, jn, jl, ilc, ikty
      LOGICAL :: ll_dmp

      REAL(wp) ::  ztra
      CHARACTER (len=25) :: charout
      REAL(wp), ALLOCATABLE, DIMENSION(:,:    ) :: zw2d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:  ) :: zw3d

      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_sms')
      !
      IF( kt == nittrc000 ) THEN
        !
        ALLOCATE( xnegtr(A2D(0),jpk) )
        !
        IF( .NOT. ln_rsttr ) THEN
            CALL p4z_che( Kbb, Kmm )                  ! initialize the chemical constants
            CALL ahini_for_at( hi, Kbb )              !  set PH at kt=nit000
            t_oce_co2_flx_cum = 0._wp
        ELSE
#if defined NEMO               
            CALL p4z_rst( nittrc000, Kbb, Kmm,  'READ' )  !* read or initialize all required fields
#else
            CALL p4z_rst_read( Kbb, Kmm )
#endif            
        ENDIF
        !
      ENDIF
      !
      ll_dmp  = ( ln_pisdmp .OR. ln_pisdmp_alk ) .AND. ( ln_rsttr .OR. kt /= 1 ) 
      ikty = nyear_len(1) * rday / NINT(rn_Dt)    ! time-step at the end of a year 
      ll_dmp = ll_dmp .AND.  MOD( kt-1, ikty ) == 0 
      !
      IF( ll_dmp )    CALL p4z_dmp( kt, Kbb, Kmm )      ! Relaxation of some tracers
      !
      rfact = rDt_trc  ! time step of PISCES
      !
      IF( ( ln_top_euler .AND. kt == nittrc000 )  .OR. ( .NOT.ln_top_euler .AND. kt <= nittrc000 + 1 ) ) THEN
         rfactr  = 1. / rfact  ! inverse of the time step
         rfact2  = rfact / REAL( nrdttrc, wp )  ! time step of the biological SMS
         rfact2r = 1. / rfact2  ! Inverse of the biological time step
         xstep = rfact2 / rday         ! Time step duration for biology relative to a day
         xfact = 1.e+3 * rfact2r
         IF(lwp) WRITE(numout,*) 
         IF(lwp) WRITE(numout,*) '    Passive Tracer  time step    rfact  = ', rfact, ' rn_Dt = ', rn_Dt
         IF(lwp) write(numout,*) '    PISCES  Biology time step    rfact2 = ', rfact2
         IF(lwp) WRITE(numout,*)
      ENDIF
      !
      IF( ll_bc )    CALL p4z_bc( kt, Kbb, Kmm, Krhs )   ! external sources of nutrients 
      !
      CALL p4z_che(     Kbb, Kmm       ) ! computation of chemical constants
      CALL p4z_int( kt, Kbb, Kmm       ) ! computation of various rates for biogeochemistry
      !
      DO jnt = 1, nrdttrc          ! Potential time splitting if requested
         !
         IF( ln_bio ) CALL p4z_bio( kt, jnt, Kbb, Kmm, Krhs )   ! Biology
         IF( ln_p2z ) THEN
            IF( ln_lys ) CALL p2z_lys( kt, jnt, Kbb, Kmm, Krhs )   ! Compute CaCO3 saturation
         ELSE
            IF( ln_lys ) CALL p4z_lys( kt, jnt, Kbb,      Krhs )   ! Compute CaCO3 saturation
         ENDIF
         IF( .NOT. lk_sed .AND. ln_sed ) CALL p4z_sed( kt, jnt, Kbb, Kmm, Krhs )   ! Surface and Bottom boundary conditions

         IF( ln_flx ) CALL p4z_flx( kt, jnt, Kbb, Kmm, Krhs )   ! Compute surface fluxes
         !
         ! Handling of the negative concentrations
         ! The biological SMS may generate negative concentrations
         ! Trends are tested at each grid cell. If a negative concentrations 
         ! is created at a grid cell, all the sources and sinks at that grid 
         ! cell are scale to avoid that negative concentration. This approach 
         ! is quite simplistic but it conserves mass.
         ! ------------------------------------------------------------------
         xnegtr(:,:,:) = 1.e0
         DO jn = jp_pcs0, jp_pcs1
            DO_3D( 0, 0, 0, 0, 1, jpk)
               IF( ( tr(ji,jj,jk,jn,Kbb) + tr(ji,jj,jk,jn,Krhs) ) < 0.e0 ) THEN
                  ztra = ABS( tr(ji,jj,jk,jn,Kbb) ) &
                    &  / ( ABS( tr(ji,jj,jk,jn,Krhs) ) + rtrn )
                  xnegtr(ji,jj,jk) = MIN( xnegtr(ji,jj,jk),  ztra )
               ENDIF
            END_3D
         END DO
         !                                ! where at least 1 tracer concentration becomes negative
         !                                ! 
         ! Concentrations are updated
         DO jn = jp_pcs0, jp_pcs1
            DO_3D( 0, 0, 0, 0, 1, jpk)
               tr(ji,jj,jk,jn,Kbb) = tr(ji,jj,jk,jn,Kbb) &
                       &           + xnegtr(ji,jj,jk) * tr(ji,jj,jk,jn,Krhs)
            END_3D
         END DO
        !
        IF(  iom_use( 'INTdtAlk' ) .OR. iom_use( 'INTdtDIC' ) .OR. iom_use( 'INTdtFer' ) .OR.  &
          &  iom_use( 'INTdtDIN' ) .OR. iom_use( 'INTdtDIP' ) .OR. iom_use( 'INTdtSil' ) )  THEN
          !
           ALLOCATE( zw3d(A2D(0),jpk), zw2d(A2D(0)) )
        ENDIF
        IF ( iom_use( 'INTdtAlk' ) ) THEN
           DO_3D( 0, 0, 0, 0, 1, jpkm1)
              zw3d(ji,jj,jk) = xnegtr(ji,jj,jk) * xfact * e3t(ji,jj,jk,Kmm) * tmask(ji,jj,jk)
           END_3D
           !
           zw2d(:,:) = 0.
           DO jk = 1, jpkm1
              DO_2D( 0, 0, 0, 0 )
                 zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk) *  tr(ji,jj,jk,jptal,Krhs) 
              END_2D
           ENDDO
           CALL iom_put( 'INTdtAlk', zw2d )
        ENDIF
          !
        IF ( iom_use( 'INTdtDIC' ) ) THEN
           zw2d(:,:) = 0.
           DO jk = 1, jpkm1
              DO_2D( 0, 0, 0, 0 )
                 zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk) *  tr(ji,jj,jk,jpdic,Krhs) 
              END_2D
           ENDDO
           CALL iom_put( 'INTdtDIC', zw2d )
        ENDIF
          !
        IF ( iom_use( 'INTdtDIN' ) ) THEN
           zw2d(:,:) = 0.
           DO jk = 1, jpkm1
              DO_2D( 0, 0, 0, 0 )
                 ztra = tr(ji,jj,jk,jpno3,Krhs) + tr(ji,jj,jk,jpnh4,Krhs)
                 zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk) * rno3 * ztra 
              END_2D
           ENDDO
           CALL iom_put( 'INTdtDIN', zw2d )
        ENDIF
          !
        IF ( iom_use( 'INTdtDIP' ) .AND. .NOT.ln_p2z ) THEN
           zw2d(:,:) = 0.
           DO jk = 1, jpkm1
              DO_2D( 0, 0, 0, 0 )
                 zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk) * po4r * tr(ji,jj,jk,jppo4,Krhs) 
              END_2D
           ENDDO
           CALL iom_put( 'INTdtDIP', zw2d )
        ENDIF
           !
        IF ( iom_use( 'INTdtFer' ) ) THEN
           zw2d(:,:) = 0.
           DO jk = 1, jpkm1
              DO_2D( 0, 0, 0, 0 )
                 zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk) *  tr(ji,jj,jk,jpfer,Krhs) 
              END_2D
           ENDDO
           CALL iom_put( 'INTdtFer', zw2d )
        ENDIF
          !
        IF ( iom_use( 'INTdtSil' ) .AND. .NOT.ln_p2z ) THEN
           zw2d(:,:) = 0.
           DO jk = 1, jpkm1
              DO_2D( 0, 0, 0, 0 )
                 zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk) *  tr(ji,jj,jk,jpsil,Krhs) 
              END_2D
           ENDDO
           CALL iom_put( 'INTdtSil', zw2d )
           !
        ENDIF 

        IF(  iom_use( 'INTdtAlk' ) .OR. iom_use( 'INTdtDIC' ) .OR. iom_use( 'INTdtFer' ) .OR.  &
          &  iom_use( 'INTdtDIN' ) .OR. iom_use( 'INTdtDIP' ) .OR. iom_use( 'INTdtSil' ) )  THEN
          DEALLOCATE( zw3d, zw2d )
        ENDIF
        !
        ! Trends are are reset to 0
        DO jn = jp_pcs0, jp_pcs1
           DO_3D( 0, 0, 0, 0, 1, jpkm1)
             tr(ji,jj,jk,jn,Kmm) = tr(ji,jj,jk,jn,Kbb) 
             tr(ji,jj,jk,jn,Krhs) = 0._wp
           END_3D
        END DO
        !
      END DO
      !
      !
      ! If ln_sediment is set to .true. then the sediment module is called
      IF( ln_sediment ) THEN 
         !
         CALL sed_model( kt, Kbb, Kmm, Krhs )     !  Main program of Sediment model
         !
      ENDIF
      !
#ifdef NEMO      
      IF( lrst_trc )  CALL p4z_rst( kt, Kbb, Kmm,  'WRITE' )           !* Write PISCES informations in restart file 
#else
      ilc = 1+iic-nit000 ! number of time step since restart
      IF( iic > nit000 ) THEN
         IF( MOD( ilc-1, nitrst ) == 0  &
#ifdef EXACT_RESTART
     &                      .OR. MOD(ilc,nitrst) == 0  &
#endif
     &                      )  THEN
#if defined PISCES
! need the CPP key to avoid compilation error
            nrecpisrst = nrecpisrst + 1
            CALL p4z_rst_wri
#endif            
       ENDIF
     ENDIF
#endif      
      !
      IF( lk_iomput )  CALL p4z_budget( kt, Kmm ) ! Budget checking

      IF( lwm .AND. kt == nittrc000    )  CALL FLUSH( numonp )         ! flush output namelist PISCES
      !
      IF( ln_timing )  CALL timing_stop('p4z_sms')
      !
   END SUBROUTINE p4z_sms


   SUBROUTINE p4z_sms_init
      !!----------------------------------------------------------------------
      !!                     ***  p4z_sms_init  ***  
      !!
      !! ** Purpose :   read the general PISCES namelist
      !!
      !! ** input   :   file 'namelist_pisces' containing the following
      !!                namelist: nampisbio
      !!----------------------------------------------------------------------
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/nampisbio/ nrdttrc, wsbio, xkmort, feratz, feratm, wsbio2, wsbio2max,    &
         &                wsbio2scale, ldocp, ldocz, lthet, no3rat3, po4rat3
         !
      NAMELIST/nampisdmp/ ln_pisdmp, ln_pisdmp_alk
#if ! defined NEMO               
      NAMELIST/nampisdbg/ ln_bio, ln_lys, ln_sed, ln_flx, ln_sink, &
         &                ln_fechem, ln_micro, ln_meso, ln_mort, &
         &                ln_prod, ln_agg, ln_rem, ln_poc, ln_diaz
      NAMELIST/nampisrst/cn_pisrst_indir, cn_pisrst_outdir, cn_pisrst_in, cn_pisrst_out 
#endif      
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_sms_init : PISCES initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF

      READ_NML_REF(numnatp,nampisbio)
      READ_NML_CFG(numnatp,nampisbio)
      IF(lwm) WRITE( numonp, nampisbio )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : nampisbio'
         WRITE(numout,*) '      frequency for the biology                 nrdttrc     =', nrdttrc
         WRITE(numout,*) '      POC sinking speed                         wsbio       =', wsbio
         WRITE(numout,*) '      half saturation constant for mortality    xkmort      =', xkmort 
         IF( ln_p5z ) THEN
            WRITE(numout,*) '      N/C in zooplankton                     no3rat3     =', no3rat3
            WRITE(numout,*) '      P/C in zooplankton                     po4rat3     =', po4rat3
         ENDIF
         WRITE(numout,*) '      Fe/C in microzooplankton                  feratz      =', feratz
         IF( .NOT. ln_p2z ) THEN
            WRITE(numout,*) '      Fe/C in mesozooplankton                   feratm      =', feratm
         ENDIF
         WRITE(numout,*) '      Big particles sinking speed               wsbio2      =', wsbio2
         WRITE(numout,*) '      Big particles maximum sinking speed       wsbio2max   =', wsbio2max
         WRITE(numout,*) '      Big particles sinking speed length scale  wsbio2scale =', wsbio2scale
         IF( ln_ligand ) THEN
            WRITE(numout,*) '      Phyto ligand production per unit doc           ldocp  =', ldocp
            WRITE(numout,*) '      Zoo ligand production per unit doc             ldocz  =', ldocz
            WRITE(numout,*) '      Proportional loss of ligands due to Fe uptake  lthet  =', lthet
         ENDIF
      ENDIF


      READ_NML_REF(numnatp,nampisdmp)
      READ_NML_CFG(numnatp,nampisdmp)
      IF(lwm) WRITE( numonp, nampisdmp )
      !
      ln_pisdmp     = .FALSE.
      ln_pisdmp_alk = .FALSE.
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist : nampisdmp --- relaxation to GLODAP'
         WRITE(numout,*) '      Relaxation of tracer to glodap mean value   ln_pisdmp =', ln_pisdmp
         WRITE(numout,*) '      Relaxation of alkalinity to glodap mean value   ln_pisdmp_alk  =', ln_pisdmp_alk
      ENDIF
      !
#if ! defined NEMO               
      READ_NML_REF(numnatp,nampisdbg)
      READ_NML_CFG(numnatp,nampisdbg)
      IF(lwm) WRITE( numonp, nampisdbg )
      !

      READ_NML_REF(numnatp,nampisrst)
      READ_NML_CFG(numnatp,nampisrst)

      cn_pisrst_in  = TRIM( cn_pisrst_indir)//'/'//TRIM(cn_pisrst_in)
      cn_pisrst_out = TRIM( cn_pisrst_outdir)//'/'//TRIM(cn_pisrst_out)

      IF (lwp) THEN
         WRITE(numout,*) ' namelist  nampisrst '
         WRITE(numout,*) '  Name of input restart file if needed = ', TRIM( cn_pisrst_in )
         WRITE(numout,*) '  Name of output restart file if needed = ', TRIM( cn_pisrst_out )
         WRITE(numout,*) ' '
      ENDIF
#endif    
      !
      ncidpisrst = -1
      nrecpisrst = 0
      !
   END SUBROUTINE p4z_sms_init

# ifdef NEMO
   SUBROUTINE p4z_rst( kt, Kbb, Kmm, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE p4z_rst  ***
      !!
      !!  ** Purpose : Read or write specific PISCES variables in restart file:
      !!
      !!  WRITE(READ) mode:
      !!       kt        : number of time step since the begining of the experiment at the
      !!                   end of the current(previous) run
      !!---------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      INTEGER         , INTENT(in) ::   Kbb, Kmm   ! time level indices
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      !!---------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN
         !
         ! Read the specific variable of PISCES
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' p4z_rst : Read specific variables from pisces model '
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
         ! 
         ! Read the pH. If not in the restart file, then it is initialized from
         ! the initial conditions
         IF( iom_varid( numrtr, 'PH', ldstop = .FALSE. ) > 0 ) THEN
            CALL iom_get( numrtr, jpdom_auto, 'PH' , hi(:,:,:)  )
         ELSE
            CALL p4z_che( Kbb, Kmm )                  ! initialize the chemical constants
            CALL ahini_for_at( hi, Kbb )
         ENDIF

         IF( ln_p2z ) THEN
            IF( iom_varid( numrtr, 'Thetanano', ldstop = .FALSE. ) > 0 ) THEN
               CALL iom_get( numrtr, jpdom_auto, 'Thetanano' , thetanano(:,:,:)  )
            ELSE
               thetanano(:,:,:) = 1.0 / 55.0
            ENDIF
         ENDIF

         IF( .NOT. ln_p2z ) THEN
            CALL iom_get( numrtr, jpdom_auto, 'Silicalim', xksi(:,:) )
   
            ! Read the Si half saturation constant and the maximum Silica concentration
            IF( iom_varid( numrtr, 'Silicamax', ldstop = .FALSE. ) > 0 ) THEN
               CALL iom_get( numrtr, jpdom_auto, 'Silicamax' , xksimax(:,:)  )
            ELSE
               xksimax(:,:) = xksi(:,:)
            ENDIF
         ENDIF

         ! Read the Fe3 consumption term by phytoplankton
         IF( iom_varid( numrtr, 'Consfe3', ldstop = .FALSE. ) > 0 ) THEN
            CALL iom_get( numrtr, jpdom_auto, 'Consfe3' , consfe3(:,:,:)  )
         ELSE
            consfe3(:,:,:) = 0._wp
         ENDIF

        ! Read mean remineralisation rate
         ll_poc_lab = .FALSE.
         IF( iom_varid( numrtr, 'remintpoc', ldstop = .FALSE. ) > 0 ) THEN
            CALL iom_get( numrtr, jpdom_auto, 'remintpoc' , remintpoc(:,:,:)  )
            IF( ln_p4z .OR. ln_p5z ) &
            CALL iom_get( numrtr, jpdom_auto, 'remintgoc' , remintgoc(:,:,:)  )
         ELSE
                                     remintpoc(:,:,:) = 0.03_wp
            IF( ln_p4z .OR. ln_p5z ) remintgoc(:,:,:) = 0.03_wp
         ENDIF

         ! Read the cumulative total flux. If not in the restart file, it is set to 0          
         IF( iom_varid( numrtr, 'tcflxcum', ldstop = .FALSE. ) > 0 ) THEN  ! cumulative total flux of carbon
            CALL iom_get( numrtr, 'tcflxcum' , t_oce_co2_flx_cum  )
         ELSE
            t_oce_co2_flx_cum = 0._wp
         ENDIF
         !
         ! PISCES size proxy
         ! Read the size of the different phytoplankton groups
         ! If not in the restart file, they are set to 1
         IF( iom_varid( numrtr, 'sizen', ldstop = .FALSE. ) > 0 ) THEN
            CALL iom_get( numrtr, jpdom_auto, 'sizen' , sizen(:,:,:)  )
               sizen(:,:,:) = MAX( 1.0, sizen(:,:,:) )
            IF( ln_p4z .OR. ln_p5z )THEN
               CALL iom_get( numrtr, jpdom_auto, 'sized' , sized(:,:,:)  )
               sized(:,:,:) = MAX( 1.0, sized(:,:,:) )
            ENDIF
            IF( ln_p5z ) THEN
               CALL iom_get( numrtr, jpdom_auto, 'sizep' , sizep(:,:,:)  )
               sizep(:,:,:) = MAX( 1.0, sizep(:,:,:) )
            ENDIF
         ELSE
                                     sizen(:,:,:) = 1.
            IF( ln_p4z .OR. ln_p5z ) sized(:,:,:) = 1.
            IF( ln_p5z )             sizep(:,:,:) = 1.
         ENDIF
        !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN
         ! write the specific variables of PISCES
         IF( kt == nitrst ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'p4z_rst : write pisces restart file  kt =', kt
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
         CALL iom_rstput( kt, nitrst, numrtw, 'PH', hi(:,:,:)           )
         CALL iom_rstput( kt, nitrst, numrtw, 'sizen', sizen(:,:,:) )  ! Size of nanophytoplankton
         CALL iom_rstput( kt, nitrst, numrtw, 'tcflxcum', t_oce_co2_flx_cum )
         IF( ln_p2z ) THEN
            CALL iom_rstput( kt, nitrst, numrtw, 'Thetanano' , thetanano(:,:,:)  )
         ENDIF

         CALL iom_rstput( kt, nitrst, numrtw, 'Consfe3', consfe3(:,:,:) ) ! Si max concentration
         CALL iom_rstput( kt, nitrst, numrtw, 'remintpoc', remintpoc(:,:,:) ) ! Mean remineralisation rate of POC
         IF ( ln_p4z .OR. ln_p5z ) THEN
            CALL iom_rstput( kt, nitrst, numrtw, 'Silicalim', xksi(:,:)    )
            CALL iom_rstput( kt, nitrst, numrtw, 'Silicamax', xksimax(:,:) )
            CALL iom_rstput( kt, nitrst, numrtw, 'remintgoc', remintgoc(:,:,:) ) ! Mean remineralisation rate of GOC
            CALL iom_rstput( kt, nitrst, numrtw, 'sized', sized(:,:,:) )  ! Size of diatoms
         ENDIF
         IF( ln_p5z ) CALL iom_rstput( kt, nitrst, numrtw, 'sizep', sizep(:,:,:) )  ! Size of picophytoplankton
      ENDIF
      !
   END SUBROUTINE p4z_rst
#endif

   SUBROUTINE p4z_dmp( kt, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                    ***  p4z_dmp  ***
      !!
      !! ** purpose  : Relaxation of the total budget of some elements
      !!               This routine avoids the model to drift far from the 
      !!               observed content in various elements
      !!               Elements that may be relaxed : Alk, P, N, Si
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in )  ::     kt            ! time step
      INTEGER, INTENT( in )  ::     Kbb, Kmm      ! time level indices
      !
      REAL(wp) ::  alkmean = 2426.     ! mean value of alkalinity ( Glodap ; for Goyet 2391. )
      REAL(wp) ::  po4mean = 2.174     ! mean value of phosphate
      REAL(wp) ::  no3mean = 31.00     ! mean value of nitrate
      REAL(wp) ::  silmean = 90.33     ! mean value of silicate
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zarea, zalksumn, zpo4sumn, zno3sumn, zsilsumn
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------

      IF(lwp)  WRITE(numout,*)
      IF(lwp)  WRITE(numout,*) ' p4z_dmp : Restoring of nutrients at time-step kt = ', kt
      IF(lwp)  WRITE(numout,*)

      IF( cn_cfg == "ORCA" .OR. cn_cfg == "orca") THEN

       ALLOCATE( zw3d(A2D(0),jpk) )  ;  zw3d(A2D(0),jpk) = 0._wp
       ! set total alkalinity, phosphate, nitrate & silicate
       zarea          = 1._wp / glob_sum( 'p4zsms', cvol(:,:,:) ) * 1e6              

       IF( ln_pisdmp ) THEN
          !      
          DO_3D( 0, 0, 0, 0, 1, jpk)
             zw3d(ji,jj,jk) = tr(ji,jj,jk,jpno3,Kmm) * cvol(ji,jj,jk)
          END_3D
          zno3sumn = glob_sum( 'p4zsms', zw3d(:,:,:) ) * zarea * rno3
          ! Correct the trn mean content of NO3
          IF(lwp) WRITE(numout,*) '       NO3N  mean : ', zno3sumn
          DO_3D( 0, 0, 0, 0, 1, jpk)
               tr(ji,jj,jk,jpno3,Kmm) = tr(ji,jj,jk,jpno3,Kmm) * no3mean / zno3sumn
          END_3D
          !
          IF ( ln_p4z .OR. ln_p5z ) THEN
             DO_3D( 0, 0, 0, 0, 1, jpk)
                zw3d(ji,jj,jk) = tr(ji,jj,jk,jppo4,Kmm) * cvol(ji,jj,jk)
             END_3D
             zpo4sumn = glob_sum( 'p4zsms', zw3d(:,:,:)  ) * zarea * po4r
             ! Correct the trn mean content of PO4
             IF(lwp) WRITE(numout,*) '       PO4N  mean : ', zpo4sumn
             DO_3D( 0, 0, 0, 0, 1, jpk)
                tr(ji,jj,jk,jppo4,Kmm) = tr(ji,jj,jk,jppo4,Kmm) * po4mean / zpo4sumn
             END_3D

             DO_3D( 0, 0, 0, 0, 1, jpk)
              zw3d(ji,jj,jk) = tr(ji,jj,jk,jpsil,Kmm) * cvol(ji,jj,jk)
             END_3D
             zsilsumn = glob_sum( 'p4zsms', zw3d(:,:,:)  ) * zarea
             ! Correct the trn mean content of SiO3
             IF(lwp) WRITE(numout,*) '       SiO3N mean : ', zsilsumn
             DO_3D( 0, 0, 0, 0, 1, jpk)
                tr(ji,jj,jk,jpsil,Kmm) = MIN( 400.e-6,tr(ji,jj,jk,jpsil,Kmm) &
                      &               * silmean / zsilsumn )
             END_3D
          ENDIF
       ENDIF   
       !
       IF( ln_pisdmp_alk ) THEN
       
          DO_3D( 0, 0, 0, 0, 1, jpk)
             zw3d(ji,jj,jk) = tr(ji,jj,jk,jptal,Kmm) * cvol(ji,jj,jk)
          END_3D
          zalksumn = glob_sum( 'p4zsms', zw3d(:,:,:) ) * zarea
 
          ! Correct the trn mean content of alkalinity
          IF(lwp) WRITE(numout,*) '       TALKN mean : ', zalksumn
          DO_3D( 0, 0, 0, 0, 1, jpk)
             tr(ji,jj,jk,jptal,Kmm) = tr(ji,jj,jk,jptal,Kmm) * alkmean / zalksumn
          END_3D
       ENDIF   

       DEALLOCATE( zw3d )  
      !
      ENDIF
      !
   END SUBROUTINE p4z_dmp


   SUBROUTINE p4z_budget( kt, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_chk_mass  ***
      !!
      !! ** Purpose :  Mass conservation check 
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      INTEGER, INTENT( in ) ::   Kmm     ! time level indices
      REAL(wp)             ::  zdenittot, znitrpottot
      REAL(wp)              :: zalkbudget, zno3budget, zsilbudget, zferbudget, zpo4budget 
      INTEGER :: ji, jj, jk
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!----------------------------------------------------------------------
      !
      IF( kt == nittrc000 ) THEN 
         xfact3 = 1.e+3 * rfact2r * rno3   ! conversion molC/l/kt ----> molN/m3/s
         l_budget =   iom_use( "pno3tot" ) .OR. iom_use( "ppo4tot" ) .OR. iom_use( "psiltot" )  &
           &     .OR. iom_use( "palktot" ) .OR. iom_use( "pfertot" ) .OR. iom_use( "tdenit"  )  &
           &     .OR. iom_use( "tnfix" ) 
        
      ENDIF

      ! Compute the budget of NO3
      IF( l_budget ) THEN
        !
        ALLOCATE( zw3d(A2D(0),jpk) )  ;  zw3d(:,:,jpk) = 0._wp
        !
        IF( ln_p2z ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jk)  =  ( tr(ji,jj,jk,jpno3,Kmm) + tr(ji,jj,jk,jpphy,Kmm) &
               &             +   tr(ji,jj,jk,jppoc,Kmm) + tr(ji,jj,jk,jpdoc,Kmm) &
               &             +   tr(ji,jj,jk,jpzoo,Kmm) ) * cvol(ji,jj,jk)
            END_3D
        ELSE IF( ln_p4z ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jk)  =  ( tr(ji,jj,jk,jpno3,Kmm) + tr(ji,jj,jk,jpnh4,Kmm) &
               &             +   tr(ji,jj,jk,jpphy,Kmm) + tr(ji,jj,jk,jpdia,Kmm) &
               &             +   tr(ji,jj,jk,jppoc,Kmm) + tr(ji,jj,jk,jpgoc,Kmm)  &
               &             +   tr(ji,jj,jk,jpdoc,Kmm)  &        
               &             +   tr(ji,jj,jk,jpzoo,Kmm) + tr(ji,jj,jk,jpmes,Kmm)  ) &
               &                 * cvol(ji,jj,jk)
            END_3D
        ELSE
            DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jk) =    ( tr(ji,jj,jk,jpno3,Kmm) + tr(ji,jj,jk,jpnh4,Kmm) &
               &             +   tr(ji,jj,jk,jpnph,Kmm)   &
               &             +   tr(ji,jj,jk,jpndi,Kmm) + tr(ji,jj,jk,jpnpi,Kmm) & 
               &             +   tr(ji,jj,jk,jppon,Kmm) + tr(ji,jj,jk,jpgon,Kmm) &
               &             +   tr(ji,jj,jk,jpdon,Kmm)   &
               &             + ( tr(ji,jj,jk,jpzoo,Kmm) + tr(ji,jj,jk,jpmes,Kmm) ) &
               &                * no3rat3 ) * cvol(ji,jj,jk)
            END_3D
        ENDIF
        !
        zno3budget = glob_sum( 'p4zsms', zw3d(:,:,:)  )  
        zno3budget = zno3budget / areatot
        CALL iom_put( "pno3tot", zno3budget )
        IF( .NOT. ln_p2z ) THEN
           !
           ! Compute the budget of PO4
           IF( ln_p4z ) THEN
               DO_3D( 0, 0, 0, 0, 1, jpk)
                  zw3d(ji,jj,jk) =    ( tr(ji,jj,jk,jppo4,Kmm)              &
                     &    + tr(ji,jj,jk,jpphy,Kmm) + tr(ji,jj,jk,jpdia,Kmm) &
                     &    + tr(ji,jj,jk,jppoc,Kmm) + tr(ji,jj,jk,jpgoc,Kmm) &
                     &    + tr(ji,jj,jk,jpdoc,Kmm)  &        
                     &    + tr(ji,jj,jk,jpzoo,Kmm) + tr(ji,jj,jk,jpmes,Kmm) ) &
                     &             * cvol(ji,jj,jk)
               END_3D
            ELSE
               DO_3D( 0, 0, 0, 0, 1, jpk)
                  zw3d(ji,jj,jk) = ( tr(ji,jj,jk,jppo4,Kmm) + tr(ji,jj,jk,jppph,Kmm) &
                     &  + tr(ji,jj,jk,jppdi,Kmm) + tr(ji,jj,jk,jpppi,Kmm) & 
                     &  + tr(ji,jj,jk,jppop,Kmm) + tr(ji,jj,jk,jpgop,Kmm) &
                     &  + tr(ji,jj,jk,jpdop,Kmm)   &
                     &  + (  tr(ji,jj,jk,jpzoo,Kmm) + tr(ji,jj,jk,jpmes,Kmm) ) &
                     &  * po4rat3 ) * cvol(ji,jj,jk)
               END_3D
            ENDIF
            !
            zpo4budget = glob_sum( 'p4zsms', zw3d(:,:,:)  )  
            zpo4budget = zpo4budget / areatot
            CALL iom_put( "ppo4tot", zpo4budget )
            !
            ! Compute the budget of SiO3
            DO_3D( 0, 0, 0, 0, 1, jpk)
               zw3d(ji,jj,jk) =  ( tr(ji,jj,jk,jpsil,Kmm) + tr(ji,jj,jk,jpgsi,Kmm) &
                       &        + tr(ji,jj,jk,jpdsi,Kmm) ) * cvol(ji,jj,jk)
            END_3D
            !
            zsilbudget = glob_sum( 'p4zsms', zw3d(:,:,:)  )  
            zsilbudget = zsilbudget / areatot
            CALL iom_put( "psiltot", zsilbudget )
        ENDIF
        !
        ! Compute the budget of Iron
        IF( ln_p2z ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk)
               zw3d(ji,jj,jk) =   ( tr(ji,jj,jk,jpfer,Kmm) &
                  &    + tr(ji,jj,jk,jpphy,Kmm) * feratz  &
                  &    + tr(ji,jj,jk,jppoc,Kmm) *feratz    &
                  &    + tr(ji,jj,jk,jpzoo,Kmm) * feratz ) &
                  &       * cvol(ji,jj,jk)
            END_3D
         ELSE IF( ln_p4z ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk)
               zw3d(ji,jj,jk) = ( tr(ji,jj,jk,jpfer,Kmm) + tr(ji,jj,jk,jpnfe,Kmm) &
                  &   + tr(ji,jj,jk,jpdfe,Kmm)   &
                  &   + tr(ji,jj,jk,jpbfe,Kmm) + tr(ji,jj,jk,jpsfe,Kmm)  &
                  &   + tr(ji,jj,jk,jpzoo,Kmm) * feratz &
                  &   + tr(ji,jj,jk,jpmes,Kmm) * feratm ) &
                  &                 * cvol(ji,jj,jk)
            END_3D
         ELSE
            DO_3D( 0, 0, 0, 0, 1, jpk)
               zw3d(ji,jj,jk) =   ( tr(ji,jj,jk,jpfer,Kmm) + tr(ji,jj,jk,jpnfe,Kmm) &
                  &    + tr(ji,jj,jk,jpdfe,Kmm) + tr(ji,jj,jk,jppfe,Kmm)  &
                  &    + tr(ji,jj,jk,jpbfe,Kmm) + tr(ji,jj,jk,jpsfe,Kmm)   &
                  &    + tr(ji,jj,jk,jpzoo,Kmm) * feratz &
                  &    + tr(ji,jj,jk,jpmes,Kmm) * feratm ) &
                  &                * cvol(ji,jj,jk)
            END_3D
         ENDIF
         !
         zferbudget = glob_sum( 'p4zsms', zw3d(:,:,:) )
         zferbudget = zferbudget / areatot
         CALL iom_put( "pfertot", zferbudget )
        !
        ! Compute the budget of total alkalinity
        IF( ln_p2z ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk)
               zw3d(ji,jj,jk) =  ( tr(ji,jj,jk,jpno3,Kmm) * rno3  &
               &                 + tr(ji,jj,jk,jptal,Kmm) ) &
                       &         * cvol(ji,jj,jk)             
            END_3D
         ELSE
            DO_3D( 0, 0, 0, 0, 1, jpk)
               zw3d(ji,jj,jk) =  ( tr(ji,jj,jk,jpno3,Kmm) * rno3 &
                 &              +  tr(ji,jj,jk,jptal,Kmm)   &
                 &              +  tr(ji,jj,jk,jpcal,Kmm) * 2. &
                 &              -  tr(ji,jj,jk,jpnh4,Kmm) * rno3 ) &
                 &           * cvol(ji,jj,jk)
            END_3D
         ENDIF
         !
         zalkbudget = glob_sum( 'p4zsms', zw3d(:,:,:)  )         !
         zalkbudget = zalkbudget / areatot
         CALL iom_put( "palktot", zalkbudget )
         !
         ! Global budget of N SMS : nitrogen fixation by the diazotrophs
         DO_3D( 0, 0, 0, 0, 1, jpkm1)
           zw3d(ji,jj,jk) =  nitrpot(ji,jj,jk) * nitrfix * cvol(ji,jj,jk)
         END_3D
         znitrpottot = glob_sum ( 'p4zsms',  zw3d) 
         CALL iom_put( "tnfix"  , znitrpottot * xfact3 )  
         !
         ! Global budget of N SMS : denitrification in the water column and in the sediment
         DO_2D( 0, 0, 0, 0 )
            zw3d(ji,jj,1) =  denitr(ji,jj,1) *  rdenit * xnegtr(ji,jj,1) * cvol(ji,jj,1) &
                &           + sdenit(ji,jj)   *  e1e2t(ji,jj) * tmask(ji,jj,1)
         END_2D
         DO_3D( 0, 0, 0, 0, 2, jpkm1)
           zw3d(ji,jj,jk) =  denitr(ji,jj,jk) *  rdenit * xnegtr(ji,jj,jk) * cvol(ji,jj,jk)
         END_3D
         zdenittot = glob_sum ( 'p4zsms',  zw3d) 
         CALL iom_put( "tdenit", zdenittot * xfact3 )  
         !
         DEALLOCATE( zw3d )
      ENDIF
      !
   END SUBROUTINE p4z_budget

#if defined PISCES
      SUBROUTINE p4z_def_rst( ncid, total_rec, ierr)  ! restart netCDF

# include "netcdf.inc"

      logical :: create_new_file
      integer :: ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim    &
      &      , r2dgrd(3),  auxil(2),  checkdims                           &
#ifdef NC4PAR
      &      , csize,cmode           &
#endif
      &      , r3dgrd(4)
#ifdef USE_CALENDAR
      CHARACTER (len=19)    :: cdate,tool_sectodat
#endif
      CHARACTER(len=20) :: cltra, cltrs, cltru

!
! Put time record index into file name. In  the case when model
! output is to be arranged into sequence of named files, the naming
! convention is as follows: 'rst_root.INDEX.[MPI_node.]nc', where
! INDEX is an integer number such that (i) it is divisible by the
! specified number of records per file; and (ii)
!
!      INDEX + record_within_the_file = total_record
!
! where, 1 =< record_within_the_file =< records_per_file, so that
! total_record changes continuously throughout the sequence of files.
!
      ierr=0
      lstr=lenstr(cn_pisrst_out)
      if (nrpfrst.gt.0) then
        lvar=total_rec - (1+mod(total_rec-1, nrpfrst))
        call insert_time_index (cn_pisrst_out, lstr, lvar, ierr)
#ifdef USE_CALENDAR
        if (nrpfrst.eq.1) then
          cdate = tool_sectodat(time)
          cn_pisrst_out=TRIM(cn_pisrst_out(1:lstr-9))//'.'//cdate(7:10)//cdate(4:5)
          cn_pisrst_out=TRIM(cn_pisrst_out)//cdate(1:2)//cdate(12:13)//cdate(15:16)
          cn_pisrst_out=TRIM(cn_pisrst_out)//cdate(18:19)//'.nc'
          lstr=lenstr(cn_pisrst_out)
        end if
#endif
        if (ierr .ne. 0) goto 99
      endif

!
! Decide whether to create a new file, or open existing one.
! Overall the whole code below is organized into 3-way switch,
!
! 10  if (create_new_file) then
!        .... create new file, save netCDF ids for all variables;
!     elseif (ncid.eq.-1) then
!        .... try to open existing file and check its dimensions
!       if (cannot be opened or rejected) then
!         create_new_file=.true.
!         goto 10
!       endif   and prepare
!        .... prepare the file for adding new data,
!        .... find and save netCDF ids for all variables
!     else
!        .... just open, no checking, all ids are assumed to be
!        .... already known (MPI single file output only).
!     endif
!
! which is designed to implement flexible opening policy:
! if ldefhis=.true., it forces creation of a new file [if the
! file already exists, it will be overwritten]; on the other hand,
! ldefhis=.false., it is assumed that the file already exists and
! an attempt to open it is made; if the attempt is successful, the
! file is prepared for appending hew data; if it fails, a new file
! is created.
!
      create_new_file = ldefhis
      IF (ncid .NE. -1) create_new_file = .false.
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      IF (mynode > 0) create_new_file = .false.
#endif
!
! Create new restart file:    Put global attributes
!======= === ======= =====    and define all variables.
!
  10  if (create_new_file) then

#ifndef NC4PAR
        ierr  = nf_create(cn_pisrst_out(1:lstr),NF_CLOBBER, ncid)
#else
        cmode = ior(nf_netcdf4,nf_classic_model)
        cmode = ior(cmode, nf_mpiio)
        csize = xi_rho*eta_rho/NNODES
        WRITE(stdout,*)'CREATE RST NC4 PARALLEL FILE'
        ierr  = nf_create_par(cn_pisrst_out(1:lstr),cmode, &
        &       MPI_COMM_WORLD,MPI_INFO_NULL,ncid)
#endif

        IF (ierr .NE. nf_noerr) THEN
           WRITE(stdout,'(/3(1x,A)/)') 'ERROR in P4Z_DEF_RST: Cannot',    &
           &             'create restart NetCDF file:', TRIM(cn_pisrst_out)
           GOTO 99                                         !--> ERROR
        ENDIF
        IF (nrpfrst == 0) total_rec = 0
!
! Put global attributes.
! --- ------ -----------
!
        CALL put_global_atts (ncid, ierr)
!
! Define dimensions of staggered fields.
! ------ ---------- -- --------- -------
!
        ierr = nf_def_dim (ncid, 'xi_rho',   xi_rho,  r2dgrd(1))
        ierr = nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr = nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))        
        ierr = nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr = nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        auxil(2)  = timedim

        r2dgrd(3) = timedim           ! Free surface
        r3dgrd(1) = r2dgrd(1)         !
        r3dgrd(2) = r2dgrd(2)         ! 3D RHO-type
        r3dgrd(4) = timedim           !
!
! Define evolving model variables:
! ------ -------- ----- ----------
!
!
! Time step number and time record numbers:
!
        ierr = nf_def_var (ncid, 'time_step', nf_int, 2, auxil,     &
        &       rstpisstep)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstpisstep,nf_collective)
#endif
        ierr = nf_put_att_text (ncid, rstpisstep, 'long_name', 48,    &
        &       'time step and record numbers from initialization')
!
! Time.
!
        lvar = lenstr(vname(1,indxTime))
        ierr = nf_def_var (ncid, vname(1,indxTime)(1:lvar),           &
        &                              NF_DOUBLE, 1, timedim, rstpistime)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstpistime,nf_collective)
#endif
        lvar = lenstr(vname(2,indxTime))
        ierr = nf_put_att_text (ncid, rstpistime, 'long_name', lvar,     &
        &                                  vname(2,indxTime)(1:lvar))
        lvar = lenstr(vname(3,indxTime))
        ierr = nf_put_att_text (ncid, rstpistime, 'units',     lvar,     &
        &                                  vname(3,indxTime)(1:lvar))
        lvar = lenstr (vname(4,indxTime))
        ierr = nf_put_att_text(ncid, rstpistime, 'field',     lvar,      &
        &                                  vname(4,indxTime)(1:lvar))

!
! Time2.
!
        lvar = lenstr(vname(1,indxTime2))
        ierr = nf_def_var (ncid, vname(1,indxTime2)(1:lvar),            &
        &                              NF_DOUBLE, 1, timedim, rstpistime2)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstpistime2,nf_collective)
#endif
        lvar = lenstr(vname(2,indxTime2))
        ierr = nf_put_att_text (ncid, rstpistime2, 'long_name', lvar,     &
        &                                  vname(2,indxTime2)(1:lvar))
        lvar = lenstr(vname(3,indxTime2))
        ierr = nf_put_att_text (ncid, rstpistime2, 'units',     lvar,     &
        &                                  vname(3,indxTime2)(1:lvar))
        lvar = lenstr (vname(4,indxTime2))
        ierr = nf_put_att_text(ncid, rstpistime2, 'field',     lvar,      &
        &                                  vname(4,indxTime2)(1:lvar))

        cltra = "PH"   ;   cltrs = "PH"   ;    cltru = "-"
        ierr = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstph)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstph,nf_collective)
#endif
        lvar = lenstr(cltrs)
        ierr = nf_put_att_text (ncid, rstph, 'long_name',    &
        &                     lvar, cltrs )
        lvar = lenstr(cltru)
        ierr = nf_put_att_text (ncid, rstph, 'units', lvar, cltru )

        cltra = "Consfe3"   ;   cltrs = "Consfe3"   ;    cltru = "molFe/L/s"
        ierr = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstfe)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstfe,nf_collective)
#endif
        lvar = lenstr(cltrs)
        ierr = nf_put_att_text (ncid, rstfe, 'long_name',    &
        &                     lvar, cltrs )
        lvar = lenstr(cltru)
        ierr = nf_put_att_text (ncid, rstfe, 'units', lvar, cltru )

        cltra = "remintpoc"   ;   cltrs = "remintpoc"   ;    cltru = "d-1"
        ierr = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstpoc)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstpoc,nf_collective)
#endif
        lvar = lenstr(cltrs)
        ierr = nf_put_att_text (ncid, rstpoc, 'long_name',    &
        &                     lvar, cltrs )
        lvar = lenstr(cltru)
        ierr = nf_put_att_text (ncid, rstpoc, 'units', lvar, cltru )

        cltra = "sizen"   ;   cltrs = "sizen"   ;    cltru = "-"
        ierr = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstszn)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstszn,nf_collective)
#endif
        lvar = lenstr(cltrs)
        ierr = nf_put_att_text (ncid, rstszn, 'long_name',    &
        &                     lvar, cltrs )
        lvar = lenstr(cltru)
        ierr = nf_put_att_text (ncid, rstszn, 'units', lvar, cltru )

        IF( ln_p2z ) THEN
           cltra = "Thetanano"   ;   cltrs = "Chl/C"   ;    cltru = "gChl/gC"
           ierr = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstthet)
#ifdef NC4PAR
           ierr = nf_var_par_access(ncid,rstthet,nf_collective)
#endif
           lvar = lenstr(cltrs)
           ierr = nf_put_att_text (ncid, rstthet, 'long_name',    &
           &                     lvar, cltrs )
           lvar = lenstr(cltru)
           ierr = nf_put_att_text (ncid, rstthet, 'units', lvar, cltru )
        ENDIF        

        IF ( ln_p4z .OR. ln_p5z ) THEN
           cltra = "sized"   ;   cltrs = "sized"   ;    cltru = "-"
           ierr = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstszd)
#ifdef NC4PAR
           ierr = nf_var_par_access(ncid,rstszd,nf_collective)
#endif
           lvar = lenstr(cltrs)
           ierr = nf_put_att_text (ncid, rstszd, 'long_name',    &
           &                     lvar, cltrs )
           lvar = lenstr(cltru)
           ierr = nf_put_att_text (ncid, rstszd, 'units', lvar, cltru )

           cltra = "Silicalim"   ;   cltrs = "xksi"   ;    cltru = "molSi/l"
           ierr = nf_def_var (ncid, cltra, NF_DOUBLE, 3, r2dgrd, rstxksi)
#ifdef NC4PAR
           ierr = nf_var_par_access(ncid,rstxksi,nf_collective)
#endif
           lvar = lenstr(cltrs)
           ierr = nf_put_att_text (ncid, rstxksi, 'long_name',    &
           &                     lvar, cltrs )
           lvar = lenstr(cltru)
           ierr = nf_put_att_text (ncid, rstxksi, 'units', lvar, cltru )

           cltra = "Silicamax"   ;   cltrs = "xksimax"   ;    cltru = "molSi/l"
           ierr = nf_def_var (ncid, cltra, NF_DOUBLE, 3, r2dgrd, rstxksim)
#ifdef NC4PAR
           ierr = nf_var_par_access(ncid,rstxksim,nf_collective)
#endif
           lvar = lenstr(cltrs)
           ierr = nf_put_att_text (ncid, rstxksim, 'long_name',    &
           &                     lvar, cltrs )
           lvar = lenstr(cltru)
           ierr = nf_put_att_text (ncid, rstxksim, 'units', lvar, cltru )
           !
           cltra = "remintgoc"   ;   cltrs = "remintgoc"   ;    cltru = "d-1"
           ierr = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstgoc)
#ifdef NC4PAR
           ierr = nf_var_par_access(ncid,rstgoc,nf_collective)
#endif
           lvar = lenstr(cltrs)
           ierr = nf_put_att_text (ncid, rstgoc, 'long_name',    &
           &                     lvar, cltrs )
           lvar = lenstr(cltru)
           ierr = nf_put_att_text (ncid, rstgoc, 'units', lvar, cltru )
           !
        ENDIF        

        IF( ln_p5z ) THEN
           cltra = "sizep"   ;   cltrs = "sizep"   ;    cltru = "-"
           ierr = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstszp)
#ifdef NC4PAR
           ierr = nf_var_par_access(ncid,rstszp,nf_collective)
#endif
           lvar = lenstr(cltrs)
           ierr = nf_put_att_text (ncid, rstszp, 'long_name',    &
           &                     lvar, cltrs )
           lvar = lenstr(cltru)
           ierr = nf_put_att_text (ncid, rstszp, 'units', lvar, cltru )
        ENDIF        
!
!
!
! Leave definition mode.                  Also initialize record
! ----- ---------- -----                  dimension size to zero.
!
        ierr = nf_enddef(ncid)
        WRITE(*,'(6x,4A,1x,A,i4)') 'P4Z_DEF_RST - Created new ',        &
        &              'netCDF file ''', TRIM(cn_pisrst_out), '''.'
!
! Open an existing file and prepare for appending data.
! ==== == ======== ==== === ======= === ========= =====
! Check consistency of the dimensions of fields from the
! file with model dimensions. Determine the current size
! of unlimited dimension and set initial record [in the
! case of MPI serialized output, at this moment the last
! time record is assumed to be **partially** written by
! MPI processes with lower rank. Thus the next write is
! expected to be into the same record rather than next
! one (except MPI-master, who initializes the record).
!
! In the case when file is rejected (whether it cannot
! be opened, or something is wrong with its dimensions, 
! create new file. 
!
      ELSEIF (ncid == -1) THEN
#ifndef NC4PAR
        ierr = nf_open (cn_pisrst_out(1:lstr), nf_write, ncid)
#else
        ierr = nf_open_par (cn_pisrst_out(1:lstr), IOR(nf_write, nf_mpiio),   &
        &     MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#endif
        IF (ierr == nf_noerr) THEN
           ierr = checkdims (ncid, cn_pisrst_out, lstr, rec)
           IF (ierr == nf_noerr) THEN
              IF (nrpfrst == 0) THEN
                 ierr = rec+1 - nrecpisrst
              ELSE
                 ierr = rec+1 - (1+mod(nrecpisrst-1, abs(nrpfrst)))
              ENDIF
              IF (ierr > 0) THEN
                 MPI_master_only write( stdout,                              &
        &                 '(/1x,A,I5,1x,A/8x,3A,I5/8x,A,I5,1x,A/)'         &
        &           ) 'P4Z_DEF_RST WARNING: Actual number of records', rec,    &
        &             'in netCDF file',  '''',  cn_pisrst_out(1:lstr),       &
        &             ''' exceeds the record number from restart data',    &
        &             rec+1-ierr,'/', total_rec,', restart is assumed.'
                 rec = rec-ierr
              ELSEIF (nrpfrst == 0) THEN
                 total_rec = rec+1           ! <-- set to the next record
#if defined MPI & !defined PARALLEL_FILES
                 IF (mynode > 0) total_rec = total_rec-1
#endif
              ENDIF
              ierr = nf_noerr
           ENDIF
        ENDIF

        IF (ierr .NE. nf_noerr) THEN
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
           IF (mynode == 0) THEN
              create_new_file = .true.
              GOTO 10
           ELSE
              WRITE(stdout,'(/1x,4A, 1x,A,I4/)')     'P4Z_DEF_RST ERROR: ',    &
              &     'Cannot open restart netCDF file ''',cn_pisrst_out(1:lstr),'''.'
              GOTO 99                                     !--> ERROR 
           ENDIF
#else
           create_new_file=.true.
           GOTO 10
#endif
        ENDIF
!
! Find netCDF IDs of evolving model variables:
! ---- ------ --- -- -------- ----- ----------
!
! Time step indices:
!
        ierr = nf_inq_varid (ncid, 'time_step', rstpisstep)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) 'time_step', cn_pisrst_out(1:lstr)
          GOTO 99                                         !--> ERROR
        ENDIF
!
! Time.
!
        lvar = lenstr(vname(1,indxTime))
        ierr = nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), rstpistime)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) vname(1,indxTime)(1:lvar), cn_pisrst_out(1:lstr)
          GOTO 99                                         !--> ERROR
        ENDIF
!
! Time2.
!
        lvar = lenstr(vname(1,indxTime2))
        ierr = nf_inq_varid (ncid, vname(1,indxTime2)(1:lvar), rstpistime2)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) vname(1,indxTime2)(1:lvar), cn_pisrst_out(1:lstr)
          GOTO 99                                         !--> ERROR
        ENDIF

       cltra="PH"
       ierr = nf_inq_varid (ncid, TRIM(cltra), rstph)
       IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) TRIM(cltra), cn_pisrst_out(1:lstr)
          GOTO 99                                       !--> ERROR
       ENDIF

       cltra="Consfe3"
       ierr = nf_inq_varid (ncid, TRIM(cltra), rstfe)
       IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) TRIM(cltra), cn_pisrst_out(1:lstr)
          GOTO 99                                       !--> ERROR
       ENDIF

       cltra="sizen"
       ierr = nf_inq_varid (ncid, TRIM(cltra), rstszn)
       IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) TRIM(cltra), cn_pisrst_out(1:lstr)
          GOTO 99                                       !--> ERROR
       ENDIF

       cltra="remintpoc"
       ierr = nf_inq_varid (ncid, TRIM(cltra), rstpoc)
       IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) TRIM(cltra), cn_pisrst_out(1:lstr)
          GOTO 99                                       !--> ERROR
       ENDIF
      

        IF( ln_p2z ) THEN
           cltra = "Thetanano" 
           ierr = nf_inq_varid (ncid, TRIM(cltra), rstthet)
           IF (ierr .NE. nf_noerr) THEN
             WRITE(stdout,1) TRIM(cltra), cn_pisrst_out(1:lstr)
             GOTO 99                                       !--> ERROR
           ENDIF
        ENDIF        

        IF ( ln_p4z .OR. ln_p5z ) THEN
           cltra = "sized"
           ierr = nf_inq_varid (ncid, TRIM(cltra), rstszd)
           IF (ierr .NE. nf_noerr) THEN
             WRITE(stdout,1) TRIM(cltra), cn_pisrst_out(1:lstr)
             GOTO 99                                       !--> ERROR
           ENDIF

           cltra = "Silicalim" 
           ierr = nf_inq_varid (ncid, TRIM(cltra), rstxksi)
           IF (ierr .NE. nf_noerr) THEN
             WRITE(stdout,1) TRIM(cltra), cn_pisrst_out(1:lstr)
             GOTO 99                                       !--> ERROR
           ENDIF

           cltra = "Silicamax" 
           ierr = nf_inq_varid (ncid, TRIM(cltra), rstxksim)
           IF (ierr .NE. nf_noerr) THEN
             WRITE(stdout,1) TRIM(cltra), cn_pisrst_out(1:lstr)
             GOTO 99                                       !--> ERROR
           ENDIF

           cltra="remintgoc"
           ierr = nf_inq_varid (ncid, TRIM(cltra), rstgoc)
           IF (ierr .NE. nf_noerr) THEN
              WRITE(stdout,1) TRIM(cltra), cn_pisrst_out(1:lstr)
              GOTO 99                                       !--> ERROR
           ENDIF
      
        ENDIF        

        IF( ln_p5z ) THEN
           cltra = "sizep"
           ierr = nf_inq_varid (ncid, TRIM(cltra), rstszp)
           IF (ierr .NE. nf_noerr) THEN
             WRITE(stdout,1) TRIM(cltra), cn_pisrst_out(1:lstr)
             GOTO 99                                       !--> ERROR
           ENDIF
        ENDIF        
!
        MPI_master_only WRITE(*,'(6x,2A,i4,1x,A,i4)')              &
        &             'P4Z_DEF_RST -- Opened ',                        &
        &             'existing restart file,  record =', rec

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      ELSE
         ierr = nf_open (cn_pisrst_out(1:lstr), nf_write, ncid)
         IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,'(/1x,4A, 1x,A,I4/)')        &
            &          'P4Z_DEF_RST ERROR: Cannot',                       &
            &          'open restart netCDF file ''', cn_pisrst_out(1:lstr), '''.'
            GOTO 99                                         !--> ERROR
         ENDIF
#endif
      ENDIF              !<-- create_new_file
   1  FORMAT(/1x,'P4Z_DEF_RST ERROR: Cannot find variable ''',        &
      &               A, ''' in netCDF file ''', A, '''.'/)
  99  RETURN                                              !--> ERROR
      END SUBROUTINE p4z_def_rst

      SUBROUTINE p4z_rst_wri      ! variables into restart
                                  ! netCDF file.
# include "netcdf.inc"

      INTEGER :: ierr, record, lstr, lvar, lenstr   &
      &  , start(2), count(2), ibuff(2), nf_fwrite
      INTEGER :: ji, jj, jk, jn
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw2d
      CHARACTER(len=20) :: cltra


#ifdef MPI
#define LOCALLM Lmmpi
#define LOCALMM Mmmpi
#else
#define LOCALLM Lm
#define LOCALMM Mm
#endif


#if defined MPI & !defined PARALLEL_FILES
      INCLUDE 'mpif.h'
      INTEGER status(MPI_STATUS_SIZE), blank
#endif

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      IF (mynode > 0) THEN
         call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,       &
         &                 1, MPI_COMM_WORLD, status, ierr)
      ENDIF
#endif
!
! Create/open restart file; write grid arrays, if so needed.
!
      CALL p4z_def_rst(ncidpisrst, nrecpisrst, ierr)
      IF (ierr .NE. nf_noerr) GOTO 99
      lstr = lenstr(cn_pisrst_out)
!                                            !!! WARNING: Here it is
! Set record within the file.                !!! assumed that global
!                                            !!! restart record index 
      nrecpisrst = max(nrecpisrst,1)                 !!! nrecrst is already
      IF (nrpfrst == 0) THEN                 !!! advanced by main.
         record = nrecpisrst
      ELSE
         record = 1+mod(nrecpisrst-1, abs(nrpfrst))
      ENDIF

!
! Write out evolving model variables:
! ----- --- -------- ----- ----------
!
! Time step number and record indices. 
!
      ibuff(1) = iic
      ibuff(2) = nrecpisrst
      start(1) = 1
      start(2) = record
      count(1) = 2 
      count(2) = 1
      ierr = nf_put_vara_int (ncidpisrst, rstpisstep, start, count, ibuff)
      IF (ierr .NE. nf_noerr) THEN
         WRITE(stdout,1) 'time_step', record, ierr      
         GOTO 99                                           !--> ERROR
      ENDIF
!
! Time.
!
      ierr = nf_put_var1_FTYPE (ncidpisrst, rstpistime, record, time)
      IF (ierr .NE. nf_noerr) THEN
         lvar = lenstr(vname(1,indxTime))
         WRITE(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
         GOTO 99                                           !--> ERROR
      ENDIF
!
      ALLOCATE( zw3d(GLOBAL_2D_ARRAY,jpk) )  ;  zw3d(:,:,:) = 0._wp
      ALLOCATE( zw2d(GLOBAL_2D_ARRAY) )      ;  zw2d(:,:) = 0._wp

      cltra="PH"   ;    zw3d(A2D(0),:) = hi(A2D(0),:)
      ierr = nf_fwrite(zw3d(START_2D_ARRAY,1), ncidpisrst,   &
      &                             rstph, record, r3dvar)
      IF (ierr .NE. nf_noerr) THEN
         WRITE(stdout,1) cltra, record, ierr
         GOTO 99                                         !--> ERROR
      ENDIF

      cltra="Consfe3"   ;    zw3d(A2D(0),:) = consfe3(A2D(0),:)
      ierr = nf_fwrite(zw3d(START_2D_ARRAY,1), ncidpisrst,   &
      &                             rstfe, record, r3dvar)
      IF (ierr .NE. nf_noerr) THEN
         WRITE(stdout,1) cltra, record, ierr
         GOTO 99                                         !--> ERROR
      ENDIF

      cltra="sizen"   ;    zw3d(A2D(0),:) = sizen(A2D(0),:)
      ierr = nf_fwrite(zw3d(START_2D_ARRAY,1), ncidpisrst,   &
      &                             rstszn, record, r3dvar)
      IF (ierr .NE. nf_noerr) THEN
         WRITE(stdout,1) cltra, record, ierr
         GOTO 99                                         !--> ERROR
      ENDIF

      cltra="remintpoc"   ;    zw3d(A2D(0),:) = remintpoc(A2D(0),:)
      ierr = nf_fwrite(zw3d(START_2D_ARRAY,1), ncidpisrst,   &
      &                             rstpoc, record, r3dvar)
      IF (ierr .NE. nf_noerr) THEN
         WRITE(stdout,1) cltra, record, ierr
         GOTO 99                                         !--> ERROR
      ENDIF
!
      IF( ln_p2z ) THEN
         cltra="Thetanano"   ;    zw3d(A2D(0),:) = thetanano(A2D(0),:)
         ierr = nf_fwrite(zw3d(START_2D_ARRAY,1), ncidpisrst,   &
         &                             rstthet, record, r3dvar)
         IF (ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) cltra, record, ierr
            GOTO 99                                         !--> ERROR
         ENDIF
      ENDIF

      IF( ln_p4z .OR. ln_p5z ) THEN
         cltra="sized"   ;    zw3d(A2D(0),:) = sized(A2D(0),:)
         ierr = nf_fwrite(zw3d(START_2D_ARRAY,1), ncidpisrst,   &
         &                             rstszd, record, r3dvar)
         IF (ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) cltra, record, ierr
            GOTO 99                                         !--> ERROR
         ENDIF

         cltra="Silicalim"   ;    zw2d(A2D(0)) = xksi(A2D(0))
         ierr = nf_fwrite(zw2d(START_2D_ARRAY), ncidpisrst,   &
         &                             rstxksi, record, r2dvar)
         IF (ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) cltra, record, ierr
            GOTO 99                                         !--> ERROR
         ENDIF

         cltra="Silicamax"   ;    zw2d(A2D(0)) = xksimax(A2D(0))
         ierr = nf_fwrite(zw2d(START_2D_ARRAY), ncidpisrst,   &
         &                             rstxksim, record, r2dvar)
         IF (ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) cltra, record, ierr
            GOTO 99                                         !--> ERROR
         ENDIF

         cltra="remintgoc"   ;    zw3d(A2D(0),:) = remintgoc(A2D(0),:)
         ierr = nf_fwrite(zw3d(START_2D_ARRAY,1), ncidpisrst,   &
         &                             rstgoc, record, r3dvar)
         IF (ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) cltra, record, ierr
            GOTO 99                                         !--> ERROR
         ENDIF
!
      ENDIF

      IF( ln_p5z ) THEN
         cltra="sizep"   ;    zw3d(A2D(0),:) = sizep(A2D(0),:)
         ierr = nf_fwrite(zw3d(START_2D_ARRAY,1), ncidpisrst,   &
         &                             rstszp, record, r3dvar)
         IF (ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) cltra, record, ierr
            GOTO 99                                         !--> ERROR
         ENDIF
      ENDIF

      DEALLOCATE(zw3d, zw2d)
!
  1   FORMAT(/1x, 'P4Z_RST_WRI ERROR while writing variable ''', A,   &
       &           ''' into restart file.', /11x, 'Time record:', &
       &               i6, 3x, 'netCDF error code', i4, 3x, A,i4) 
      GOTO 100 
  99  may_day_flag=3
 100  CONTINUE

!
! Synchronize restart netCDF file to disk to allow other
! processes to access data immediately after it is written.
!
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      ierr = nf_close (ncidpisrst)
      IF (nrpfrst > 0 .AND. record >= nrpfrst) ncidpisrst = -1
#else
      IF (nrpfrst > 0 .AND. record >= nrpfrst) THEN
        ierr = nf_close (ncidpisrst)
        ncidpisrst = -1
      ELSE
        ierr = nf_sync(ncidpisrst)
      ENDIF
#endif
      IF (ierr == nf_noerr) THEN
         MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')    & 
         &            'P4Z_RST_WRI -- wrote ',                          &
         &            'restart fields into time record =', record, '/',  &
         &             nrecpisrst  
      ELSE
         MPI_master_only  write(stdout,'(/1x,2A/)')     & 
         &             'P4Z_RST_WRI ERROR: Cannot ',        &
         &             'synchronize/close restart netCDF file.'
         may_day_flag = 3
      ENDIF

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      IF (mynode < NNODES-1) THEN
         CALL MPI_Send (blank, 1, MPI_INTEGER, mynode+1, 1, MPI_COMM_WORLD,  ierr)
      ENDIF
#endif
      RETURN
      END SUBROUTINE p4z_rst_wri 

                              ! Read initial conditions for the
      SUBROUTINE p4z_rst_read( Kbb, Kmm ) ! primitive variables from NetCDF
                              ! initialization file.

        INTEGER, INTENT( in ) ::   Kbb, Kmm  ! time level indices
!======================================================
!
!======================================================

# include "netcdf.inc"

      real(wp) :: time_scale
      integer  :: itrc
      integer  :: ji, jj, jk, jn
      integer  :: ncid, indx, varid,  ierr, lstr, lvar, latt, lenstr,    &
      &        start(2), count(2), ibuff(2), nf_fread, checkdims, max_rec
      character :: units*180
      logical   :: llexist
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw2d
      CHARACTER(len=20) :: cltra


#ifdef MPI
#define LOCALLM Lmmpi
#define LOCALMM Mmmpi
#else
#define LOCALLM Lm
#define LOCALMM Mm
#endif

!
! Open initial conditions netCDF file for reading. Check that all
! spatial dimensions in that file are consistent with the model
! arrays, determine how many time records are available in the file
! and set record from which the dada will be read.
!
! The record is set as follows: (1) if only one time record is
! available in the file, then that record is used REGARDLESS of
! value of nrrec supplied from the parameter file; (2) if the
! file has multiple records and nrrec is positive, then nrrec is
! available record is used.
!
      if (may_day_flag.ne.0) return      !-->  EXIT
      lstr = lenstr(cn_pisrst_in)
      ierr = nf_open(TRIM(cn_pisrst_in), nf_nowrite, ncid)
      if (ierr .eq. nf_noerr) then
        llexist = .true.      
        ierr=checkdims (ncid, cn_pisrst_in, lstr, max_rec)
        if (ierr .eq. nf_noerr) then
          if (max_rec.gt.0) then
            if (nrrec.gt.0) then
              if (nrrec.le.max_rec) then
                indx=nrrec
              else
       MPI_master_only write(stdout,'(/1x,2A,I4,1x,A/12x,A,I4,1x,3A/)')  &
     &           'P4Z_RST_READ ERROR :: requested restart time ', &
     &           'record',  nrrec, 'exceeds number',  'of records', &
     &                       max_rec,  'available in netCDF file ''',  &
     &                                        TRIM(cn_pisrst_in), '''.'

                goto 99         !--> ERROR
              endif
            else
              indx=max_rec
            endif
          else
            indx=1
          endif
        else
          goto 99
        endif
      else
        llexist = .false.      
        MPI_master_only write(stdout,'(/1x,4A/12x,A/)')  &
     &               'P4Z_RST_READ WARNING :: ',    &
     &               'Cannot open netCDF file ''',   TRIM(cn_pisrst_in),  &
     &                                  '''.',     nf_strerror(ierr)
!        goto 99                                           !--> ERROR
      endif
!
! Read in evolving model variables:
! ---- -- -------- ----- ----------
!
! Time: find netCDF id, read value, read attribute 'units'
! and set starting time index and time clock in days.
!
      IF( llexist ) THEN
      lvar = lenstr(vname(1,indxTime))
      ierr = nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), varid)
      IF (ierr == nf_noerr) THEN
        ierr = nf_get_var1_FTYPE (ncid, varid, indx, time)
        IF (ierr == nf_noerr) THEN
          ierr = nf_get_att_text(ncid, varid, 'units', units)
          IF (ierr == nf_noerr) THEN
            latt = lenstr(units)
            IF (units(1:6) == 'second') THEN
               time_scale = 1.
            ELSEIF (units(1:3) == 'day') THEN
              time_scale = day2sec
            ELSE
              WRITE (stdout,'(/1x,4A/8x,3A/)') 'P4Z_RST_READ ',      &
       &              'ERROR: unknown units of for variable ''',      &
       &               vname(1,indxTime)(1:lvar), '''',               &
       &              'in netCDF file ''', TRIM(cn_pisrst_in),'''.'
              GOTO 99                                    !--> ERROR
            ENDIF
          ELSE
            WRITE (stdout,'(/1x,2A/8x,5A/)') 'P4Z_RST_READ ERROR: ',   &
       &             'cannot read attribute ''units'' for variable',  &
       &             '''', vname(1,indxTime)(1:lvar),                 &
       &             ''' in netCDF file ''',  TRIM(cn_pisrst_in), '''.'
            GOTO 99                                       !--> ERROR
          ENDIF
        ELSE
          MPI_master_only write(stdout,2) vname(1,indxTime)(1:lvar)  &
          &                                , indx, TRIM(cn_pisrst_in)
          GOTO 99                                         !--> ERROR
        ENDIF
      ELSE
        MPI_master_only write(stdout,1) vname(1,indxTime)(1:lvar), TRIM(cn_pisrst_in)
        GOTO 99                                           !--> ERROR
      ENDIF

!      time = time*time_scale
!      tdays = time*sec2day

      ierr = nf_inq_varid (ncid, 'time_step', varid)
      IF (ierr == nf_noerr) THEN
         start(1) = 1
         start(2) = indx
         count(1) = 2
         count(2) = 1
         ierr = nf_get_vara_int (ncid, varid, start, count, ibuff)
         IF (ierr == nf_noerr) THEN
!            ntstart = ibuff(1)
            nrecpisrst = ibuff(2)

            MPI_master_only WRITE(stdout,                            &
            &     '(6x,A,G12.4,A,I2,A,I6,A,I3,A)')              &
            &     'P4Z_RST_READ: Restarted from day =', tdays, ' rec =',   &
            &      indx, '(', ntstart, ',', nrecpisrst, ').'

         ELSE
            MPI_master_only write(stdout,'(/1x,2A/)')                     &
            &                            'P4Z_RST_READ ERROR: Cannot ',    &
            &                            'read time and record indices.'
            GOTO 99                                         !--> ERROR
         ENDIF
      ELSE
!         ntstart = 1
         nrecpisrst = 0
         MPI_master_only WRITE(stdout,'(6x,2A,G12.4,1x,A,I4)')      &
         &          'P4Z_RST_READ -- ',                              &
         &          'Processing data for time =', tdays, 'record =', indx
      ENDIF
!      IF (ntstart < 1) ntstart = 1
!      ntimes = ntstart+ntimes-1
!
! Tracer variables.
!
      ALLOCATE( zw3d(GLOBAL_2D_ARRAY,jpk), zw2d(GLOBAL_2D_ARRAY) )

      cltra ="PH"   
      ierr = nf_inq_varid (ncid, cltra, varid)
      IF(ierr == nf_noerr) THEN
         ierr = nf_fread (zw3d(START_2D_ARRAY,1), ncid, varid, indx, r3dvar)
         IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_pisrst_in)
            GOTO 99                                       !--> ERROR
         ENDIF
         MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_pisrst_in)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            hi(ji,jj,jk) = zw3d(ji,jj,jk)
         END_3D
      ELSE
         MPI_master_only WRITE(stdout,4) cltra, TRIM(cn_pisrst_in)
         CALL p4z_che( Kbb, Kmm )                  ! initialize the chemical constants
         CALL ahini_for_at( hi, Kbb )
      ENDIF

      cltra ="Consfe3"   
      ierr = nf_inq_varid (ncid, cltra, varid)
      IF(ierr == nf_noerr) THEN
         ierr = nf_fread (zw3d(START_2D_ARRAY,1), ncid, varid, indx, r3dvar)
         IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_pisrst_in)
            GOTO 99                                       !--> ERROR
         ENDIF
         MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_pisrst_in)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            consfe3(ji,jj,jk) = zw3d(ji,jj,jk)
         END_3D
      ELSE
         MPI_master_only WRITE(stdout,4) cltra, TRIM(cn_pisrst_in)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            consfe3(ji,jj,jk) = 0.
         END_3D
      ENDIF

      cltra ="sizen"   
      ierr = nf_inq_varid (ncid, cltra, varid)
      IF(ierr == nf_noerr) THEN
         ierr = nf_fread (zw3d(START_2D_ARRAY,1), ncid, varid, indx, r3dvar)
         IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_pisrst_in)
            GOTO 99                                       !--> ERROR
         ENDIF
         MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_pisrst_in)
         DO_3D( 0, 0, 0, 0, 1, jpk)
             sizen(ji,jj,jk) = zw3d(ji,jj,jk)
         END_3D
      ELSE
         MPI_master_only WRITE(stdout,4) cltra, TRIM(cn_pisrst_in)
         DO_3D( 0, 0, 0, 0, 1, jpk)
             sizen(ji,jj,jk) = 1.
         END_3D
      ENDIF

      cltra ="remintpoc"   
      ierr = nf_inq_varid (ncid, cltra, varid)
      IF(ierr == nf_noerr) THEN
         ierr = nf_fread (zw3d(START_2D_ARRAY,1), ncid, varid, indx, r3dvar)
         IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_pisrst_in)
            GOTO 99                                       !--> ERROR
         ENDIF
         MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_pisrst_in)
         DO_3D( 0, 0, 0, 0, 1, jpk)
             remintpoc(ji,jj,jk) = zw3d(ji,jj,jk)
         END_3D
      ELSE
         MPI_master_only WRITE(stdout,4) cltra, TRIM(cn_pisrst_in)
         DO_3D( 0, 0, 0, 0, 1, jpk)
             remintpoc(ji,jj,jk) = 0.03
         END_3D
      ENDIF

      IF( ln_p2z ) THEN
        cltra ="Thetanano"   
        ierr = nf_inq_varid (ncid, cltra, varid)
        IF(ierr == nf_noerr) THEN
           ierr = nf_fread (zw3d(START_2D_ARRAY,1), ncid, varid, indx, r3dvar)
           IF (ierr .NE. nf_noerr) THEN
              MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_pisrst_in)
              GOTO 99                                       !--> ERROR
           ENDIF
           MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_pisrst_in)
           DO_3D( 0, 0, 0, 0, 1, jpk)
              thetanano(ji,jj,jk) = zw3d(ji,jj,jk)
           END_3D
        ELSE
           MPI_master_only WRITE(stdout,4) cltra, TRIM(cn_pisrst_in)
           DO_3D( 0, 0, 0, 0, 1, jpk)
              thetanano(ji,jj,jk) = 1.0 / 55.0
           END_3D
        ENDIF
      ENDIF  

     IF ( ln_p4z .OR. ln_p5z ) THEN
        cltra ="Silicalim"   
        ierr = nf_inq_varid (ncid, cltra, varid)
        IF(ierr == nf_noerr) THEN
           ierr = nf_fread (zw2d(START_2D_ARRAY), ncid, varid, indx, r2dvar)
           IF (ierr .NE. nf_noerr) THEN
              MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_pisrst_in)
              GOTO 99                                       !--> ERROR
           ENDIF
           MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_pisrst_in)
           DO_2D( 0, 0, 0, 0 )
              xksi(ji,jj) = zw2d(ji,jj)
           END_2D
        ELSE
           MPI_master_only WRITE(stdout,4) cltra, TRIM(cn_pisrst_in)
           DO_2D( 0, 0, 0, 0 )
              xksi(ji,jj) = 2.e-6
           END_2D
        ENDIF

        cltra ="Silicamax"   
        ierr = nf_inq_varid (ncid, cltra, varid)
        IF(ierr == nf_noerr) THEN
           ierr = nf_fread (zw2d(START_2D_ARRAY), ncid, varid, indx, r2dvar)
           IF (ierr .NE. nf_noerr) THEN
              MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_pisrst_in)
              GOTO 99                                       !--> ERROR
           ENDIF
           MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_pisrst_in)
           DO_2D( 0, 0, 0, 0 )
              xksimax(ji,jj) = zw2d(ji,jj)
           END_2D
        ELSE
           MPI_master_only WRITE(stdout,4) cltra, TRIM(cn_pisrst_in)
           DO_2D( 0, 0, 0, 0 )
              xksimax(ji,jj) = xksi(ji,jj)
           END_2D
        ENDIF

        cltra ="sized"   
        ierr = nf_inq_varid (ncid, cltra, varid)
        IF(ierr == nf_noerr) THEN
           ierr = nf_fread (zw3d(START_2D_ARRAY,1), ncid, varid, indx, r3dvar)
           IF (ierr .NE. nf_noerr) THEN
              MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_pisrst_in)
              GOTO 99                                       !--> ERROR
           ENDIF
           MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_pisrst_in)
           DO_3D( 0, 0, 0, 0, 1, jpk)
              sized(ji,jj,jk) = zw3d(ji,jj,jk)
           END_3D
        ELSE
           MPI_master_only WRITE(stdout,4) cltra, TRIM(cn_pisrst_in)
           DO_3D( 0, 0, 0, 0, 1, jpk)
              sized(ji,jj,jk) = 1.
           END_3D
        ENDIF

        cltra ="remintgoc"   
        ierr = nf_inq_varid (ncid, cltra, varid)
        IF(ierr == nf_noerr) THEN
           ierr = nf_fread (zw3d(START_2D_ARRAY,1), ncid, varid, indx, r3dvar)
           IF (ierr .NE. nf_noerr) THEN
              MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_pisrst_in)
              GOTO 99                                       !--> ERROR
           ENDIF
           MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_pisrst_in)
           DO_3D( 0, 0, 0, 0, 1, jpk)
              remintgoc(ji,jj,jk) = zw3d(ji,jj,jk)
           END_3D
        ELSE
           MPI_master_only WRITE(stdout,4) cltra, TRIM(cn_pisrst_in)
           DO_3D( 0, 0, 0, 0, 1, jpk)
             remintgoc(ji,jj,jk) = 0.03
           END_3D
        ENDIF

     ENDIF

     IF( ln_p5z ) THEN
        cltra ="sizep"   
        ierr = nf_inq_varid (ncid, cltra, varid)
        IF(ierr == nf_noerr) THEN
           ierr = nf_fread (zw3d(START_2D_ARRAY,1), ncid, varid, indx, r3dvar)
           IF (ierr .NE. nf_noerr) THEN
              MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_pisrst_in)
              GOTO 99                                       !--> ERROR
           ENDIF
           MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_pisrst_in)
           DO_3D( 0, 0, 0, 0, 1, jpk)
              sizep(ji,jj,jk) = zw3d(ji,jj,jk)
           END_3D
        ELSE
           MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_pisrst_in)
           DO_3D( 0, 0, 0, 0, 1, jpk)
              sizep(ji,jj,jk) = 1.
           END_3D
        ENDIF
     ENDIF

     DEALLOCATE( zw3d, zw2d )
!
!  Close input NetCDF file.
      ierr = nf_close(ncid)

     ELSE  ! llexits = .false. restart file doesn't exist. Initalisation of PISCES variable with analytical value
         ! PH
         cltra ="PH"   
         MPI_master_only WRITE(stdout,5) cltra
         CALL p4z_che( Kbb, Kmm )                  ! initialize the chemical constants
         CALL ahini_for_at( hi, Kbb )
         ! Consfe3
         cltra ="Consfe3"   
         MPI_master_only WRITE(stdout,5) cltra
         consfe3(:,:,:) = 0.
         ! sizen
         cltra ="sizen"   
         MPI_master_only WRITE(stdout,5) cltra
         sizen(:,:,:) = 1.
         ! remintpoc
         cltra ="remintpoc"   
         MPI_master_only WRITE(stdout,5) cltra
         remintpoc(:,:,:) = 0.03
         ! Thetanano
         IF( ln_p2z ) THEN
            cltra ="Thetanano"   
            MPI_master_only WRITE(stdout,5) cltra
            thetanano(:,:,:) = 1.0 / 55.0
         ENDIF
         ! Xksi/Xksimax/sized
         IF( ln_p4z .OR. ln_p5z ) THEN
            cltra ="xksi"   
            MPI_master_only WRITE(stdout,5) cltra
            xksi(:,:) = 2.e-6
            cltra ="xksimax"   
            MPI_master_only WRITE(stdout,5) cltra
            xksimax(:,:) = xksi(:,:)
            cltra ="sized"   
            MPI_master_only WRITE(stdout,5) cltra
            sized(:,:,:) = 1.
            ! remintgoc
            cltra ="remintgoc"   
            MPI_master_only WRITE(stdout,5) cltra
            remintgoc(:,:,:) = 0.03
         ENDIF
         ! sizep
         IF( ln_p5z ) THEN
            cltra ="sized"   
            MPI_master_only WRITE(stdout,5) cltra
            sizep(:,:,:) = 1.
         ENDIF
      ENDIF



  1   FORMAT(/1x,'P4Z_RST_READ - unable to find variable:',    1x,A,    &
      &                            /15x,'in input NetCDF file:',1x,A/)
  2   FORMAT(/1x,'P4Z_RST_READ - error while reading variable:',1x, A,  &
      &    2x,'at time record =',i4/15x,'in input NetCDF file:',1x,A/)
  3   format(/1x,'P4Z_RST_READ : Read variable ',    1x,A,    &
     &                            /15x,'in input NetCDF file:',1x,A/)
  4   FORMAT(/1x,'P4Z_RST_READ - Unable to find variable:',    1x,A,    &
      &                            /15x,'in input NetCDF file:',1x,A,  &
      &    1x,'-> set constant value'/)
  5   format(/1x,'P4Z_RST_READ : Initialize variable ',    1x,A/)
      RETURN
  99  may_day_flag = 2
      RETURN
      END SUBROUTINE p4z_rst_read
#endif

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sms                   ! Empty routine
   END SUBROUTINE p4z_sms
#endif

   !!======================================================================
END MODULE p4zsms 

