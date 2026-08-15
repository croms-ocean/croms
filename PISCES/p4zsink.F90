#include "cppdefs.h"

MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  PISCES  vertical flux of particulate matter due to 
   !!        gravitational sinking
   !!        This module is the same for both PISCES and PISCES-QUOTA
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!             4.0  !  2019     (O. Aumont) an external subroutine is called
   !!                                          to compute the impact of sinking
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE trcsink         !  General routine to compute sedimentation
   USE prtctl          !  print control for debugging
   USE iom             !  I/O manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcini_pisces.F90
   PUBLIC   p4z_sink_alloc   ! called in trcini_pisces.F90

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sinkpocb  !: POC sinking fluxes at bottom 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sinkcalb  !: CaCO3 sinking fluxes at bottom
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sinksilb  !: BSi sinking fluxes at bottom
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sinkponb  !: POC sinking fluxes at bottom 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sinkpopb  !: POC sinking fluxes at bottom 
   INTEGER,  PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: nlev100

   REAL(wp) :: xfact
   LOGICAL  :: l_dia_sink, l_diag

   !! * Substitutions
#  include "ocean2pisces.h90"   
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsink.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sink ( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking. 
      !!
      !! ** Method  : - An external advection subroutine is called to compute
      !!                the impact of sinking on the particles. The tracers
      !!                concentrations are updated in this subroutine which
      !!                is mandatory to deal with negative concentrations
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER, INTENT(in) :: Kbb, Kmm, Krhs  ! time level indices
      INTEGER  ::   ji, jj, jk, ikb, ik100
      CHARACTER (len=25) :: charout
      REAL(wp) :: zmax, zfact
      REAL(wp), DIMENSION(A2D(0),jpk+1) :: zsinking, zsinking2
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_sink')

      IF( kt == nittrc000 )  THEN
         l_dia_sink = iom_use( "EPC100" ) .OR. iom_use( "EPFE100" ) .OR. iom_use( "EPCAL100" ) .OR. iom_use( "EPSI100" )  &
           &    .OR.  iom_use( "EXPC"   ) .OR. iom_use( "EXPFE"   ) .OR. iom_use( "EXPCAL"   ) .OR. iom_use( "EXPSI"   )  &
           &    .OR.  iom_use( "tcexp"  ) 
         !
         xfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
      ENDIF
      l_diag = l_dia_sink .AND. knt == nrdttrc

      !
      nlev100(:,:) = jpk        !  last level where depth less than 100 m
      DO_2D( 0, 0, 0, 0 )
         DO jk = jpkm1, 1, -1
            IF( gdept(ji,jj,jk,Kmm) > 100. )  nlev100(ji,jj) = jk - 1
         END DO
      END_2D

      ! Initialize the local arrays
      zsinking(:,:,:) = 0._wp ; zsinking2(:,:,:) = 0._wp

      ! Sinking speeds of big detritus is increased with depth as shown
      ! by data and from the coagulation theory. This is controled by
      ! wsbio2max and wsbio2scale. If wsbio2max is set to wsbio2, then
      ! sinking speed is constant with depth.
      ! CaCO3 and bSi are supposed to sink at the big particles speed 
      ! due to their high density
      ! ---------------------------------------------------------------
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         zmax  = MAX( heup_01(ji,jj), hmld(ji,jj) )
         zfact = MAX( 0., gdepw(ji,jj,jk+1,Kmm) - zmax ) / wsbio2scale
         wsbio4(ji,jj,jk) = wsbio2 + MAX(0., ( wsbio2max - wsbio2 )) * zfact
      END_3D
      wsbio4(:,:,jpk) = wsbio2

      ! Sinking speed of the small particles is always constant
      ! except in PISCES reduced
      IF( ln_p2z ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpkm1)
            zmax  = MAX( heup_01(ji,jj), hmld(ji,jj) )
            zfact = MAX( 0., gdepw(ji,jj,jk+1,Kmm) - zmax ) &
              &      * tgfunc(ji,jj,jk) * ( 0.03 / wsbio2 - 0.015/ wsbio)
            zfact = MIN(3.0, EXP(-zfact) )
            wsbio3(ji,jj,jk) = ( wsbio + wsbio2 * ( sizen(ji,jj,1) - 1.0 ) * 0.05 * zfact )    &
               &               / ( 1.0 + ( sizen(ji,jj,1) - 1.0 ) * 0.05 * zfact )
         END_3D
         wsbio3(:,:,jpk) = wsbio
      ELSE
         wsbio3(:,:,:) = wsbio
      ENDIF

      ! Compute the sedimentation term using trc_sink for all the sinking particles
      ! ---------------------------------------------------------------------------
      CALL trc_sink( kt, Kbb, Kmm, wsbio3, zsinking , jppoc, rfact2 )
      IF(  .NOT. ln_p2z ) THEN
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking2, jpgoc, rfact2 )
      ENDIF
      DO_2D( 0, 0, 0, 0 )
         ikb = mbkt(ji,jj)+1
         sinkpocb(ji,jj) = zsinking(ji,jj,ikb) + zsinking2(ji,jj,ikb)
      END_2D
      !
      IF( l_diag ) THEN
         ALLOCATE( zw3d(GLOBAL_2D_ARRAY,jpk) )  ;  zw3d(:,:,:) = 0._wp     
         ALLOCATE( zw2d(GLOBAL_2D_ARRAY)     )  ;  zw2d(:,:) = 0._wp
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jkR) = ( zsinking(ji,jj,jk) + zsinking2(ji,jj,jk) ) &
                   &        * xfact * tmask(ji,jj,jk)
         END_3D
         DO_2D( 0, 0, 0, 0 )
            ik100 = nlev100(ji,jj)
            zw2d(ji,jj) = ( zsinking(ji,jj,ik100) + zsinking2(ji,jj,ik100) ) &
                   &        * xfact * tmask(ji,jj,ik100)
         END_2D
         CALL iom_put( "EPC100",  zw2d )  ! Export of carbon at 100m
         CALL iom_put( "EXPC"  ,  zw3d )             ! Export of carbon in the water column
         DO_2D( 0, 0, 0, 0 )
            zw2d(ji,jj) = zw2d(ji,jj) * e1e2t(ji,jj)
         END_2D 
         t_oce_co2_exp  = glob_sum( 'p4zsink',  zw2d )
         CALL iom_put( "tcexp", t_oce_co2_exp )      ! Total cabon exort
      ENDIF
#if defined key_trc_diaadd
      DO_2D( 0, 0, 0, 0 )
         ik100 = nlev100(ji,jj)
         trc2d(ji,jj,jp_sinkco2) = ( zsinking(ji,jj,ik100) + zsinking2(ji,jj,ik100) ) &
                            &        * xfact * tmask(ji,jj,ik100)  
      END_2D
#endif      
      IF(  .NOT. ln_p2z ) THEN
         !
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking2, jpcal, rfact2 )
         DO_2D( 0, 0, 0, 0 )
            ikb = mbkt(ji,jj)+1
            sinkcalb(ji,jj) = zsinking2(ji,jj,ikb)
         END_2D
         IF( l_diag ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk)
               zw3d(ji,jj,jkR) =  zsinking2(ji,jj,jk) * xfact * tmask(ji,jj,jk)
            END_3D
            DO_2D( 0, 0, 0, 0 )
               ik100 = nlev100(ji,jj)
               zw2d(ji,jj) = zsinking2(ji,jj,ik100) * xfact * tmask(ji,jj,ik100)
            END_2D
            CALL iom_put( "EPCAL100",  zw2d )  ! Export of calcite at 100m
            CALL iom_put( "EXPCAL"  ,  zw3d )  ! Export of calcite in the water column
         ENDIF
#if defined key_trc_diaadd
         DO_2D( 0, 0, 0, 0 )
            ik100 = nlev100(ji,jj)
            trc2d(ji,jj,jp_sinkcal) = zsinking2(ji,jj,ik100) * xfact * tmask(ji,jj,ik100)  
         END_2D
#endif      
         !
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking2, jpgsi, rfact2 )
         DO_2D( 0, 0, 0, 0 )
            ikb = mbkt(ji,jj)+1
            sinksilb(ji,jj) = zsinking2(ji,jj,ikb)
         END_2D
         IF( l_diag ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk)
               zw3d(ji,jj,jkR) =  zsinking2(ji,jj,jk) * xfact * tmask(ji,jj,jk)
            END_3D
            DO_2D( 0, 0, 0, 0 )
               ik100 = nlev100(ji,jj)
               zw2d(ji,jj) = zsinking2(ji,jj,ik100) * xfact * tmask(ji,jj,ik100)
            END_2D
            CALL iom_put( "EPSI100",  zw2d )  ! Export of Silicate at 100m
            CALL iom_put( "EXPSI"  ,  zw3d )             ! Export of Silicate in the water column
         ENDIF
#if defined key_trc_diaadd
         DO_2D( 0, 0, 0, 0 )
            ik100 = nlev100(ji,jj)
            trc2d(ji,jj,jp_sinksil) = zsinking2(ji,jj,ik100) * xfact * tmask(ji,jj,ik100)  
         END_2D
#endif      
         !
         CALL trc_sink( kt, Kbb, Kmm, wsbio3, zsinking, jpsfe, rfact2 )
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking2, jpbfe, rfact2 )
         IF( l_diag ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk)
              zw3d(ji,jj,jkR) = ( zsinking(ji,jj,jk) + zsinking2(ji,jj,jk) ) &
                      &        * xfact * tmask(ji,jj,jk)
            END_3D
            DO_2D( 0, 0, 0, 0 )
               ik100 = nlev100(ji,jj)
               zw2d(ji,jj) = ( zsinking(ji,jj,ik100) + zsinking2(ji,jj,ik100) ) &
                      &        * xfact * tmask(ji,jj,ik100)
            END_2D
            CALL iom_put( "EPFE100",  zw2d )  ! Export of iron at 100m
            CALL iom_put( "EXPFE"  ,  zw3d )             ! Export of iron in the water column
         ENDIF
#if defined key_trc_diaadd
         DO_2D( 0, 0, 0, 0 )
            ik100 = nlev100(ji,jj)
            trc2d(ji,jj,jp_sinkfer) = ( zsinking(ji,jj,ik100) + zsinking2(ji,jj,ik100) ) &
                           &        * xfact * tmask(ji,jj,ik100)  
         END_2D
#endif      
      ENDIF
      !
      ! PISCES-QUOTA part
      IF( ln_p5z ) THEN
         !
         CALL trc_sink( kt, Kbb, Kmm, wsbio3, zsinking, jppon, rfact2 )
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking2, jpgon, rfact2 )
         DO_2D( 0, 0, 0, 0 )
            ikb = mbkt(ji,jj)+1
            sinkponb(ji,jj) = zsinking(ji,jj,ikb) + zsinking2(ji,jj,ikb)
         END_2D
         !
         CALL trc_sink( kt, Kbb, Kmm, wsbio3, zsinking, jppop, rfact2 )
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking2, jpgop, rfact2 )
         DO_2D( 0, 0, 0, 0 )
            ikb = mbkt(ji,jj)+1
            sinkpopb(ji,jj) = zsinking(ji,jj,ikb) + zsinking2(ji,jj,ikb)
         END_2D
         !
      ENDIF
      !
      IF( l_diag )  DEALLOCATE( zw3d, zw2d )
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
      !   CALL prt_ctl(tab4d_1=tr(:,:,:,:,Kbb), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink


   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!
      !! ** Purpose :   Initialization of sinking parameters
      !!
      !! ** Method  :   
      !!
      !! ** input   :   
      !!----------------------------------------------------------------------
      !
      wsbio3(:,:,:) = wsbio
      wsbio4(:,:,:) = wsbio2
      sinkpocb(:,:) = 0.
      sinkcalb(:,:) = 0.
      IF( .NOT. ln_p2z )  sinksilb(:,:) = 0.
      IF( ln_p5z ) THEN
         sinkponb(:,:) = 0.
         sinkpopb(:,:) = 0.
      ENDIF

   END SUBROUTINE p4z_sink_init

   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(4)
      !!----------------------------------------------------------------------
      !
      ierr(:) = 0
      !
      ALLOCATE( sinkpocb(A2D(0)), sinkcalb(A2D(0)), STAT=ierr(1) )
      ALLOCATE( nlev100(A2D(0)), STAT=ierr(2) )
      !
      IF( .NOT. ln_p2z ) THEN
         ALLOCATE( sinksilb(A2D(0)), STAT=ierr(3) )
         !
         IF( ln_p5z ) &
             &   ALLOCATE( sinkponb(A2D(0)), sinkpopb(A2D(0)), STAT=ierr(4) )
      ENDIF
      !
      p4z_sink_alloc = MAXVAL( ierr )
      IF( p4z_sink_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_sink_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_sink_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sink                    ! Empty routine
   END SUBROUTINE p4z_sink
#endif
   
   !!======================================================================
END MODULE p4zsink
