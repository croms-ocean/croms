! $Id: wrt_rst.F 1571 2014-07-01 12:38:05Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"

MODULE sedwri

#if defined key_sediment

   USE sed
   USE sedinorg
   USE lib_mpp         ! distribued memory computing library
   USE iom
#ifdef AGRIF
      USE param, ONLY : Lmmpi,Mmmpi
#endif

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC sed_wri         ! routine called by opa.F90

   !! * Substitutions
#  include "ocean2pisces.h90"
#  include "do_loop_substitute.h90"

CONTAINS   

   SUBROUTINE sed_wri
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_wri  ***
      !!
      !! ** Purpose :  output of sediment passive tracer
      !!
      !!   History :
      !!        !  06-07  (C. Ethe)  original
      !!----------------------------------------------------------------------
      !
      INTEGER  :: ji, jj, jk, jn
      INTEGER  :: it
      CHARACTER(len = 20)  ::  cltra
      REAL(wp), DIMENSION(jpoce, jptrased+1) :: zflx
      REAL(wp), DIMENSION(A2D(0), jpksed)   :: trcsedi
      REAL(wp), DIMENSION(A2D(0)) :: flxsedi2d
      REAL(wp), DIMENSION(A2D(0), jpksed) :: zw3d
      REAL(wp), DIMENSION(jpoce,jpksed) :: zdta
      !
      REAL(wp), DIMENSION(GLOBAL_2D_ARRAY) :: zw2d_out
      REAL(wp), DIMENSION(GLOBAL_2D_ARRAY,jpksed) :: zw3d_out
      !!-------------------------------------------------------------------

      ! 1.  Initilisations
      ! -----------------------------------------------------------------
      IF( ln_timing )  CALL timing_start('sed_wri')
!
      !
      zw3d_out(:,:,:) = 0._wp 
      zw2d_out(:,:)   = 0._wp 


      ! 2.  Back to 2D geometry
      ! -----------------------------------------------------------------
      ! Calculation of fluxes mol/cm2/s
      DO jn = 1, jpwat
         DO ji = 1, jpoce
            zflx(ji,jn) = ( pwcp(ji,1,jn) - pwcp_dta(ji,jn) ) * ( 1.e-3 * dzkbot(ji) ) / dtsed
         ENDDO
      ENDDO

      ! Calculation of fluxes g/cm2/s
      ! Calculation of accumulation rate per dt
      zflx(:,jptrased+1) = 0.0
      DO jn = 1, jpsol
         DO ji = 1, jpoce
            zflx(ji,jpwat+jn) = ( tosed(ji,jn) - fromsed(ji,jn) ) / dtsed
            zflx(ji,jptrased+1) = zflx(ji,jptrased+1) + ( tosed(ji,jn) - fromsed(ji,jn) ) / ( dtsed * por1(jpksed) * dens_sol(jn) )
         ENDDO
      ENDDO
      
     !
      ! Start writing data
      ! ---------------------
     DO jn = 1, jptrased
         cltra = TRIM( sedtrcd(jn) ) ! short title for 3D diagnostic
         IF ( iom_use( cltra ) ) THEN
            IF ( jn <= jpsol ) THEN
               DO jk = 1, jpksed
                  trcsedi(:,:,jk) = UNPACK( solcp(:,jk,jn), sedmask == 1.0, 0.0 )
               END DO
            ELSE
               DO jk = 1, jpksed
                  trcsedi(:,:,jk) = UNPACK( pwcp(:,jk,jn-jpsol)*1E6, sedmask == 1.0, 0.0 )
               END DO
            ENDIF
            zw3d_out(A2D(0),:) = trcsedi(A2D(0),:)
            CALL iom_put( cltra, zw3d_out )
         ENDIF
      END DO

      DO jn = 1, jptrased+1
         cltra = TRIM( seddia2d(jn) ) ! short title for 2D diagnostic
         IF ( iom_use( cltra ) ) THEN
            flxsedi2d(:,:) = UNPACK( zflx(:,jn), sedmask == 1.0, 0.0 )
            zw2d_out(A2D(0)) = flxsedi2d(A2D(0))
            CALL iom_put( cltra, zw2d_out )
         ENDIF
      END DO

      IF ( iom_use( "dzdep" ) ) THEN
         zflx(:,1) = dzdep(:) / dtsed
         flxsedi2d(:,:) = UNPACK( zflx(:,1), sedmask == 1.0, 0.0 )
         zw2d_out(A2D(0)) = flxsedi2d(A2D(0))
         CALL iom_put( "dzdep", zw2d_out )
      ENDIF

      IF ( iom_use( "Rstepros" ) ) THEN
         flxsedi2d(:,:) = UNPACK( rstepros(:), sedmask == 1.0, 0.0 )
         zw2d_out(A2D(0)) = flxsedi2d(A2D(0))
         CALL iom_put( "Rstepros", zw2d_out )
      ENDIF

      IF ( iom_use( "SaturCO3" ) ) THEN
         DO jk = 1, jpksed
            DO ji = 1, jpoce
               saturco3(ji,jk) = (1.0 - co3por(ji,jk) /  co3sat(ji) )
            END DO
            zw3d(:,:,jk) = UNPACK( saturco3(:,jk), sedmask == 1.0, 0.0)
         END DO
         zw3d_out(A2D(0),:) = zw3d(A2D(0),:)
         CALL iom_put( "SaturCO3", zw3d_out )
      ENDIF
      IF ( iom_use( "SedCO3por" ) ) THEN
         DO jk = 1, jpksed
            zw3d(:,:,jk) = UNPACK( co3por(:,jk), sedmask == 1.0, 0.0)
         END DO
         zw3d_out(A2D(0),:) = zw3d(A2D(0),:)
         CALL iom_put( "SedCO3por", zw3d_out )
      ENDIF
      IF ( iom_use( "SedpH" ) ) THEN
         DO jk = 1, jpksed
            DO ji = 1, jpoce
               zdta(ji,jk) = -LOG10( hipor(ji,jk) / ( densSW(ji) + rtrn ) + rtrn )
            END DO
            zw3d(:,:,jk) = UNPACK( zdta(:,jk), sedmask == 1.0, 0.0)
         END DO
         zw3d_out(A2D(0),:) = zw3d(A2D(0),:)
         CALL iom_put( "SedpH", zw3d_out )
      ENDIF

      IF( ln_timing )  CALL timing_stop('sed_wri')

   END SUBROUTINE sed_wri
#else

CONTAINS

   SUBROUTINE sed_wri
   END SUBROUTINE sed_wri

#endif

END MODULE sedwri
