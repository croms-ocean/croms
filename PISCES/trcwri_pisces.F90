#include "cppdefs.h"

MODULE trcwri_pisces
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    PISCES :   Output of PISCES tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined XIOS && defined key_pisces 
   !!----------------------------------------------------------------------
   !! trc_wri_pisces   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE sms_pisces  ! PISCES variables
   USE trc         ! passive tracers common variables 
   USE sedwri
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_pisces 

   !! * Substitutions
#include "ocean2pisces.h90"
#include "do_loop_substitute.h90"


CONTAINS

   SUBROUTINE trc_wri_pisces
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      CHARACTER (len=20)   :: cltra
      INTEGER              :: ji, jj, jk, jn
      REAL(wp), DIMENSION(GLOBAL_2D_ARRAY,jpk) :: ztra
      !!---------------------------------------------------------------------
 
      IF( ln_sediment )  CALL sed_wri

      ztra(:,:,:) = 0.
      DO jn = jp_pcs0, jp_pcs1
         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         DO_3D( 0, 0, 0, 0, 1, jpk)
            ztra(ji,jj,N+1-jk) = MAX( 0., tr(ji,jj,jk,jn,nnew) )
         END_3D   
         CALL iom_put( cltra, ztra )
      END DO
      !
   END SUBROUTINE trc_wri_pisces
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_pisces
CONTAINS
   SUBROUTINE trc_wri_pisces                     ! Empty routine  
   END SUBROUTINE trc_wri_pisces
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcwri_pisces.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_pisces
