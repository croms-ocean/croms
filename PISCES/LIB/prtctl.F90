#include "cppdefs.h"

MODULE prtctl

#if defined key_pisces

   USE trc

   IMPLICIT NONE
   PUBLIC

   PUBLIC prt_ctl_ini, prt_ctl_info, prt_ctl   

   !! * Substitutions
#  include "ocean2pisces.h90"
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p5zprod.F90 15459 2021-10-29 08:19:18Z cetlod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

      SUBROUTINE prt_ctl_ini

 !     ALLOCATE( tra_ctl(jptra) )
 !     tra_ctl(:) = 0.e0           ! Initialization to zero

   END SUBROUTINE prt_ctl_ini

   SUBROUTINE prt_ctl_info (clinfo, cdcomp )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE prt_ctl_trc_info  ***
      !!
      !! ** Purpose : - print information without any computation
      !!----------------------------------------------------------------------
      CHARACTER(len=*),           INTENT(in) ::   clinfo
      CHARACTER(len=3), OPTIONAL, INTENT(in) ::   cdcomp   ! only 'top' is accepted
      CHARACTER(len=3) :: clcomp

      !!
      IF( PRESENT(cdcomp) ) THEN   ;   clcomp = cdcomp
      ELSE                         ;   clcomp = 'top'
      ENDIF
      !
   END SUBROUTINE prt_ctl_info

   SUBROUTINE prt_ctl( ptab, pmask, clinfo )

      REAL(wp), DIMENSION(:,:,:,:), INTENT(in) :: ptab
      CHARACTER(len=*),  INTENT(in)           :: clinfo   ! information about the tab3d array
      REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL ::   pmask

      INTEGER :: ji, jj, jk, jn
      REAL(wp), DIMENSION(A2D(0),jpk,jptra)        :: ztab
      REAL(wp)  :: zsum  

      DO jn = 1, SIZE( ptab, 4 )
         DO jk = 1, SIZE( ptab, 3 )
            DO jj = 1, SIZE( ptab, 2 )
               DO ji = 1, SIZE( ptab, 1 )
                  ztab(ji,jj,jk,jn) = ptab(ji,jj,jk,jn) * pmask(ji,jj,jk)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO jn = 1, jptra
         zsum   = SUM( ztab(:,:,:,jn) )
         IF( lk_mpp ) CALL mpp_sum( zsum )      ! min over the global domain
         IF( lwp ) WRITE(numout,FMT="(3x,a10,' : ',D23.16)") TRIM(ctrcnm(jn)), zsum
      END DO

   END SUBROUTINE prt_ctl      
# endif
! Empty module
END MODULE prtctl

