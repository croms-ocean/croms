#include "cppdefs.h"

MODULE sed_oce
   !!======================================================================
   !!                        ***  sed  ***
   !! Sediment :   set sediment global variables
   !!======================================================================

   !! History :
   !!        !  06-12  (C. Ethe)  Orignal
   !!----------------------------------------------------------------------
#if defined key_sediment
   USE par_sed
   USE trc, ONLY : profsed

   !!* Substitution
#  include "ocean2pisces.h90"

   IMPLICIT NONE
   PUBLIC

   PUBLIC sed_oce_alloc

   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  profsedw       !: depth of middle of each layer


   !! $Id: sed.F90 7646 2017-02-06 09:25:03Z timgraham $
CONTAINS

   INTEGER FUNCTION sed_oce_alloc()
      !!-------------------------------------------------------------------
      !!                    *** ROUTINE sed_alloc ***
      !!-------------------------------------------------------------------
      USE lib_mpp, ONLY: ctl_stop
      !!-------------------------------------------------------------------
      !
      ALLOCATE( profsed(jpksed) , profsedw(jpksed) , STAT=sed_oce_alloc )

      IF( sed_oce_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sed_oce_alloc: failed to allocate arrays' )
      !
   END FUNCTION sed_oce_alloc

#endif

END MODULE sed_oce
