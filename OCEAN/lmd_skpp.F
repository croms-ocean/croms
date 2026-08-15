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

#ifdef LMD_SKPP
# if defined LMD_SKPP2005
#  include "lmd_skpp2005.F"
# else
#  include "lmd_skpp1994.F"
# endif
#else
!========================
      subroutine lmd_skpp_empty
      return
      end
!========================
#endif
