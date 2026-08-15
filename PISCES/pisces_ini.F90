#include "cppdefs.h"

MODULE pisces_ini
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   USE oce_trc
   USE trcini_pisces


   IMPLICIT NONE
   PRIVATE

   PUBLIC   pisces_ini_tile   ! called by trcini.F90 module

   !!* Substitution
#  include "ocean2pisces.h90"

CONTAINS

    SUBROUTINE pisces_ini_tile( tile )

       INTEGER :: tile, Nbb

#include "compute_tile_bounds.h"
      
       Nbb = nnew 
       CALL ocean_2_pisces( nit000,Istr,Iend,Jstr,Jend )
       CALL trc_ini_pisces( Nbb )

    END SUBROUTINE pisces_ini_tile
 
#else
   !!----------------------------------------------------------------------
   !!   Dummy module                            No PISCES biochemical model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE pisces_ini_tile
   END SUBROUTINE pisces_ini_tile
#endif

   !!======================================================================
END MODULE pisces_ini
