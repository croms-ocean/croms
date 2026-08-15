#include "cppdefs.h"

#if defined OBSTRUCTION

MODULE plug_OBSTRUCTIONS
   ! interface between croco and obstruction module

   USE module_OBSTRUCTIONS

   USE init_OBSTRUCTIONS, ONLY: obst_init
   USE OBSTRUCTIONS, ONLY: obst_update

   IMPLICIT NONE

   PRIVATE

   PUBLIC obst_update_main
   PUBLIC obst_init_main

contains
!
!-----------------------------------------------------------------------
   SUBROUTINE obst_update_main(tile)

      INTEGER :: tile
# include "ocean2d.h"
# include "compute_tile_bounds.h"

      CALL obst_update(Istr, Iend, Jstr, Jend)

   END SUBROUTINE

!-----------------------------------------------------------------------
   SUBROUTINE obst_init_main()

      REAL :: h0fond
# include "ocean2d.h"

# ifdef WET_DRY
      h0fond = D_wetdry
# else
      h0fond = 0.
# endif

      CALL obst_init(h0fond)
   END SUBROUTINE

!-----------------------------------------------------------------------
END MODULE plug_OBSTRUCTIONS

#else

MODULE plug_OBSTRUCTIONS
END MODULE plug_OBSTRUCTIONS

#endif /* OBSTRUCTIONS */
