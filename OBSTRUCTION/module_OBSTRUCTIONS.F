#include "cppdefs.h"

!!! give access to dimensions, variables, parameters from hydro model to OBSTRUCTIONS

#if defined OBSTRUCTION

      MODULE module_OBSTRUCTIONS

   ! pass from fix (F) to free format (F90)

#include "param.h"
#include "ncscrum.h"
#include "scalars.h"
#include "forces.h"
#include "grid.h"
#include "ocean3d.h"
#include "mixing.h"
#include "coupling.h"

      END MODULE module_OBSTRUCTIONS
#else
      MODULE module_OBSTRUCTIONS_empty
      END MODULE module_OBSTRUCTIONS_empty
      
#endif /* OBSTRUCTION */