#include "cppdefs.h"

!!! give access to dimensions, variables, parameters from hydro model to substance

#if defined SUBSTANCE

      module module_substance

   ! pass from fix (F) to free format (F90)

#include "param.h"
#include "ncscrum.h"
#include "scalars.h"
#include "forces.h"
#include "grid.h"
#include "sources.h"
#include "ocean2d.h"
#include "ocean3d.h"  
      ! ocean3d needed for Hz

      end module module_substance
#else
      module module_substance
      end module module_substance
      
#endif /* SUBSTANCE */
