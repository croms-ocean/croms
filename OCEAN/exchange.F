#include "cppdefs.h"
#ifdef MPI
# define LOCALLM Lmmpi
# define LOCALMM Mmmpi
#else
# define LOCALLM Lm
# define LOCALMM Mm
#endif
#if defined BD_EXCHANGE

# define exchange_2d_tile exchange_r2d_tile
# define exchange_2d_3pts_tile exchange_r2d_3pts_tile
# define exchange_2d_1pts_tile exchange_r2d_1pts_tile
# define JSTART JstrR
# define ISTART IstrR
# include "exchange_2d_tile.h"
# undef ISTART
# undef JSTART
# undef exchange_2d_tile
# undef exchange_2d_3pts_tile
# undef exchange_2d_1pts_tile

# define exchange_2d_tile exchange_u2d_tile
# define exchange_2d_3pts_tile exchange_u2d_3pts_tile
# define exchange_2d_1pts_tile exchange_u2d_1pts_tile
# define JSTART JstrR
# define ISTART Istr
# include "exchange_2d_tile.h"
# undef ISTART
# undef JSTART
# undef exchange_2d_tile
# undef exchange_2d_3pts_tile
# undef exchange_2d_1pts_tile

# define exchange_2d_tile exchange_v2d_tile
# define exchange_2d_3pts_tile exchange_v2d_3pts_tile
# define exchange_2d_1pts_tile exchange_v2d_1pts_tile
# define JSTART Jstr
# define ISTART IstrR
# include "exchange_2d_tile.h"
# undef ISTART
# undef JSTART
# undef exchange_2d_tile
# undef exchange_2d_3pts_tile
# undef exchange_2d_1pts_tile

# define exchange_2d_tile exchange_p2d_tile
# define exchange_2d_3pts_tile exchange_p2d_3pts_tile
# define exchange_2d_1pts_tile exchange_p2d_1pts_tile
# define JSTART Jstr
# define ISTART Istr
# include "exchange_2d_tile.h"
# undef ISTART
# undef JSTART
# undef exchange_2d_tile
# undef exchange_2d_3pts_tile
# undef exchange_2d_1pts_tile

# ifdef SOLVE3D
#  define KSTART 1
#  define exchange_3d_tile exchange_r3d_tile
#  define exchange_3d_3pts_tile exchange_r3d_3pts_tile
#  define JSTART JstrR
#  define ISTART IstrR
#  include "exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef exchange_3d_tile
#  undef exchange_3d_3pts_tile

#  define exchange_3d_tile exchange_u3d_tile
#  define exchange_3d_3pts_tile exchange_u3d_3pts_tile
#  define JSTART JstrR
#  define ISTART Istr
#  include "exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef exchange_3d_tile
#  undef exchange_3d_3pts_tile

#  define exchange_3d_tile exchange_v3d_tile
#  define exchange_3d_3pts_tile exchange_v3d_3pts_tile
#  define JSTART Jstr
#  define ISTART IstrR
#  include "exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef exchange_3d_tile
#  undef exchange_3d_3pts_tile

#  define exchange_3d_tile exchange_p3d_tile
#  define exchange_3d_3pts_tile exchange_p3d_3pts_tile
#  define JSTART Jstr
#  define ISTART Istr
#  include "exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef exchange_3d_tile
#  undef exchange_3d_3pts_tile
#  undef KSTART

#  define KSTART 0
#  define exchange_3d_tile exchange_w3d_tile
#  define exchange_3d_3pts_tile exchange_w3d_3pts_tile
#  define JSTART JstrR
#  define ISTART IstrR
#  include "exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef exchange_3d_tile
#  undef exchange_3d_3pts_tile
#  undef KSTART



#  ifdef M3FAST_SEDLAYERS
#   define MP_M3FAST_SEDLAYERS
#   define KSTART (-N_sl+1)
#   define exchange_3d_tile exchange_r3d_sedlay_tile
#   define exchange_3d_3pts_tile exchange_r3d_3pts_sedlay_tile
#   define JSTART JstrR
#   define ISTART IstrR
#   include "exchange_3d_tile.h"
#   undef ISTART
#   undef JSTART
#   undef exchange_3d_tile
#   undef exchange_3d_3pts_tile

#   define exchange_3d_tile exchange_u3d_sedlay_tile
#   define exchange_3d_3pts_tile exchange_u3d_sedlay_3pts_tile
#   define JSTART JstrR
#   define ISTART Istr
#   include "exchange_3d_tile.h"
#   undef ISTART
#   undef JSTART
#   undef exchange_3d_tile
#   undef exchange_3d_3pts_tile

#   define exchange_3d_tile exchange_v3d_sedlay_tile
#   define exchange_3d_3pts_tile exchange_v3d_sedlay_3pts_tile
#   define JSTART Jstr
#   define ISTART IstrR
#   include "exchange_3d_tile.h"
#   undef ISTART
#   undef JSTART
#   undef exchange_3d_tile
#   undef exchange_3d_3pts_tile

#   define exchange_3d_tile exchange_p3d_sedlay_tile
#   define exchange_3d_3pts_tile exchange_p3d_sedlay_3pts_tile
#   define JSTART Jstr
#   define ISTART Istr
#   include "exchange_3d_tile.h"
#   undef ISTART
#   undef JSTART
#   undef exchange_3d_tile
#   undef exchange_3d_3pts_tile

#   define exchange_3d_tile exchange_w3d_sedlay_tile
#   define exchange_3d_3pts_tile exchange_w3d_sedlay_3pts_tile
#   define JSTART JstrR
#   define ISTART IstrR
#   include "exchange_3d_tile.h"
#   undef ISTART
#   undef JSTART
#   undef exchange_3d_tile
#   undef exchange_3d_3pts_tile

#   undef KSTART

#   define KSTART (-N_sl)
#   define exchange_3d_tile exchange_r3d_sedlayb_tile
#   define exchange_3d_3pts_tile exchange_r3d_3pts_sedlayb_tile
#   define JSTART JstrR
#   define ISTART IstrR
#   include "exchange_3d_tile.h"
#   undef ISTART
#   undef JSTART
#   undef exchange_3d_tile
#   undef exchange_3d_3pts_tile

#   define exchange_3d_tile exchange_u3d_sedlayb_tile
#   define exchange_3d_3pts_tile exchange_u3d_sedlayb_3pts_tile
#   define JSTART JstrR
#   define ISTART Istr
#   include "exchange_3d_tile.h"
#   undef ISTART
#   undef JSTART
#   undef exchange_3d_tile
#   undef exchange_3d_3pts_tile

#   define exchange_3d_tile exchange_v3d_sedlayb_tile
#   define exchange_3d_3pts_tile exchange_v3d_sedlayb_3pts_tile
#   define JSTART Jstr
#   define ISTART IstrR
#   include "exchange_3d_tile.h"
#   undef ISTART
#   undef JSTART
#   undef exchange_3d_tile
#   undef exchange_3d_3pts_tile

#   define exchange_3d_tile exchange_p3d_sedlayb_tile
#   define exchange_3d_3pts_tile exchange_p3d_sedlayb_3pts_tile
#   define JSTART Jstr
#   define ISTART Istr
#   include "exchange_3d_tile.h"
#   undef ISTART
#   undef JSTART
#   undef exchange_3d_tile
#   undef exchange_3d_3pts_tile

#   define exchange_3d_tile exchange_w3d_sedlayb_tile
#   define exchange_3d_3pts_tile exchange_w3d_sedlayb_3pts_tile
#   define JSTART JstrR
#   define ISTART IstrR
#   include "exchange_3d_tile.h"
#   undef ISTART
#   undef JSTART
#   undef exchange_3d_tile
#   undef exchange_3d_3pts_tile

#   undef KSTART
#   undef MP_M3FAST_SEDLAYERS
# endif
!! <---- M3FAST_SEDLAYERS

# endif
#else
      subroutine exchange_empty
      end
#endif




