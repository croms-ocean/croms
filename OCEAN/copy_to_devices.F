! $Id: copy_to_devices.F 1615 2014-12-17 13:27:07Z rblod $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!

#include "cppdefs.h"
#if defined OPENACC
       subroutine copy_to_devices(tile)
       integer tile
# include "param.h"
# include "compute_tile_bounds.h"
       call copy_to_devices_tile(Istr,Iend,Jstr,Jend)
       return
       end

       subroutine copy_to_devices_tile(Istr,Iend,Jstr,Jend)
       implicit none

# include "param.h"
       integer Istr,Iend,Jstr,Jend
# include "scalars.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "grid.h"
# include "coupling.h"
# include "climat.h"
# if defined PSOURCE || defined PSOURCE_NCFILE
#  include "sources.h"
# endif
# ifdef MPI
#  include "mpi_buffers.h"
# endif
# if defined NBQ
#  include "nbq.h"
# endif

!$acc enter data if(compute_on_device) copyin(h
!$acc$ ,hinv,zeta
!$acc& ,pm_u,pn_u,pm,pn,pn_v,pm_v
!$acc& ,om_v,on_u,on_r,om_r,fomn,pnom_p,pmon_p,om_p,on_p
!$acc& ,on_v,om_u
!$acc& ,pnom_r,pmon_r, pnom_v, pmon_u
!$acc& ,Zob
!$acc& ,f
# ifdef SOLVE3D
!$acc&, z_w,z_r,Hz_bak,Hz,u,v,Zt_avg1
!$acc&, cs_r,sc_r,cs_w,sc_w
# endif
# ifdef MASKING
!$acc& ,rmask, pmask, umask, vmask, pmask, pmask2
# endif
# if defined CURVGRID && defined UV_ADV
!$acc& ,dmde,dndx
# endif
# if defined CURVGRID
!$acc& ,angler
# endif
# ifdef SPHERICAL
!$acc&, latr, lonr, latu, lonu, latv, lonv
# else
!$acc& ,xp,xr,yp,yr
# endif
# ifdef MPI
!$acc&, buf_snd4,buf_snd2,buf_rev4,buf_rev2
!$acc&, buf_snd1,buf_snd3,buf_rev1,buf_rev3
!$acc&, ibuf_sndN, ibuf_revN,
!$acc&     ibuf_sndS, ibuf_revS,
!$acc&     jbuf_sndW, jbuf_sndE,
!$acc&     jbuf_revW, jbuf_revE
#  ifdef SOLVE3D
!$acc&, buf_snd4_3D, buf_snd2_3D,
!$acc&                  buf_rev4_3D, buf_rev2_3D,
!$acc&                  buf_snd1_3D, buf_snd3_3D,
!$acc&                  buf_rev1_3D, buf_rev3_3D
!$acc&, ibuf_sndN_3D, ibuf_revN_3D,
!$acc&     ibuf_sndS_3D, ibuf_revS_3D,
!$acc&     jbuf_sndW_3D, jbuf_sndE_3D,
!$acc&     jbuf_revW_3D, jbuf_revE_3D
#  endif
# endif
# if defined NBQ
!$acc&  , rho_nbq, rho_grd
!$acc&  , wz
#  if defined NBQ_MASS
!$acc&  , Hzr
#  endif
# endif
# if defined ZCLIMATOLOGY || defined AGRIF
!$acc& ,ssh
# endif
# if defined TRACERS && (defined TCLIMATOLOGY || (defined AGRIF && !defined T_FRC_BRY))
!$acc& ,tclm
# endif
# if defined M2CLIMATOLOGY || (defined AGRIF && !defined M2_FRC_BRY)
!$acc&, ubclm, vbclm
# endif
# if defined SOLVE3D && (defined M3CLIMATOLOGY || \
                        (defined AGRIF && !defined M3_FRC_BRY))
!$acc&, uclm, vclm
# endif
# if defined ZNUDGING
!$acc&, Znudgcof
# endif
# if defined M2NUDGING
!$acc&, M2nudgcof
# endif
# if defined M3NUDGING
!$acc&, M3nudgcof
# endif
# if defined TNUDGING
!$acc&, Tnudgcof
# endif
# ifdef ZONAL_NUDGING
!$acc& ,zetazon,ubzon,vbzon,uzon,vzon
!$acc& ,sshzon,ubclmzon,vbclmzon,uclmzon,vclmzon
#  ifdef TRACERS
!$acc& ,tzon,tclmzon
#  endif
# endif
# if defined PSOURCE || defined PSOURCE_NCFILE
!$acc& ,Isrc,Jsrc
#  ifdef MPI
!$acc& ,Isrc_mpi,Jsrc_mpi
#  endif
# endif
# ifdef WET_DRY
!$acc& ,Dcrit, wetdry,rmask_wet
# endif

!$acc& )

       return
       end

       subroutine copy_from_devices(tile)
       integer tile
# include "param.h"
# include "compute_tile_bounds.h"
       call copy_from_devices_tile(Istr,Iend,Jstr,Jend)
       return
       end

       subroutine copy_from_devices_tile(Istr,Iend,Jstr,Jend)
       implicit none
# include "param.h"
       integer Istr,Iend,Jstr,Jend
# include "scalars.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "grid.h"
# include "coupling.h"

# include "compute_auxiliary_bounds.h"

!$acc exit data if(compute_on_device)
!$acc& copyout (f,h
# ifdef SOLVE3D
!$acc&, hinv(IstrR:IendR,JstrR:JendR)
!$acc&, z_w(IstrR:IendR,JstrR:JendR,:)
!$acc& ,z_r(IstrR:IendR,JstrR:JendR,:)
!$acc& ,Hz_bak(IstrR:IendR,JstrR:JendR,:)
!$acc& ,Hz(IstrR:IendR,JstrR:JendR,:)
!$acc& ,Zt_avg1(IstrR:IendR,JstrR:JendR)
# endif
!$acc& ,pm,pn
# ifdef SOLVE3D
!$acc&, cs_r,sc_r,cs_w,sc_w
# endif
# ifdef MASKING
!$acc& ,rmask, pmask, umask, vmask, pmask, pmask2
# endif
# if defined CURVGRID && defined UV_ADV
!$acc& ,dmde,dndx
# endif
# if defined CURVGRID
!$acc& ,angler
# endif
# ifdef SPHERICAL
!$acc&, latr, lonr, latu, lonu, latv, lonv
# else
!$acc& ,xp,xr,yp,yr
# endif
# ifdef WET_DRY
!$acc& , wetdry, Dcrit
# endif
!$acc& )

       return
       end

       subroutine copy_to_devices_2(tile)
       integer tile
# include "param.h"
# include "compute_tile_bounds.h"
       call copy_to_devices_2_tile(Istr,Iend,Jstr,Jend)
       return
       end

       subroutine copy_to_devices_2_tile(Istr,Iend,Jstr,Jend)
       implicit none
# include "param.h"
       integer Istr,Iend,Jstr,Jend
# include "scalars.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "grid.h"
# include "coupling.h"
# include "private_scratch.h"
# include "mixing.h"
# include "forces.h"
# include "work.h"
# include "ncscrum.h"
# include "averages.h"
# include "lmd_kpp.h"
# include "climat.h"
# include "boundary.h"
# include "sources.h"
# if defined TIDES
#  include "tides.h"
# endif
# if defined WKB_WWAVE
#  include "wkb_wwave.h"
# endif
# if defined BBL
#  include "bbl.h"
# endif
# if defined NBQ
#  include "nbq.h"
# endif

# include "copy_to_devices.h"

       return
       end
#else
       subroutine copy_to_devices_empty
       return
       end
#endif
