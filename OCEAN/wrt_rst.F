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
                              ! Write model prognostic
      subroutine wrt_rst      ! variables into restart
                              ! netCDF file.

#if defined MUSTANG
      USE comMUSTANG, only: rstMust, vname_rstMust, ksmi, ksma, cv_sed, dzs
      USE dredging, ONLY : l_dredging, dredg_hsed_init
#endif

!     implicit none
#include "param.h"
#include "scalars.h"
#include "ncscrum.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "forces.h"
#include "netcdf.inc"
#include "mpi_cpl.h"
#if defined GLS_MIXING || defined LMD_MIXING
# include "mixing.h"
#endif
#if defined EXACT_RESTART
# include "coupling.h"
# if defined M3FAST
#  include "nbq.h"
# endif
#endif
#ifdef SEDIMENT
# include "sediment.h"
#endif
#ifdef MORPHODYN
# include "grid.h"
#endif
#ifdef BBL
# include "bbl.h"
#endif
#ifdef ABL1D
# include "abl1d.h"
#endif
      integer ierr, record, lstr, lvar, lenstr
     &  , start(2), count(2), nf_fwrite
#if defined OUTPUTS_SURFACE && ! defined XIOS
     &   , ibuff(6)
#else
     &   , ibuff(4)
#endif
#ifdef SOLVE3D
     &      , itrc
#endif
#ifdef SEDIMENT
     &      , indxWrk
#endif
#if defined MPI & !defined PARALLEL_FILES
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE), blank
#endif

#if defined MUSTANG
      real workrsed(GLOBAL_2D_ARRAY,ksdmax)
      real work2d(GLOBAL_2D_ARRAY)
#endif

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, MPI_COMM_WORLD, status, ierr)
      endif
#endif

!
! Create/open restart file; write grid arrays, if so needed.
!
      call def_rst (ncidrst, nrecrst, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(rstname)
!                                            !!! WARNING: Here it is
! Set record within the file.                !!! assumed that global
!                                            !!! restart record index
      nrecrst=max(nrecrst,1)                 !!! nrecrst is already
      if (nrpfrst.eq.0) then                 !!! advanced by main.
        record=nrecrst
      else
        record=1+mod(nrecrst-1, abs(nrpfrst))
      endif

!#define CR
!
! Write out evolving model variables:
! ----- --- -------- ----- ----------
!
! Time step number and record indices.
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
#ifdef AVERAGES
      ibuff(4)=nrecavg
#else
      ibuff(4)=0
#endif
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (ncidrst, rstTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        MPI_master_only write(stdout,1) 'time_step', record, ierr
     &                   MYID
        goto 99                                           !--> ERROR
      endif

#if defined OUTPUTS_SURFACE && ! defined XIOS
      count(1)=6
      ibuff(5)=nrecsurf
# ifdef AVERAGES
      ibuff(6)=nrecsurf_avg
# else
      ibuff(6)=0
# endif
      start(1)=1
      start(2)=record
      count(2)=1
      ierr=nf_put_vara_int (ncidrst, rstTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        MPI_master_only write(stdout,1) 'time_step_surf', record, ierr
     &                   MYID
      endif
#endif


!
! Time.
!
      ! scrum_time
      ierr=nf_put_var1_FTYPE (ncidrst, rstTime, record, time
#ifdef USE_CALENDAR
     &  - origin_date_in_sec
#endif
     &  )
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

      ! time
      ierr=nf_put_var1_FTYPE (ncidrst, rstTime2, record, time
#ifdef USE_CALENDAR
     &  - origin_date_in_sec
#endif
     &  )
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime2))
        write(stdout,1) vname(1,indxTime2)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

!
! Free-surface.
!
!$acc update if(compute_on_device) host(zeta(:,:,fast_indx_out))
      ierr=nf_fwrite(zeta(START_2D_ARRAY,fast_indx_out), ncidrst,
     &                                      rstZ, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxZ))
        write(stdout,1) vname(1,indxZ)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

!
! 2D momentum component in XI-direction.
!
!$acc update if(compute_on_device) host(ubar(:,:,fast_indx_out))
      ierr=nf_fwrite(ubar(START_2D_ARRAY,fast_indx_out), ncidrst,
     &                                     rstUb, record, u2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxUb))
        write(stdout,1) vname(1,indxUb)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! 2D momentum component in ETA-direction.
!
!$acc update if(compute_on_device) host(vbar(:,:,fast_indx_out))
      ierr=nf_fwrite(vbar(START_2D_ARRAY,fast_indx_out), ncidrst,
     &                                     rstVb, record, v2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxVb))
        write(stdout,1) vname(1,indxVb)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

#ifdef ABL1D
!
! 3D U ABL
!
      ierr=nf_fwrite(u_abl(START_2D_ARRAY,1,nt_a), ncidrst, rstAblU,
     &                                        record, abl3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_abl1d(1,indxabl_u_abl))
        write(stdout,1) vname_abl1d(1,indxabl_u_abl)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! 3D V ABL
!
      ierr=nf_fwrite(v_abl(START_2D_ARRAY,1,nt_a), ncidrst, rstAblV,
     &                                        record, abl3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_abl1d(1,indxabl_v_abl))
        write(stdout,1) vname_abl1d(1,indxabl_v_abl)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! 3D T ABL
!
      ierr=nf_fwrite(tq_abl(START_2D_ARRAY,1,nt_a,jp_ta), ncidrst, rstAblT,
     &                                        record, abl3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_abl1d(1,indxabl_t_abl))
        write(stdout,1) vname_abl1d(1,indxabl_t_abl)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! 3D Q ABL
!
      ierr=nf_fwrite(tq_abl(START_2D_ARRAY,1,nt_a,jp_qa), ncidrst, rstAblQ,
     &                                        record, abl3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_abl1d(1,indxabl_q_abl))
        write(stdout,1) vname_abl1d(1,indxabl_q_abl)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif                  
!
! 3D TKE ABL
!
      ierr=nf_fwrite(tke_abl(START_2D_ARRAY,1,nt_a), ncidrst, rstAblTke,
     &                                        record, abl3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_abl1d(1,indxabl_tke_abl))
        write(stdout,1) vname_abl1d(1,indxabl_tke_abl)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! 3D AVM ABL
!
      ierr=nf_fwrite(avm_abl(START_2D_ARRAY,1,nt_a), ncidrst, rstAblAvm,
     &                                        record, abl3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_abl1d(1,indxabl_avm_abl))
        write(stdout,1) vname_abl1d(1,indxabl_avm_abl)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! 3D AVT ABL
!
      ierr=nf_fwrite(avt_abl(START_2D_ARRAY,1,nt_a), ncidrst, rstAblAvt,
     &                                        record, abl3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_abl1d(1,indxabl_avt_abl))
        write(stdout,1) vname_abl1d(1,indxabl_avt_abl)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif     
!
! 3D MXLD ABL
!
      ierr=nf_fwrite(mxld_abl(START_2D_ARRAY,1,nt_a), ncidrst, rstAblMxld,
     &                                        record, abl3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_abl1d(1,indxabl_mxld_abl))
        write(stdout,1) vname_abl1d(1,indxabl_mxld_abl)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! 3D MXLM ABL
!
      ierr=nf_fwrite(mxlm_abl(START_2D_ARRAY,1,nt_a), ncidrst, rstAblMxlm,
     &                                        record, abl3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_abl1d(1,indxabl_mxlm_abl))
        write(stdout,1) vname_abl1d(1,indxabl_mxlm_abl)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif         
!
! 2D ABLH
!
      ierr=nf_fwrite(ablh(START_2D_ARRAY,nt_a), ncidrst, rstAblAblh,
     &                                        record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_abl1d(1,indxabl_ablh_abl))
        write(stdout,1) vname_abl1d(1,indxabl_ablh_abl)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif               
#endif

#ifdef SOLVE3D
!
! 3D momentum component in XI-direction.
!
!$acc update if(compute_on_device) host(u(:,:,:,nstp))
      ierr=nf_fwrite(u(START_2D_ARRAY,1,nstp), ncidrst, rstU,
     &                                        record, u3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxU))
        write(stdout,1) vname(1,indxU)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! 3D momentum component in ETA-direction.
!
!$acc update if(compute_on_device) host(v(:,:,:,nstp))
      ierr=nf_fwrite(v(START_2D_ARRAY,1,nstp), ncidrst, rstV,
     &                                        record, v3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxV))
        write(stdout,1) vname(1,indxV)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Tracer variables.
!
# ifdef TRACERS
!$acc update if(compute_on_device) host(t(:,:,:,nstp,:))
      do itrc=1,NT
# ifdef DIAGNOSTICS_TS_MLD_CRIT
        if (itrc.ne.iCRT2 .and. itrc.ne.iCRT3
     &  .and. itrc.ne.iCRT4) then
# endif
        ierr=nf_fwrite(t(START_2D_ARRAY,1,nstp,itrc), ncidrst,
     &                             rstT(itrc), record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxV+itrc))
          write(stdout,1) vname(1,indxV+itrc)(1:lvar), record,
     &                    ierr MYID
          goto 99                                         !--> ERROR
        endif
# ifdef DIAGNOSTICS_TS_MLD_CRIT
       endif
# endif
      enddo
# endif

# if defined LMD_MIXING
!
!     Write out depth of planetary boundary layer (m).
!
#  ifdef LMD_SKPP
#   ifdef LMD_SKPP2005
!$acc update if(compute_on_device) host(hbls(:,:,3-nstp))
        ierr=nf_fwrite(hbls(START_2D_ARRAY,3-nstp), ncidrst,
     &                             rstHbl, record, r2dvar)
!#   else
!        ierr=nf_fwrite(hbl(START_2D_ARRAY), ncidrst,
!     &                             rstHbl, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxHbl))
          write(stdout,1) vname(1,indxHbl)(1:lvar), record,
     &         ierr MYID
          goto 99               !--> ERROR
        endif
#   endif /* LMD_SKPP2005  */
#  endif /* LMD_SKPP  */

#  if defined LMD_BKPP
#   if defined LMD_BKPP2005
!
!     Write out depth of bottom planetary boundary layer (m).
!
!$acc update if(compute_on_device) host(hbbl(:,:))
        ierr=nf_fwrite(hbbl(START_2D_ARRAY), ncidrst,
     &                             rstHbbl, record, r2dvar)
        if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxHbbl))
            write(stdout,1) vname(1,indxHbbl)(1:lvar), record,
     &                    ierr MYID
            goto 99                                         !--> ERROR
        endif
#   endif /* LMD_BKPP2005  */
#  endif /* LMD_BKPP  */
# endif /* LMD_MIXING */

# if defined GLS_MIXING
!
! Turbulent kinetic energy.
!
!$acc update if(compute_on_device) host(trb(:,:,nstp,:))
        ierr=nf_fwrite(trb(START_2D_ARRAY,0,nstp,itke), ncidrst,
     &                             rstTke, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxTke))
          write(stdout,1) vname(1,indxTke)(1:lvar), record,
     &                    ierr MYID
          goto 99                                         !--> ERROR
        endif
!
! Generic length scale.
!
        ierr=nf_fwrite(trb(START_2D_ARRAY,0,nstp,igls), ncidrst,
     &                             rstGls, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxGls))
          write(stdout,1) vname(1,indxGls)(1:lvar), record,
     &                    ierr MYID
          goto 99                                         !--> ERROR
        endif
!
! Vertical viscosity coefficient.
!
!$acc update if(compute_on_device) host(Akv(:,:,:))
        ierr=nf_fwrite(Akv(START_2D_ARRAY,0), ncidrst,
     &                             rstAkv, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkv))
          write(stdout,1) vname(1,indxAkv)(1:lvar), record,
     &                    ierr MYID
          goto 99                                         !--> ERROR
        endif
#  ifdef TEMPERATURE
!
! Write out vertical diffusion coefficient for potential temperature.
!
!$acc update if(compute_on_device) host(Akt(:,:,:,itemp))
        ierr=nf_fwrite(Akt(START_2D_ARRAY,0,itemp), ncidrst,
     &                             rstAkt, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkt))
          write(stdout,1) vname(1,indxAkt)(1:lvar), record,
     &                    ierr MYID
          goto 99                                         !--> ERROR
        endif
#  endif /* TEMPERATURE  */
#  ifdef SALINITY
!
! Write out vertical diffusion coefficient for potential temperature.
!
!$acc update if(compute_on_device) host(Akt(:,:,:,isalt))
        ierr=nf_fwrite(Akt(START_2D_ARRAY,0,isalt), ncidrst,
     &                             rstAks, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAks))
          write(stdout,1) vname(1,indxAks)(1:lvar), record,
     &                    ierr MYID
          goto 99                                         !--> ERROR
        endif
#  endif /* SALINITY  */
# endif  /* GLS_MIXING  */

# if defined LMD_MIXING
#  if defined M3FAST
!
! Write out bottom stress u component
!
!$acc update if(compute_on_device) host(bustr(:,:))
        ierr=nf_fwrite(bustr(START_2D_ARRAY), ncidrst,
     &                             rstBustr, record, u2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxBustr))
          write(stdout,1) vname(1,indxBustr)(1:lvar), record,
     &                    ierr MYID
          goto 99                                         !--> ERROR
        endif
!
! Write out bottom stress v component
!
!$acc update if(compute_on_device) host(bvstr(:,:))
        ierr=nf_fwrite(bvstr(START_2D_ARRAY), ncidrst,
     &                             rstBvstr, record, v2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxBvstr))
          write(stdout,1) vname(1,indxBvstr)(1:lvar), record,
     &                    ierr MYID
          goto 99                                         !--> ERROR
        endif
#  endif
# endif

# ifdef EXACT_RESTART
!
! Write out forcing for barotropic equations in XI-direction.
!
        ierr=nf_fwrite(rufrc_bak(START_2D_ARRAY,3-nstp), ncidrst,
     &                             rstrufrc, record, u2dvar)

        if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxrufrc))
            write(stdout,1) vname(1,indxrufrc)(1:lvar), record,
     &                    ierr MYID
            goto 99                                         !--> ERROR
        endif
!
! Write out forcing for barotropic equations in XI-direction.
!
        ierr=nf_fwrite(rvfrc_bak(START_2D_ARRAY,3-nstp), ncidrst,
     &                             rstrvfrc, record, v2dvar)

        if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxrvfrc))
            write(stdout,1) vname(1,indxrvfrc)(1:lvar), record,
     &                    ierr MYID
            goto 99                                         !--> ERROR
          endif

#  ifdef M3FAST
!
! 3D rhs for M3FAST in XI-direction.
!
      ierr=nf_fwrite(ru_nbq(START_2D_ARRAY,1),
     &                               ncidrst, rstru_nbq,
     &                               record, u3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxru_nbq))
        write(stdout,1) vname(1,indxru_nbq)(1:lvar),
     &                  record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! 3D rhs for M3FAST in ETA-direction.
!
       ierr=nf_fwrite(rv_nbq(START_2D_ARRAY,1),
     &                               ncidrst, rstrv_nbq,
     &                               record, v3dvar)
       if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxrv_nbq))
         write(stdout,1) vname(1,indxrv_nbq)(1:lvar),
     &                  record, ierr
     &                  MYID
         goto 99                                           !--> ERROR
       endif

!
! 3D rhs for M3FAST in XI-direction.
!
      ierr=nf_fwrite(ru_nbq_avg2(START_2D_ARRAY,1),
     &                               ncidrst, rstru_nbq_avg2,
     &                               record, u3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxru_nbq_avg2))
        write(stdout,1) vname(1,indxru_nbq_avg2)(1:lvar),
     &                  record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! 3D rhs for M3FAST in ETA-direction.
!
       ierr=nf_fwrite(rv_nbq_avg2(START_2D_ARRAY,1),
     &                               ncidrst, rstrv_nbq_avg2,
     &                               record, v3dvar)
       if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxrv_nbq_avg2))
         write(stdout,1) vname(1,indxrv_nbq_avg2)(1:lvar),
     &                  record, ierr
     &                  MYID
         goto 99                                           !--> ERROR
       endif

!
! 3D rhs for M3FAST in XI-direction.
!
      ierr=nf_fwrite(qdmu_nbq(START_2D_ARRAY,1),
     &                               ncidrst, rstqdmu_nbq,
     &                               record, u3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxqdmu_nbq))
        write(stdout,1) vname(1,indxqdmu_nbq)(1:lvar),
     &                  record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! 3D rhs for M3FAST in ETA-direction.
!
       ierr=nf_fwrite(qdmv_nbq(START_2D_ARRAY,1),
     &                               ncidrst, rstqdmv_nbq,
     &                               record, v3dvar)
       if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxqdmv_nbq))
         write(stdout,1) vname(1,indxqdmv_nbq)(1:lvar),
     &                  record, ierr
     &                  MYID
         goto 99                                           !--> ERROR
       endif

#  endif  /* M3FAST */

#  ifdef TS_MIX_ISO_FILT
!
! density gradients for use in t3dmix
!

      ierr=nf_fwrite(dRdx(START_2D_ARRAY,1),
     &                               ncidrst, rstdRdx,
     &                               record, u3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxdRdx))
        write(stdout,1) vname(1,indxdRdx)(1:lvar),
     &                  record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

       ierr=nf_fwrite(dRde(START_2D_ARRAY,1),
     &                               ncidrst, rstdRde,
     &                               record, v3dvar)
       if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxdRde))
         write(stdout,1) vname(1,indxdRde)(1:lvar),
     &                  record, ierr
     &                  MYID
         goto 99                                           !--> ERROR
       endif

#  endif /* TS_MIX_ISO_FILT */
# endif /* EXACT_RESTART */

# ifdef SEDIMENT
!
! Write out sediment bed layer thickness, porosity, volume
! fraction of size class in sediment bed (2+2*NST b3dgrd variables)
!
      ierr=nf_fwrite(bed_thick, ncidrst, rstSed(1),
     &                                       record, b3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxBTHK))
        write(stdout,1) vname(1,indxBTHK)(1:lvar), record, ierr
     &                  MYID
        goto 99                                         !--> ERROR
      endif
CR      write(*,*) 'wrt_rst:  Bed_thick' MYID
      ierr=nf_fwrite(bed_poros, ncidrst, rstSed(2),
     &                                       record, b3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxBPOR))
        write(stdout,1) vname(1,indxBPOR)(1:lvar), record, ierr
     &                  MYID
        goto 99                                         !--> ERROR
      endif
CR      write(*,*) 'wrt_rst:  Bed_poros' MYID

      do itrc=1,NST
        indxWrk=indxBFRA(1)+itrc-1
        ierr=nf_fwrite(bed_frac(START_2D_ARRAY,1,itrc),
     &             ncidrst, rstSed(itrc+2), record, b3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxWrk))
          write(stdout,1) vname(1,indxWrk)(1:lvar), record,
     &                ierr
     &                MYID
          goto 99                                         !--> ERROR
        endif
      enddo
CR      write(*,*) 'wrt_rst:  Bed_frac' MYID

# endif /* SEDIMENT */

# ifdef MUSTANG

! ksmi, ksma
      work2d=ksmi
      indx=1
      ierr=nf_fwrite(work2d, ncidrst, rstMust(1),
     &                                       record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_rstMust(1,indx))
        write(stdout,1) vname_rstMust(1,indx)(1:lvar), record, ierr
     &                  MYID
        goto 99                                         !--> ERROR
      endif

      work2d=ksma
      indx=indx+1
      ierr=nf_fwrite(work2d, ncidrst, rstMust(2),
     &                                       record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_rstMust(1,indx))
        write(stdout,1) vname_rstMust(1,indx)(1:lvar), record, ierr
     &                  MYID
        goto 99                                         !--> ERROR
      endif

! dzs
      do k=1,ksdmax
        workrsed(:,:,k)=dzs(k,:,:)
      enddo
      indx=indx+1
      ierr=nf_fwrite(workrsed, ncidrst, rstMust(3),
     &                                       record, b3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_rstMust(1,indx))
        write(stdout,1) vname_rstMust(1,indx)(1:lvar), record, ierr
     &                  MYID
        goto 99                                         !--> ERROR
      endif

! tempsed
      do k=1,ksdmax
        workrsed(:,:,k)=cv_sed(-1,k,:,:)
      enddo
      indx=indx+1
      ierr=nf_fwrite(workrsed, ncidrst, rstMust(4),
     &                                       record, b3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_rstMust(1,indx))
        write(stdout,1) vname_rstMust(1,indx)(1:lvar), record, ierr
     &                  MYID
        goto 99                                         !--> ERROR
      endif

! salsed
      do k=1,ksdmax
         workrsed(:,:,k)=cv_sed(0,k,:,:)
      enddo
       indx=indx+1

      ierr=nf_fwrite(workrsed, ncidrst, rstMust(5),
     &                                       record, b3dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname_rstMust(1,indx))
        write(stdout,1) vname_rstMust(1,indx)(1:lvar), record, ierr
     &                  MYID
        goto 99                                         !--> ERROR
      endif

! cvsed
       do itrc=1,ntrc_subs
         indx=indx+1           ! 3D: traceurSed
         do k=1,ksdmax
           workrsed(:,:,k)=cv_sed(itrc,k,:,:)
         enddo

        ierr=nf_fwrite(workrsed, ncidrst, rstMust(indx),
     &                                       record, b3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname_rstMust(1,indx))
          write(stdout,1) vname_rstMust(1,indx)(1:lvar), record, ierr
     &                  MYID
          goto 99                                         !--> ERROR
        endif
       enddo

! dredging
       if (l_dredging) then
          work2d = dredg_hsed_init
          indx = indx + 1
          ierr = nf_fwrite(work2d, ncidrst, rstMust(indx),
     &                                       record, r2dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname_rstMust(1,indx))
            write(stdout,1) vname_rstMust(1,indx)(1:lvar), record, ierr
     &                  MYID
            goto 99                                         !--> ERROR
          endif
       endif


# endif /* MUSTANG */

#endif /* SOLVE3D */

#ifdef MORPHODYN
!
!     Write out time evolving bathymetry
!
      ierr=nf_fwrite(h, ncidrst, rstHm, record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxHm))
        write(stdout,1) vname(1,indxHm)(1:lvar), record, ierr
     &                    MYID
        goto 99                                         !--> ERROR
      endif
CR      write(*,*) 'wrt_rst:  Hm' MYID
#endif

#ifdef BBL
      ierr=nf_fwrite(Hripple, ncidrst, rstBBL(1), record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxHrip))
        write(stdout,1) vname(1,indxHrip)(1:lvar), record, ierr
     &                    MYID
        goto 99                                         !--> ERROR
      endif
CR      write(*,*) 'wrt_rst:  Hripple' MYID
      ierr=nf_fwrite(Lripple, ncidrst, rstBBL(2), record, r2dvar)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxLrip))
        write(stdout,1) vname(1,indxLrip)(1:lvar), record, ierr
     &                  MYID
        goto 99                                         !--> ERROR
      endif
CR      write(*,*) 'wrt_rst:  Lripple' MYID
#endif /* BBL */

  1   format(/1x, 'WRT_RST ERROR while writing variable ''', A,
     &           ''' into restart file.', /11x, 'Time record:',
     &               i6, 3x, 'netCDF error code', i4, 3x, A,i4)
      goto 100
  99  may_day_flag=3
 100  continue

!
! Synchronize restart netCDF file to disk to allow other
! processes to access data immediately after it is written.
!
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      ierr=nf_close (ncidrst)
      if (nrpfrst.gt.0 .and. record.ge.nrpfrst) ncidrst=-1
#else
      if (nrpfrst.gt.0 .and. record.ge.nrpfrst) then
        ierr=nf_close (ncidrst)
        ncidrst=-1
      else
        ierr=nf_sync(ncidrst)
      endif
#endif
      if (ierr .eq. nf_noerr) then
      MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_RST -- wrote ',
     &            'restart fields into time record =', record, '/',
     &             nrecrst  MYID
      else
      MPI_master_only  write(stdout,'(/1x,2A/)')
     &             'WRT_RST ERROR: Cannot ',
     &             'synchronize/close restart netCDF file.'
        may_day_flag=3
      endif

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      if (mynode .lt. NNODES-1) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, MPI_COMM_WORLD,  ierr)
      endif
#endif
      return
      end
