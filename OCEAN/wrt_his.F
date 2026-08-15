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
                               ! Writes requested model
      subroutine wrt_his       ! fields at requested levels
                               ! into history netCDF file.
#ifdef SUBSTANCE
      USE comsubstance, ONLY : cvfix_wat, nvpc
#endif
#ifdef MUSTANG
      USE sed_MUSTANG,  ONLY : sed_MUSTANG_outres
      USE comsubstance, ONLY : ws_part
      USE comMUSTANG,   ONLY : vname_Must, hisMust, outMust,
     &                         nk_nivsed_out, var3D_cvsed, var3D_dzs,
     &                         var3D_TEMP, var3D_SAL, var3D_poro,
     &                         var2D_hsed, ksma, tauskin, tauskin_c, tauskin_w,
     &                         var2D_toce, var2D_flx_s2w, var2D_flx_w2s,
     &                         var2D_pephm_fcor, var2D_flx_bx, var2D_flx_by, 
     &                         var2D_bil_bedload, var2D_fsusp, 
     &                         var2D_frmudsup, var2D_dzs_ksmax, 
     &                         var2D_theoric_active_layer, 
     &                         var2D_tero_noncoh, var2D_tero_coh, 
     &                         var2D_pct_iter_noncoh, var2D_pct_iter_coh, 
     &                         var2D_niter_ero, 
     &                         var2D_flx_s2w_coh, var2D_flx_s2w_noncoh, 
     &                         var2D_flx_w2s_coh, var2D_flx_w2s_noncoh, 
     &                         var2D_flx_bx_int, var2D_flx_by_int, 
     &                         var2D_bil_bedload_int, z0sed,
     &                         l_dyn_insed, var3Dksed_dtsdzs, 
     &                         var3Dksed_loadograv, var3Dksed_permeab, 
     &                         var3Dksed_sed_rate, var3Dksed_sigmadjge, 
     &                         var3Dksed_sigmapsg, var3Dksed_hinder, 
     &                         var3Dksed_stateconsol
# ifdef key_sand2D
      USE comsubstance, ONLY : l_subs2D, nv_grav, nv_sand, l_outsandrouse
      USE comMUSTANG,   ONLY : sum_tmp, rouse2D
# endif
#endif

#ifdef OBSTRUCTION
      USE com_OBSTRUCTIONS, ONLY : obst_nbvar,
     &     vname_obst, hisObst, l_obstout_pos,
     &     l_obstout_height_f, l_obstout_height_e, 
     &     l_obstout_dens_f, l_obstout_dens_e, 
     &     l_obstout_width_f, l_obstout_width_e, 
     &     l_obstout_thick_f, l_obstout_thick_e,
     &     l_obstout_theta, l_obstout_frac_xy, 
     &     l_obstout_frac_z, l_obstout_fuzvz, 
     &     l_obstout_a2d, l_obstout_a3d, l_obstout_s2d, 
     &     l_obstout_s3d, l_obstout_drag, l_obstout_tau,
     &     obst_position, 
     &     obst_height_inst, obst_height, 
     &     obst_dens_inst, obst_dens3d, 
     &     obst_width_inst, obst_width3d, 
     &     obst_thick_inst, obst_thick3d,
     &     obst_theta3d, obst_fracxy, obst_fracz3d,
     &     obst_drag3d, obst_a2d, obst_a3d,
     &     obst_s2d, obst_s3d,
     &     obst_tau, obst_fuz, obst_fvz
#endif

!     implicit none
      integer ierr,  record, lstr, lvar, lenstr, type
     &   , start(2), count(2), nf_fwrite, cff
#if defined OUTPUTS_SURFACE && ! defined XIOS
     &   , ibuff(6)
#else
     &   , ibuff(4)
#endif
      real eps
      parameter (eps=1.D-20)
      real stf_cff
      parameter(stf_cff=86400./0.01)
#ifdef SOLVE3D
      integer tile, itrc,i,j,k
#endif
#if defined MPI & !defined PARALLEL_FILES
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE), blank
#endif
#include "param.h"
#include "scalars.h"
#include "ncscrum.h"
#include "forces.h"
#include "grid.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "mixing.h"
#include "mpi_cpl.h"
#ifdef SEDIMENT
# include "sediment.h"
      integer indxWrk
#endif
#ifdef BBL
# include "bbl.h"
#endif
#ifdef WKB_WWAVE
# include "wkb_wwave.h"
#endif
#ifdef NBQ
# include "nbq.h"
#endif
#include "work.h"
#include "netcdf.inc"

#ifdef MUSTANG
      real     workrsed(GLOBAL_2D_ARRAY,nk_nivsed_out),
     &         workrsed1(GLOBAL_2D_ARRAY,nk_nivsed_out)
      real     workrsed_nv(GLOBAL_2D_ARRAY),
     &         workrsed1_nv(GLOBAL_2D_ARRAY)
#endif

#ifdef MASKING
# define SWITCH *
#else
# define SWITCH !
#endif

#if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, MPI_COMM_WORLD, status, ierr)
      endif
#endif
!
! Make sure nrechis is larger than 0 before creating file
!
      nrechis=max(1,nrechis) 
!
! Create/open history file; write grid arrays, if so needed.
!
      call def_his (ncidhis, nrechis, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(hisname)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been
!                                          !!! started, it is assumed
      nrechis=max(1,nrechis)               !!! that global history
      if (nrpfhis.eq.0) then               !!! record index is
        record=nrechis                     !!! advanced by main.
      else
        record=1+mod(nrechis-1, nrpfhis)
      endif


#undef CR
CR      write(*,*) 'wrt_his: Entry  ' MYID
!
!     Write out evolving model variables:
!     ----- --- -------- ----- ----------
!
!     Time step number and record numbers.
!
      type=filetype_his
!
      if (iic .eq. 0) then
        ibuff(1)=iic+ntstart
      else
        ibuff(1)=iic
      endif
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
      ierr=nf_put_vara_int (ncidhis, hisTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
         write(*,*) nf_strerror(ierr)
         MPI_master_only write(stdout,1) 'time_step', record, ierr
     &        MYID
         goto 99                !--> ERROR
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
      ierr=nf_put_vara_int (ncidhis, hisTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        MPI_master_only write(stdout,1)
     &        'time_step_surf', record, ierr
     &        MYID
      endif
#endif

CR      write(*,*) 'wrt_his:  time  ' MYID
!
!     Time
!
      ierr=nf_put_var1_FTYPE (ncidhis, hisTime, record, time
#ifdef USE_CALENDAR
     &  - origin_date_in_sec
#endif
     &  )
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        MPI_master_only write(stdout,1)
     &        vname(1,indxTime)(1:lvar), record, ierr
     &        MYID
        goto 99                !--> ERROR
      endif

!
!     Time2
!
      ierr=nf_put_var1_FTYPE (ncidhis, hisTime2, record, time
#ifdef USE_CALENDAR
     &  - origin_date_in_sec
#endif
     &  )
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime2))
        MPI_master_only write(stdout,1)
     &        vname(1,indxTime2)(1:lvar), record, ierr
     &        MYID
        goto 99                !--> ERROR
      endif
CR      write(*,*) 'wrt_his:  time  ' MYID
!
!     Barotropic mode variables: free-surface and 2D momentum
!     components in XI-,ETA-directions.
!
      if (wrthis(indxZ)) then
!$acc update if(compute_on_device) host(zeta(:,:,fast_indx_out))
        work2d=zeta(:,:,fast_indx_out)
        call fillvalue2d(work2d,ncidhis,hisZ,vname(1,indxZ),
     &        record,r2dvar,type)
      endif
!
      if (wrthis(indxUb)) then
!$acc update if(compute_on_device) host(ubar(:,:,fast_indx_out))
        work2d(:,:)=ubar(:,:,fast_indx_out)
        call fillvalue2d(work2d,ncidhis,hisUb,vname(1,indxUb),
     &        record,u2dvar,type)
CR      write(*,*) 'wrt_his: wrtUBAR' MYID
      endif
!
      if (wrthis(indxVb)) then
!$acc update if(compute_on_device) host(vbar(:,:,fast_indx_out))
        work2d(:,:)=vbar(:,:,fast_indx_out)
        call fillvalue2d(work2d,ncidhis,hisVb,vname(1,indxVb),
     &        record,v2dvar,type)
CR           write(*,*) 'wrt_his: wrtVBAR' MYID
      endif

#ifdef MORPHODYN
!
!     Write out time evolving bathymetry
!
      if (wrthis(indxHm)) then
!$acc update if(compute_on_device) host(h(:,:))
        work2d(:,:)=h(:,:)
        call fillvalue2d(work2d,ncidhis,hisHm,vname(1,indxHm),
     &        record,r2dvar,type)
CR      write(*,*) 'wrt_his: wrtHm' MYID
      endif
#endif

!
!     Write out kinematic bottom stress (N/m2).
!
      if (wrthis(indxBostr)) then
!$acc update if(compute_on_device) host(bustr(:,:))
!$acc update if(compute_on_device) host(bvstr(:,:))
        do j=0,Mm
          do i=0,Lm
            work2d(i,j)=0.5*sqrt((bustr(i,j)+bustr(i+1,j))**2
     &                          +(bvstr(i,j)+bvstr(i,j+1))**2)
     &                                                  *rho0
          enddo
        enddo
        call fillvalue2d(work2d,ncidhis,hisBostr,vname(1,indxBostr),
     &        record,r2dvar,type)
CR      write(*,*) 'wrt_his:  bostr' MYID
      endif

!
!     Write out kinematic U-bottom stress (N/m2).
!
      if (wrthis(indxBustr)) then
!$acc update if(compute_on_device) host(bustr(:,:))
         work2d=bustr*rho0
         call fillvalue2d(work2d,ncidhis,hisBustr,vname(1,indxBustr),
     &        record,u2dvar,type)
CR      write(*,*) 'wrt_his:  bustr' MYID
      endif



!
!     Write out kinematic V-bottom stress (N/m2).
!
      if (wrthis(indxBvstr)) then
!$acc update if(compute_on_device) host(bvstr(:,:))
         work2d=bvstr*rho0
         call fillvalue2d(work2d,ncidhis,hisBvstr,vname(1,indxBvstr),
     &        record,v2dvar,type)
CR      write(*,*) 'wrt_his:  bvstr' MYID
      endif



!
!--   Atmospheric forcing : no mask, no fill value
!
!
!     Write out kinematic surface stress (N/m2).
!
      if (wrthis(indxWstr)) then
!$acc update if(compute_on_device) host(sustr(:,:))
!$acc update if(compute_on_device) host(svstr(:,:))
        do j=1,Mm
          do i=1,Lm
#ifdef OA_COUPLING
            work2d2(i,j)=smstr(i,j)*rho0
#elif defined STRESS_AT_RHO_POINTS
            work2d2(i,j)=sqrt(sustr(i,j)**2+svstr(i,j)**2)
     &                                                       *rho0
#else
            work2d2(i,j)=0.5*sqrt((sustr(i,j)+sustr(i+1,j))**2
     &                           +(svstr(i,j)+svstr(i,j+1))**2)
     &                                                       *rho0
#endif
          enddo
        enddo
!        write(*,*)'WSTR : work2d2(:,5)=',work2d2(:,5)
        ierr=nf_fwrite(work2d2, ncidhis, hisWstr, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxWstr))
          MPI_master_only write(stdout,1)
     &        vname(1,indxWstr)(1:lvar), record, ierr
     &           MYID
          goto 99             !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  Wstr' MYID
      endif
!
!     Write out kinematic U surface stress component (N/m2).
!
      if (wrthis(indxUWstr)) then
!$acc update if(compute_on_device) host(sustr(:,:))
        work2d=sustr*rho0
#ifdef STRESS_AT_RHO_POINTS
        ierr=nf_fwrite(work2d, ncidhis, hisUWstr, record, r2dvar)
#else
        ierr=nf_fwrite(work2d, ncidhis, hisUWstr, record, u2dvar)
#endif
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxUWstr))
          MPI_master_only write(stdout,1)
     &        vname(1,indxUWstr)(1:lvar), record, ierr
     &           MYID
          goto 99             !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  UWstr' MYID
      endif
!
!     Write out kinematic V surface stress component (N/m2).
!
      if (wrthis(indxVWstr)) then
!$acc update if(compute_on_device) host(svstr(:,:))
        work2d=svstr*rho0
#ifdef STRESS_AT_RHO_POINTS
         ierr=nf_fwrite(work2d,ncidhis,hisVWstr,record,r2dvar)
#else
         ierr=nf_fwrite(work2d,ncidhis,hisVWstr,record,v2dvar)
#endif
         if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxVWstr))
         MPI_master_only write(stdout,1)
     &        vname(1,indxVWstr)(1:lvar), record, ierr
     &           MYID
          goto 99             !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  VWstr' MYID
      endif

#ifdef WAVE_IO
      if (wrthis(indxHRM)) then
!$acc update if(compute_on_device) host(whrm(:,:))
        ierr=nf_fwrite (whrm(START_2D_ARRAY), ncidhis,
     &                                   hisWAVE(1), record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxHRM))
         MPI_master_only write(stdout,1)
     &        vname(1,indxHRM)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtHrm' MYID
      endif
      if (wrthis(indxFRQ)) then
!$acc update if(compute_on_device) host(wfrq(:,:))
        ierr=nf_fwrite (wfrq(START_2D_ARRAY), ncidhis,
     &                                   hisWAVE(2), record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxFRQ))
         MPI_master_only write(stdout,1)
     &        vname(1,indxFRQ)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtFrq' MYID
      endif
# ifdef WKB_WWAVE
      if (wrthis(indxWAC)) then
!$acc update if(compute_on_device) host(wac(:,:,wstp))
        ierr=nf_fwrite (wac(START_2D_ARRAY,wstp), ncidhis,
     &                                   hisWAVE(3), record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxWAC))
         MPI_master_only write(stdout,1)
     &        vname(1,indxWAC)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtWac' MYID
      endif
# endif
      if (wrthis(indxWKX)) then
!$acc update if(compute_on_device) host(wwkx(:,:))
        ierr=nf_fwrite (wwkx(START_2D_ARRAY), ncidhis,
     &                                   hisWAVE(4), record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxWKX))
         MPI_master_only write(stdout,1)
     &        vname(1,indxWKX)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtWkx' MYID
      endif
      if (wrthis(indxWKE)) then
!$acc update if(compute_on_device) host(wwke(:,:))
        ierr=nf_fwrite (wwke(START_2D_ARRAY), ncidhis,
     &                                   hisWAVE(5), record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxWKE))
         MPI_master_only write(stdout,1)
     &        vname(1,indxWKE)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtWke' MYID
      endif
      if (wrthis(indxEPB)) then
!$acc update if(compute_on_device) host(wepb0(:,:))
        ierr=nf_fwrite (wepb0(START_2D_ARRAY), ncidhis,
     &                                   hisWAVE(6), record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxEPB))
         MPI_master_only write(stdout,1)
     &        vname(1,indxEPB)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtEpb' MYID
      endif
      if (wrthis(indxEPD)) then
!$acc update if(compute_on_device) host(wepd(:,:))
        ierr=nf_fwrite (wepd(START_2D_ARRAY), ncidhis,
     &                                   hisWAVE(7), record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxEPD))
         MPI_master_only write(stdout,1)
     &        vname(1,indxEPD)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtEpd' MYID
      endif
# ifdef WAVE_ROLLER
#  ifdef WKB_WWAVE
      if (wrthis(indxWAR)) then
!$acc update if(compute_on_device) host(war(:,:,wstp))
        ierr=nf_fwrite (war(START_2D_ARRAY,wstp), ncidhis,
     &                                   hisWAVE(8), record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxWAR))
         MPI_master_only write(stdout,1)
     &        vname(1,indxWAR)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtWar' MYID
      endif
#  endif
      if (wrthis(indxEPR)) then
!$acc update if(compute_on_device) host(wepr(:,:))
        ierr=nf_fwrite (wepr(START_2D_ARRAY), ncidhis,
     &                                   hisWAVE(9), record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxEPR))
         MPI_master_only write(stdout,1)
     &        vname(1,indxEPR)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtEpr' MYID
      endif
# endif
#endif /* WKB_WWAVE || OW_COUPLING || WAVE_OFFLINE */

#ifdef MRL_WCI
!
!  Add output variables linked to MRL_WCI
!
      if (wrthis(indxSUP)) then
!$acc update if(compute_on_device) host(sup(:,:))
        ierr=nf_fwrite (sup, ncidhis, hisSUP, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxSUP))
         MPI_master_only write(stdout,1)
     &        vname(1,indxSUP)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtSup' MYID
      endif
!Ustokes 2D
      if (wrthis(indxUST2D)) then
!$acc update if(compute_on_device) host(ust2d(:,:))
        ierr=nf_fwrite (ust2d, ncidhis, hisUST2D, record, u2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxUST2D))
         MPI_master_only write(stdout,1)
     &        vname(1,indxUST2D)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtUst2D' MYID
      endif
!Vstokes 2D
      if (wrthis(indxVST2D)) then
!$acc update if(compute_on_device) host(vst2d(:,:))
        ierr=nf_fwrite (vst2d, ncidhis, hisVST2D, record, v2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxVST2D))
         MPI_master_only write(stdout,1)
     &        vname(1,indxVST2D)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtVst2D' MYID
      endif
# ifdef SOLVE3D
!Ustokes
      if (wrthis(indxUST)) then
!$acc update if(compute_on_device) host(ust(:,:,:))
        ierr=nf_fwrite (ust, ncidhis, hisUST,  record,  u3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxUST))
         MPI_master_only write(stdout,1)
     &        vname(1,indxUST)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtUst' MYID
      endif
!Vstokes
      if (wrthis(indxVST)) then
!$acc update if(compute_on_device) host(vst(:,:,:))
        ierr=nf_fwrite (vst, ncidhis, hisVST,  record,  v3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxVST))
         MPI_master_only write(stdout,1)
     &        vname(1,indxVST)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtVst' MYID
      endif
!Wstokes
      if (wrthis(indxWST)) then
        do tile=0,NSUB_X*NSUB_E-1
          call wstokes (tile)
        enddo
!$acc update if(compute_on_device) host(wst(:,:,:))
        ierr=nf_fwrite (wst, ncidhis, hisWST, record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxWST))
         MPI_master_only write(stdout,1)
     &        vname(1,indxWST)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtWst' MYID
      endif
!Mixing coeff
      if (wrthis(indxAkb)) then
!$acc update if(compute_on_device) host(Akb(:,:,:))
        ierr=nf_fwrite (Akb, ncidhis, hisAkb, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkb))
         MPI_master_only write(stdout,1)
     &        vname(1,indxAkb)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtAkb' MYID
      endif
      if (wrthis(indxAkw)) then
!$acc update if(compute_on_device) host(Akw(:,:,:))
        ierr=nf_fwrite (Akw, ncidhis, hisAkw, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkw))
         MPI_master_only write(stdout,1)
     &        vname(1,indxAkw)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtAkw' MYID
      endif
!Vortex force
      if (wrthis(indxKVF)) then
!$acc update if(compute_on_device) host(kvf(:,:,:))
        ierr=nf_fwrite (kvf, ncidhis, hisKVF, record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxKVF))
         MPI_master_only write(stdout,1)
     &        vname(1,indxKVF)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtKvf' MYID
      endif
!Bernouilli head
      if (wrthis(indxCALP)) then
!$acc update if(compute_on_device) host(calP(:,:))
        ierr=nf_fwrite (calP, ncidhis, hisCALP, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxCALP))
         MPI_master_only write(stdout,1)
     &        vname(1,indxCALP)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtCalp' MYID
      endif
      if (wrthis(indxKAPS)) then
!$acc update if(compute_on_device) host(Kapsrf(:,:))
        ierr=nf_fwrite (Kapsrf, ncidhis, hisKAPS, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxKAPS))
         MPI_master_only write(stdout,1)
     &        vname(1,indxKAPS)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtKaps' MYID
      endif
# endif  /* SOLVE3D */
#endif   /* MRL_WCI */
!
#ifdef SOLVE3D
!
!     3D momentum components in XI- and ETA-directions.
!
      if (wrthis(indxU)) then
!$acc update if(compute_on_device) host(u(:,:,:,nstp))
        workr=u(:,:,:,nstp)
        call fillvalue3d(workr,ncidhis,hisU,vname(1,indxU),
     &        record,u3dvar,type)
CR      write(*,*) 'wrt_his:  wrtU  ' MYID
      endif
!
      if (wrthis(indxV)) then
!$acc update if(compute_on_device) host(v(:,:,:,nstp))
        workr=v(:,:,:,nstp)
        call fillvalue3d(workr,ncidhis,hisV,vname(1,indxV),
     &        record,v3dvar,type)
CR      write(*,*) 'wrt_his:  wrtV  ' MYID
      endif
!
!     Tracer variables.
!
# ifdef TRACERS
      do itrc=1,NT
# ifdef DIAGNOSTICS_TS_MLD_CRIT
         if (itrc.ne.iCRT2 .and. itrc.ne.iCRT3
     &   .and. itrc.ne.iCRT4) then
# endif
         if (wrthis(indxV+itrc)) then
!$acc update if(compute_on_device) host(t(:,:,:,nstp,itrc))
            workr=t(:,:,:,nstp,itrc)
#  ifdef MASKING
            if (iic.eq.0) then
!     We do that because the mask has been applied if iic ==0 !
               do k=1,N
                  workr(:,:,k)=workr(:,:,k)
#   ifdef WET_DRY
     &                         SWITCH rmask_wet(:,:)
#   else
     &                         SWITCH rmask(:,:)
#   endif
               enddo
            endif
#  endif
            call fillvalue3d(workr,ncidhis,hisT(itrc),
     &           vname(1,indxV+itrc),record,r3dvar,type)
CR      write(*,*) 'wrt_his:  wrtT  ' MYID

      endif
# ifdef DIAGNOSTICS_TS_MLD_CRIT
      endif
# endif
      enddo
# endif /* TRACERS */
!
!     Density anomaly.
!
      if (wrthis(indxR)) then
!$acc update if(compute_on_device) host(rho(:,:,:))
        workr=rho+rho0-1000.
        call fillvalue3d(workr,ncidhis,hisR,vname(1,indxR),
     &        record,r3dvar,type)
CR      write(*,*) 'wrt_his:  wrtRHO' MYID
      endif
# if defined ANA_VMIX || defined GLS_MIXING \
  || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP
!
!     Brunt vaisala frequency.
!
      if (wrthis(indxbvf)) then
!$acc update if(compute_on_device) host(bvf(:,:,:))
         work=bvf
         call fillvalue3d_w(work,ncidhis,hisbvf,vname(1,indxbvf),
     &           record,w3dvar,type)
CR      write(*,*) 'wrt_his:  wrtbvf' MYID
      endif
# endif

!
!     S-coordinate omega vertical velocity (m/s).
!
      if (wrthis(indxO)) then
!$acc update if(compute_on_device) host(we(:,:,:))
# ifdef VADV_ADAPT_IMP
!$acc update if(compute_on_device) host(wi(:,:,:))
# endif
        do k=0,N
          do j=0,Mm+1
            do i=0,Lm+1
              work(i,j,k)= ( We(i,j,k)
# ifdef VADV_ADAPT_IMP
     &                   +   Wi(i,j,k)
# endif
     &                              ) *pm(i,j)*pn(i,j)
!#ifdef NBQ_MASS
!     &                            /rho_nbq_avg1(i,j,k)
!#endif
            enddo
          enddo
        enddo
        call fillvalue3d_w(work,ncidhis,hisO,vname(1,indxO),
     &        record,w3dvar,type)
CR      write(*,*) 'wrt_his:  wrtO  ' MYID
      endif
!
!     Write out true vertical velocity (m/s).
!
      if (wrthis(indxW)) then
# ifdef NBQ
!$acc update if(compute_on_device) host(wz(:,:,:))
        do k=0,N
          do j=0,Mm+1
            do i=0,Lm+1
              work(i,j,k)=wz(i,j,k,nstp)
            enddo
          enddo
        enddo
        call fillvalue3d(work,ncidhis,hisW,
     &      vname(1,indxW),record,w3dvar,type)
# else
        do tile=0,NSUB_X*NSUB_E-1
          call Wvlcty (tile, workr)
        enddo
!$acc update if(compute_on_device) host(workr(:,:,:))        
        call fillvalue3d(workr,ncidhis,hisW,
     &      vname(1,indxW),record,r3dvar,type)
# endif

CR      write(*,*) 'wrt_his:  wrtW  ' MYID
      endif
!
# ifdef VIS_COEF_3D
!
!     Write out Horizontal viscosity coefficient.
!
      if (wrthis(indxVisc)) then
!$acc update if(compute_on_device) host(visc3d_r(:,:,:))
        workr=visc3d_r
        call fillvalue3d(workr,ncidhis,hisVisc,vname(1,indxVisc),
     &        record,r3dvar,type)
CR      write(*,*) 'wrt_his:  wrtVisc' MYID
      endif
# endif
!
# ifdef DIF_COEF_3D
!
!     Write out Horizontal Diffusivity coefficient.
!
      if (wrthis(indxDiff)) then
#  ifdef TS_DIF2
#   ifdef TEMPERATURE
!$acc update if(compute_on_device) host(diff2(:,:itemp))
#   endif
#   ifdef TS_DIF_SMAGO
!$acc update if(compute_on_device) host(diff3d_r(:,:,:))
#   endif
#  elif defined TS_DIF4
!$acc update if(compute_on_device) host(diff4(:,:,itemp))
#   ifdef TS_DIF_SMAGO
!$acc update if(compute_on_device) host(diff3d_r(:,:,:))
#   endif
!$acc update if(compute_on_device) host(diff3d_u(:,:,:))
!$acc update if(compute_on_device) host(diff3d_v(:,:,:))
#  endif
 
        do k=1,N
          do j=1,Mm
            do i=1,Lm
              workr(i,j,k)=
#  ifdef TS_DIF2
#   ifdef TEMPERATURE
     &                 diff2(i,j,itemp)
#   endif
#   ifdef TS_DIF_SMAGO
     &              +diff3d_r(i,j,k)
#   endif
#  elif defined TS_DIF4
     &               diff4(i,j,itemp)
#   ifdef TS_DIF_SMAGO
     &              +diff3d_r(i,j,k)*om_r(i,j)*on_r(i,j)
#   endif
     &          +0.25*(diff3d_u(i,j,k)+diff3d_u(i+1,j,k)
     &                +diff3d_v(i,j,k)+diff3d_v(i,j+1,k))
#  endif
            enddo
          enddo
        enddo
        call fillvalue3d(workr,ncidhis,hisDiff,vname(1,indxDiff),
     &        record,r3dvar,type)
      endif
CR      write(*,*) 'wrt_his:  wrtDiff' MYID
# endif
!
!     Write out vertical viscosity coefficient.
!
      if (wrthis(indxAkv)) then
!$acc update if(compute_on_device) host(Akv(:,:,:))
        work=Akv
        call fillvalue3d_w(work,ncidhis,hisAkv,vname(1,indxAkv),
     &        record,w3dvar,type)
CR      write(*,*) 'wrt_his:  wrtAkv' MYID
      endif
!
!     Write out vertical diffusion coefficient for potential temperature.
!
# ifdef TEMPERATURE
      if (wrthis(indxAkt)) then
!$acc update if(compute_on_device) host(Akt(:,:,:,itemp))
        work=Akt(:,:,:,itemp)
        call fillvalue3d_w(work,ncidhis,hisAkt,vname(1,indxAkt),
     &        record,w3dvar,type)
CR      write(*,*) 'wrt_his:  wrtAkt' MYID
      endif
# endif
# ifdef SALINITY
!
!     Write out vertical diffusion coefficient for salinity.
!
      if (wrthis(indxAks)) then
!$acc update if(compute_on_device) host(Akt(:,:,:,isalt))
        work=Akt(:,:,:,isalt)
        call fillvalue3d_w(work,ncidhis,hisAks,vname(1,indxAks),
     &        record,w3dvar,type)
CR      write(*,*) 'wrt_his:  wrtAks' MYID
      endif
# endif

# if defined LMD_SKPP || defined GLS_MIXING
!
!     Write out depth of planetary boundary layer (m).
!
      if (wrthis(indxHbl)) then
#  ifdef LMD_SKPP2005
!$acc update if(compute_on_device) host(hbls(:,:))
        work2d=hbls(:,:,nstp)
#  else
!$acc update if(compute_on_device) host(hbl(:,:))
        work2d=hbl
#  endif
        call fillvalue2d(work2d,ncidhis,hisHbl,vname(1,indxHbl),
     &        record,r2dvar,type)
CR      write(*,*) 'wrt_his:  wrtHBL' MYID
      endif
# endif
# ifdef LMD_BKPP
!
!     Write out depth of bottom planetary boundary layer (m).
!
      if (wrthis(indxHbbl)) then
!$acc update if(compute_on_device) host(hbbl(:,:))
        work2d=hbbl
        call fillvalue2d(work2d,ncidhis,hisHbbl,vname(1,indxHbbl),
     &        record,r2dvar,type)
CR      write(*,*) 'wrt_his:  wrtHBBL' MYID
      endif
# endif
# ifdef GLS_MIXING
!
!     Write out turbulent kinetic energy.
!
      if (wrthis(indxTke)) then
!$acc update if(compute_on_device) host(trb(:,:,:,nstp,itke))
        work=trb(:,:,:,nstp,itke)
        call fillvalue3d_w(work,ncidhis,hisTke,vname(1,indxTke),
     &        record,w3dvar,type)
CR      write(*,*) 'wrt_his:  wrtTke' MYID
      endif
!
!     Write out generic length scale
!
      if (wrthis(indxGls)) then
!$acc update if(compute_on_device) host(trb(:,:,:,nstp,igls))
        work=trb(:,:,:,nstp,igls)
        call fillvalue3d_w(work,ncidhis,hisGls,vname(1,indxGls),
     &        record,w3dvar,type)
CR      write(*,*) 'wrt_his:  wrtGls' MYID
      endif
!
!     Write out vertical mixing length scale
!
      if (wrthis(indxLsc)) then
!$acc update if(compute_on_device) host(Lscale(:,:,:))
        work=Lscale
        call fillvalue3d_w(work,ncidhis,hisLsc,vname(1,indxLsc),
     &        record,w3dvar,type)
CR      write(*,*) 'wrt_his:  wrtLsc' MYID
      endif
# endif

!
!     Write out total surface heat flux
!
# if defined TEMPERATURE
      if (wrthis(indxShflx)) then
!$acc update if(compute_on_device) host(stflx(:,:,itemp))
        work2d=stflx(:,:,itemp)*rho0*Cp
     &          SWITCH rmask
        ierr=nf_fwrite(work2d, ncidhis, hisShflx, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxShflx))
          MPI_master_only write(stdout,1)
     &        vname(1,indxShflx)(1:lvar), record, ierr
     &           MYID
          goto 99             !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  stflx(:,:,itemp)' MYID
      endif
# endif

# ifdef BHFLUX
!
!     Write out bottom heat flux
!
      if (wrthis(indxBhflx)) then
!$acc update if(compute_on_device) host(btflx(:,:,itemp))
        work2d=btflx(:,:,itemp)*rho0*Cp
     &       SWITCH rmask
        ierr=nf_fwrite(work2d, ncidhis, hisBhflx, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxBhflx))
          MPI_master_only write(stdout,1)
     &         vname(1,indxBhflx)(1:lvar), record, ierr
     &         MYID
            goto 99             !--> ERROR
          endif
CR          write(*,*) 'wrt_his:  btflx(:,:,itemp)' MYID
        endif
# endif

# ifdef SALINITY
!     Write out surface salt flux
!
!       --> convert [psu.m.s-1] => [cm/days] :
!       stf_cff= 86400/0.01 = 1/stf_scale
!       After forcing reading (get_stflx.F), multiply by 0.01/86400.
!       Then in get_vbc.F, mult. by salinity t(i,j,N,nstp,isalt)
!
      if (wrthis(indxSwflx)) then
!$acc update if(compute_on_device) host(stflx(:,:,isalt))
        do j=0,Mm+1
          do i=0,Lm+1
            work2d(i,j)=stf_cff*stflx(i,j,isalt)/
     &              ( max(eps,t(i,j,N,nstp,isalt)))
     &              SWITCH rmask(i,j)
          enddo
        enddo
        ierr=nf_fwrite(work2d, ncidhis, hisSwflx, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxSwflx))
           MPI_master_only write(stdout,1)
     &        vname(1,indxSwflx)(1:lvar), record, ierr
     &           MYID

            goto 99             !--> ERROR
         endif
CR         write(*,*) 'wrt_his:  stflx(:,:,isalt)' MYID
       endif
#  ifdef BWFLUX
!
!     Write out bottom salt flux
!

!       --> convert [psu.m.s-1] => [cm/days] :
!       btf_cff= 86400/0.01 = 1/btf_scale
!       After forcing reading (get_btflx.F), multiply by 0.01/86400.
!       Then in get_vbc.F, mult. by salinity t(i,j,1,nstp,isalt)
!
      if (wrthis(indxBwflx)) then
!$acc update if(compute_on_device) host(btflx(:,:,isalt))
         do j=0,Mm+1
            do i=0,Lm+1
               work2d(i,j)=btf_cff*btflx(i,j,isalt)/
     &              ( max(eps,t(i,j,1,nstp,isalt)))
     &              SWITCH rmask(i,j)
            enddo
        enddo
         ierr=nf_fwrite(work2d, ncidhis, hisBwflx, record, r2dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxBwflx))
           MPI_master_only write(stdout,1)
     &        vname(1,indxBwflx)(1:lvar), record, ierr
     &           MYID
            goto 99             !--> ERROR
         endif
      write(*,*) 'wrt_his:  btflx(:,:,isalt)' MYID
      endif
#  endif
# endif
!
!     Write out surface heat flux component (W/m2)
!
      if (wrthis(indxShflx_rsw)) then
# if defined BULK_FLUX && defined TEMPERATURE
!$acc update if(compute_on_device) host(shflx_rsw(:,:))
         work2d=shflx_rsw*rho0*Cp
         ierr=nf_fwrite(work2d, ncidhis, hisShflx_rsw, record,
     &        r2dvar)
# else
!$acc update if(compute_on_device) host(srflx(:,:))
        work2d=srflx*rho0*Cp
        ierr=nf_fwrite(work2d, ncidhis, hisShflx_rsw, record, r2dvar)
# endif
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxShflx_rsw))
          MPI_master_only write(stdout,1)
     &        vname(1,indxShflx_rsw)(1:lvar),
     &           record, ierr
     &           MYID
          goto 99             !--> ERROR
        endif
      endif
# if defined BULK_FLUX && defined TEMPERATURE
      if (wrthis(indxShflx_rlw)) then
!$acc update if(compute_on_device) host(shflx_rlw(:,:))
        work2d=shflx_rlw*rho0*Cp
        ierr=nf_fwrite(work2d, ncidhis, hisShflx_rlw, record,
     &        r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxShflx_rlw))
           MPI_master_only write(stdout,1)
     &        vname(1,indxShflx_rlw)(1:lvar),
     &           record, ierr
     &           MYID
          goto 99             !--> ERROR
        endif
      endif

      if (wrthis(indxShflx_lat)) then
!$acc update if(compute_on_device) host(shflx_lat(:,:))
        work2d=shflx_lat*rho0*Cp
        ierr=nf_fwrite(work2d, ncidhis, hisShflx_lat, record,
     &        r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxShflx_lat))
           MPI_master_only write(stdout,1)
     &        vname(1,indxShflx_lat)(1:lvar),
     &           record, ierr
     &           MYID
          goto 99             !--> ERROR
        endif
      endif

      if (wrthis(indxShflx_sen)) then
!$acc update if(compute_on_device) host(shflx_sen(:,:))
        work2d=shflx_sen*rho0*Cp
        ierr=nf_fwrite(work2d, ncidhis, hisShflx_sen, record,
     &        r2dvar)
         if (ierr .ne. nf_noerr) then
           lvar=lenstr(vname(1,indxShflx_sen))
           MPI_master_only write(stdout,1)
     &        vname(1,indxShflx_sen)(1:lvar),
     &           record, ierr
     &           MYID
            goto 99             !--> ERROR
         endif
      endif
# endif
# if defined SST_SKIN && defined TEMPERATURE
!
!     Write out skin temperature (degC)
!
      if (wrthis(indxT)) then
!$acc update if(compute_on_device) host(sst_skin(:,:))
        work2d=sst_skin
        call fillvalue2d(work2d,ncidhis,hisSST_skin,
     &        vname(1,indxSST_skin), record,r2dvar,type)
      endif
# endif
!
# if defined BIOLOGY && !defined PISCES
!
!     Write out depth of the euphotic layer (m).
!
      if (wrthis(indxHel)) then
!$acc update if(compute_on_device) host(hel(:,:))
        work2d=hel
        call fillvalue2d(work2d,ncidhis,hisHel,
     &        vname(1,indxHel), record,r2dvar,type)
      endif
!
!     Write out depth of the euphotic layer (m).
!
#  if (defined BIO_NChlPZD ||  defined BIO_N2ChlPZD2)
      if (wrthis(indxChC)) then
!$acc update if(compute_on_device) host(theta(:,:))
        workr=theta
        call fillvalue3d(workr,ncidhis,hisChC,
     &        vname(1,indxChC),record,r3dvar,type)
      endif
#   ifdef OXYGEN
      if (wrthis(indxU10)) then
!$acc update if(compute_on_device) host(u10(:,:))
        work2d=u10
        call fillvalue2d(work2d,ncidhis,hisU10,
     &        vname(1,indxU10),record,r2dvar,type)
      endif
!
      if (wrthis(indxKvO2)) then
!$acc update if(compute_on_device) host(Kv_O2(:,:))
        work2d=Kv_O2
        call fillvalue2d(work2d,ncidhis,hisKvO2,
     &       vname(1,indxKvO2), record,r2dvar,type)
      endif
!
      if (wrthis(indxO2sat)) then
!$acc update if(compute_on_device) host(O2satu(:,:))
        work2d=O2satu
        call fillvalue2d(work2d,ncidhis,hisO2sat,
     &        vname(1,indxO2sat),record,r2dvar,type)
      endif
#   endif /* OXYGEN */
#  elif defined BIO_BioEBUS
      if (wrthis(indxAOU)) then
!$acc update if(compute_on_device) host(AOU(:,:))
        workr=AOU
        call fillvalue3d(workr,ncidhis,hisAOU,
     &       vname(1,indxAOU), record,r3dvar,type)
      endif
      if (wrthis(indxWIND10)) then
!$acc update if(compute_on_device) host(wind10(:,:))
        work2d=wind10
        call fillvalue2d(work2d,ncidhis,hiswind10,
     &        vname(1,indxWIND10),record,r2dvar,type)
      endif
#  endif
# endif /* BIOLOGY */
!
# ifdef SEDIMENT
!
!     Write out sediment bed layer thickness, porosity, volume
!     fraction of size class in sediment bed (2+2*NST b3dgrd variables)
!
      if (wrthis(indxATHK)) then
!$acc update if(compute_on_device) host(bot_thick(:,:))
        work2d=bot_thick
        call fillvalue2d(work2d,ncidhis,hisSed(1),
     &        vname(1,indxATHK),record,r2dvar,type)
      endif

      if (wrthis(indxBTHK)) then
!$acc update if(compute_on_device) host(bed_thick(:,:,:))
        worksed_bed=bed_thick
        call fillvalue3d(worksed_bed,ncidhis,hisSed(2),
     &        vname(1,indxBTHK),record,b3dvar,type)
      endif

      if (wrthis(indxBPOR)) then
!$acc update if(compute_on_device) host(bed_poros(:,:,:))
        worksed_bed=bed_poros
        call fillvalue3d(worksed_bed,ncidhis,hisSed(3),
     &        vname(1,indxBPOR),record,b3dvar,type)
      endif

      do itrc=1,NST
        indxWrk=indxBFRA(1)+itrc-1
        if (wrthis(indxWrk)) then
!$acc update if(compute_on_device) host(bed_frac(:,:,:,itrc))
          worksed_frac=bed_frac(:,:,:,itrc)
          call fillvalue3d(worksed_frac,ncidhis,hisSed(itrc+3),
     &           vname(1,indxWrk),record,b3dvar,type)
         endif
      enddo

#  ifdef SUSPLOAD
      do itrc=1,NST
        indxWrk=indxDFLX(1)+itrc-1
        if (wrthis(indxWrk)) then
!$acc update if(compute_on_device) host(settling_flux(:,:,itrc))
          work2d=settling_flux(:,:,itrc)/dt
          call fillvalue2d(work2d,ncidhis,hisSed(itrc+3+NST),
     &           vname(1,indxWrk),record,r2dvar,type)
        endif
      enddo
      do itrc=1,NST
        indxWrk=indxEFLX(1)+itrc-1
        if (wrthis(indxWrk)) then
!$acc update if(compute_on_device) host(ero_flux(:,:,itrc))
          work2d=ero_flux(:,:,itrc)/dt
          call fillvalue2d(work2d,ncidhis,hisSed(itrc+3+2*NST),
     &           vname(1,indxWrk),record,r2dvar,type)
        endif
      enddo
#  endif
#  ifdef BEDLOAD
      do itrc=1,NST
        indxWrk=indxBDLU(1)+itrc-1
        if (wrthis(indxWrk)) then
!$acc update if(compute_on_device) host(bedldu(:,:,itrc))
          work2d=bedldu(:,:,itrc)
#   ifdef SUSPLOAD
          call fillvalue2d(work2d,ncidhis,hisSed(itrc+3+3*NST),
     &           vname(1,indxWrk),record,r2dvar,type)
#   else
          call fillvalue2d(work2d,ncidhis,hisSed(itrc+3+NST),
     &           vname(1,indxWrk),record,r2dvar,type)
#   endif
        endif
      enddo

      do itrc=1,NST
        indxWrk=indxBDLV(1)+itrc-1
        if (wrthis(indxWrk)) then
!$acc update if(compute_on_device) host(bedldv(:,:,itrc))
          work2d=bedldv(:,:,itrc)
#   ifdef SUSPLOAD
          call fillvalue2d(work2d,ncidhis,hisSed(itrc+3+4*NST),
     &           vname(1,indxWrk),record,r2dvar,type)
#   else
          call fillvalue2d(work2d,ncidhis,hisSed(itrc+3+2*NST),
     &           vname(1,indxWrk),record,r2dvar,type)
#   endif
        endif
      enddo
#  endif

#  if defined MIXED_BED || defined COHESIVE_BED
      if (wrthis(indxBTCR)) then
!$acc update if(compute_on_device) host(bed_bctr(:,:,:))
        worksed_bed=bed_bctr
        call fillvalue3d(worksed_bed,ncidhis,hisSed(ubound(hisSed)),
     &                   vname(1,indxBTCR), record,b3dvar,type)
      endif
#  endif

# endif /* SEDIMENT */
!
# ifdef BBL
      if (wrthis(indxAbed)) then
!$acc update if(compute_on_device) host(Abed(:,:))
        work2d=Abed
        call fillvalue2d(work2d,ncidhis,hisBBL(1),vname(1,indxAbed),
     &        record,r2dvar,type)
CR        write(*,*) 'wrt_his:  Abed' MYID
      endif

      if (wrthis(indxHrip)) then
!$acc update if(compute_on_device) host(Hripple(:,:))
        work2d=Hripple
        call fillvalue2d(work2d,ncidhis,hisBBL(2),vname(1,indxHrip),
     &        record,r2dvar,type)
CR        write(*,*) 'wrt_his:  Hripple' MYID
      endif

      if (wrthis(indxLrip)) then
!$acc update if(compute_on_device) host(Lripple(:,:))
        work2d=Lripple
        call fillvalue2d(work2d,ncidhis,hisBBL(3),vname(1,indxLrip),
     &        record,r2dvar,type)
CR        write(*,*) 'wrt_his:  Lripple' MYID
      endif

      if (wrthis(indxZbnot)) then
!$acc update if(compute_on_device) host(Zbnot(:,:))
        work2d=Zbnot
        call fillvalue2d(work2d,ncidhis,hisBBL(4),vname(1,indxZbnot),
     &        record,r2dvar,type)
CR        write(*,*) 'wrt_his:  Zbnot' MYID
      endif

      if (wrthis(indxZbapp)) then
!$acc update if(compute_on_device) host(Zbapp(:,:))
         work2d=Zbapp
         call fillvalue2d(work2d,ncidhis,hisBBL(5),vname(1,indxZbapp),
     &        record,r2dvar,type)
CR        write(*,*) 'wrt_his:  Zbapp' MYID
      endif

      if (wrthis(indxBostrw)) then
!$acc update if(compute_on_device) host(bustrw(:,:))
!$acc update if(compute_on_device) host(bvstrw(:,:))
        do j=0,Mm
          do i=0,Lm
            work2d(i,j)=0.5*sqrt((bustrw(i,j)+bustrw(i+1,j))**2
     &                          +(bvstrw(i,j)+bvstrw(i,j+1))**2)
     &                                                     *rho0
          enddo
        enddo
        call fillvalue2d(work2d,ncidhis,hisBBL(6),vname(1,indxBostrw),
     &        record,r2dvar,type)
CR        write(*,*) 'wrt_his:  bostrw' MYID
      endif
!
# endif /* BBL */

!
!     Tracer variables.
!
      do itrc=1,NT
# ifdef DIAGNOSTICS_TS_MLD_CRIT
         if (itrc.ne.iCRT2 .and. itrc.ne.iCRT3
     &   .and. itrc.ne.iCRT4) then
# endif
        if (wrthis(indxT+itrc-1)) then
          workr=t(:,:,:,nstp,itrc)
# if defined MUSTANG & defined key_sand2D
          if (itrc .ge. itsubs1+nv_grav .and.
     &        itrc .le. itsubs1+nv_grav+nv_sand-1) then ! isand1,isand2
            if (l_subs2D(itrc-(itsubs1+nv_grav)+1)
     &         .and. l_outsandrouse(itrc-(itsubs1+nv_grav)+1)) then
              do j=0,Mm+1
                do i=0,Lm+1
                  do k=1,N
                    if ( sum_tmp(itrc-itsubs1-nv_grav+1,i,j) .gt. 0.) then
                      workr(i,j,k)=t(i,j,1,nstp,itrc)*
     &                             (z_w(i,j,1)-z_w(i,j,0))*
     &                           ( (z_w(i,j,N)-z_r(i,j,k))
     &                            /(z_r(i,j,k)-z_w(i,j,0))
     &                           )**rouse2D(itrc-itsubs1-nv_grav+1,i,j)
                      workr(i,j,k)=workr(i,j,k)/
     &                             sum_tmp(itrc-itsubs1-nv_grav+1,i,j)
                    else
                      workr(i,j,k)=0.
                    endif
                  enddo
                enddo
              enddo
            endif
          endif
# endif
          call fillvalue3d(workr,ncidhis,hisT(itrc),
     &                     vname(1,indxT+itrc-1),record,r3dvar,type)
CR      write(*,*) 'wrt_his:  wrtT  ' MYID
        endif
# ifdef DIAGNOSTICS_TS_MLD_CRIT
         endif
# endif
      enddo

!
!    MUSTANG output
!
# ifdef SUBSTANCE
      do itrc=1,ntfix
        if (wrthis(indxT+NT+itrc-1)) then
          workr=cvfix_wat(:,:,:,itrc)
          call fillvalue3d(workr,ncidhis,hisT(NT+itrc),
     &                     vname(1,indxT+NT+itrc-1),record,r3dvar,type)
CR      write(*,*) 'wrt_his: cvfix_wat   ' MYID
        endif
      enddo
# endif
# ifdef MUSTANG
      CALL sed_MUSTANG_outres(0, Lm+1, 0, Mm+1, ntrc_subs, h)
      
      indx=1
      if (outMust(indx)) then
        work2d=ksma(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his ksma'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        work2d=var2D_hsed(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_hsed'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        work2d=tauskin(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his tauskin'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        work2d=tauskin_c(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his tauskin_c'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        work2d=tauskin_w(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his tauskin_w'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        do k=1,nk_nivsed_out
          workrsed(:,:,k)=var3D_dzs(k,:,:)
        enddo
        call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3D_dzs'   MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        do k=1,nk_nivsed_out
          workrsed(:,:,k)=var3D_poro(k,:,:)
        enddo
        call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3D_poro'   MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        do k=1,nk_nivsed_out
          workrsed(:,:,k)=var3D_TEMP(k,:,:)
        enddo
        call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3D_TEMP'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        do k=1,nk_nivsed_out
          workrsed(:,:,k)=var3D_SAL(k,:,:)
        enddo
        call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3D_SAL'  MYID
      endif
      do itrc=1,ntrc_subs
        indx = indx + 1
        if (outMust(indx)) then
          do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3D_cvsed(k,:,:,itrc)
          enddo
          workrsed1=workrsed(:,:,:)
          call fillvalue3d_must(workrsed1,ncidhis,hisMust(indx),
     &                     vname_Must(1,indx),record,b3dvar,type)
CR        write(*,*) 'wrt_his var3D_cvsed'  MYID
        endif
      enddo
      do itrc=1,ntrc_subs
            indx = indx + 1
            if (outMust(indx)) then
              workr=ws_part(:,:,:,itemp+ntrc_salt+itrc)
              call fillvalue3d(workr,ncidhis,hisMust(indx),
     &                     vname_Must(1,indx),record,r3dvar,type)
CR        write(*,*) 'wrt_his ws_part'  MYID
            endif
          enddo

      do itrc=1,nvpc
            indx=indx+1
            if (outMust(indx)) then
            work2d=var2D_toce(itrc,:,:)
            call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_toce'  MYID
            endif
            indx=indx+1
            if (outMust(indx)) then
            work2d=var2D_flx_s2w(itrc,:,:)
            call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_flx_s2w'  MYID
            endif      
            indx=indx+1
            if (outMust(indx)) then
            work2d=var2D_flx_w2s(itrc,:,:)
            call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_flx_w2s'  MYID
            endif
#ifdef  key_MUSTANG_V2
            indx=indx+1
            if (outMust(indx)) then
            work2d=var2D_pephm_fcor(itrc,:,:)
            call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_pephm_fcor'  MYID
            endif
#ifdef key_MUSTANG_bedload
            indx=indx+1
            if (outMust(indx)) then
            work2d=var2D_flx_bx(itrc,:,:)
            call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_flx_bx'  MYID
            endif
            indx=indx+1
            if (outMust(indx)) then
            work2d=var2D_flx_by(itrc,:,:)
            call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_flx_by'  MYID
            endif
            indx=indx+1
            if (outMust(indx)) then
            work2d=var2D_bil_bedload(itrc,:,:)
            call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_bil_bedload'  MYID
            endif
            indx=indx+1
            if (outMust(indx)) then
            work2d=var2D_fsusp(itrc,:,:)
            call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_fsusp'  MYID
            endif
# endif /* key_MUSTANG_bedload */
# endif /* key_MUSTANG_V2 */
      enddo

      indx=indx+1
      if (outMust(indx)) then
        work2d=var2D_frmudsup(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_frmudsup'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        work2d=var2D_dzs_ksmax(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_dzs_ksmax'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        work2d=z0sed(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his z0sed'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        work2d=zob(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his zob'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        work2d=var2D_flx_s2w_noncoh(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_flx_s2w_noncoh'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        work2d=var2D_flx_w2s_noncoh(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_flx_w2s_noncoh'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        work2d=var2D_flx_s2w_coh(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_flx_s2w_ncoh'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
        work2d=var2D_flx_w2s_coh(:,:)
        call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_flx_w2s_coh'  MYID
      endif
#ifdef  key_MUSTANG_V2
      indx=indx+1
      if (outMust(indx)) then
      work2d=var2D_theoric_active_layer(:,:)
      call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_theoric_active_layer'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
      work2d=var2D_tero_noncoh(:,:)
      call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_tero_noncoh'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
      work2d=var2D_tero_coh(:,:)
      call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_tero_coh'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
      work2d=var2D_pct_iter_noncoh(:,:)
      call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_pct_iter_noncoh'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
      work2d=var2D_pct_iter_coh(:,:)
      call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_pct_iter_coh'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
      work2d=var2D_niter_ero(:,:)
      call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_niter_ero'  MYID
      endif
#ifdef key_MUSTANG_bedload
      indx=indx+1
      if (outMust(indx)) then
      work2d=var2D_flx_bx_int(:,:)
      call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_flx_bx_int'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
      work2d=var2D_flx_by_int(:,:)
      call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_flx_by_int'  MYID
      endif
      indx=indx+1
      if (outMust(indx)) then
      work2d=var2D_bil_bedload_int(:,:)
      call fillvalue2d(work2d,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,r2dvar,type)
CR      write(*,*) 'wrt_his var2D_bil_bedload_int'  MYID
      endif

# endif /* key_MUSTANG_bedload */
# endif /* key_MUSTANG_V2 */

      IF (l_dyn_insed) THEN
            indx=indx+1
            if (outMust(indx)) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3Dksed_loadograv(k,:,:)
            enddo
            call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3Dksed_loadograv'   MYID
            endif
            indx=indx+1
            if (outMust(indx)) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3Dksed_sigmadjge(k,:,:)
            enddo
            call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3Dksed_sigmadjge'   MYID
            endif
            indx=indx+1
            if (outMust(indx)) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3Dksed_sigmapsg(k,:,:)
            enddo
            call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3Dksed_sigmapsg'   MYID
            endif
            indx=indx+1
            if (outMust(indx)) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3Dksed_stateconsol(k,:,:)
            enddo
            call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3Dksed_stateconsol'   MYID
            endif
            indx=indx+1
            if (outMust(indx)) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3Dksed_permeab(k,:,:)
            enddo
            call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3Dksed_permeab'   MYID
            endif
            indx=indx+1
            if (outMust(indx)) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3Dksed_hinder(k,:,:)
            enddo
            call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3Dksed_hinder'   MYID
            endif
            indx=indx+1
            if (outMust(indx)) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3Dksed_sed_rate(k,:,:)
            enddo
            call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3Dksed_sed_rate'   MYID
            endif
            indx=indx+1
            if (outMust(indx)) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3Dksed_dtsdzs(k,:,:)
            enddo
            call fillvalue3d_must(workrsed,ncidhis,hisMust(indx),
     &                   vname_Must(1,indx),record,b3dvar,type)
CR      write(*,*) 'wrt_his var3Dksed_dtsdzs'   MYID
            endif
      ENDIF ! l_dyn_insed

# endif /* MUSTANG */

# ifdef OBSTRUCTION
        do itrc=1,obst_nbvar
            if (l_obstout_pos) then
                indx=itrc
                work2d(:,:)=obst_position(itrc,:,:)
                call fillvalue2d(work2d,ncidhis,hisObst(itrc),
     &                     vname_obst(1,indx),record,r2dvar,type)
            endif

            if (l_obstout_height_f) then
                indx=obst_nbvar+itrc
                work2d(:,:)=obst_height_inst(itrc,:,:)
                call fillvalue2d(work2d,ncidhis,
     &                     hisObst(obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r2dvar,type)
            endif

            if (l_obstout_height_e) then
                indx=2*obst_nbvar+itrc
                work2d(:,:)=obst_height(itrc,:,:)
                call fillvalue2d(work2d,ncidhis,
     &                     hisObst(2*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r2dvar,type)
            endif

            if (l_obstout_dens_f) then
                indx=3*obst_nbvar+itrc
                work2d(:,:)=obst_dens_inst(itrc,:,:)
                call fillvalue2d(work2d,ncidhis,
     &                     hisObst(3*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r2dvar,type)
            endif

            if (l_obstout_dens_e) then
                indx=4*obst_nbvar+itrc
                do k=1,N
                    workr(:,:,k)=obst_dens3d(itrc,k,:,:)
                enddo
                call fillvalue3d(workr,ncidhis,
     &                     hisObst(4*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif

            if (l_obstout_width_f) then
                indx=5*obst_nbvar+itrc
                work2d(:,:)=obst_width_inst(itrc,:,:)
                call fillvalue2d(work2d,ncidhis,
     &                     hisObst(5*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r2dvar,type)
            endif

            if (l_obstout_width_e) then
                indx=6*obst_nbvar+itrc
                do k=1,N
                    workr(:,:,k)=obst_width3d(itrc,k,:,:)
                enddo
                call fillvalue3d(workr,ncidhis,
     &                     hisObst(6*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif

            if (l_obstout_thick_f) then
               indx=7*obst_nbvar+itrc
               work2d(:,:)=obst_thick_inst(itrc,:,:)
               call fillvalue2d(work2d,ncidhis,
     &                     hisObst(7*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r2dvar,type)
            endif

            if (l_obstout_thick_e) then
               indx=8*obst_nbvar+itrc
               do k=1,N
                   workr(:,:,k)=obst_thick3d(itrc,k,:,:)
               enddo
               call fillvalue3d(workr,ncidhis,
     &                     hisObst(8*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif

            if (l_obstout_theta) then
               indx=9*obst_nbvar+itrc
               do k=1,N
                   workr(:,:,k)=obst_theta3d(itrc,k,:,:)*180.0/3.1415
               enddo
               call fillvalue3d(workr,ncidhis,
     &                     hisObst(9*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif

            if (l_obstout_frac_xy) then
               indx=10*obst_nbvar+itrc
               work2d(:,:)=obst_fracxy(itrc,:,:)
               call fillvalue2d(work2d,ncidhis,
     &                     hisObst(10*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r2dvar,type)
            endif

            if (l_obstout_frac_z) then
               indx=11*obst_nbvar+itrc
               do k=1,N
                   workr(:,:,k)=obst_fracz3d(itrc,k,:,:)
               enddo
               call fillvalue3d(workr,ncidhis,
     &                     hisObst(11*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif

            if (l_obstout_drag) then
               indx=12*obst_nbvar+itrc
               do k=1,N
                   workr(:,:,k)=obst_drag3d(itrc,k,:,:)
               enddo
               call fillvalue3d(workr,ncidhis,
     &                     hisObst(12*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif

            if (l_obstout_a2d) then
               indx=13*obst_nbvar+itrc
               work2d(:,:)=obst_a2d(itrc,:,:)
               call fillvalue2d(work2d,ncidhis,
     &                     hisObst(13*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r2dvar,type)
            endif

            if (l_obstout_a3d) then
               indx=14*obst_nbvar+itrc
               do k=1,N
                   workr(:,:,k)=obst_a3d(itrc,k,:,:)
               enddo
               call fillvalue3d(workr,ncidhis,
     &                     hisObst(14*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif

            if (l_obstout_s2d) then
               indx=15*obst_nbvar+itrc
               work2d(:,:)=obst_s2d(itrc,:,:)
               call fillvalue2d(work2d,ncidhis,
     &                     hisObst(15*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r2dvar,type)
            endif

            if (l_obstout_s3d) then
               indx=16*obst_nbvar+itrc
               do k=1,N
                   workr(:,:,k)=obst_s3d(itrc,k,:,:)
               enddo
               call fillvalue3d(workr,ncidhis,
     &                     hisObst(16*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif

        enddo

        do itrc=1,3

            if (l_obstout_a2d) then
               indx=17*obst_nbvar+itrc
               work2d(:,:)=obst_a2d(obst_nbvar+itrc,:,:)
               call fillvalue2d(work2d,ncidhis,
     &                     hisObst(17*obst_nbvar+itrc),
     &                     vname_obst(1,indx),record,r2dvar,type)
            endif

            if (l_obstout_a3d) then
               indx=17*obst_nbvar+3+itrc
               do k=1,N
                   workr(:,:,k)=obst_a3d(obst_nbvar+itrc,k,:,:)
               enddo
               call fillvalue3d(workr,ncidhis,
     &                     hisObst(17*obst_nbvar+3+itrc),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif

            if (l_obstout_s2d) then
               indx=17*obst_nbvar+6+itrc
               work2d(:,:)=obst_s2d(obst_nbvar+itrc,:,:)
               call fillvalue2d(work2d,ncidhis,
     &                     hisObst(17*obst_nbvar+6+itrc),
     &                     vname_obst(1,indx),record,r2dvar,type)
            endif

            if (l_obstout_s3d) then
               indx=17*obst_nbvar+9+itrc
               do k=1,N
                   workr(:,:,k)=obst_s3d(obst_nbvar+itrc,k,:,:)
               enddo
               call fillvalue3d(workr,ncidhis,
     &                     hisObst(17*obst_nbvar+9+itrc),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif
         enddo

            if (l_obstout_tau) then
               indx=17*obst_nbvar+13
               workr(:,:,:)=obst_tau(:,:,:)
               call fillvalue3d(workr,ncidhis,
     &                     hisObst(17*obst_nbvar+13),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif

            if (l_obstout_fuzvz) then
               indx=17*obst_nbvar+14
               workr(:,:,:)=obst_fuz(:,:,:)
               call fillvalue3d(workr,ncidhis,
     &                     hisObst(17*obst_nbvar+14),
     &                     vname_obst(1,indx),record,r3dvar,type)

               indx=17*obst_nbvar+15
               workr(:,:,:)=obst_fvz(:,:,:)
               call fillvalue3d(workr,ncidhis,
     &                     hisObst(17*obst_nbvar+15),
     &                     vname_obst(1,indx),record,r3dvar,type)
            endif


#endif /* OBSTRUCTION */

#endif /* SOLVE3D */

 1    format(/1x,'WRT_HIS ERROR while writing variable ''', A,
     &     ''' into history file.', /11x, 'Time record:',
     &     I6,3x,'netCDF error code',i4,3x,a,i4)
      goto 100
 99   may_day_flag=3
 100  continue

!
!     Synchronize netCDF file to disk to allow other processes
!     to access data immediately after it is written.
!
#if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      ierr=nf_close (ncidhis)
      if (nrpfhis.gt.0 .and. record.ge.nrpfhis) ncidhis=-1
#else
      if (nrpfhis.gt.0 .and. record.ge.nrpfhis) then
        ierr=nf_close (ncidhis)
        ncidhis=-1
      else
        ierr=nf_sync(ncidhis)
      endif
#endif
      if (ierr .eq. nf_noerr) then
       MPI_master_only write(stdout,'(6x,A,2(A,I5,1x),A,I5)')
     &       'WRT_HIS -- wrote ',
     &       'history fields into time record =', record, '/',
     &        nrechis  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)')
     &        'WRT_HIS ERROR: Cannot ',
     &        'synchronize/close history netCDF file.'
        may_day_flag=3
      endif

#if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode .lt. NNODES-1) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &        1, MPI_COMM_WORLD,  ierr)
      endif
#endif
      return
      end
