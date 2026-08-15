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
#ifdef ABL1D
# ifndef AVRH
      subroutine def_abl_his (ncid, total_rec, ierr)
# else
      subroutine def_abl_avg (ncid, total_rec, ierr)
# endif
!
! Create/open averages/history netCDF file. In the case when a new
! netCDF file is created, define all variables, their dimensions and
! attributes. In the case when a previously existing netCDF file is
! to be opened for addition of new data, verify that all dimensions
! of the file are consistent with the present model configuration
! and all necessary variables exist. Save netCDF IDs for all needed
! variables. Also determine size of the unlimited dimension.
!
! The difference between def_abl_his and def_abl_avg is as follows: they
! have different netCDF file name (ablhisname/ablavgname); netCDF file ID
! (passed as argument); time record index (ablhisindx/ablavgindx); array
! of switches which variables to write (wrtabl/wrtabl_avg); and different
! sets of netCDF variable IDs (ablhisTime...anlhisHbl/ablavgTime...ablavgHbl);
! and the first attribute of each variable, long_name, has prefix
! 'averaged'. Because most of the code is identical for both
! routines, the second one is generated from the first entirely
! by CPP.
!
      implicit none
      logical create_new_file
      integer ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),  u2dgrd(3), v2dgrd(3), auxil(2), checkdims
     &      , abl3dgrd(4), abl2dgrd(3)
# ifdef NC4PAR
     &      , csize,cmode
# endif

# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "netcdf.inc"
# ifdef NC4PAR
#  include "mpi_cpl.h"
      include 'mpif.h'
# endif
      character*70 text

# ifndef AVRH
#  define ncname ablname
#  define rec_per_file nrpfablhis
#  define wrt wrtabl
#  define vidTime ablhisTime
#  define vidTime2 ablhisTime2
#  define vidTstep ablhisTstep
#  define vidabl_pu_dta ablhis_pu_dta
#  define vidabl_pv_dta ablhis_pv_dta
#  define vidabl_pt_dta ablhis_pt_dta
#  define vidabl_pq_dta ablhis_pq_dta
#  define vidabl_pgu_dta ablhis_pgu_dta
#  define vidabl_pgv_dta ablhis_pgv_dta
#  define vidabl_u_abl ablhis_u_abl
#  define vidabl_v_abl ablhis_v_abl
#  define vidabl_t_abl ablhis_t_abl
#  define vidabl_q_abl ablhis_q_abl
#  define vidabl_tke_abl ablhis_tke_abl
#  define vidabl_mxlm_abl ablhis_mxlm_abl
#  define vidabl_mxld_abl ablhis_mxld_abl
#  define vidabl_avm_abl ablhis_avm_abl
#  define vidabl_avt_abl ablhis_avt_abl
#  define vidabl_ablh_abl ablhis_ablh_abl
#  define vidabl_zr_abl ablhis_zr_abl
#  define vidabl_zw_abl ablhis_zw_abl
#  define vidabl_Hzr_abl ablhis_Hzr_abl
#  define vidabl_Hzw_abl ablhis_Hzw_abl
# else
#  define ncname ablname_avg
#  define rec_per_file nrpfablavg
#  define wrt wrtabl_avg
#  define vidTime ablavgTime
#  define vidTime2 ablavgTime2
#  define vidTstep ablavgTstep
#  define vidabl_pu_dta ablavg_pu_dta
#  define vidabl_pv_dta ablavg_pv_dta
#  define vidabl_pt_dta ablavg_pt_dta
#  define vidabl_pq_dta ablavg_pq_dta
#  define vidabl_pgu_dta ablavg_pgu_dta
#  define vidabl_pgv_dta ablavg_pgv_dta
#  define vidabl_u_abl ablavg_u_abl
#  define vidabl_v_abl ablavg_v_abl
#  define vidabl_t_abl ablavg_t_abl
#  define vidabl_q_abl ablavg_q_abl
#  define vidabl_tke_abl ablavg_tke_abl
#  define vidabl_mxlm_abl ablavg_mxlm_abl
#  define vidabl_mxld_abl ablavg_mxld_abl
#  define vidabl_avm_abl ablavg_avm_abl
#  define vidabl_avt_abl ablavg_avt_abl
#  define vidabl_ablh_abl ablavg_ablh_abl
#  define vidabl_zr_abl ablavg_zr_abl
#  define vidabl_zw_abl ablavg_zw_abl
#  define vidabl_Hzr_abl ablavg_Hzr_abl
#  define vidabl_Hzw_abl ablavg_Hzw_abl
# endif
!
! Put time record index into file name. In  the case when model
! output is to be arranged into sequence of named files, the naming
! convention is as follows: 'his_root.INDEX.[MPI_node.]nc', where
! INDEX is an integer number such that (i) it is divisible by the
! specified number of records per file; and (ii)
!
!      INDEX + record_within_the_file = total_record
!
! where, 1 =< record_within_the_file =< records_per_file, so that
! total_record changes continuously throughout the sequence of files.
!
      ierr=0
      lstr=lenstr(ncname)
      if (rec_per_file.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, rec_per_file))
        call insert_time_index (ncname, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
!
! Decide whether to create a new file, or open existing one.
! Overall the whole code below is organized into 3-way switch,
!
! 10  if (create_new_file) then
!        .... create new file, save netCDF ids for all variables;
!     elseif (ncid.eq.-1) then
!        .... try to open existing file and check its dimensions
!       if (cannot be opened or rejected) then
!         create_new_file=.true.
!         goto 10
!       endif   and prepare
!        .... prepare the file for adding new data,
!        .... find and save netCDF ids for all variables
!     else
!        .... just open, no checking, all ids are assumed to be
!        .... already known (MPI single file output only).
!     endif
!
! which is designed to implement flexible opening policy:
! if ldefhis=.true., it forces creation of a new file [if the
! file already exists, it will be overwritten]; on the other hand,
! ldefhis=.false., it is assumed that the file already exists and
! an attempt to open it is made; if the attempt is successful, the
! file is prepared for appending hew data; if it fails, a new file
! is created.
!

# ifdef DO_NOT_OVERWRITE
      if (lvar.eq.total_rec-1) then
        create_new_file=ldefablhis
      else
        create_new_file=.false.
      endif
# else
      create_new_file=ldefablhis
# endif

      if (ncid.ne.-1) create_new_file=.false.
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode.gt.0) create_new_file=.false.
# endif
!
!=================================================================
!
! Create new history/averages file:    Put global attributes
! ====== === ======= ======== =====    and define all variables.
!
!=================================================================
!
  10  if (create_new_file) then
# ifndef NC4PAR
        ierr=nf_create(ncname(1:lstr),NF_CLOBBER, ncid)
# else
        cmode = ior(nf_netcdf4,nf_classic_model)
        cmode = ior(cmode, nf_mpiio)
        csize=xi_rho*eta_rho/NNODES
        MPI_master_only write(stdout,*)
     &       'CREATE HIS/AVG NC4 PARALLEL FILE'
        ierr=nf_create_par(ncname(1:lstr),cmode,
     &        MPI_COMM_WORLD,MPI_INFO_NULL,ncid)
# endif
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,'(/3(1x,A)/)')
     &           'ERROR in def_his/avg:',
     &           'Cannot create netCDF file:', ncname(1:lstr)
          goto 99                                         !--> ERROR
        endif
        if (rec_per_file.eq.0) total_rec=0
!
! Put global attributes.
! --- ------ -----------
!
        call put_global_atts (ncid, ierr)
!
! Define dimensions of staggered fields.
! ------ ---------- -- --------- -------
!
        ierr=nf_def_dim (ncid, 'xi_rho' , xi_rho      , r2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho', eta_rho     , r2dgrd(2))
        ierr=nf_def_dim (ncid, 'N_abl'  , N_abl       , abl3dgrd(3))
        ierr=nf_def_dim (ncid, 'time'   , nf_unlimited, timedim)

        ierr=nf_def_dim (ncid, 'auxil'  , 4           , auxil(1))

        auxil(2)   = timedim

        abl3dgrd(1)=r2dgrd(1)         !
        abl3dgrd(2)=r2dgrd(2)         ! 3D ABL-type
        abl3dgrd(4)=timedim           !

        abl2dgrd(1)=r2dgrd(1)         !
        abl2dgrd(2)=r2dgrd(2)         ! 2D ABL-type
        abl2dgrd(3)=timedim           !

!
! Define evolving model variables.
! ------ -------- ----- ----------
!
!-----------------------------------------------------------------------
! Time step number and time record indices:
!-----------------------------------------------------------------------
!
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 vidTstep)
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,vidTstep,nf_collective)
# endif
        ierr=nf_put_att_text (ncid, vidTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
!
!-----------------------------------------------------------------------
! Time.
!-----------------------------------------------------------------------
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, vidTime)

# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,vidTime,nf_collective)
# endif
# ifndef AVRH
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', lvar,
     &                                vname(2,indxTime)(1:lvar))
# else
        text='averaged '/ /vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', lvar,
     &                                             text(1:lvar))
# endif
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
        lvar=lenstr(vname(4,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'field',  lvar,
     &                                vname(4,indxTime)(1:lvar))

        call nf_add_attribute(ncid, vidTime, vname(:,indxTime), 5,
     &       NF_FOUT, ierr)

!
!-----------------------------------------------------------------------
! Time2.
!-----------------------------------------------------------------------
!
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_def_var (ncid, vname(1,indxTime2)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, vidTime2)
# ifdef NC4PAR
        ierr=nf_var_par_access(ncid,vidTime2,nf_collective)
# endif
# ifndef AVRH
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name', lvar,
     &                                vname(2,indxTime2)(1:lvar))
# else
        text='averaged '/ /vname(2,indxTime2)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime2, 'long_name', lvar,
     &                                             text(1:lvar))
# endif
        lvar=lenstr(vname(3,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'units',  lvar,
     &                                vname(3,indxTime2)(1:lvar))
        lvar=lenstr(vname(4,indxTime2))
        ierr=nf_put_att_text (ncid, vidTime2, 'field',  lvar,
     &                                vname(4,indxTime2)(1:lvar))

        call nf_add_attribute(ncid, vidTime2, vname(:,indxTime2), 5,
     &       NF_FOUT, ierr)

!=======================================================================
!
! Define primary variables
!
!=======================================================================
!
!-----------------------------------------------------------------------
! U-wind
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pu_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pu_dta))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_pu_dta)(1:lvar),
     &        NF_FOUT, 4, abl3dgrd, vidabl_pu_dta)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_pu_dta,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_pu_dta))
          ierr=nf_put_att_text (ncid, vidabl_pu_dta, 'long_name',
     &          lvar, vname_abl1d(2,indxabl_pu_dta)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_pu_dta)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_pu_dta, 'long_name',
     &         lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_pu_dta))
          ierr=nf_put_att_text (ncid, vidabl_pu_dta, 'units',
     &          lvar, vname_abl1d(3,indxabl_pu_dta)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_pu_dta))
          ierr=nf_put_att_text (ncid, vidabl_pu_dta, 'field',
     &         lvar, vname_abl1d(4,indxabl_pu_dta)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_pu_dta,
     &         vname_abl1d(:,indxabl_pu_dta), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! V-wind
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pv_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pv_dta))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_pv_dta)(1:lvar),
     &        NF_FOUT, 4, abl3dgrd, vidabl_pv_dta)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_pv_dta,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_pv_dta))
          ierr=nf_put_att_text (ncid, vidabl_pv_dta, 'long_name',
     &        lvar, vname_abl1d(2,indxabl_pv_dta)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_pv_dta)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_pv_dta, 'long_name',
     &         lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_pv_dta))
          ierr=nf_put_att_text (ncid, vidabl_pv_dta, 'units',
     &          lvar, vname_abl1d(3,indxabl_pv_dta)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_pv_dta))
          ierr=nf_put_att_text (ncid, vidabl_pv_dta, 'field',
     &          lvar, vname_abl1d(4,indxabl_pv_dta)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_pv_dta,
     &          vname_abl1d(:,indxabl_pv_dta), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! Air temperature
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pt_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pt_dta))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_pt_dta)(1:lvar),
     &          NF_FOUT, 4, abl3dgrd, vidabl_pt_dta)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_pt_dta,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_pt_dta))
          ierr=nf_put_att_text (ncid, vidabl_pt_dta, 'long_name',
     &            lvar, vname_abl1d(2,indxabl_pt_dta)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_pt_dta)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_pt_dta, 'long_name',
     &             lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_pt_dta))
          ierr=nf_put_att_text (ncid, vidabl_pt_dta, 'units',
     &              lvar, vname_abl1d(3,indxabl_pt_dta)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_pt_dta))
          ierr=nf_put_att_text (ncid, vidabl_pt_dta, 'field',
     &              lvar, vname_abl1d(4,indxabl_pt_dta)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_pt_dta,
     &              vname_abl1d(:,indxabl_pt_dta), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! Air specific humidity
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pq_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pq_dta))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_pq_dta)(1:lvar),
     &           NF_FOUT, 4, abl3dgrd, vidabl_pq_dta)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_pq_dta,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_pq_dta))
          ierr=nf_put_att_text (ncid, vidabl_pq_dta, 'long_name',
     &           lvar, vname_abl1d(2,indxabl_pq_dta)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_pq_dta)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_pq_dta, 'long_name',
     &            lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_pq_dta))
          ierr=nf_put_att_text (ncid, vidabl_pq_dta, 'units',
     &            lvar, vname_abl1d(3,indxabl_pq_dta)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_pq_dta))
          ierr=nf_put_att_text (ncid, vidabl_pq_dta, 'field',
     &            lvar, vname_abl1d(4,indxabl_pq_dta)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_pq_dta,
     &            vname_abl1d(:,indxabl_pq_dta), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! U gradient
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pgu_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pgu_dta))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_pgu_dta)(1:lvar),
     &          NF_FOUT, 4, abl3dgrd, vidabl_pgu_dta)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_pgu_dta,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_pgu_dta))
          ierr=nf_put_att_text (ncid, vidabl_pgu_dta, 'long_name',
     &           lvar, vname_abl1d(2,indxabl_pgu_dta)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_pgu_dta)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_pgu_dta, 'long_name',
     &           lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_pgu_dta))
          ierr=nf_put_att_text (ncid, vidabl_pgu_dta, 'units',
     &           lvar, vname_abl1d(3,indxabl_pgu_dta)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_pgu_dta))
          ierr=nf_put_att_text (ncid, vidabl_pgu_dta, 'field',
     &            lvar, vname_abl1d(4,indxabl_pgu_dta)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_pgu_dta,
     &            vname_abl1d(:,indxabl_pgu_dta), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! V gradient
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pgv_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pgv_dta))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_pgv_dta)(1:lvar),
     &          NF_FOUT, 4, abl3dgrd, vidabl_pgv_dta)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_pgv_dta,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_pgv_dta))
          ierr=nf_put_att_text (ncid, vidabl_pgv_dta, 'long_name',
     &          lvar, vname_abl1d(2,indxabl_pgv_dta)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_pgv_dta)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_pgv_dta, 'long_name',
     &          lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_pgv_dta))
          ierr=nf_put_att_text (ncid, vidabl_pgv_dta, 'units',
     &          lvar, vname_abl1d(3,indxabl_pgv_dta)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_pgv_dta))
          ierr=nf_put_att_text (ncid, vidabl_pgv_dta, 'field',
     &          lvar, vname_abl1d(4,indxabl_pgv_dta)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_pgv_dta,
     &          vname_abl1d(:,indxabl_pgv_dta), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! U abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_u_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_u_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_u_abl)(1:lvar),
     &           NF_FOUT, 4, abl3dgrd, vidabl_u_abl)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_u_abl,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_u_abl))
          ierr=nf_put_att_text (ncid, vidabl_u_abl, 'long_name',
     &          lvar, vname_abl1d(2,indxabl_u_abl)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_u_abl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_u_abl, 'long_name',
     &          lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_u_abl))
          ierr=nf_put_att_text (ncid, vidabl_u_abl, 'units',
     &          lvar, vname_abl1d(3,indxabl_u_abl)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_u_abl))
          ierr=nf_put_att_text (ncid, vidabl_u_abl, 'field',
     &          lvar, vname_abl1d(4,indxabl_u_abl)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_u_abl,
     &          vname_abl1d(:,indxabl_u_abl), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! v-abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_v_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_v_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_v_abl)(1:lvar),
     &          NF_FOUT, 4, abl3dgrd, vidabl_v_abl)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_v_abl,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_v_abl))
          ierr=nf_put_att_text (ncid, vidabl_v_abl, 'long_name',
     &           lvar, vname_abl1d(2,indxabl_v_abl)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_v_abl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_v_abl, 'long_name',
     &           lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_v_abl))
          ierr=nf_put_att_text (ncid, vidabl_v_abl, 'units',
     &           lvar, vname_abl1d(3,indxabl_v_abl)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_v_abl))
          ierr=nf_put_att_text (ncid, vidabl_v_abl, 'field',
     &           lvar, vname_abl1d(4,indxabl_v_abl)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_v_abl,
     &           vname_abl1d(:,indxabl_v_abl), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! t_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_t_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_t_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_t_abl)(1:lvar),
     &         NF_FOUT, 4, abl3dgrd, vidabl_t_abl)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_t_abl,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_t_abl))
          ierr=nf_put_att_text (ncid, vidabl_t_abl, 'long_name',
     &          lvar, vname_abl1d(2,indxabl_t_abl)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_t_abl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_t_abl, 'long_name',
     &           lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_t_abl))
          ierr=nf_put_att_text (ncid, vidabl_t_abl, 'units',
     &           lvar, vname_abl1d(3,indxabl_t_abl)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_t_abl))
          ierr=nf_put_att_text (ncid, vidabl_t_abl, 'field',
     &           lvar, vname_abl1d(4,indxabl_t_abl)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_t_abl,
     &           vname_abl1d(:,indxabl_t_abl), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! q_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_q_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_q_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_q_abl)(1:lvar),
     &          NF_FOUT, 4, abl3dgrd, vidabl_q_abl)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_q_abl,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_q_abl))
          ierr=nf_put_att_text (ncid, vidabl_q_abl, 'long_name',
     &           lvar, vname_abl1d(2,indxabl_q_abl)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_q_abl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_q_abl, 'long_name',
     &            lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_q_abl))
          ierr=nf_put_att_text (ncid, vidabl_q_abl, 'units',
     &            lvar, vname_abl1d(3,indxabl_q_abl)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_q_abl))
          ierr=nf_put_att_text (ncid, vidabl_q_abl, 'field',
     &            lvar, vname_abl1d(4,indxabl_q_abl)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_q_abl,
     &            vname_abl1d(:,indxabl_q_abl), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! tke_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_tke_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_tke_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_tke_abl)(1:lvar),
     &          NF_FOUT, 4, abl3dgrd, vidabl_tke_abl)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_tke_abl,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_tke_abl))
          ierr=nf_put_att_text (ncid, vidabl_tke_abl, 'long_name',
     &          lvar, vname_abl1d(2,indxabl_tke_abl)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_tke_abl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_tke_abl, 'long_name',
     &          lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_tke_abl))
          ierr=nf_put_att_text (ncid, vidabl_tke_abl, 'units',
     &           lvar, vname_abl1d(3,indxabl_tke_abl)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_tke_abl))
          ierr=nf_put_att_text (ncid, vidabl_tke_abl, 'field',
     &            lvar, vname_abl1d(4,indxabl_tke_abl)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_tke_abl,
     &            vname_abl1d(:,indxabl_tke_abl), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! mxlm_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_mxlm_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_mxlm_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_mxlm_abl)(1:lvar),
     &          NF_FOUT, 4, abl3dgrd, vidabl_mxlm_abl)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_mxlm_abl,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_mxlm_abl))
          ierr=nf_put_att_text (ncid, vidabl_mxlm_abl, 'long_name',
     &           lvar, vname_abl1d(2,indxabl_mxlm_abl)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_mxlm_abl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_mxlm_abl, 'long_name',
     &            lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_mxlm_abl))
          ierr=nf_put_att_text (ncid, vidabl_mxlm_abl, 'units',
     &            lvar, vname_abl1d(3,indxabl_mxlm_abl)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_mxlm_abl))
          ierr=nf_put_att_text (ncid, vidabl_mxlm_abl, 'field',
     &             lvar, vname_abl1d(4,indxabl_mxlm_abl)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_mxlm_abl,
     &            vname_abl1d(:,indxabl_mxlm_abl), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! mxld_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_mxld_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_mxld_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_mxld_abl)(1:lvar),
     &           NF_FOUT, 4, abl3dgrd, vidabl_mxld_abl)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_mxld_abl,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_mxld_abl))
          ierr=nf_put_att_text (ncid, vidabl_mxld_abl, 'long_name',
     &           lvar, vname_abl1d(2,indxabl_mxld_abl)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_mxld_abl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_mxld_abl, 'long_name',
     &           lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_mxld_abl))
          ierr=nf_put_att_text (ncid, vidabl_mxld_abl, 'units',
     &            lvar, vname_abl1d(3,indxabl_mxld_abl)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_mxld_abl))
          ierr=nf_put_att_text (ncid, vidabl_mxld_abl, 'field',
     &            lvar, vname_abl1d(4,indxabl_mxld_abl)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_mxld_abl,
     &            vname_abl1d(:,indxabl_mxld_abl), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! avm_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_avm_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_avm_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_avm_abl)(1:lvar),
     &          NF_FOUT, 4, abl3dgrd, vidabl_avm_abl)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_avm_abl,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_avm_abl))
          ierr=nf_put_att_text (ncid, vidabl_avm_abl, 'long_name',
     &           lvar, vname_abl1d(2,indxabl_avm_abl)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_avm_abl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_avm_abl, 'long_name',
     &           lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_avm_abl))
          ierr=nf_put_att_text (ncid, vidabl_avm_abl, 'units',
     &            lvar, vname_abl1d(3,indxabl_avm_abl)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_avm_abl))
          ierr=nf_put_att_text (ncid, vidabl_avm_abl, 'field',
     &            lvar, vname_abl1d(4,indxabl_avm_abl)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_avm_abl,
     &            vname_abl1d(:,indxabl_avm_abl), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! avt_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_avt_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_avt_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_avt_abl)(1:lvar),
     &          NF_FOUT, 4, abl3dgrd, vidabl_avt_abl)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_avt_abl,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_avt_abl))
          ierr=nf_put_att_text (ncid, vidabl_avt_abl, 'long_name',
     &           lvar, vname_abl1d(2,indxabl_avt_abl)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_avt_abl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_avt_abl, 'long_name',
     &           lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_avt_abl))
          ierr=nf_put_att_text (ncid, vidabl_avt_abl, 'units',
     &            lvar, vname_abl1d(3,indxabl_avt_abl)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_avt_abl))
          ierr=nf_put_att_text (ncid, vidabl_avt_abl, 'field',
     &            lvar, vname_abl1d(4,indxabl_avt_abl)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_avt_abl,
     &            vname_abl1d(:,indxabl_avt_abl), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! ablh
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_ablh_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_ablh_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_ablh_abl)(1:lvar),
     &           NF_FOUT, 3, abl2dgrd, vidabl_ablh_abl)
# ifdef NC4PAR
          ierr=nf_var_par_access(ncid,vidabl_ablh_abl,nf_collective)
# endif
# ifndef AVRH
          lvar=lenstr(vname_abl1d(2,indxabl_ablh_abl))
          ierr=nf_put_att_text (ncid, vidabl_ablh_abl, 'long_name',
     &            lvar, vname_abl1d(2,indxabl_ablh_abl)(1:lvar))
# else
          text='averaged '/ /vname_abl1d(2,indxabl_ablh_abl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidabl_ablh_abl, 'long_name',
     &            lvar, text(1:lvar))
# endif
          lvar=lenstr(vname_abl1d(3,indxabl_ablh_abl))
          ierr=nf_put_att_text (ncid, vidabl_ablh_abl, 'units',
     &            lvar, vname_abl1d(3,indxabl_ablh_abl)(1:lvar))
          lvar=lenstr(vname_abl1d(4,indxabl_ablh_abl))
          ierr=nf_put_att_text (ncid, vidabl_ablh_abl, 'field',
     &            lvar, vname_abl1d(4,indxabl_ablh_abl)(1:lvar))
!
          call nf_add_attribute(ncid,vidabl_ablh_abl,
     &            vname_abl1d(:,indxabl_ablh_abl), 5, NF_FOUT, ierr)
        endif
!
!-----------------------------------------------------------------------
! zr_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_zr_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_zr_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_zr_abl)(1:lvar),
     &            NF_FOUT, 1, abl3dgrd(3), vidabl_zr_abl)
          ierr=nf_put_att_text (ncid, vidabl_zr_abl, 'long_name',
     &           lvar, vname_abl1d(2,indxabl_zr_abl)(1:lvar))
        endif
!
!-----------------------------------------------------------------------
! zw_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_zw_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_zw_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_zw_abl)(1:lvar),
     &            NF_FOUT, 1, abl3dgrd(3), vidabl_zw_abl)
          ierr=nf_put_att_text (ncid, vidabl_zw_abl, 'long_name',
     &            lvar, vname_abl1d(2,indxabl_zw_abl)(1:lvar))
        endif
!
!-----------------------------------------------------------------------
! Hzr_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_Hzr_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_Hzr_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_Hzr_abl)(1:lvar),
     &          NF_FOUT, 1, abl3dgrd(3), vidabl_Hzr_abl)
          ierr=nf_put_att_text (ncid, vidabl_Hzr_abl, 'long_name',
     &          lvar, vname_abl1d(2,indxabl_Hzr_abl)(1:lvar))
        endif
!
!-----------------------------------------------------------------------
! Hzw_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_Hzw_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_Hzw_abl))
          ierr=nf_def_var (ncid, vname_abl1d(1,indxabl_Hzw_abl)(1:lvar),
     &          NF_FOUT, 1, abl3dgrd(3), vidabl_Hzw_abl)
          ierr=nf_put_att_text (ncid, vidabl_Hzw_abl, 'long_name',
     &          lvar, vname_abl1d(2,indxabl_Hzw_abl)(1:lvar))
        endif
!
!-----------------------------------------------------------------------
! Leave definition mode.
!-----------------------------------------------------------------------
!
        ierr=nf_enddef(ncid)
        MPI_master_only write(stdout,'(6x,4A,1x,A,i4)')
     &                'DEF_ABL_HIS/AVG - Created ',
     &                'new netCDF file ''', ncname(1:lstr), '''.'
     &                 MYID
!
!===================================================================
! Open an existing file and prepare for appending data.
! ==== == ======== ==== === ======= === ========= =====
! Inquire about the dimensions and variables. Check for
! consistency with model dimensions. In the case when file
! is rejected (whether it cannot be opened, or something
! is wrong with its dimensions) create a new file.
!
! After that verify that all necessary variables are already
! defined there and find their netCDF IDs.
!===================================================================
!
      elseif (ncid.eq.-1) then
# ifndef NC4PAR
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
# else
        ierr=nf_open_par (ncname(1:lstr), IOR(nf_write, nf_mpiio),
     &     MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
# endif
        if (ierr. eq. nf_noerr) then
          ierr=checkdims (ncid, ncname, lstr, rec)
          if (ierr .eq. nf_noerr) then
            if (rec_per_file.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, rec_per_file))
            endif
            if (ierr.gt.0) then
              MPI_master_only write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &            ) 'DEF_ABL_HIS/AVG WARNING: Actual number of records',
     &               rec,  'in netCDF file',  '''',  ncname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (rec_per_file.eq.0) then
              total_rec=rec+1           ! <-- set to the next record
# if defined MPI & !defined PARALLEL_FILES
              if (mynode.gt.0) total_rec=total_rec-1
# endif
            endif
            ierr=nf_noerr
          endif
        endif
!
        if (ierr. ne. nf_noerr) then
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
            MPI_master_only write(stdout,'(/1x,4A,2x,A,I4/)')
     &                  'DEF_ABL_HIS/AVG ERROR: ',
     &                  'Cannot open file ''', ncname(1:lstr), '''.'
     &                   MYID
            goto 99                                       !--> ERROR
          endif
# else
          create_new_file=.true.
          goto 10
# endif
        endif
!
! Find netCDF IDs of evolving model variables:
! ---- ------ --- -- -------- ----- ----------
!
!-----------------------------------------------------------------------
! Time step indices:
!-----------------------------------------------------------------------
!
        ierr=nf_inq_varid (ncid, 'time_step', vidTstep)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', ncname(1:lstr)
          goto 99                                         !--> ERROR
        endif
!
!-----------------------------------------------------------------------
! Time.
!-----------------------------------------------------------------------
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),vidTime)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), ncname(1:lstr)
          goto 99                                         !--> ERROR
        endif
!
!-----------------------------------------------------------------------
! Time2.
!-----------------------------------------------------------------------
!
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_inq_varid (ncid,vname(1,indxTime2)(1:lvar),vidTime2)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime2)(1:lvar), ncname(1:lstr)
          goto 99                                         !--> ERROR
        endif
!
!-----------------------------------------------------------------------
! U-wind
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pu_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pu_dta))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_pu_dta)(1:lvar),vidabl_pu_dta)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_pu_dta)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! V-wind
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pv_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pv_dta))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_pv_dta)(1:lvar),vidabl_pv_dta)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_pv_dta)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! Air temperature
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pt_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pt_dta))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_pt_dta)(1:lvar),vidabl_pt_dta)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_pt_dta)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! Air specific humidity
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pq_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pq_dta))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_pq_dta)(1:lvar),vidabl_pq_dta)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_pq_dta)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! U-gradient
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pgu_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pgu_dta))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_pgu_dta)(1:lvar),vidabl_pgu_dta)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_pgu_dta)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! V-gradient
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_pgv_dta)) then
          lvar=lenstr(vname_abl1d(1,indxabl_pgv_dta))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_pgv_dta)(1:lvar),vidabl_pgv_dta)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_pgv_dta)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! u_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_u_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_u_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_u_abl)(1:lvar),vidabl_u_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_u_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! v_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_v_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_v_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_v_abl)(1:lvar),vidabl_v_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_v_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! t_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_t_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_t_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_t_abl)(1:lvar),vidabl_t_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_t_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! q_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_q_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_q_abl))
          ierr=nf_inq_varid (ncid, 
     &          vname_abl1d(1,indxabl_q_abl)(1:lvar),vidabl_q_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_q_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! tke_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_tke_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_tke_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_tke_abl)(1:lvar),vidabl_tke_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_tke_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! mxlm_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_mxlm_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_mxlm_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_mxlm_abl)(1:lvar),vidabl_mxlm_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_mxlm_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! mxld_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_mxld_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_mxld_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_mxld_abl)(1:lvar),vidabl_mxld_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_mxld_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! avm_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_avm_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_avm_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_avm_abl)(1:lvar),vidabl_avm_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_avm_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! avt_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_avt_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_avt_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_avt_abl)(1:lvar),vidabl_avt_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_avt_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! ablh
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_ablh_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_ablh_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_ablh_abl)(1:lvar),vidabl_ablh_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_ablh_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! zr_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_zr_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_zr_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_zr_abl)(1:lvar),vidabl_zr_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_zr_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! zw_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_zw_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_zw_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_zw_abl)(1:lvar),vidabl_zw_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_zw_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! Hzr_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_Hzr_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_Hzr_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_Hzr_abl)(1:lvar),vidabl_Hzr_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_Hzr_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
! Hzw_abl
!-----------------------------------------------------------------------
!
        if (wrt(indxabl_Hzw_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_Hzw_abl))
          ierr=nf_inq_varid (ncid, 
     &         vname_abl1d(1,indxabl_Hzw_abl)(1:lvar),vidabl_Hzw_abl)
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname_abl1d(1,indxabl_Hzw_abl)(1:lvar),
     &                                ncname(1:lstr)
            goto 99                                     !--> ERROR
          endif
        endif
!
!-----------------------------------------------------------------------
!
      MPI_master_only write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'DEF_ABL_HIS/AVG -- Opened ',
     &                     'existing file  from record =', rec
     &                      MYID

# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      else
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,'(/1x,4A,2x,A,I4/)')
     &                'DEF_ABL_HIS/AVG ERROR: ',
     &                'Cannot open file ''', ncname(1:lstr), '''.'
     &                 MYID
          goto 99                                         !--> ERROR
        endif
# endif
      endif             !<--  create_new_file

      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_ABL_HIS/AVG ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr
      endif
   1  format(/1x,'DEF_ABL_HIS/AVG ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)

  99  return
      end

# undef ncname
# undef rec_per_file
# undef wrt
# undef vidTime
# undef vidTime2
# undef vidTstep
# undef vidabl_pu_dta
# undef vidabl_pv_dta
# undef vidabl_pt_dta
# undef vidabl_pq_dta
# undef vidabl_pgu_dta
# undef vidabl_pgv_dta
# undef vidabl_u_abl
# undef vidabl_v_abl
# undef vidabl_t_abl
# undef vidabl_q_abl
# undef vidabl_tke_abl
# undef vidabl_mxlm_abl
# undef vidabl_mxld_abl
# undef vidabl_avm_abl
# undef vidabl_avt_abl
# undef vidabl_ablh_abl
# undef vidabl_zr_abl
# undef vidabl_zw_abl
# undef vidabl_Hzr_abl
# undef vidabl_Hzw_abl

!======================================================================
# ifdef AVERAGES
#  ifndef AVRH
#   define AVRH
#   include "def_abl_his.F"
#  endif
# endif
!======================================================================
#else
# ifndef AVRH
      subroutine def_abl_his_empty
# else
      subroutine def_abl_avg_empty
# endif
      return
      end
#endif /* !ABL1D */
