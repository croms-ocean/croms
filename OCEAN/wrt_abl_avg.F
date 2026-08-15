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
#if defined ABL1D & defined AVERAGES
                               ! Writes requested model
      subroutine wrt_abl_avg   ! fields at requested levels
                               ! into averages netCDF file.
!     implicit none
      integer ierr,  record, lstr, lvar, lenstr, type
     &   , start(2), count(2), nf_fwrite, cff,ik, varid
     &   , ibuff(4)
      real eps
      parameter (eps=1.D-20)
      real stf_cff
      parameter(stf_cff=86400./0.01)
# if defined MPI & !defined PARALLEL_FILES
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE), blank
# endif
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "averages.h"
# include "forces.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "mixing.h"
# include "mpi_cpl.h"
# include "work.h"
# include "netcdf.inc"
# include "abl1d.h"

# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif

# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, MPI_COMM_WORLD, status, ierr)
      endif
# endif
!
! Create/open averages file; write grid arrays, if so needed.
!
      call def_abl_avg (ncidablavg, nrecablavg, ierr)

      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(ablname_avg)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been
!                                          !!! started, it is assumed
      nrecablavg=max(1,nrecablavg)         !!! that global averages
      if (nrpfablavg.eq.0) then            !!! record index is
        record=nrecablavg                  !!! advanced by main.
      else
        record=1+mod(nrecablavg-1, nrpfablavg)
      endif


# undef CR
CR      write(*,*) 'wrt_avg: Entry  ' MYID
!
!     Write out evolving model variables:
!     ----- --- -------- ----- ----------
!
!-----------------------------------------------------------------------
!     Time step number and record numbers.
!-----------------------------------------------------------------------
!
      type=filetype_abl_avg
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrecablhis
      ibuff(4)=nrecablavg

      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1

      ierr=nf_put_vara_int (ncidablavg, ablavgTstep, start,
     &                      count, ibuff)
      if (ierr .ne. nf_noerr) then
         MPI_master_only write(stdout,1) 'time_step', record, ierr
     &        MYID
         goto 99                !--> ERROR
      endif

CR      write(*,*) 'wrt_avg:  time  ' MYID
!-----------------------------------------------------------------------
!     Time
!-----------------------------------------------------------------------
      ierr=nf_put_var1_FTYPE (ncidablavg, ablavgTime, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        MPI_master_only write(stdout,1)
     &        vname(1,indxTime)(1:lvar), record, ierr
     &        MYID
        goto 99                !--> ERROR
      endif

!-----------------------------------------------------------------------
!     Time2
!-----------------------------------------------------------------------
      ierr=nf_put_var1_FTYPE (ncidablavg, ablavgTime2, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime2))
        MPI_master_only write(stdout,1)
     &        vname(1,indxTime2)(1:lvar), record, ierr
     &        MYID
        goto 99                !--> ERROR
      endif
!
! ============================================================
!     Large scale ABL1D variables
! ============================================================

!-----------------------------------------------------------------------
!    pu_dta
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_pu_dta)) then
        work3dabl=pu_dta_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_pu_dta,
     &                   vname_abl1d(1,indxabl_pu_dta),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    pv_dta
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_pv_dta)) then
        work3dabl=pv_dta_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_pv_dta,
     &                   vname_abl1d(1,indxabl_pv_dta),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    pt_dta
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_pt_dta)) then
        work3dabl=pt_dta_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_pt_dta,
     &                   vname_abl1d(1,indxabl_pt_dta),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    pq_dta
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_pq_dta)) then
        work3dabl=pq_dta_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_pq_dta,
     &                   vname_abl1d(1,indxabl_pq_dta),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    pgu_dta
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_pgu_dta)) then
        work3dabl=pgu_dta_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_pgu_dta,
     &                   vname_abl1d(1,indxabl_pgu_dta),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    pgv_dta
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_pgv_dta)) then
        work3dabl=pgv_dta_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_pgv_dta,
     &                   vname_abl1d(1,indxabl_pgv_dta),
     &                   record,abl3dvar,type)
      endif
!
! ============================================================
!     ABL1D variables
! ============================================================

!-----------------------------------------------------------------------
!    u_abl
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_u_abl)) then
        work3dabl=u_abl_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_u_abl,
     &                   vname_abl1d(1,indxabl_u_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    v_abl
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_v_abl)) then
        work3dabl=v_abl_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_v_abl,
     &                   vname_abl1d(1,indxabl_v_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    t_abl
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_t_abl)) then
        work3dabl=t_abl_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_t_abl,
     &                   vname_abl1d(1,indxabl_t_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    q_abl
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_q_abl)) then
        work3dabl=q_abl_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_q_abl,
     &                   vname_abl1d(1,indxabl_q_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    tke_abl
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_tke_abl)) then
        work3dabl=tke_abl_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_tke_abl,
     &                   vname_abl1d(1,indxabl_tke_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    mxlm_abl
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_mxlm_abl)) then
        work3dabl=mxlm_abl_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_mxlm_abl,
     &                   vname_abl1d(1,indxabl_mxlm_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    mxld_abl
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_mxld_abl)) then
        work3dabl=mxld_abl_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_mxld_abl,
     &                   vname_abl1d(1,indxabl_mxld_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    avm_abl
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_avm_abl)) then
        work3dabl=avm_abl_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_avm_abl,
     &                   vname_abl1d(1,indxabl_avm_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    avt_abl
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_avt_abl)) then
        work3dabl=avt_abl_avg(:,:,:)
        call fillvalue3d(work3dabl,ncidablavg,ablavg_avt_abl,
     &                   vname_abl1d(1,indxabl_avt_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    ablh
!-----------------------------------------------------------------------
      if (wrtabl_avg(indxabl_ablh_abl)) then
        work2d=ablh_avg(:,:)
        call fillvalue2d(work2d,ncidablavg,ablavg_ablh_abl,
     &                   vname_abl1d(1,indxabl_ablh_abl),
     &                   record,r2dvar,type)
      endif
! ============================================================
!     ABL1D vertical grid
! ============================================================
      if (iic.eq.0 .or. nrpfablavg.gt.0) then
!-----------------------------------------------------------------------
!    zr_abl
!-----------------------------------------------------------------------
        if (wrtabl_avg(indxabl_zr_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_zr_abl))
          ierr=nf_inq_varid (ncidablavg, 
     &         vname_abl1d(1,indxabl_zr_abl)(1:lvar), varid)
          ierr=nf_put_vara_FTYPE(ncidablavg,varid,1,N_abl,zr_abl)
        endif
!-----------------------------------------------------------------------
!    zw_abl
!-----------------------------------------------------------------------
        if (wrtabl_avg(indxabl_zw_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_zw_abl))
          ierr=nf_inq_varid (ncidablavg, 
     &         vname_abl1d(1,indxabl_zw_abl)(1:lvar), varid)
          ierr=nf_put_vara_FTYPE(ncidablavg,varid,1,N_abl,zw_abl)
        endif
!-----------------------------------------------------------------------
!    Hzr_abl
!-----------------------------------------------------------------------
        if (wrtabl_avg(indxabl_Hzr_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_Hzr_abl))
          ierr=nf_inq_varid (ncidablavg, 
     &         vname_abl1d(1,indxabl_Hzr_abl)(1:lvar), varid)
          ierr=nf_put_vara_FTYPE(ncidablavg,varid,1,N_abl,Hzr_abl)
        endif
!-----------------------------------------------------------------------
!    Hzw_abl
!-----------------------------------------------------------------------
        if (wrtabl_avg(indxabl_Hzw_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_Hzw_abl))
          ierr=nf_inq_varid (ncidablavg, 
     &         vname_abl1d(1,indxabl_Hzw_abl)(1:lvar), varid)
          ierr=nf_put_vara_FTYPE(ncidablavg,varid,1,N_abl,Hzw_abl)
        endif
!
      endif
!
! ============================================================

 1    format(/1x,'WRT_ABL_AVG ERROR while writing variable ''', A,
     &     ''' into averages file.', /11x, 'Time record:',
     &     I6,3x,'netCDF error code',i4,3x,a,i4)
      goto 100
 99   may_day_flag=3
 100  continue

!
!     Synchronize netCDF file to disk to allow other processes
!     to access data immediately after it is written.
!
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      ierr=nf_close (ncidablavg)
      if (nrpfablavg.gt.0 .and. record.ge.nrpfablavg) ncidablavg=-1
# else
      if (nrpfablavg.gt.0 .and. record.ge.nrpfablavg) then
        ierr=nf_close (ncidablavg)
        ncidablavg=-1
      else
        ierr=nf_sync(ncidablavg)
      endif
# endif
      if (ierr .eq. nf_noerr) then
       MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &       'WRT_ABL_AVG -- wrote ',
     &       'averages fields into time record =', record, '/',
     &        nrecablavg  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)')
     &        'WRT_ABL_avg ERROR: Cannot ',
     &        'synchronize/close avgtory netCDF file.'
        may_day_flag=3
      endif

# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode .lt. NNODES-1) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &        1, MPI_COMM_WORLD,  ierr)
      endif
# endif
      return
      end
!======================================================================
#else
      subroutine wrt_abl_avg_empty
      return
      end
#endif /* ABL1D and AVERAGES */
