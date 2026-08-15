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
                               ! Writes requested model
      subroutine wrt_abl_his   ! fields at requested levels
                               ! into history netCDF file.
!     implicit none
      integer ierr,  record, lstr, lvar, lenstr, type
     &   , start(2), count(2), nf_fwrite, cff,ik, varid
# if defined OUTPUTS_SURFACE && ! defined XIOS
     &   , ibuff(6)
# else
     &   , ibuff(4)
# endif
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
! Create/open history file; write grid arrays, if so needed.
!
      call def_abl_his (ncidablhis, nrecablhis, ierr)

      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(ablname)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been
!                                          !!! started, it is assumed
      nrecablhis=max(1,nrecablhis)         !!! that global history
      if (nrpfablhis.eq.0) then            !!! record index is
        record=nrecablhis                  !!! advanced by main.
      else
        record=1+mod(nrecablhis-1, nrpfablhis)
      endif


# undef CR
CR      write(*,*) 'wrt_his: Entry  ' MYID
!
!     Write out evolving model variables:
!     ----- --- -------- ----- ----------
!
!-----------------------------------------------------------------------
!     Time step number and record numbers.
!-----------------------------------------------------------------------
!
      type=filetype_abl
!
      if (iic .eq. 0) then
        ibuff(1)=iic+ntstart
      else
        ibuff(1)=iic
      endif
      ibuff(2)=nrecrst
      ibuff(3)=nrecablhis
# ifdef AVERAGES
      ibuff(4)=nrecablavg
# else
      ibuff(4)=0
# endif
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1


      ierr=nf_put_vara_int (ncidablhis, ablhisTstep, start,
     &                      count, ibuff)
      if (ierr .ne. nf_noerr) then
         MPI_master_only write(stdout,1) 'time_step', record, ierr
     &        MYID
         goto 99                !--> ERROR
      endif

CR      write(*,*) 'wrt_his:  time  ' MYID
!-----------------------------------------------------------------------
!     Time
!-----------------------------------------------------------------------
      ierr=nf_put_var1_FTYPE (ncidablhis, ablhisTime, record, time)
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
      ierr=nf_put_var1_FTYPE (ncidablhis, ablhisTime2, record, time)
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
      if (wrtabl(indxabl_pu_dta)) then
        work3dabl=pu_dta(:,:,:)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_pu_dta,
     &                   vname_abl1d(1,indxabl_pu_dta),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    pv_dta
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_pv_dta)) then
        work3dabl=pv_dta(:,:,:)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_pv_dta,
     &                   vname_abl1d(1,indxabl_pv_dta),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    pt_dta
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_pt_dta)) then
        work3dabl=pt_dta(:,:,:)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_pt_dta,
     &                   vname_abl1d(1,indxabl_pt_dta),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    pq_dta
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_pq_dta)) then
        work3dabl=pq_dta(:,:,:)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_pq_dta,
     &                   vname_abl1d(1,indxabl_pq_dta),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    pgu_dta
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_pgu_dta)) then
        work3dabl=pgu_dta(:,:,:)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_pgu_dta,
     &                   vname_abl1d(1,indxabl_pgu_dta),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    pgv_dta
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_pgv_dta)) then
        work3dabl=pgv_dta(:,:,:)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_pgv_dta,
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
      if (wrtabl(indxabl_u_abl)) then
        work3dabl=u_abl(:,:,:,nt_a)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_u_abl,
     &                   vname_abl1d(1,indxabl_u_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    v_abl
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_v_abl)) then
        work3dabl=v_abl(:,:,:,nt_a)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_v_abl,
     &                   vname_abl1d(1,indxabl_v_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    t_abl
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_t_abl)) then
        work3dabl=tq_abl(:,:,:,nt_a,jp_ta)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_t_abl,
     &                   vname_abl1d(1,indxabl_t_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    q_abl
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_q_abl)) then
        work3dabl=tq_abl(:,:,:,nt_a,jp_qa)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_q_abl,
     &                   vname_abl1d(1,indxabl_q_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    tke_abl
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_tke_abl)) then
        work3dabl=tke_abl(:,:,:,nt_a)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_tke_abl,
     &                   vname_abl1d(1,indxabl_tke_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    mxlm_abl
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_mxlm_abl)) then
        work3dabl=mxlm_abl(:,:,:,nt_a)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_mxlm_abl,
     &                   vname_abl1d(1,indxabl_mxlm_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    mxld_abl
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_mxld_abl)) then
        work3dabl=mxld_abl(:,:,:,nt_a)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_mxld_abl,
     &                   vname_abl1d(1,indxabl_mxld_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    avm_abl
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_avm_abl)) then
        work3dabl=avm_abl(:,:,:,nt_a)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_avm_abl,
     &                   vname_abl1d(1,indxabl_avm_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    avt_abl
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_avt_abl)) then
        work3dabl=avt_abl(:,:,:,nt_a)
        call fillvalue3d(work3dabl,ncidablhis,ablhis_avt_abl,
     &                   vname_abl1d(1,indxabl_avt_abl),
     &                   record,abl3dvar,type)
      endif
!-----------------------------------------------------------------------
!    ablh
!-----------------------------------------------------------------------
      if (wrtabl(indxabl_ablh_abl)) then
        work2d=ablh(:,:,nt_a)
        call fillvalue2d(work2d,ncidablhis,ablhis_ablh_abl,
     &                   vname_abl1d(1,indxabl_ablh_abl),
     &                   record,r2dvar,type)
      endif
! ============================================================
!     ABL1D vertical grid
! ============================================================
      if (iic.eq.0 .or. nrpfablhis.gt.0) then
!-----------------------------------------------------------------------
!    zr_abl
!-----------------------------------------------------------------------
        if (wrtabl(indxabl_zr_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_zr_abl))
          ierr=nf_inq_varid (ncidablhis, 
     &         vname_abl1d(1,indxabl_zr_abl)(1:lvar), varid)
          ierr=nf_put_vara_FTYPE(ncidablhis,varid,1,N_abl,zr_abl)
        endif
!-----------------------------------------------------------------------
!    zw_abl
!-----------------------------------------------------------------------
        if (wrtabl(indxabl_zw_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_zw_abl))
          ierr=nf_inq_varid (ncidablhis, 
     &         vname_abl1d(1,indxabl_zw_abl)(1:lvar), varid)
          ierr=nf_put_vara_FTYPE(ncidablhis,varid,1,N_abl,zw_abl)
        endif
!-----------------------------------------------------------------------
!    Hzr_abl
!-----------------------------------------------------------------------
        if (wrtabl(indxabl_Hzr_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_Hzr_abl))
          ierr=nf_inq_varid (ncidablhis, 
     &         vname_abl1d(1,indxabl_Hzr_abl)(1:lvar), varid)
          ierr=nf_put_vara_FTYPE(ncidablhis,varid,1,N_abl,Hzr_abl)
        endif
!-----------------------------------------------------------------------
!    Hzw_abl
!-----------------------------------------------------------------------
        if (wrtabl(indxabl_Hzw_abl)) then
          lvar=lenstr(vname_abl1d(1,indxabl_Hzw_abl))
          ierr=nf_inq_varid (ncidablhis, 
     &         vname_abl1d(1,indxabl_Hzw_abl)(1:lvar), varid)
          ierr=nf_put_vara_FTYPE(ncidablhis,varid,1,N_abl,Hzw_abl)
        endif
!
      endif
!
! ============================================================

 1    format(/1x,'WRT_ABL_HIS ERROR while writing variable ''', A,
     &     ''' into history file.', /11x, 'Time record:',
     &     I6,3x,'netCDF error code',i4,3x,a,i4)
      goto 100
 99   may_day_flag=3
 100  continue

!
!     Synchronize netCDF file to disk to allow other processes
!     to access data immediately after it is written.
!
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      ierr=nf_close (ncidablhis)
      if (nrpfablhis.gt.0 .and. record.ge.nrpfablhis) ncidablhis=-1
# else
      if (nrpfablhis.gt.0 .and. record.ge.nrpfablhis) then
        ierr=nf_close (ncidablhis)
        ncidablhis=-1
      else
        ierr=nf_sync(ncidablhis)
      endif
# endif
      if (ierr .eq. nf_noerr) then
       MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &       'WRT_ABL_HIS -- wrote ',
     &       'history fields into time record =', record, '/',
     &        nrecablhis  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)')
     &        'WRT_ABL_HIS ERROR: Cannot ',
     &        'synchronize/close history netCDF file.'
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
      subroutine wrt_abl_his_empty
      return
      end
#endif /* !ABL1D */
