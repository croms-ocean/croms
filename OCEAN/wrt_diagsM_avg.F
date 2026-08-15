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
#if defined DIAGNOSTICS_UV && defined AVERAGES
!
!---------------------------------------------------------------
! Write time-averaged diags flux fields into averages netCDF file
! Writes requested model fields at requested levels into
! diagnostics netCDF file.
!---------------------------------------------------------------
!
      subroutine wrt_diagsM_avg

      implicit none
      integer ierr,  record, lstr, lvar, lenstr
     &   , start(2), count(2), ibuff(4), nf_fwrite, type
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
# include "diagnostics.h"
# include "mpi_cpl.h"
# ifdef SEDIMENT
#  include "sediment.h"
# endif
# ifdef BBL
#  include "bbl.h"
# endif
# ifdef SOLVE3D
      integer tile, itrc,i,j,k,ivar
#  ifdef SEDIMENT
     &      , indxWrk
#  endif
#  include "work.h"
# endif
# include "netcdf.inc"
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, MPI_COMM_WORLD, status, ierr)
      endif
# endif
!
!
! Create/open diagnostic file; write grid arrays, if so needed.
!
      call def_diagsM_avg (nciddiaM_avg, nrecdiaM_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(dianameM_avg)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been
!                                          !!! started, it is assumed
      nrecdiaM_avg=max(nrecdiaM_avg,1)     !!! that global history
      if (nrpfdiaM_avg.eq.0) then          !!! record index is
        record=nrecdiaM_avg                !!! advanced by main.
      else
        record=1+mod(nrecdiaM_avg-1, nrpfdiaM_avg)
      endif
!
! Write out evolving model variables:
! ----- --- -------- ----- ----------
!
! Time step number and record numbers.
!
      type=filetype_diaM_avg
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecdiaM_avg
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiaM_avg, diaTstepM_avg,
     &                                             start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step_avg', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged diag Time
!
      ierr=nf_put_var1_FTYPE (nciddiaM_avg, diaTimeM_avg, record,
     &                                              timediaM_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiaM_avg, diaTime2M_avg, record,
     &                                              time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime2))
        write(stdout,1) vname(1,indxTime2)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

!---------------------------------------------------------------
! Momentum diagnostic variables.
!---------------------------------------------------------------
!
!
      do itrc=1,2
        if (wrtdiaM_avg(itrc)) then
!

# if defined TENDENCY
          ivar=r3dvar
# else
          if (itrc.eq.1) then
            ivar=u3dvar
          else
            ivar=v3dvar
          endif
# endif

!
! indxMXadv
!
           workr=MXadv_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMXadv_avg(itrc),
     &                      vname(1,indxMXadv+itrc-1),
     &                      record,ivar,type)
!
! indxMYadv
!
           workr=MYadv_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMYadv_avg(itrc),
     &                      vname(1,indxMYadv+itrc-1),
     &                      record,ivar,type)
!
! indxMVadv
!
           workr=MVadv_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMVadv_avg(itrc),
     &                      vname(1,indxMVadv+itrc-1),
     &                      record,ivar,type)
!
! indxMCor
!
           workr=MCor_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMCor_avg(itrc),
     &                      vname(1,indxMCor+itrc-1),
     &                      record,ivar,type)
!
! indxMPrsgrd
!
           workr=MPrsgrd_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMPrsgrd_avg(itrc),
     &                      vname(1,indxMPrsgrd+itrc-1),
     &                      record,ivar,type)
!
! indxMHmix
!
           workr=MHmix_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMHmix_avg(itrc),
     &                      vname(1,indxMHmix+itrc-1),
     &                      record,ivar,type)
!
! indxMHdiff
!
           workr=MHdiff_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMHdiff_avg(itrc),
     &                      vname(1,indxMHdiff+itrc-1),
     &                      record,ivar,type)
!
! indxMVmix
!
           workr=MVmix_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMVmix_avg(itrc),
     &                      vname(1,indxMVmix+itrc-1),
     &                      record,ivar,type)
!
! indxMVmix2
!
           workr=MVmix2_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMVmix2_avg(itrc),
     &                      vname(1,indxMVmix2+itrc-1),
     &                      record,ivar,type)
# if defined DIAGNOSTICS_BARO
!
! indxMBaro
!
           workr=MBaro_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMBaro_avg(itrc),
     &                      vname(1,indxMBaro+itrc-1),
     &                      record,ivar,type)
# endif
# if defined M3FAST
!
! indxMfast
!
           workr=Mfast_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMfast_avg(itrc),
     &                      vname(1,indxMfast+itrc-1),
     &                      record,ivar,type)
# endif
# ifdef MRL_WCI
!
! indxMvf
!
           workr=Mvf_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMvf_avg(itrc),
     &                      vname(1,indxMvf+itrc-1),
     &                      record,ivar,type)
!
! indxMbrk
!
           workr=Mbrk_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMbrk_avg(itrc),
     &                      vname(1,indxMbrk+itrc-1),
     &                      record,ivar,type)
!
! indxMStCo
!
           workr=MStCo_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMStCo_avg(itrc),
     &                      vname(1,indxMStCo+itrc-1),
     &                      record,ivar,type)
!
! indxMVvf
!
           workr=MVvf_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMVvf_avg(itrc),
     &                      vname(1,indxMVvf+itrc-1),
     &                      record,ivar,type)
!
! indxMPrscrt
!
           workr=MPrscrt_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMPrscrt_avg(itrc),
     &                      vname(1,indxMPrscrt+itrc-1),
     &                      record,ivar,type)
!
! indxMsbk
!
           workr=Msbk_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMsbk_avg(itrc),
     &                      vname(1,indxMsbk+itrc-1),
     &                      record,ivar,type)
!
! indxMbwf
!
           workr=Mbwf_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMbwf_avg(itrc),
     &                      vname(1,indxMbwf+itrc-1),
     &                      record,ivar,type)
!
! indxMfrc
!
           workr=Mfrc_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMfrc_avg(itrc),
     &                      vname(1,indxMfrc+itrc-1),
     &                      record,ivar,type)
# endif
!
! indxMrate
!
           workr=Mrate_avg(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM_avg,diaMrate_avg(itrc),
     &                      vname(1,indxMrate+itrc-1),
     &                      record,ivar,type)
!

        endif
      enddo

  1   format(/1x,'WRT_DIAG_AVG ERROR while writing variable ''', A,
     &          ''' into diag_avg file.', /11x, 'Time record:',
     &                   I6,3x,'netCDF error code',i4,3x,a,i4)
      goto 100
  99  may_day_flag=3
 100  continue

!
!---------------------------------------------------------------
! Synchronize netCDF file to disk to allow other processes
! to access data immediately after it is written.
!---------------------------------------------------------------
!
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      ierr=nf_close (nciddiaM_avg)
      if (nrpfdiaM_avg.gt.0 .and. record.ge.nrpfdiaM_avg) nciddiaM_avg=-1
# else
      if (nrpfdiaM_avg.gt.0 .and. record.ge.nrpfdiaM_avg) then
        ierr=nf_close (nciddiaM_avg)
        nciddiaM_avg=-1
      else
        ierr=nf_sync(nciddiaM_avg)
      endif
# endif


      if (ierr .eq. nf_noerr) then
       MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_DIAGM_AVG -- wrote ',
     &            'diagM_avg fields into time record =', record, '/',
     &             nrecdiaM_avg  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)')
     &             'WRT_DIAGM_AVG ERROR: Cannot ',
     &             'synchronize/close diag_avg netCDF file.'
        may_day_flag=3
      endif

# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode .lt. NNODES-1) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, MPI_COMM_WORLD,  ierr)
      endif
# endif
      return
      end

#else
      subroutine wrt_diagsM_avg_empty
      end
#endif /* (DIAGNOSTICS_UV) && AVERAGES */
