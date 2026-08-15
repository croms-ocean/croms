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
#if defined DIAGNOSTICS_UV
!
!---------------------------------------------------------------
! Write diagnostics fields at requested levels into diagnostics
! netCDF file.

!---------------------------------------------------------------
!
      subroutine wrt_diagsM

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
      integer tile,itrc,i,j,k,ivar
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
! Create/open diagnostic file; write grid arrays, if so needed.
!
      call def_diagsM (nciddiaM, nrecdiaM, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(dianameM)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been
!
                                           !!! started, it is assumed
      if (iic.eq.0) nrecdiaM=nrecdiaM+1      !!! that global history
      if (nrpfdiaM.eq.0) then               !!! record index is
        record=nrecdiaM                     !!! advanced by main.
      else
        record=1+mod(nrecdiaM-1, nrpfdiaM)
      endif
!
!---------------------------------------------------------------
! Write out evolving model variables:
!---------------------------------------------------------------
!
! Time step number and record numbers.
!
      type=filetype_diaM
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecdiaM
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiaM, diaTstepM, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiaM, diaTimeM, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time2
!
      ierr=nf_put_var1_FTYPE (nciddiaM, diaTime2M, record, time)
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
        if (wrtdiaM(itrc)) then
          if (itrc.eq.1) then
            ivar=u3dvar
          else
            ivar=v3dvar
          endif
!
! indxMXadv
!
           workr=MXadv(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMXadv(itrc),
     &                      vname(1,indxMXadv+itrc-1),
     &                      record,ivar,type)
!
! indxMYadv
!
           workr=MYadv(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMYadv(itrc),
     &                      vname(1,indxMYadv+itrc-1),
     &                      record,ivar,type)
!
! indxMVadv
!
           workr=MVadv(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMVadv(itrc),
     &                      vname(1,indxMVadv+itrc-1),
     &                      record,ivar,type)
!
! indxMCor
!
           workr=MCor(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMCor(itrc),
     &                      vname(1,indxMCor+itrc-1),
     &                      record,ivar,type)
!
! indxMPrsgrd
!
           workr=MPrsgrd(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMPrsgrd(itrc),
     &                      vname(1,indxMPrsgrd+itrc-1),
     &                      record,ivar,type)
!
! indxMHmix
!
           workr=MHmix(:,:,:,itrc,nstp)
           call fillvalue3d(workr,nciddiaM,diaMHmix(itrc),
     &                      vname(1,indxMHmix+itrc-1),
     &                      record,ivar,type)
!
! indxMHdiff
!
           workr=MHdiff(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMHdiff(itrc),
     &                      vname(1,indxMHdiff+itrc-1),
     &                      record,ivar,type)
!
! indxMVmix
!
           workr=MVmix(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMVmix(itrc),
     &                      vname(1,indxMVmix+itrc-1),
     &                      record,ivar,type)
!
! indxMVmix2
!
           workr=MVmix2(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMVmix2(itrc),
     &                      vname(1,indxMVmix2+itrc-1),
     &                      record,ivar,type)
!
! indxMBaro
!
# if defined DIAGNOSTICS_BARO
           workr=MBaro(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMBaro(itrc),
     &                      vname(1,indxMBaro+itrc-1),
     &                      record,ivar,type)
# endif
!
! indxMfast
!
# if defined M3FAST
           workr=Mfast(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMfast(itrc),
     &                      vname(1,indxMfast+itrc-1),
     &                      record,ivar,type)
# endif
# if defined MRL_WCI
!
! indxMvf
!
           workr=Mvf(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMvf(itrc),
     &                      vname(1,indxMvf+itrc-1),
     &                      record,ivar,type)
!
! indxMbrk
!
           workr=Mbrk(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMbrk(itrc),
     &                      vname(1,indxMbrk+itrc-1),
     &                      record,ivar,type)
!
! indxMStCo
!
           workr=MStCo(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMStCo(itrc),
     &                      vname(1,indxMStCo+itrc-1),
     &                      record,ivar,type)
!
! indxMVvf
!
           workr=MVvf(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMVvf(itrc),
     &                      vname(1,indxMVvf+itrc-1),
     &                      record,ivar,type)
!
! indxMPrscrt
!
           workr=MPrscrt(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMPrscrt(itrc),
     &                      vname(1,indxMPrscrt+itrc-1),
     &                      record,ivar,type)
!
! indxMsbk
!
           workr=Msbk(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMsbk(itrc),
     &                      vname(1,indxMsbk+itrc-1),
     &                      record,ivar,type)
!
! indxMbwf
!
           workr=Mbwf(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMbwf(itrc),
     &                      vname(1,indxMbwf+itrc-1),
     &                      record,ivar,type)
!
! indxMfrc
!
           workr=Mfrc(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMfrc(itrc),
     &                      vname(1,indxMfrc+itrc-1),
     &                      record,ivar,type)
# endif
!
! indxMrate
!
           workr=Mrate(:,:,:,itrc)
           call fillvalue3d(workr,nciddiaM,diaMrate(itrc),
     &                      vname(1,indxMrate+itrc-1),
     &                      record,ivar,type)
!
        endif
      enddo
!

  1   format(/1x,'WRT_DIAGM ERROR while writing variable ''', A,
     &          ''' into diag file.', /11x, 'Time record:',
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
      ierr=nf_close (nciddiaM)
      if (nrpfdiaM.gt.0 .and. record.ge.nrpfdiaM) nciddiaM=-1
# else
      if (nrpfdiaM.gt.0 .and. record.ge.nrpfdiaM) then
        ierr=nf_close (nciddiaM)
        nciddiaM=-1
      else
        ierr=nf_sync(nciddiaM)
      endif
# endif


      if (ierr .eq. nf_noerr) then
        MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_DIAGM -- wrote',
     &            ' diag fields into time record =', record, '/',
     &             nrecdiaM  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)')
     &             'WRT_DIAGM ERROR: Cannot ',
     &             'synchronize/close diag netCDF file.'
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
      subroutine wrt_diagsM_empty
      end
#endif /* (DIAGNOSTICS_UV) */
