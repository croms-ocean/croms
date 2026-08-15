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
#if defined DIAGNOSTICS_TS
!
!---------------------------------------------------------------
! Write diagnostics fields at requested levels into diagnostics
! netCDF file.
!---------------------------------------------------------------
!
      subroutine wrt_diags

      implicit none
      integer ierr,  record, lstr, lvar, lenstr, type
     &   , start(2), count(2), ibuff(4), nf_fwrite
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
# include "mpi_cpl.h"
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, MPI_COMM_WORLD, status, ierr)
      endif
# endif
!
! Create/open diagnostic file; write grid arrays, if so needed.
!
      call def_diags (nciddia, nrecdia, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(dianame)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been
!
                                           !!! started, it is assumed
      if (iic.eq.0) nrecdia=nrecdia+1      !!! that global history
      if (nrpfdia.eq.0) then               !!! record index is
        record=nrecdia                     !!! advanced by main.
      else
        record=1+mod(nrecdia-1, nrpfdia)
      endif
!
!---------------------------------------------------------------
! Write out evolving model variables:
!---------------------------------------------------------------
!
! Time step number and record numbers.
!
      type=filetype_dia
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecdia
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddia, diaTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddia, diaTime, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time2
!
      ierr=nf_put_var1_FTYPE (nciddia, diaTime2, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime2))
        write(stdout,1) vname(1,indxTime2)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

!
!---------------------------------------------------------------
! Tracer diagnostic variables.
!---------------------------------------------------------------
!
!
      do itrc=1,NT
        if (wrtdia3D(itrc)) then
!
! indxTXadv
!
           workr=TXadv(:,:,:,itrc)
           call fillvalue3d(workr,nciddia,diaTXadv(itrc),
     &                      vname(1,indxTXadv+itrc-1),
     &                      record,r3dvar,type)
!
! indxTYadv
!
           workr=TYadv(:,:,:,itrc)
           call fillvalue3d(workr,nciddia,diaTYadv(itrc),
     &                      vname(1,indxTYadv+itrc-1),
     &                      record,r3dvar,type)
!
! indxTVadv
!
           workr=TVadv(:,:,:,itrc)
           call fillvalue3d(workr,nciddia,diaTVadv(itrc),
     &                      vname(1,indxTVadv+itrc-1),
     &                      record,r3dvar,type)
!
! indxTHmix
!
           workr=THmix(:,:,:,itrc)
           call fillvalue3d(workr,nciddia,diaTHmix(itrc),
     &                      vname(1,indxTHmix+itrc-1),
     &                      record,r3dvar,type)
!
! indxTVmix
!
           workr=TVmix(:,:,:,itrc)
           call fillvalue3d(workr,nciddia,diaTVmix(itrc),
     &                      vname(1,indxTVmix+itrc-1),
     &                      record,r3dvar,type)
# ifdef DIAGNOSTICS_TSVAR
!
! indxTVmixt
!
           workr=TVmixt(:,:,:,itrc)
           call fillvalue3d(workr,nciddia,diaTVmixt(itrc),
     &                      vname(1,indxTVmixt+itrc-1),
     &                      record,r3dvar,type)
# endif
!
! indxTForc
!
           workr=TForc(:,:,:,itrc)
           call fillvalue3d(workr,nciddia,diaTForc(itrc),
     &                      vname(1,indxTForc+itrc-1),
     &                      record,r3dvar,type)
!
! indxTrate
!
           workr=Trate(:,:,:,itrc)
           call fillvalue3d(workr,nciddia,diaTrate(itrc),
     &                      vname(1,indxTrate+itrc-1),
     &                      record,r3dvar,type)
!
        endif
!
# ifdef DIAGNOSTICS_TS_MLD
!
        if (wrtdia2D(itrc)) then
!
!---------------------------------------------------------------
! Tracer diagnostic variables averaged over the MLD
!---------------------------------------------------------------
!
! indxTXadv_mld
!
           work2d=TXadv_mld(:,:,itrc)
           call fillvalue2d(work2d,nciddia,diaTXadv_mld(itrc),
     &                      vname(1,indxTXadv_mld+itrc-1),
     &                      record,r2dvar,type)
!
! indxTYadv_mld
!
           work2d=TYadv_mld(:,:,itrc)
           call fillvalue2d(work2d,nciddia,diaTYadv_mld(itrc),
     &                      vname(1,indxTYadv_mld+itrc-1),
     &                      record,r2dvar,type)
!
! indxTVadv_mld
!
           work2d=TVadv_mld(:,:,itrc)
           call fillvalue2d(work2d,nciddia,diaTVadv_mld(itrc),
     &                      vname(1,indxTVadv_mld+itrc-1),
     &                      record,r2dvar,type)
!
! indxTHmix_mld
!
           work2d=THmix_mld(:,:,itrc)
           call fillvalue2d(work2d,nciddia,diaTHmix_mld(itrc),
     &                      vname(1,indxTHmix_mld+itrc-1),
     &                      record,r2dvar,type)
!
! indxTVmix_mld
!
           work2d=TVmix_mld(:,:,itrc)
           call fillvalue2d(work2d,nciddia,diaTVmix_mld(itrc),
     &                      vname(1,indxTVmix_mld+itrc-1),
     &                      record,r2dvar,type)
!
! indxTForc_mld
!
           work2d=TForc_mld(:,:,itrc)
           call fillvalue2d(work2d,nciddia,diaTForc_mld(itrc),
     &                      vname(1,indxTForc_mld+itrc-1),
     &                      record,r2dvar,type)
!
! indxTrate_mld
!
           work2d=Trate_mld(:,:,itrc)
           call fillvalue2d(work2d,nciddia,diaTrate_mld(itrc),
     &                      vname(1,indxTrate_mld+itrc-1),
     &                      record,r2dvar,type)
!
! indxTentr_mld
!
           work2d=Tentr_mld(:,:,itrc)
           call fillvalue2d(work2d,nciddia,diaTentr_mld(itrc),
     &                      vname(1,indxTentr_mld+itrc-1),
     &                      record,r2dvar,type)

!
! indxTaver_mld
!
           work2d=Taver_mld(:,:,itrc)
           call fillvalue2d(work2d,nciddia,diaTaver_mld(itrc),
     &                      vname(1,indxTaver_mld+itrc-1),
     &                      record,r2dvar,type)


# ifdef DIAGNOSTICS_TS_MLD_CRIT
       if (itrc.eq.itemp .or. itrc.eq.iCRT2 .or. itrc.eq.iCRT3
     & .or. itrc.eq.iCRT4) then
!
! indxTcrit_mld
!
           work2d=Tcrit_mld(:,:,itrc)
           call fillvalue2d(work2d,nciddia,diaTcrit_mld(itrc),
     &                      vname(1,indxTcrit_mld+itrc-1),
     &                      record,r2dvar,type)
       endif
# endif
        endif
!
# endif /*DIAGNOSTICS_TS_MLD*/
!
      enddo
!

  1   format(/1x,'WRT_DIAG ERROR while writing variable ''', A,
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
      ierr=nf_close (nciddia)
      if (nrpfdia.gt.0 .and. record.ge.nrpfdia) nciddia=-1
# else
      if (nrpfdia.gt.0 .and. record.ge.nrpfdia) then
        ierr=nf_close (nciddia)
        nciddia=-1
      else
        ierr=nf_sync(nciddia)
      endif
# endif


      if (ierr .eq. nf_noerr) then
        MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_DIAG -- wrote ',
     &            'diag fields into time record =', record, '/',
     &             nrecdia  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)')
     &             'WRT_DIAG ERROR: Cannot ',
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
      subroutine wrt_diag_empty
      end
#endif /*DIAGNOSTICS_TS */
