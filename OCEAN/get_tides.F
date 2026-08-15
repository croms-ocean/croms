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
                              ! Read gridded tidal parameters
      subroutine get_tides    ! from forcing netCDF file.
                              ! Scale them.
#if defined SSH_TIDES || defined UV_TIDES
# ifdef TIDES_MAS
      USE MOD_TIDES_MAS
# endif
#ifdef OPENACC
      use exchange_device_host
#endif
      implicit none

# ifndef TIDES_MAS
#  include "param.h"
# endif
# include "scalars.h"
# include "ncscrum.h"
# include "tides.h"
# include "grid.h"
      real cff
      integer i,j, itide, ierr, lstr,lvar, lenstr, nf_fread,
     &        status, varid,ind
# ifdef TIDES_MAS
      integer :: nn,kount,k,l
      integer :: xt
# endif

# ifdef MPI
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif
# include "netcdf.inc"

# ifdef TIDES_MAS
! initialisation du tableau des constituants harmoniques
!-------------------------------------------------
      call tide_data

      l_presence(:)=.false.

! ACTIVATION DE TOUTES LES ONDES
      !l_presence(nu(1:Ntides))=.true.
! ACTIVATION M2
      ! on prend que M2
!      l_presence(im2)=.true.


      rr0(53)=-1000
      rr1(54)=32
# endif
!
! Inquire about the contents of forcing netCDF file:
! variables and dimensions. Check for consistency.
!
      if (may_day_flag.ne.0) return      !-->  EXIT
      lstr=lenstr(frcname)
!
! If not opened yet, open forcing netCDF file for reading. Find and
! save IDs for relevant variables, determine whether momentum stress
! components exist as fields or scalars.
!
      if (ncidfrc.eq.-1) then
        ierr=nf_open (frcname(1:lstr), nf_nowrite, ncidfrc)
        if (ierr .ne. nf_noerr) goto 4                 !--> ERROR
      endif

# ifdef TIDES_MAS
!test
        status=nf_inq_varid(ncidfrc,'tide_period',varid)
        if (status.ne.NF_NOERR) then
          write(6,3) 'id_tide_per',frcname(1:lstr)
          STOP
        endif
        status=nf_get_var(ncidfrc,varid,Tperiod)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_per',frcname(1:lstr)
          STOP
        endif
# endif

      do itide=1,Ntides
# ifdef TIDES_MAS
       ! search for the corresponding tide
       do l=1 ,nmarmax
        if (abs(Tperiod(itide)-360./frequence(l)) .le. 1e-4  ) then
          nu(itide)=l
          ind=itide
          l_presence(nu(itide))=.true.
          exit
       endif
      end do
# else
      ind=itide
# endif
# ifndef TIDES_MAS
!
!  Tidal Period.
!
        status=nf_inq_varid(ncidfrc,'tide_period',varid)
        status=nf_get_var1_FTYPE(ncidfrc,varid,itide,Tperiod(itide))
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_period',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif
        Tperiod(itide)=Tperiod(itide)*3600.
#ifdef OPENACC
!$acc update device( Tperiod(itide) ) 
#endif
# endif

# ifdef SSH_TIDES
!
!  Tidal elevation amplitude and phase.
!
        status=nf_inq_varid(ncidfrc,'tide_Eamp',varid)
        status=nf_fread(SSH_Tamp(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
#ifdef OPENACC
          call update_device(SSH_Tamp(:,:,itide))
#endif
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Eamp',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif

        status=nf_inq_varid(ncidfrc,'tide_Ephase',varid)
        status=nf_fread(SSH_Tphase(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Ephase',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif
#  ifndef TIDES_MAS
        do j=0,LOCALMM+1
          do i=0,LOCALLM+1
            SSH_Tphase(i,j,itide)=SSH_Tphase(i,j,itide)*deg2rad
          enddo
        enddo
#  else
        do j=0,LOCALMM+1
          do i=0,LOCALLM+1
           if (SSH_Tphase(i,j,itide)<0.0)
     &        SSH_Tphase(i,j,itide)=SSH_Tphase(i,j,itide)+360.
          enddo
        enddo
#  endif
#ifdef OPENACC
          call update_device(SSH_Tphase(:,:,itide))
#endif     

# endif
# ifdef UV_TIDES
!
!  Tidal currents angle, phase, major and minor ellipse axis.
!
        status=nf_inq_varid(ncidfrc,'tide_Cangle',varid)
        status=nf_fread(UV_Tangle(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Cangle',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif

        status=nf_inq_varid(ncidfrc,'tide_Cphase',varid)
        status=nf_fread(UV_Tphase(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Cphase',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif

        status=nf_inq_varid(ncidfrc,'tide_Cmax',varid)
        status=nf_fread(UV_Tmajor(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
#ifdef OPENACC
          call update_device(UV_Tmajor(:,:,itide))
#endif
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Cmax',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif

        status=nf_inq_varid(ncidfrc,'tide_Cmin',varid)
        status=nf_fread(UV_Tminor(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
#ifdef OPENACC
          call update_device(UV_Tminor(:,:,itide))
#endif
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Cmin',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif
        do j=0,LOCALMM+1
          do i=0,LOCALLM+1
            UV_Tangle(i,j,itide)=UV_Tangle(i,j,itide)*deg2rad
            UV_Tphase(i,j,itide)=UV_Tphase(i,j,itide)*deg2rad
          enddo
        enddo
#ifdef OPENACC
          call update_device(UV_Tangle(:,:,itide))
          call update_device(UV_Tphase(:,:,itide))
#endif
# endif
# ifdef POT_TIDES
!
!  Tidal potential amplitude and phase.
!
        status=nf_inq_varid(ncidfrc,'tide_Pamp',varid)
        status=nf_fread(POT_Tamp(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,itide,r2dvar)
#ifdef OPENACC
          call update_device(POT_Tamp(:,:,itide))
#endif
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Pamp',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif
        status=nf_inq_varid(ncidfrc,'tide_Pphase',varid)
        status=nf_fread(POT_Tphase(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,itide,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Pphase',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif
        do j=0,LOCALMM+1
          do i=0,LOCALLM+1
            POT_Tphase(i,j,itide)=POT_Tphase(i,j,itide)*deg2rad
          enddo
        enddo
#ifdef OPENACC
          call update_device(POT_Tphase(:,:,itide))
#endif
# endif

      enddo !Ntides
# ifdef TIDES_MAS
      Tperiod=Tperiod*3600.0
# endif


      return

!
! Sort out error messages: The following portion of the code is
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_TIDES - unable to find forcing variable: ',A,
     &                          /,14x,'in forcing netCDF file: ',A)
      goto 99
  4   write(stdout,5) frcname(1:lstr)
  5   format(/,' GET_TIDES - unable to open forcing netCDF file:',
     &                                                        1x,A)
  6   format(/,' GET_TIDES - error while reading variable: ',A,2x,
     &                                      ' at TIME index = ',i4)

 99   may_day_flag=2

#endif /* SSH_TIDES || UV_TIDES */
      return
      end
