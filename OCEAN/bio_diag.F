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
#ifdef BIOLOGY
      subroutine bio_diag (tile)
      implicit none
# include "param.h"
      integer tile
# include "compute_tile_bounds.h"
      call bio_diag_tile  (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine bio_diag_tile (Istr,Iend,Jstr,Jend)
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"
# include "mpi_cpl.h"
      integer Istr,Iend,Jstr,Jend, i,j,k,itrc, iocheck,
     &                             trd,omp_get_thread_num
      integer nbtrc
# if defined PISCES
      parameter(nbtrc = 9)
# else
      parameter(nbtrc = 4)
# endif
      integer my_itrc(1:nbtrc)
!
!
      real cff1, NSUB
      real*QUAD cff, my_sum(0:nbtrc)
# ifdef MPI
      include 'mpif.h'
      integer size, step, status(MPI_STATUS_SIZE), ierr
      real*QUAD buff(0:nbtrc)
      common /xyz_bio/ buff
# endif

# ifdef BIO_NChlPZD
      my_itrc(1)=iNO3_
      my_itrc(2)=iPhy1
      my_itrc(3)=iZoo1
      my_itrc(4)=iDet1
# elif defined BIO_N2ChlPZD2
      my_itrc(1)=iNO3_
      my_itrc(2)=iPhy1
      my_itrc(3)=iZoo1
      my_itrc(4)=iDet1
# elif defined BIO_BioEBUS
      my_itrc(1)=iNO3_
      my_itrc(2)=iPhy2
      my_itrc(3)=iZoo2
      my_itrc(4)=iDet2
# elif defined PISCES
      my_itrc(1)=iDIC_
      my_itrc(2)=iNO3_
      my_itrc(3)=iOXY_
      my_itrc(4)=iFer_
      my_itrc(5)=iTAL_
      my_itrc(6)=iDOC_
      my_itrc(7)=iPOC_
      my_itrc(8)=iPHY_
      my_itrc(9)=iZOO_
# endif

      if (mod(iic-1,ninfo).eq.0) then
        do itrc=0,nbtrc
          my_sum(itrc)=QuadZero          ! <-- Reset local sums
        enddo

        do j=Jstr,Jend
         do i=Istr,Iend
# ifdef MASKING
           cff1=rmask(i,j)/(pm(i,j)*pn(i,j))  ! <-- grid box area
# else
           cff1=1./(pm(i,j)*pn(i,j))
# endif
           do k=1,N
            cff=cff1*Hz(i,j,k)       ! <-- volume of grid box(i,j,k)
            my_sum(0)=my_sum(0)+cff  ! <-- accumulate volume
             do itrc=1,nbtrc
              my_sum(itrc)=my_sum(itrc)+cff*t(i,j,k,nstp,my_itrc(itrc))
             enddo
           enddo
         enddo
        enddo

!
! Perform global summation: whoever gets first to the critical region
! resets global sums before global summation starts; after the global
! summation is completed, thread, which is the last one to enter the
! critical region, finalizes the computation of diagnostics and
! prints them out.
!
        if (SINGLE_TILE_MODE) then
          NSUB=1
        else
          NSUB=NSUB_X*NSUB_E
        endif

C$OMP CRITICAL (bio_cr_rgn)
          if (bio_count.eq.0) then
            do itrc=0,nbtrc
              global_sum(itrc)=QuadZero   ! <-- Reset global sums
            enddo
          endif
                                       !-->  Perform global summation
          do itrc=0,nbtrc
            global_sum(itrc)=global_sum(itrc)+my_sum(itrc)
          enddo

          bio_count=bio_count+1             ! This counter identifies
          if (bio_count.eq.NSUB) then      ! the last thread, whoever
            bio_count=0                    ! it is, not always master.
# ifdef MPI
            if (NNODES.gt.1) then         ! Perform global summation
              size=NNODES                 ! among MPI processes
   10          step=(size+1)/2
                if (mynode.ge.step .and. mynode.lt.size) then
                  do itrc=0,nbtrc
                    buff(itrc)=global_sum(itrc)  ! This is MPI_Reduce
                  enddo
                  call MPI_Send (buff, nbtrc+1, MPI_DOUBLE_PRECISION,
     &                 mynode-step, 17, MPI_COMM_WORLD,      ierr)
                elseif (mynode .lt. size-step) then
                  call MPI_Recv (buff, nbtrc+1, MPI_DOUBLE_PRECISION,
     &                 mynode+step, 17, MPI_COMM_WORLD, status, ierr)
                  do itrc=0,nbtrc
                    global_sum(itrc)=global_sum(itrc)+buff(itrc)
                  enddo
                endif
               size=step
              if (size.gt.1) goto 10
            endif
            if (mynode.eq.0) then
# endif
            trd=omp_get_thread_num()

            cff=1./global_sum(0)
            do itrc=1,nbtrc
              global_sum(itrc)=cff*global_sum(itrc)
            enddo
            if (first_time.eq.0) then
              first_time=1
# ifdef BIO_NChlPZD
              write(stdout,1) 'STEP','time[DAYS]','       NO3',
     &                   '     PHYTO','       ZOO','    DETRIT',
     &                   'trd'
   1          format(1x,A4,3x,A10,1x,A10,4x,A10,4x,A10,4x,A10,3x,A3)
# elif defined BIO_N2ChlPZD2
              write(stdout,1) 'STEP','time[DAYS]','       NO3',
     &                   '    PHY2 ','    ZOO2 ', '  DET2',
     &                   'trd'
   1          format(1x,A4,3x,A10,1x,A10,4x,A10,4x,A10,4x,A10,3x,A3)
# elif defined BIO_BioEBUS
              write(stdout,1) 'STEP','time[DAYS]','       NO3',
     &                   '    PHY2 ','    ZOO2 ', '  DET2',
     &                   'trd'
   1          format(1x,A4,3x,A10,1x,A10,4x,A10,4x,A10,4x,A10,3x,A3)
# elif defined PISCES
              write(stdout,1) 'STEP','time[DAYS]','       DIC',
     &                   '      NO3','     OXY ', '   FER',
     &                   '      ALK','     DOC ', '   POC',
     &                   '      PHY','     ZOO ', 'trd'
   1          format(1x,A4,3x,A10,1x,A10,4x,A10,4x,A10,4x,A10,4x,
     &               A10,4x,A10,4x,A10,4x,A10,4x,A10,4x,3x,A30)
# endif
            endif
# if ! defined PISCES
            write(stdout,2) iic-1,tdays,
     &                      global_sum(1),
     &                      global_sum(2),
     &                      global_sum(3),
     &                      global_sum(4)
     &                    , trd
   2        format(I6, F12.5, 4(1PE14.7), I3)
# else
            write(stdout,2) iic-1,tdays,
     &                      global_sum(1),
     &                      global_sum(2),
     &                      global_sum(3),
     &                      global_sum(4),
     &                      global_sum(5),
     &                      global_sum(6),
     &                      global_sum(7),
     &                      global_sum(8),
     &                      global_sum(9),
     &                      trd
   2        format(I6, F12.5, 9(1PE14.7), I3)
# endif
# ifdef MPI
            endif    ! <-- mynode.eq.0
# endif
          endif
C$OMP END CRITICAL (bio_cr_rgn)
      endif
      return
      end
#else
      subroutine bio_diag_empty
      return
      end
#endif /* BIOLOGY */

