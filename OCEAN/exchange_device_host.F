! $Id: exchange_device_host.F 1615 2014-12-17 13:27:07Z rblod $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!

#include "cppdefs.h"
#if defined OPENACC
          module exchange_device_host

          interface update_device
            module procedure update_device_1d,
     &                       update_device_2d,
     &                       update_device_3d
          end interface

          interface update_host
            module procedure update_host_2D,
     &                       update_host_3D
          end interface

          contains

!
! update_device
! this procedure updates the values of the array tab
! from the host to the (possibly multiple) device(s)
!
          subroutine update_device_1D(tab)
          real,dimension(:) :: tab
# include "param.h"
# include "scalars.h"
          integer tile, subs, trd

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
# ifdef MPI
!$acc set device_num(my_acc_device)
# else
!$acc set device_num(tile)
# endif
!$acc update if(compute_on_device) device(tab)
      enddo
          end subroutine update_device_1D

          subroutine update_device_2D(tab)
          real,dimension(:,:) :: tab
# include "param.h"
# include "scalars.h"
          integer tile, subs, trd

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
# ifdef MPI
!$acc set device_num(my_acc_device)
# else
!$acc set device_num(tile)
# endif
!$acc update if(compute_on_device) device(tab)
      enddo
          end subroutine update_device_2D

          subroutine update_device_3D(tab)
          real,dimension(:,:,:) :: tab
# include "param.h"
# include "scalars.h"
          integer tile, subs, trd

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
# ifdef MPI
!$acc set device_num(my_acc_device)
# else
!$acc set device_num(tile)
# endif
!$acc update if(compute_on_device) device(tab)
      enddo
          end subroutine update_device_3D
! 
! update_host
! this procedure updates the values of the array tab
! from the device to the host
! possibly from multiple devices
!
          subroutine update_host_2D(tab)
      implicit none
# include "param.h"
# include "scalars.h"
          real,dimension(GLOBAL_2D_ARRAY) :: tab
         
          integer tile, subs, trd

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call do_update_host_2D(tile,tab)
      enddo
          end subroutine update_host_2D
          
          subroutine do_update_host_2D(tile,tab)
          implicit none
          integer tile
# include "param.h"
          real,dimension(GLOBAL_2D_ARRAY) :: tab
C$    integer omp_get_thread_num
# include "compute_tile_bounds.h"

          call do_update_host_2D_tile(Istr, Iend, Jstr, Jend, tab)

          return
          end

          subroutine do_update_host_2D_tile(Istr, Iend, Jstr, Jend, tab)
          implicit none
# include "param.h"
# include "scalars.h"
          integer Istr, Iend, Jstr, Jend
          real,dimension(GLOBAL_2D_ARRAY) :: tab
# include "compute_auxiliary_bounds.h"

!$acc update if(compute_on_device) host(tab(IstrR:IendR,JstrR:JendR))

          return
          end

          subroutine update_host_3D(tab)
      implicit none
# include "param.h"
# include "scalars.h"
          real,dimension(GLOBAL_2D_ARRAY,N) :: tab
         
          integer tile, subs, trd

C$OMP PARALLEL DO PRIVATE(tile)
      do tile=0,NSUB_X*NSUB_E-1
        call do_update_host_3D(tile,tab)
      enddo
          end subroutine update_host_3D
          
          subroutine do_update_host_3D(tile,tab)
          implicit none
          integer tile
# include "param.h"
          real,dimension(GLOBAL_2D_ARRAY,N) :: tab
C$    integer omp_get_thread_num
# include "compute_tile_bounds.h"

          call do_update_host_3D_tile(Istr, Iend, Jstr, Jend, tab)

          return
          end

          subroutine do_update_host_3D_tile(Istr, Iend, Jstr, Jend, tab)
          implicit none
# include "param.h"
# include "scalars.h"
          integer Istr, Iend, Jstr, Jend
          real,dimension(GLOBAL_2D_ARRAY,N) :: tab
# include "compute_auxiliary_bounds.h"

!$acc update if(compute_on_device)
!$acc&       host(tab(IstrR:IendR,JstrR:JendR,1:N))

          return
          end
          end module exchange_device_host
#else
          module exchange_device_host_empty

          end module exchange_device_host_empty
#endif
