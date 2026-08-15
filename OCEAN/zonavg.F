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

#ifdef ZONAL_NUDGING

      subroutine zonavg_2d (istr,iend,jstr,jend,A,Azon)
      implicit none
      integer istr,iend,jstr,jend

# include "param.h"
# include "scalars.h"

      real    A(GLOBAL_2D_ARRAY), Azon(0:Mm+1)
      real*8 Ad(GLOBAL_2D_ARRAY), Azond(0:Mm+1)
!$acc declare create(Ad, Azond)
      integer i,j

!-- for parallelization (mean of A in xi direction)
# ifdef MPI
#  include "mpi_cpl.h"
      include 'mpif.h'
      integer ierr, lenode
      real*8 allsumd(0:Mm+1,NNODES)
!$acc declare create( allsumd )      
# endif

!$acc kernels if(compute_on_device) default(present)

      Ad = A
!$acc end kernels

!$acc kernels if(compute_on_device) default(present)
!$acc loop vector reduction(+:Azond)
      do j=jstr-1,jend+1
         Azond(j)=sum(A(istr:iend,j))
      enddo
!$acc end kernels

# ifdef MPI
!$acc host_data if(compute_on_device)
!$acc&   use_device (Azond,allsumd )
      call MPI_ALLGATHER(Azond,(Mm+2),MPI_DOUBLE_PRECISION,
     &                   allsumd,(Mm+2),MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD,ierr)
!$acc end host_data

!$acc kernels if(compute_on_device) default(present)
      Azond=0.
!$acc end kernels

!$acc kernels if(compute_on_device) default(present)
!$acc loop gang reduction(+:Azond)
      do i=1,NP_XI
         lenode=jj*NP_XI + i
!$acc loop vector reduction(+:Azond)
         do j=jstr-1,jend+1
            Azond(j)=Azond(j)+allsumd(j,lenode)
         enddo
      enddo
!$acc end kernels
# endif /* MPI */

!$acc kernels if(compute_on_device) default(present)
      Azon=nint(Azond*1.d6,kind=8)/1.d6/float(LLm)
!$acc end kernels


!     rewind (80+mynode)
!     write(80+mynode) Azon
!     call crash

      return
      end

!===================================================================

      subroutine zonavg_3d (istr,iend,jstr,jend,A,Azon)
      implicit none
      integer istr,iend,jstr,jend

# include "param.h"
# include "scalars.h"

      real    A(GLOBAL_2D_ARRAY,N), Azon(0:Mm+1,N)
      real*8 Ad(GLOBAL_2D_ARRAY,N), Azond(0:Mm+1,N)
!$acc declare create(Ad, Azond)
      integer i,j,k

!-- for parallelization (mean of A in xi direction)
# ifdef MPI
#  include "mpi_cpl.h"
      include 'mpif.h'
      integer ierr, lenode
      real*8 allsumd(0:Mm+1,N,NNODES)
!$acc declare create( allsumd )      
# endif

!$acc kernels if(compute_on_device) default(present)
      Ad = A
!$acc end kernels
      
!$acc kernels if(compute_on_device) default(present)
      do k=1,N
         do j=jstr-1,jend+1
         Azond(j,k)=sum(A(istr:iend,j,k))
         enddo
      enddo
!$acc end kernels

# ifdef MPI
#  ifdef OPENACC
!$acc host_data if(compute_on_device)
!$acc&       use_device(Azond,allsumd)           
#  endif
      call MPI_ALLGATHER(Azond,(Mm+2)*N,MPI_DOUBLE_PRECISION,
     &                   allsumd,(Mm+2)*N,MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD,ierr)
!$acc end host_data

!$acc kernels if(compute_on_device) default(present)
      Azond=0.
!$acc end kernels      
!$acc kernels if(compute_on_device) default(present)
      do i=1,NP_XI
         lenode=jj*NP_XI + i
!!$acc loop gang reduction(+:Azond)
            do j=jstr-1,jend+1
!!$acc loop vector reduction(+:Azond)
            do k=1,N
               Azond(j,k)=Azond(j,k)+allsumd(j,k,lenode)
            enddo
         enddo
      enddo
!$acc end kernels

      
# endif /* MPI */

!$acc kernels if(compute_on_device) default(present)
       Azon=nint(Azond*1.d6,kind=8)/1.d6/float(LLm)
!$acc end kernels

!     rewind (80+mynode)
!     write(80+mynode) Azon
!     call crash

#else
      subroutine zonavg_empty
#endif /* ZONAL_NUDGING */
      return
      end



