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
#if defined CVTK_DEBUG || defined CVTK_DEBUG_ADVANCED

          module debug
# if defined OPENACC
          use exchange_device_host
# endif
!$AGRIF_DO_NOT_TREAT
          integer nbprocs_in
          integer,dimension(:),allocatable :: debug_procs_units
!$AGRIF_END_DO_NOT_TREAT
# if defined CVTK_DEBUG_PERFRST
#  define DEBUG_MPI_MASTER_ONLY
!
!     define DEBUG_MPI_MASTER_ONLY : do the restartability test
!     in case of a mpi run, but it only test the master node 0
          
          integer start_debug,start_debug0,start_debug_read
          parameter (start_debug0=3)
          parameter (start_debug=start_debug0)
          parameter (start_debug_read=start_debug)
# endif        
          logical, parameter :: non_strict = .false.
          real, parameter :: max_tolerance = 1.D-10
          real :: tolerance
          contains
!
!======================================================================
!                    debug_ini
!======================================================================
!
          subroutine debug_ini()
          implicit none
# include "param.h"
# include "scalars.h"
          character(len=80) :: filename
          character(len=80) :: tmp_string
          integer :: nb, new_unit
          
          filename="check_file"
# if !defined CVTK_DEBUG_WRITE
          open(107,file="debug_infos",form="formatted",
     &  status='OLD')
            read(107,*)nbprocs_in
          close(107)
          allocate(debug_procs_units(0:nbprocs_in-1))
          do nb = 0,nbprocs_in-1
            write(filename,'(a11,i0,a1,i0)')
     &    "check_file_",nbprocs_in,'_',nb
!RB : quickfix before updation conv to deal with newunit
!$AGRIF_DO_NOT_TREAT
            open(newunit=new_unit,
     &           file=filename,form='unformatted',
     &           status='OLD')
!$AGRIF_END_DO_NOT_TREAT
            debug_procs_units(nb)=new_unit
          enddo
# else
#  if defined MPI && ! defined DEBUG_MPI_MASTER_ONLY
          if (mynode == 0) then
          open(107,file="debug_infos",form="formatted")
            write(107,*)NNODES
          close(107)
          endif
          print *,'MYNODE = ',mynode
          write(filename,'(a11,i0,a1,i0)')
     &    "check_file_",NNODES,'_',mynode
#  else
          open(107,file="debug_infos",form="formatted")
            write(107,*)1
          close(107)
          write(filename,'(a14)')"check_file_1_0"
#  endif
          open(107,file=filename,form='unformatted',
     &            status='UNKNOWN')
# endif
 
          if (non_strict) then
            tolerance = max_tolerance
          else
            tolerance = 0.
          endif

          return

          end subroutine debug_ini
!
!======================================================================
!                    check_tab2d
!======================================================================
!
          subroutine check_tab2d(tab,comment,typevar,ondevice,correct)
      implicit none
# include "param.h"
# include "scalars.h"
# include "ocean2d.h"
# ifdef MPI
#  include "mpi_cpl.h"
      include 'mpif.h'
# endif
      real,dimension(GLOBAL_2D_ARRAY) :: tab
      integer Lmseq,Mmseq
      integer lb(2),ub(2)
      real,dimension(:,:),allocatable :: tabread
      character*(*) :: comment,typevar
      logical, optional :: ondevice, correct
      integer i0,j0,i,j,i1,j1,i2,j2
# if !defined CVTK_DEBUG_WRITE
      integer i1_r,i2_r,j1_r,j2_r,nb
# endif
# if !defined CVTK_DEBUG_WRITE && defined OPENACC
      real,dimension(GLOBAL_2D_ARRAY) :: temp_host
      logical :: ondevice_l
# else
      logical, parameter :: ondevice_l=.FALSE.
# endif
      logical :: correct_l
# ifdef OPENMP
#  undef mynode
      integer mynode
# endif
      integer iseq1, iseq2, jseq1, jseq2
      logical mystop
# ifdef MPI
      logical globstop
      integer ierr
# endif
      real max_abs
      integer imax,jmax
      character*200 :: comment_bug
          
# if defined OPENMP && !defined MPI
      mynode=0
# endif

!# if !defined MPI && !defined OPENMP
# if defined CVTK_DEBUG_WRITE
#  include "debug_indices.h"

#  if defined CVTK_DEBUG_PERFRST
      if (iic.gt.start_debug
#   if defined DEBUG_MPI_MASTER_ONLY && defined MPI
     &    .and.mynode.eq.0
#   endif
     &    ) then
      write(*,*)'iic.gt.start_debug check_tab2d write'
#  endif /* CVTK_DEBUG_PERFRST */
      
#  if defined MPI
      i1 = i1 + iminmpi - 1
      i2 = i2 + iminmpi - 1
      j1 = j1 + jminmpi - 1
      j2 = j2 + jminmpi - 1
#  endif
      write(107)i1,i2,j1,j2
      write(107)lbound(tab)
#  if defined MPI
     &  + (/iminmpi-1,jminmpi-1/)
#  endif
      write(107)ubound(tab)
#  if defined MPI
     &  + (/iminmpi-1,jminmpi-1/)
#  endif
      write(107)tab

#  if defined CVTK_DEBUG_PERFRST
      endif
#  endif

# else /* defined CVTK_DEBUG_WRITE   */
#  ifdef OPENACC
      ondevice_l = .false.
      if (present(ondevice)) ondevice_l = ondevice
      if (ondevice_l) then
        temp_host=tab
        call update_host(tab)
      endif
#  endif

      correct_l = .false.
      if (present(correct)) correct_l = correct
#  if defined CVTK_DEBUG_PERFRST      
      if (iic.gt.start_debug_read
#   if defined DEBUG_MPI_MASTER_ONLY && defined MPI
     &    .and.mynode.eq.0
#   endif
     &    ) then
        !write(*,*)'iic=',iic
        write(*,*)'iic.gt.start_debug_read check_tab2d read'
#  endif /* CVTK_DEBUG_PERFRST */
        
      if (non_strict) max_abs=0.

      do nb = 0, nbprocs_in-1
      read(debug_procs_units(nb))i1_r,i2_r,j1_r,j2_r
      read(debug_procs_units(nb))lb
      read(debug_procs_units(nb))ub
      allocate(tabread(lb(1):ub(1),lb(2):ub(2)))
      read(debug_procs_units(nb))tabread

#  include "debug_indices.h"
      
      iseq1 = i1_r
      iseq2 = i2_r
      jseq1 = j1_r
      jseq2 = j2_r


      mystop =.false.
      DebugExit : do j=max(j1,jseq1-jminmpi+1),min(j2,jseq2-jminmpi+1)
                  j0 = j+jminmpi-1
                  do i=max(i1,iseq1-iminmpi+1),min(i2,iseq2-iminmpi+1)
                    i0 = i+iminmpi-1
                    if (non_strict) then
                        if (abs(tabread(i0,j0)-tab(i,j))>max_abs) then
                                imax=i;jmax=j
                                max_abs=abs(tabread(i0,j0)-tab(i,j))
                        endif
                    endif
                    if (abs(tabread(i0,j0)-tab(i,j))>tolerance) then
            write(comment_bug,'(A,A,2x,5i4,3e20.12)')'BUGBIN = ',
     &                 comment,mynode,i0,j0,i,j,
     &                 tabread(i0,j0),tab(i,j),
     &                 abs(tabread(i0,j0)-tab(i,j))
!                       print *,'iif = ',iif,iic
#  ifdef AGRIF
                       print *,'GRID# ',Agrif_CFixed()
#  endif                    
                       if (.not.correct_l) then
                        write(*,'(A)')comment_bug
                       mystop=.true.
                       exit DebugExit
                      else
                        comment_bug = TRIM(comment_bug) // '-- CORRECTED'
                        write(*,'(A)')comment_bug
                        tab(i,j) = tabread(i0,j0)
                      endif
                     endif
                   enddo
                 enddo DebugExit
#  if defined MPI && ! defined DEBUG_MPI_MASTER_ONLY
      call MPI_allreduce(mystop,globstop,1,
     &                   MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
      mystop=globstop
#  endif
      if (mystop) then
#  if defined MPI && ! defined DEBUG_MPI_MASTER_ONLY
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Finalize (ierr)
#  endif
        stop                            !-->  EXIT
      endif
            
      deallocate(tabread)
      enddo /* do nb = 0, nbprocs_in-1 */
#  if defined CVTK_DEBUG_PERFRST      
      endif
#  endif
# endif /* !defined CVTK_DEBUG_WRITE */      

# ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
# endif

# if !defined CVTK_DEBUG_WRITE
#  if defined CVTK_DEBUG_PERFRST      
      if (iic.gt.start_debug_read) then
#  endif
      if (non_strict.AND.(max_abs > 0.)) then
        MPI_master_only print *,'CHECK ',comment,' PASSED with MAX_ABS = ',
     & max_abs
#  ifdef AGRIF
     &                       ,' ON GRID ',Agrif_CFixed()
#  endif
     &, 'AT POINT ',imax,jmax
      else
      MPI_master_only print *,'CHECK ',comment,' PASSED'
#  ifdef AGRIF
     &                       ,' ON GRID ',Agrif_CFixed()
#  endif
      endif

#  if defined CVTK_DEBUG_PERFRST      
      endif
#  endif
#  ifdef OPENACC
      if (ondevice_l) tab=temp_host
#  endif
# endif /*  !defined CVTK_DEBUG_WRITE     */

          return
          end subroutine check_tab2d
!
!======================================================================
!                    check_tab3d
!======================================================================
!
          subroutine check_tab3d(tab,comment,typevar,ondevice,correct)
      implicit none
# include "param.h"
# include "scalars.h"  
# include "ocean2d.h" 
# ifdef MPI
#  include "mpi_cpl.h"
      include 'mpif.h'
# endif   
# ifdef MPI
      logical globstop
      integer ierr
# endif 
      real,dimension(GLOBAL_2D_ARRAY,N) :: tab
      integer Lmseq,Mmseq
      integer lb(3),ub(3)
      real,dimension(:,:,:),allocatable :: tabread
      character*(*) :: comment,typevar
      logical, optional :: ondevice, correct
      integer i0,j0,k,i,j,i1,i2,j1,j2
# if !defined CVTK_DEBUG_WRITE
      integer i1_r,i2_r,j1_r,j2_r,nb
# endif
# if !defined CVTK_DEBUG_WRITE && defined OPENACC
      real,dimension(GLOBAL_2D_ARRAY,N) :: temp_host
      logical :: ondevice_l 
# else
      logical, parameter :: ondevice_l=.FALSE.
# endif
      character*80 :: comment_k
# ifdef OPENMP
      integer mynode
# endif
      integer iseq1, iseq2, jseq1, jseq2 
      logical mystop     
      real :: max_abs
      integer :: imax,jmax,kmax
      real :: diff
      logical :: correct_l

# ifdef OPENACC
      ondevice_l = .false.
      if (present(ondevice)) ondevice_l = ondevice
# endif

      correct_l = .false.
      if (present(correct)) correct_l = correct

      if (typevar == 'w') then
        do k=1,N+1
          write(comment_k,'(a,a,i0)')TRIM(comment),' layer : ',k-1
          call check_tab2d(tab(START_2D_ARRAY,k),
     &                     TRIM(comment_k),'r',ondevice_l,correct_l)
        enddo
        return
      elseif (typevar == 'wint') then
        do k=1,N+1
          write(comment_k,'(a,a,i0)')TRIM(comment),' layer : ',k-1
          call check_tab2d(tab(START_2D_ARRAY,k),
     &                         TRIM(comment_k),'rint',ondevice_l,correct_l)
        enddo
        return
      endif 

# if defined OPENMP && !defined MPI
          mynode=0 
# endif

# if defined CVTK_DEBUG_WRITE
#  include "debug_indices.h"
          
#  if defined CVTK_DEBUG_PERFRST      
          if (iic.gt.start_debug
#   if defined DEBUG_MPI_MASTER_ONLY && defined MPI
     &    .and.mynode.eq.0
#   endif
     &         ) then         
            write(*,*)'iic.gt.start_debug check_tab3d write'
#  endif /* defined CVTK_DEBUG_PERFRST   */
            
#  if defined MPI
      i1 = i1 + iminmpi - 1
      i2 = i2 + iminmpi - 1
      j1 = j1 + jminmpi - 1
      j2 = j2 + jminmpi - 1
#  endif
      write(107)i1,i2,j1,j2
      write(107)lbound(tab)
#  if defined MPI
     &  + (/iminmpi-1,jminmpi-1,0/)
#  endif
      write(107)ubound(tab)
#  if defined MPI
     &  + (/iminmpi-1,jminmpi-1,0/)
#  endif
      write(107)tab
#  if defined CVTK_DEBUG_PERFRST      
      endif
#  endif
      
# else /* defined CVTK_DEBUG_WRITE */
      
#  if defined CVTK_DEBUG_PERFRST      
      if (iic.gt.start_debug_read
#   if defined DEBUG_MPI_MASTER_ONLY && defined MPI
     &    .and.mynode.eq.0
#   endif
     &     ) then
        !write(*,*)'iic ',iic
        write(*,*)'iic.gt.start_debug_read check_tab3d read'
#  endif /* defined CVTK_DEBUG_PERFRST  */
      
#  ifdef OPENACC
      if (ondevice_l) then
        temp_host=tab
        call update_host(tab)
      endif
#  endif

      if (non_strict) max_abs=0.

      do nb = 0, nbprocs_in-1
      read(debug_procs_units(nb))i1_r,i2_r,j1_r,j2_r
      read(debug_procs_units(nb))lb
      read(debug_procs_units(nb))ub
      allocate(tabread(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
      read(debug_procs_units(nb))tabread

#  include "debug_indices.h"

      iseq1 = i1_r
      iseq2 = i2_r
      jseq1 = j1_r
      jseq2 = j2_r
      
      mystop =.false.
      DebugExit : do k=1,N
                do j=max(j1,jseq1-jminmpi+1),min(j2,jseq2-jminmpi+1)
                 j0 = j+jminmpi-1
                do i=max(i1,iseq1-iminmpi+1),min(i2,iseq2-iminmpi+1)
                    i0 = i+iminmpi-1
                     if (non_strict) then
 !                       if (abs(tabread(i0,j0,k))>1.D-8) then
 !                         diff =
 !    &             abs(tabread(i0,j0,k)-tab(i,j,k))/abs(tabread(i0,j0,k))
 !                       else
                          diff = abs(tabread(i0,j0,k)-tab(i,j,k))
 !                       endif
                      if (diff>max_abs) then
                              imax=i;jmax=j;kmax=k
                              max_abs=diff
                      endif
                     else
                        diff=abs(tabread(i0,j0,k)-tab(i,j,k))
                     endif
                     if (diff>tolerance)
     &                              then
                       write(*,'(A,A,2x,6i4,3e20.12)')'BUGBIN = ',
     &                 comment,
     &                 mynode,i0,j0,i,j,k,
     &                 tabread(i0,j0,k),tab(i,j,k),
     &                 diff
!                      print *,'iif = ',iif,iic
#  ifdef AGRIF
                       print *,'GRID# ',Agrif_CFixed()
#  endif    
                       if (.not.correct_l) then
                        mystop=.true.
                       exit DebugExit
                      else
                        print *,'CORRECTED'
                        tab(i,j,k) = tabread(i0,j0,k)
                      endif
                      endif
                     enddo
                    enddo
                   enddo DebugExit
      
#  if defined MPI && ! defined DEBUG_MPI_MASTER_ONLY
      call MPI_allreduce(mystop,globstop,1,
     &                   MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
      mystop=globstop
#  endif
      if (mystop) then
#  if defined MPI && ! defined DEBUG_MPI_MASTER_ONLY
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Finalize (ierr)
#  endif
        stop                            !-->  EXIT
      endif
           
      deallocate(tabread)
      enddo /* do nb = 0, nbprocs_in-1 */
#  if defined CVTK_DEBUG_PERFRST      
      endif
#  endif     
# endif /* !defined CVTK_DEBUG_WRITE */

# ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
# endif

# if !defined CVTK_DEBUG_WRITE
#  if defined CVTK_DEBUG_PERFRST        
      if (iic.gt.start_debug_read
#   if defined DEBUG_MPI_MASTER_ONLY && defined MPI
     &    .and.mynode.eq.0
#   endif
     &     ) then
#  endif
      if (non_strict.AND.(max_abs > 0.)) then
        MPI_master_only print *,'CHECK ',comment,' PASSED with MAX_ABS = ',
     & max_abs
#  ifdef AGRIF
     &                       ,' ON GRID ',Agrif_CFixed()
#  endif
     &, 'AT POINT ',imax,jmax,kmax
      else
        MPI_master_only print *,'CHECK ',comment,' PASSED'
#  ifdef AGRIF
     &                       ,' ON GRID ',Agrif_CFixed()
#  endif
      endif
#  if defined CVTK_DEBUG_PERFRST      
      endif
#  endif
#  ifdef OPENACC
      if (ondevice_l) tab=temp_host
#  endif
# endif /* !defined CVTK_DEBUG_WRITE  */
          return
          end subroutine check_tab3d    

      end module debug
          
#else /*  CVTK_DEBUG */

      module debug
      end module debug

#endif                
