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

module tools_string
   ! set of utils to manage string in Fortran

   implicit none

   ! default
   private
   public to_uppercase

contains

   function to_uppercase(str) result(upper_str)
      ! Return uppercase of a string
      implicit none
      character(len=*), intent(in) :: str
      character(len=len(str)) :: upper_str
      integer :: i, ich

      do i = 1, len(str)
         ich = iachar(str(i:i))
         if (ich >= iachar('a') .and. ich <= iachar('z')) then
            upper_str(i:i) = achar(ich - 32)
         else
            upper_str(i:i) = str(i:i)
         end if
      end do
   end function to_uppercase

end module tools_string
