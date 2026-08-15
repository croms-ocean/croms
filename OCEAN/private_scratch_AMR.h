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
      real A2d(N2d,NiSA,0:NPP-1), A3d(N3d,4,0:NPP-1)
     &  ,A3dHz(N3dHz,4,0:NPP-1)	
# ifdef SEDIMENT
      integer B2d(N2d,0:NPP-1)
# endif

      common /private_scratch/ A2d,A3d,A3dHz
#  ifdef SEDIMENT
      common /private_scratch_bis/ B2d
#  endif


