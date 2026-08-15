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
#if defined LMD_SKPP || defined LMD_BKPP || defined GLS_MIXING
      integer Jwtype(GLOBAL_2D_ARRAY)
  !    real,parameter,dimension(5)::lmd_mu1=(/0.35,0.6,1.0,1.5,1.4/)
  !    lmd_mu1(1)=0.35    !  Define reciprocal of the absorption
  !    lmd_mu1(2)=0.6     !  coefficient for each of two solar
  !    lmd_mu1(3)=1.0     !  wavelength bands as a function 
  !    lmd_mu1(4)=1.5     !  of water type (Ref: Paulson and
  !    lmd_mu1(5)=1.4     !  Simpson, 1977).
  !    real,parameter,dimension(5)::lmd_mu2=(/23.,20.,17.,14.,7.9/)

  !    real,parameter,dimension(5)::lmd_r1=(/0.58,0.62,0.67,0.77,0.78/)
  !    lmd_r1(1)=0.58    !  Define fraction of the total radiance 
  !    lmd_r1(2)=0.62    !  for wavelength band 1 as a function of 
  !    lmd_r1(3)=0.67    !  Jerlov water type. The fraction for
  !    lmd_r1(4)=0.77    !  wavelength band 2 is lmd_r2=1-lmd_r1.
  !    lmd_r1(5)=0.78
       common/lmd_kpp_coeff/Jwtype
#endif
