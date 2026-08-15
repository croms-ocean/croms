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
#if defined  VIS_COEF_3D || (!defined SOLVE3D && defined UV_VIS_SMAGO)

      subroutine hvisc_coef (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call hvisc_coef_tile (Istr,Iend,Jstr,Jend, A3d(1,1,trd))
      return
      end
!
      subroutine hvisc_coef_tile (Istr,Iend,Jstr,Jend, defrate)
!
!==================================================================
!
! Compute horizontal viscosity
!
! Patrick Marchesiello and Pierrick Penven, IRD 2007
!-----------------------------------------------------------------
! NN_HMIX extension: optional neural-network replacement for
! Smagorinsky viscosity.
!==================================================================
!
!!!---------------------------------------------------------------
# ifdef NN_HMIX
      use nn_hmix_coef_mod, only: nn_hmix_predict, nn_hmix_buffer,
       &                                               alpha_hmix
# endif
!!!---------------------------------------------------------------
!
      implicit none
# include "param.h"
# include "grid.h"
# include "mixing.h"
# include "ocean3d.h"
# include "ocean2d.h"
# include "scalars.h"
      integer Istr,Iend,Jstr,Jend, i,j,k,itrc,
     &        imin,imax,jmin,jmax, kp
      real horcon, cff,cff0,cff1,cff2,cff3,cff4,cff5
# ifdef UV_VIS_SMAGO_3D
      real Ric, Pr
      parameter (horcon=0.01, Ric=0.25, Pr=0.3)
# else
      parameter (horcon=0.1)
# endif
# ifdef SWASH
      real dZdt,dZdt0,Bbr
# endif
      real defrate(PRIVATE_2D_SCRATCH_ARRAY,1:N)
!
!!!---------------------------------------------------------------
# ifdef NN_HMIX
      real kap_nn, nu_nn, cfl_kap_max
      real dudx2_ij, dvdy2_ij, ri_ij
      integer my_rank_hv
#  ifdef MPI
      include 'mpif.h'
      integer mpi_err_hv
#  endif
# endif
!!!---------------------------------------------------------------
!
# include "compute_auxiliary_bounds.h"
!
!!!---------------------------------------------------------------
# ifdef NN_HMIX
#  ifdef MPI
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank_hv, mpi_err_hv)
#  else
      my_rank_hv = 0
#  endif
# endif
!!!---------------------------------------------------------------
!
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        imin=istr
      else
        imin=istr-1
      endif
      if (EASTERN_EDGE) then
        imax=iend
      else
        imax=iend+1
      endif
# else
      imin=istr-1
      imax=iend+1
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        jmin=jstr
      else
        jmin=jstr-1
      endif
      if (NORTHERN_EDGE) then
        jmax=jend
      else
        jmax=jend+1
      endif
# else
      jmin=jstr-1
      jmax=jend+1
# endif
!
# if defined UV_VIS2 && defined GLS_MIXING_3D
!---------------------------------------------------------------
!  GLS horizontal viscosity at RHO points
!
      do j=jmin,jmax
        do i=imin,imax
          cff0=om_r(i,j)*on_r(i,j)
          cff1=0.01*cff0/dt
          do k=1,N
            cff2=0.5*(Akv(i,j,k)+Akv(i,j,k-1))
            visc3d_r(i,j,k)=min(cff1,cff2)
          enddo
        enddo
      enddo

# elif defined UV_VIS2 && defined UV_VIS_SMAGO
!---------------------------------------------------------------
! Smagorinsky viscosity at RHO points
!
!    A = Cs L^2 D
!
!  with deformation rate D:
!
!    D = sqrt[ 2 DijDij ]  -- Dij = 0.5*(dui/dxj + duj/dxi) --
!      = sqrt[ 2 dudx^2 + 2 dvdy^2 + (dvdx+dudy)^2         --> 2D
!            + 2 dwdz^2 + (dudz+dwdx)^2 + (dvdz+dwdy)^2 ]  --> 3D
!
!  L Length-scale: horizontal Lh = sqrt(dx*dy) --> visc3d
!                  vertical   Lv = Hz          --> Akv
!
!  For LES problems, Lilly's buoyancy correction is applied based
!  on the ratio of gradient and critical Ri (mixing goes to 0 for
!  Rig > Ric)
!---------------------------------------------------------------
!
!  Horizontal deformation rate (squared & halved) at RHO point:
!
!     dudx^2 + dvdy^2 + 0.5(dvdx+dudy)^2
!
      do k=1,N
        do j=jmin,jmax
          do i=imin,imax
            defrate(i,j,k)=
     &                    ((u(i+1,j,k,nrhs)-u(i,j,k,nrhs))*pm(i,j))**2
     &                   +((v(i,j+1,k,nrhs)-v(i,j,k,nrhs))*pn(i,j))**2
     &               +0.125*( pn(i,j)*
     &                        (u(i,j+1,k,nrhs)+u(i+1,j+1,k,nrhs)
     &                        -u(i,j-1,k,nrhs)-u(i+1,j-1,k,nrhs))
     &                       +pm(i,j)*
     &                        (v(i+1,j,k,nrhs)+v(i+1,j+1,k,nrhs)
     &                        -v(i-1,j,k,nrhs)-v(i-1,j+1,k,nrhs)) )**2
          enddo
        enddo
      enddo
!
#  ifdef UV_VIS_SMAGO_3D
!
!  Add vertical deformation rate
!
!     dwdz^2 + 0.5*(dudz+dwdx)^2 + 0.5*(dvdz+dwdy)^2
!
      do k=2,N-1
        do j=jmin,jmax
          do i=imin,imax
            defrate(i,j,k)=defrate(i,j,k)
#   ifdef NBQ
     &           +((wz(i,j,k,nrhs)-wz(i,j,k-1,nrhs))/Hzr(i,j,k))**2
#   endif
     &               +0.125*( (u(i,j,k+1,nrhs)+u(i+1,j,k+1,nrhs)
     &                        -u(i,j,k-1,nrhs)-u(i+1,j,k-1,nrhs))
     &                               /(z_r(i,j,k+1)-z_r(i,j,k-1))
#   ifdef NBQ
     &                        +pm(i,j)*
     &                        (wz(i+1,j,k,nrhs)+wz(i+1,j,k+1,nrhs)
     &                        -wz(i-1,j,k,nrhs)-wz(i-1,j,k+1,nrhs))
#   endif
     &                      )**2
     &               +0.125*( (v(i,j,k+1,nrhs)+v(i,j+1,k+1,nrhs)
     &                        -v(i,j,k-1,nrhs)-v(i,j+1,k-1,nrhs))
     &                               /(z_r(i,j,k+1)-z_r(i,j,k-1))
#   ifdef NBQ
     &                        +pn(i,j)*
     &                        (wz(i,j+1,k,nrhs)+wz(i,j+1,k+1,nrhs)
     &                        -wz(i,j-1,k,nrhs)-wz(i,j-1,k+1,nrhs))
#   endif
     &                      )**2
          enddo
        enddo
      enddo
      do j=jmin,jmax
        do i=imin,imax
          defrate(i,j,N)=max(1.5*defrate(i,j,N-1)-0.5*defrate(i,j,N-2),
     &                       0.0)
        enddo
      enddo
#  endif
!
!  Compute lateral viscosity at RHO points
!
      do j=jmin,jmax
        do i=imin,imax
          cff=horcon*om_r(i,j)*on_r(i,j)
          cfl_kap_max = om_r(i,j)*on_r(i,j) / (2.*dt)
          do k=1,N
            cff1=2.*defrate(i,j,k)               ! D^2
#  ifdef UV_VIS_SMAGO_3D
!  Buoyancy correction (Lilly, 1962)
            cff2=sqrt(max(0.,1.-(bvf(i,j,k)/max(cff1,1.e-12))/Ric))
            cff3=cff*sqrt(cff1)*cff2              ! Smago reference nu
#  else
            cff3=cff*sqrt(cff1)
#  endif
!!!---------------------------------------------------------------
# ifdef NN_HMIX
            ! du/dx and dv/dy squared (diagonal strain components)
            dudx2_ij = ((u(i+1,j,k,nrhs)-u(i,j,k,nrhs))*pm(i,j))**2
            dvdy2_ij = ((v(i,j+1,k,nrhs)-v(i,j,k,nrhs))*pn(i,j))**2
            ! Richardson number from buoyancy and deformation rate
#  ifdef UV_VIS_SMAGO_3D
            ri_ij = bvf(i,j,k) / max(cff1, 1.e-12)
#  else
            ri_ij = 0.0
#  endif
            call nn_hmix_predict(
     &        dudx2_ij, dvdy2_ij,
     &        defrate(i,j,k),          ! D^2/2 -- consistent with hdiff
     &        bvf(i,j,k), ri_ij, Hz(i,j,k),
     &        om_r(i,j)*on_r(i,j),
     &        cff3, cff3,
     &        kap_nn, nu_nn,
     &        cfl_kap_max, my_rank_hv, iic,
     &        Istr, Iend, Jstr, Jend, N)
!!!---------------------------------------------------------------
#  ifdef NN_HMIX_VIS
            visc3d_r(i,j,k) = nu_nn
#  else
            visc3d_r(i,j,k) = cff3
#  endif
!!!---------------------------------------------------------------
# else
            visc3d_r(i,j,k) = cff3
# endif  /* NN_HMIX */
!!!---------------------------------------------------------------
          enddo
        enddo
      enddo

#  ifdef UV_VIS_SMAGO_3D
!
!  Compute vertical Smagorinsky viscosity at W points
!
      do k=1,N
        kp=min(k+1,N)
        do j=jmin,jmax
          do i=imin,imax
            Akv(i,j,k)=0.5*(visc3d_r(i,j,k)+visc3d_r(i,j,kp))
     &                                       *pm(i,j)*pn(i,j)
     &                                       *Hz(i,j,k)*Hz(i,j,k)
#   if !defined LMD_MIXING && !defined GLS_MIXING
     &                                              + Akv_bak
#   endif
            Akt(i,j,k,itemp)=Pr*Akv(i,j,k)
#   ifdef SALINITY
            Akt(i,j,k,isalt)=Akt(i,j,k,itemp)
#   endif
          enddo
        enddo
      enddo
#  endif /* UV_VIS_SMAGO_3D */
!
!------------------------------------------------------------
! Boundary conditions
!------------------------------------------------------------
!
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do k=1,N
          do j=JstrV-1,Jend
           visc3d_r(Istr-1,j,k)=visc3d_r(Istr,j,k)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do k=1,N
          do j=JstrV-1,Jend
            visc3d_r(Iend+1,j,k)=visc3d_r(Iend,j,k)
          enddo
        enddo
      endif
#  endif /* !EW_PERIODIC */
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=IstrU-1,Iend
            visc3d_r(i,Jstr-1,k)=visc3d_r(i,Jstr,k)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=IstrU-1,Iend
            visc3d_r(i,Jend+1,k)=visc3d_r(i,Jend,k)
          enddo
        enddo
      endif
#  endif /* !NS_PERIODIC */
!
! Corners
!
#  if !defined EW_PERIODIC && !defined NS_PERIODIC

      if (SOUTHERN_EDGE .and. WESTERN_EDGE) then
        do k=1,N
          visc3d_r(Istr-1,Jstr-1,k)=0.5*
     &                          ( visc3d_r(Istr,Jstr-1,k)
     &                           +visc3d_r(Istr-1,Jstr,k))
        enddo
      endif
      if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
        do k=1,N
          visc3d_r(Iend+1,Jstr-1,k)=0.5*
     &                          (visc3d_r(Iend,Jstr-1,k)
     &                          +visc3d_r(Iend+1,Jstr,k))
        enddo
      endif
      if (NORTHERN_EDGE .and. WESTERN_EDGE) then
        do k=1,N
          visc3d_r(Istr-1,Jend+1,k)=0.5*
     &                          ( visc3d_r(Istr,Jend+1,k)
     &                           +visc3d_r(Istr-1,Jend,k))
        enddo
      endif
      if (NORTHERN_EDGE .and. EASTERN_EDGE) then
        do k=1,N
          visc3d_r(Iend+1,Jend+1,k)=0.5*
     &                          ( visc3d_r(Iend,Jend+1,k)
     &                           +visc3d_r(Iend+1,Jend,k))
        enddo
      endif
#  endif /* !EW_PERIODIC && !NS_PERIODIC */
!
#  if defined BD_EXCHANGE
!      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,visc3d_r)
#  endif
!--------------------------
!  Viscosity at PSI points
!--------------------------
      do k=1,N
        do j=Jstr,JendR
          do i=Istr,IendR
            visc3d_p(i,j,k)=0.25*
     &                  ( visc3d_r(i,j  ,k)+visc3d_r(i-1,j  ,k)
     &                   +visc3d_r(i,j-1,k)+visc3d_r(i-1,j-1,k))
          enddo
        enddo
      enddo
!
#  if defined BD_EXCHANGE
      call exchange_p3d_tile(Istr,Iend,Jstr,Jend,visc3d_p)
#  endif

#  ifdef UV_VIS_SMAGO_3D
#   define k0 1
#   ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          do k=k0,N
            Akv(istr-1,j,k)=Akv(istr,j,k)
            Akt(istr-1,j,k,itemp)=Akt(istr,j,k,itemp)
#    ifdef SALINITY
            Akt(istr-1,j,k,isalt)=Akt(istr,j,k,isalt)
#    endif
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          do k=k0,N
            Akv(iend+1,j,k)=Akv(iend,j,k)
            Akt(iend+1,j,k,itemp)=Akt(iend,j,k,itemp)
#    ifdef SALINITY
            Akt(iend+1,j,k,isalt)=Akt(iend,j,k,isalt)
#    endif
          enddo
        enddo
      endif
#   endif
#   ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          do k=k0,N
            Akv(i,jstr-1,k)=Akv(i,jstr,k)
            Akt(i,jstr-1,k,itemp)=Akt(i,jstr,k,itemp)
#    ifdef SALINITY
            Akt(i,jstr-1,k,isalt)=Akt(i,jstr,k,isalt)
#    endif
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          do k=k0,N
            Akv(i,jend+1,k)=Akv(i,jend,k)
            Akt(i,jend+1,k,itemp)=Akt(i,jend,k,itemp)
#    ifdef SALINITY
            Akt(i,jend+1,k,isalt)=Akt(i,jend,k,isalt)
#    endif
          enddo
        enddo
      endif
#    ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=k0,N
          Akv(istr-1,jstr-1,k)=Akv(istr,jstr,k)
          Akt(istr-1,jstr-1,k,itemp)=Akt(istr,jstr,k,itemp)
#     ifdef SALINITY
          Akt(istr-1,jstr-1,k,isalt)=Akt(istr,jstr,k,isalt)
#     endif
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=k0,N
          Akv(istr-1,jend+1,k)=Akv(istr,jend,k)
          Akt(istr-1,jend+1,k,itemp)=Akt(istr,jend,k,itemp)
#     ifdef SALINITY
          Akt(istr-1,jend+1,k,isalt)=Akt(istr,jend,k,isalt)
#     endif
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=k0,N
          Akv(iend+1,jstr-1,k)=Akv(iend,jstr,k)
          Akt(iend+1,jstr-1,k,itemp)=Akt(iend,jstr,k,itemp)
#     ifdef SALINITY
          Akt(iend+1,jstr-1,k,isalt)=Akt(iend,jstr,k,isalt)
#     endif
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=k0,N
          Akv(iend+1,jend+1,k)=Akv(iend,jend,k)
          Akt(iend+1,jend+1,k,itemp)=Akt(iend,jend,k,itemp)
#     ifdef SALINITY
          Akt(iend+1,jend+1,k,isalt)=Akt(iend,jend,k,isalt)
#     endif
        enddo
      endif
#    endif
#   endif
#  endif /* UV_VIS_SMAGO_3D */
# endif /* UV_VIS_SMAGO  */

!!!---------------------------------------------------------------
# ifdef NN_HMIX
      call nn_hmix_buffer(my_rank_hv, iic)
# endif
!!!---------------------------------------------------------------
#else
      subroutine hvisc_coef_empty
#endif
      return
      end

!
!===================================================================
!                          DIFFUSIVITY
!===================================================================
!
#ifdef DIF_COEF_3D

      subroutine hdiff_coef (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call hdiff_coef_tile (Istr,Iend,Jstr,Jend,A2d(1,1,trd))
      return
      end
!
      subroutine hdiff_coef_tile (Istr,Iend,Jstr,Jend,grd_scale)
!
!==================================================================
!
! Compute horizontal diffusivity coefficients
!
! Patrick Marchesiello and Pierrick Penven, IRD 2007
!------------------------------------------------------------------
! NN_HMIX :: kappa_T from NN replaces Smagorinsky (diff3d_r)
!
!==================================================================
!
!!!---------------------------------------------------------------
# ifdef NN_HMIX
      use nn_hmix_coef_mod, only: nn_hmix_predict, nn_hmix_buffer,
     &                            alpha_hmix
# endif
!!!---------------------------------------------------------------
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k,itrc
      real grd_scale(PRIVATE_2D_SCRATCH_ARRAY)
      real cff, surf, defrate, horcon
# if defined TS_DIF_SMAGO && (defined UV_VIS_SMAGO_3D ||\
                              defined GLS_MIXING_3D)
      real Pr
      parameter (horcon=0.05, Pr=0.3)  ! --> LES SGS model
# else
#  ifdef TS_HADV_C6
      parameter (horcon=1./20.)
#  else
      parameter (horcon=1./12.)
#  endif
# endif
# include "param.h"
# include "grid.h"
# include "mixing.h"
# include "ocean3d.h"
# include "scalars.h"
!
!!!---------------------------------------------------------------
# ifdef NN_HMIX
      real kap_nn, nu_nn, cfl_kap_max, smago_ref
      real dudx2_ij, dvdy2_ij, defrate2_ij
      integer my_rank_hd
#  ifdef MPI
      include 'mpif.h'
      integer mpi_err_hd
#  endif
# endif
!!!---------------------------------------------------------------
!
# include "compute_auxiliary_bounds.h"
!
!!!---------------------------------------------------------------
# ifdef NN_HMIX
#  ifdef MPI
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank_hd, mpi_err_hd)
#  else
      my_rank_hd = 0
#  endif
# endif
!!!---------------------------------------------------------------
!
!$acc kernels if(compute_on_device) default(present)

# if defined GLS_MIXING_3D
!---------------------------------------------------------------
!  GLS horizontal diffusivity at RHO points
!
      do k=1,N
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            diff3d_r(i,j,k)=Pr*visc3d_r(i,j,k)
          enddo
        enddo
      enddo

# elif defined TS_DIF_SMAGO
!---------------------------------------------------------------
!  Smagorinsky diffusion coefficient
!
!  For Laplacian diffusion, the POM formulation is used:
!    A = CA*DX*DY*DEFRATE
!        with DEFRATE=sqrt[du/dx^2 +dvdy^2 + 0.5(dvdx+dudy)^2]
!        DEFRATE is the deformation rate
!
!---------------------------------------------------------------
!
      do k=1,N
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
#  ifndef UV_VIS_SMAGO
            ! Compute deformation rate components explicitly
            dudx2_ij = ((u(i+1,j,k,nrhs)-u(i,j,k,nrhs))*pm(i,j))**2
            dvdy2_ij = ((v(i,j+1,k,nrhs)-v(i,j,k,nrhs))*pn(i,j))**2
            defrate  = sqrt( dudx2_ij + dvdy2_ij
     &               +0.5*(0.25*pn(i,j)*(
     &                     u(i,j+1,k,nrhs)+u(i+1,j+1,k,nrhs)
     &                    -u(i,j-1,k,nrhs)-u(i+1,j-1,k,nrhs))
     &                   + 0.25*pm(i,j)*(
     &                     v(i+1,j,k,nrhs)+v(i+1,j+1,k,nrhs)
     &                    -v(i-1,j,k,nrhs)-v(i-1,j+1,k,nrhs)) )**2)
            surf      = om_r(i,j)*on_r(i,j)
            smago_ref = horcon*surf*defrate
            ! defrate2_ij = D^2/2 to match hvisc_coef defrate convention
            defrate2_ij = 0.5 * defrate * defrate
!!!-------------------TODO-sm-loop-------------------------------
# ifdef NN_HMIX
            cfl_kap_max = surf / (2.*dt)
            call nn_hmix_predict(
     &        dudx2_ij, dvdy2_ij,
     &        defrate2_ij,
     &        bvf(i,j,k), 0.0, Hz(i,j,k), surf,
     &        smago_ref, smago_ref,
     &        kap_nn, nu_nn,
     &        cfl_kap_max, my_rank_hd, iic,
     &        Istr, Iend, Jstr, Jend, N)
            diff3d_r(i,j,k)=kap_nn
!!!---------------------------------------------------------------
# else
            diff3d_r(i,j,k)=smago_ref
# endif  /* NN_HMIX */
!!!---------------------------------------------------------------
#  elif defined UV_VIS_SMAGO_3D
            diff3d_r(i,j,k)=Pr*visc3d_r(i,j,k)
#  else
            diff3d_r(i,j,k)=visc3d_r(i,j,k)
#  endif
          enddo
        enddo
      enddo
# endif /* GLS || SMAGO */
!
! Boundary conditions
!
# if defined GLS_MIXING_3D || defined TS_DIF_SMAGO
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do k=1,N
          do j=JstrV-1,Jend
           diff3d_r(Istr-1,j,k)=diff3d_r(Istr,j,k)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do k=1,N
          do j=JstrV-1,Jend
            diff3d_r(Iend+1,j,k)=diff3d_r(Iend,j,k)
          enddo
        enddo
      endif
#  endif /* !EW_PERIODIC */
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=IstrU-1,Iend
            diff3d_r(i,Jstr-1,k)=diff3d_r(i,Jstr,k)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=IstrU-1,Iend
            diff3d_r(i,Jend+1,k)=diff3d_r(i,Jend,k)
          enddo
        enddo
      endif
#  endif /* !NS_PERIODIC */
# endif

# if defined TS_HADV_RSUP3 || defined TS_HADV_RSUP5

!--------------------------------------------------------------------
! Tracer mixing coefficient for rotated split 3rd/5th order upstream
!--------------------------------------------------------------------
!
      do k=1,N
        do j=Jstr-1,Jend+1
          do i=IstrU-1,Iend+1
            diff3d_u(i,j,k)=horcon*abs(u(i,j,k,nrhs))*om_u(i,j)**3
          enddo
        enddo
        do j=JstrV-1,Jend+1
          do i=Istr-1,Iend+1
            diff3d_v(i,j,k)=horcon*abs(v(i,j,k,nrhs))*on_v(i,j)**3
          enddo
        enddo
      enddo
# endif /* TS_HADV_RSUP3 */
!$acc end kernels

!!!---------------------------------------------------------------
# ifdef NN_HMIX
      call nn_hmix_buffer(my_rank_hd, iic)
# endif
!!!---------------------------------------------------------------

# if defined BD_EXCHANGE
#  ifdef TS_DIF_SMAGO
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,diff3d_r)
#  endif
#  if defined TS_HADV_RSUP3 || defined TS_HADV_RSUP5
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,diff3d_u)
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,diff3d_v)
#  endif
# endif

#else
      subroutine hdiff_coef_empty
#endif
      return
      end
