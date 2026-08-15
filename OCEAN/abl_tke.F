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
#ifdef ABL1D
!======================================================================
      subroutine abl_tke (tile)
      implicit none
# include "param.h"
# include "scalars.h"
# include "private_scratch.h"
      integer tile, trd, omp_get_thread_num
# include "abl1d.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call abl_tke_tile(Istr,Iend,Jstr,Jend,T3d(1,1,trd),T2d(1,1,trd),
     &                         T2d(1,2,trd),T2d(1,3,trd),T2d(1,4,trd),
     &                         T2d(1,5,trd),T1d(1,trd)              )
      return
      end
!======================================================================
      subroutine abl_tke_tile(Istr,Iend,Jstr,Jend,bn2,sh2,
     &                           elem_a,elem_b,elem_c,CF,ikbl)
      implicit none
# include "param.h"
# include "scalars.h"
# include "grid.h"
# include "abl1d.h"
      integer Istr,Iend,Jstr,Jend,i,j,k
      integer imin,imax,jmin,jmax
      real    cff,cff1,cff2
      real    dU,dV,iRic,esrf,etop,zxx
      real    shear,buoy
      real sh2(PRIVATE_1D_SCRATCH_ARRAY,N_abl)
      real bn2(PRIVATE_2D_SCRATCH_ARRAY,N_abl)
      real elem_a(PRIVATE_1D_SCRATCH_ARRAY,N_abl)
      real elem_b(PRIVATE_1D_SCRATCH_ARRAY,N_abl)
      real elem_c(PRIVATE_1D_SCRATCH_ARRAY,N_abl)
      real CF    (PRIVATE_1D_SCRATCH_ARRAY,N_abl)
      integer ikbl (PRIVATE_1D_SCRATCH_ARRAY)
!
!======================================================================
! Extended range  (needed for subgrid scale closures)
!======================================================================
# ifdef EW_PERIODIC
        imin=Istr-2
        imax=Iend+2
# else
      if (WESTERN_EDGE) then
        imin=Istr-1
      else
        imin=Istr-2
      endif
      if (EASTERN_EDGE) then
        imax=Iend+1
      else
        imax=Iend+2
      endif
# endif
# ifdef NS_PERIODIC
        jmin=Jstr-2
        jmax=Jend+2
# else
      if (SOUTHERN_EDGE) then
        jmin=Jstr-1
      else
        jmin=Jstr-2
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend+1
      else
        jmax=Jend+2
      endif
# endif
!======================================================================
      do j=jmin,jmax  ! outer loop
        ! Compute vertical shear
        do k=2,N_abl-1
          do i=imin,imax
            cff      = 1.0 / Hzw_abl(k)**2
            dU = cff* Avm_abl(i,j,k,nt_n)*(u_abl(i,j,k+1,nt_n)
     &                                    -u_abl(i,j,k,nt_n) )**2
            dV = cff* Avm_abl(i,j,k,nt_n)*(v_abl(i,j,k+1,nt_n)
     &                                    -v_abl(i,j,k,nt_n) )**2
            sh2(i,k) = dU+dV   !<-- sh2 = Km ( ( du/dz )^2 + ( dv/dz )^2 )
          enddo
        enddo
        ! Compute brunt-vaisala frequency
        do k=2,N_abl-1
          do i=imin,imax
            cff  = grav * itvref / Hzw_abl(k)
            cff1 =  tq_abl(i,j,k+1,nt_n, jp_ta) - tq_abl(i,j,k  ,nt_n,jp_ta)
            cff2 =  tq_abl(i,j,k+1,nt_n, jp_ta) * tq_abl(i,j,k+1,nt_n,jp_qa)
     &        - tq_abl(i,j,k,nt_n,jp_ta) * tq_abl(i,j,k,nt_n,jp_qa)
            bn2(i,j,k) = cff * ( cff1 + rctv0 * cff2 )
          enddo
        enddo
        ! Terms for the tridiagonal problem
        do k=2,N_abl-1
          do i=imin,imax
            shear    = sh2(i,k)                    ! zsh2 is already multiplied by Avm_abl at this point
            sh2(i,k) = sh2(i,k) / Avm_abl(i,j,k,nt_n)   ! reformulate zsh2 as a 'true' vertical shear for PBLH computation
            buoy     = - avt_abl(i,j,k,nt_n)*bn2(i,j,k)
            !
            elem_a(i,k) = -0.5*dt*Sch* ( Avm_abl(i,j,k,nt_n)
     &                      +Avm_abl(i,j,k-1,nt_n) ) / Hzr_abl(k  ) ! lower-diagonal
            elem_c(i,k) = -0.5*dt*Sch* ( Avm_abl(i,j,k,nt_n)
     &                      +Avm_abl(i,j,k+1,nt_n) ) / Hzr_abl(k+1) ! upper-diagonal
            !
            if( (buoy + shear) .gt. 0.) then    ! Patankar trick to avoid negative values of TKE
              elem_b(i,k) = Hzw_abl(k)-elem_a(i,k)-elem_c(i,k)+Hzw_abl(k)*dt
     &             * Ceps * SQRT(tke_abl(i,j,k,nt_n)) / mxld_abl(i,j,k,nt_n)     ! diagonal
              tke_abl(i,j,k,nt_a) = Hzw_abl(k)*(tke_abl(i,j,k,nt_n)
     &                                             +dt*(buoy+shear))             ! right-hand-side
            else
              elem_b(i,k) = Hzw_abl(k)-elem_a(i,k)-elem_c(i,k )+Hzw_abl(k)*dt
     &             * Ceps * SQRT(tke_abl(i,j,k,nt_n)) / mxld_abl(i,j,k,nt_n)     ! diagonal
     &                     - Hzw_abl(k) * dt * buoy
              tke_abl(i,j,k,nt_a)  = Hzw_abl(k)*(tke_abl(i,j,k,nt_n)
     &                                                    +dt*shear)             ! right-hand-side
            endif
          enddo
        enddo

        do i=imin,imax
          esrf   =  MAX( esfc * ustar2(i,j), tke_min )
          etop   = tke_min
          tke_abl(i,j,1,nt_a) = esrf
          elem_a(i,1) = 0.
          elem_c(i,1) = 0.
          elem_b(i,1) = 1.
          !++ Top Dirichlet B.C.
          elem_a(i,N_abl) = 0.
          elem_c(i,N_abl) = 0.
          elem_b(i,N_abl) = 1.
          tke_abl(i,j,N_abl,nt_a) = etop
          !++
          bn2(i,j,    1) = bn2( i,j,      2)
          sh2(i,      1) = sh2( i,        2)
          bn2(i,j,N_abl) = bn2( i,j,N_abl-1)
          sh2(i,  N_abl) = sh2( i  ,N_abl-1)
        enddo
        !! Matrix inversion
        !! ----------------------------------------------------------
        do i=imin,imax
           cff     =  1. / elem_b(i,1)
           CF(i,1) = - cff * elem_c(i,1)
           tke_abl(i,j,1,nt_a) = cff*tke_abl(i,j,1,nt_a)
        enddo

        do k=2,N_abl
          do i=imin,imax
            cff=1./(elem_b(i,k)+elem_a(i,k)*CF(i,k-1))
            CF(i,k)=-cff*elem_c(i,k)
            tke_abl(i,j,k,nt_a)=cff *( tke_abl(i,j,k  ,nt_a)
     &                - elem_a(i,k) *  tke_abl(i,j,k-1,nt_a) )
           enddo
        enddo

        do k=N_abl-1,1,-1
          do i=imin,imax
            tke_abl(i,j,k,nt_a)=tke_abl(i,j,k  ,nt_a)
     &              + CF(i,k) * tke_abl(i,j,k+1,nt_a)
          enddo
        enddo

        tke_abl(imin:imax,j,1:N_abl,nt_a) =
     &           MAX( tke_abl(imin:imax,j,1:N_abl,nt_a), tke_min )

     !! Diagnose PBL height
     !! ----------------------------------------------------------
        CF  (:,1      ) = 0.
        ikbl(imin:imax) = 0
        iRic            = 1. / Ric_abl
        do k=2,N_abl
          do i=imin,imax
            cff  = Hzw_abl( k-1 )
            cff1 = cff / ( cff + epssfc * ablh (i,j,nt_n) )
            cff  = Hzw_abl( k )
            cff2 = cff / ( cff + epssfc * ablh (i,j,nt_n) )
            cff  = f(i,j)*f(i,j)
            CF(i,k) = CF(i,k-1) + 0.5 * Hzr_abl(k)*(
     &           cff2 * (sh2(i,k  )-iRic*bn2(i,j,k  )-Cek*cff )
     &         + cff1 * (sh2(i,k-1)-iRic*bn2(i,j,k-1)-Cek*cff ))
            if( ikbl(i)==0 .and. CF(i,k).lt.0.) ikbl(i)=k
          enddo
        enddo
        ! finalize the computation of the PBL height
        do i=imin,imax
          k=ikbl(i)
          if(k>2) then ! linear interpolation to get subgrid value of pblh
            ablh(i,j,nt_a) =  ( zw_abl(k-1)*CF(i,k)-zw_abl(k)*CF(i,k-1)
     &                          ) / (  CF(i,k)-          CF(i,k-1) )
          elseif(k==2) then
            ablh(i,j,nt_a) = zw_abl(2   )
          else
            ablh(i,j,nt_a) = zw_abl(N_abl)
          endif
        enddo
      enddo !<-- Outer loop
!======================================================================
!! Diagnostic mixing length computation
!! ----------------------------------------------------------
      zxx = -(LOG(Ceps)-3.*LOG(Cm)+4.*LOG(vkarmn))/LOG(16.)
# define lup elem_a
# define ldw elem_b
      do j=jmin,jmax  ! outer loop
        !
        do i=imin,imax
          mxld_abl (i,j,    1,nt_a) = mxl_min
          mxld_abl (i,j,N_abl,nt_a) = mxl_min
          mxlm_abl (i,j,    1,nt_a) = mxl_min
          mxlm_abl (i,j,N_abl,nt_a) = mxl_min
          ldw(i,    1) = z0_abl(i,j) * Lsfc
          lup(i,N_abl) = mxl_min
        enddo
        !
        do k=2,N_abl-1
          do i=imin,imax
            buoy     = MAX( bn2(i,j,k), 1.E-12 )
            mxlm_abl(i,j,k,nt_a) = MAX( mxl_min,
     &                             SQRT( 2.*tke_abl(i,j,k,nt_a)/buoy))
          enddo
        enddo
        ! Limit mxl
        do k=N_abl-1,1,-1
          do i=imin,imax
            lup(i,k) = MIN(lup(i,k+1)+(zw_abl(k+1)-zw_abl(k)),
     &                                                  mxlm_abl(i,j,k,nt_a))
          enddo
        enddo
                  !
        do k=2,N_abl
          do i=imin,imax
            ldw(i,k) = MIN(ldw(i,k-1)+(zw_abl(k)-zw_abl(k-1)),
     &                                                  mxlm_abl(i,j,k,nt_a))
          enddo
        enddo
        !
        do k=1,N_abl
          do i=imin,imax
            cff = SQRT( ldw(i,k)*lup(i,k) )
            mxlm_abl(i,j,k,nt_a) = MIN ( ldw(i,k),lup(i,k) )
            mxld_abl(i,j,k,nt_a) = MAX ( cff,mxl_min )
          enddo
        enddo
       !
      enddo  !<- outer loop
# undef lup
# undef ldw
!======================================================================
!! Finalize the computation of turbulent visc./diff.
!! ----------------------------------------------------------
      do j=jmin,jmax  ! outer loop
        do k=1,N_abl
          do i=imin,imax
            cff = MAX(phimax,Ric_abl*mxlm_abl(i,j,k,nt_a)*mxld_abl(i,j,k,nt_a)
     &          * MAX(bn2(i,j,k),1.E-12) / tke_abl(i,j,k,nt_a) )
            cff2 =  1. / ( 1. + cff )   !<-- phi_z(z)
            cff  = mxlm_abl(i,j,k,nt_a)*SQRT(tke_abl(i,j,k,nt_a))
            avm_abl(i,j,k,nt_a) = MAX( Cm * cff        , avm_bak   )
            avt_abl(i,j,k,nt_a) = MAX( Ct * cff * cff2 , avt_bak   )
          enddo
        enddo
      enddo
      return
      end
!======================================================================
#else
      subroutine abl_tke_empty
      return
      end
#endif /* !ABL1D */
