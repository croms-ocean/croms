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
#ifndef CHILD_SPG
      subroutine t3dmix (tile)
      implicit none
      integer tile, itrc, trd, omp_get_thread_num
# include "param.h"
# include "ncscrum.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      do itrc=1,NTRA_T3DMIX

# ifdef AGRIF
        if (AGRIF_Root()) then
            call t3dmix_tile (istr,iend,jstr,jend, itrc, 
     &               A2d(1,1,trd),A2d(1,2,trd), A2d(1,3,trd),
     &                            A2d(1,4,trd), A2d(1,5,trd))

        else
            call t3dmix_child_tile (istr,iend,jstr,jend, itrc,
     &               A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &                             A2d(1,4,trd), A2d(1,5,trd))

        endif
# else 
       call t3dmix_tile(istr,iend,jstr,jend, itrc, A2d(1,1,trd),
     &                               A2d(1,2,trd), A2d(1,3,trd),
     &                               A2d(1,4,trd), A2d(1,5,trd))
# endif   /* AGRIF */

       enddo
       return
       end
!
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!
!PARENT
! 
      subroutine t3dmix_tile (istr,iend,jstr,jend, itrc,
     &                    FX,FE, LapT, diff3u,diff3v)
!
#else
!
! CHILD
!
      subroutine t3dmix_child_tile(istr,iend,jstr,jend, itrc,
     &                     FX,FE, LapT, diff3u,diff3v)
!
#endif /* CHILD_SPG */
!
!---------------------------------------------------------------------
! ******************************Common Code***************************
!---------------------------------------------------------------------
!!   
      implicit none
#include "param.h"
#include "ncscrum.h"
      integer itrc, istr,iend,jstr,jend, i,j,k, kmld,
     &        imin,imax,jmin,jmax
      integer kblt(PRIVATE_2D_SCRATCH_ARRAY)
      real    FX(PRIVATE_2D_SCRATCH_ARRAY),     cff,
     &        FE(PRIVATE_2D_SCRATCH_ARRAY),     cff1,
     &        LapT(PRIVATE_2D_SCRATCH_ARRAY),   cff2,
     &        diff3u(PRIVATE_2D_SCRATCH_ARRAY),
     &        diff3v(PRIVATE_2D_SCRATCH_ARRAY)
#include "grid.h"
#include "ocean3d.h"
#include "mixing.h"
#ifdef CLIMAT_TS_MIXH
# include "climat.h"
#endif
#include "scalars.h"
#ifdef DIAGNOSTICS_TS
# include "diagnostics.h"
#endif
#ifdef DIAGNOSTICS_PV
# include "diags_pv.h"
#endif
!
#ifdef AGRIF
# include "zoom.h"
#endif
!
#include "compute_auxiliary_bounds.h"
!
#  if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_CRIT
      real    nearest_diff
      integer indx_depth
      real rho0i, zN2_c
      real hmlp(Istr:Iend,Jstr:Jend)
# endif
#ifdef CHILD_SPG
# define TCLM tsponge
#else
# define TCLM tclm
#endif
!
#ifdef MASKING
# define SWITCH *
#else
# define SWITCH !
#endif
!
#ifndef EW_PERIODIC
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
#else
      imin=istr-1
      imax=iend+1
#endif
#ifndef NS_PERIODIC
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
#else
      jmin=jstr-1
      jmax=jend+1
#endif
!$acc kernels if(compute_on_device) default(present)
      do k=1,N

!++=============================================================
!++ Compute total diffusivity according to model configuration
!++=============================================================
        do j=jmin,jmax
          do i=imin,imax+1
            diff3u(i,j)= 
#ifdef TS_DIF2
     &                      0.5*(diff2(i,j,itrc)+diff2(i-1,j,itrc))
# if defined DIF_COEF_3D && defined TS_DIF_SMAGO
     &                     +0.5*(diff3d_r(i,j,k)+diff3d_r(i-1,j,k))
# endif
#elif defined TS_DIF4
     &                     +sqrt(
     &                      0.5*(diff4(i,j,itrc)+diff4(i-1,j,itrc))
# ifdef DIF_COEF_3D
#  ifdef TS_DIF_SMAGO
     &                     +0.5*(diff3d_r(i,j,k)+diff3d_r(i-1,j,k))
     &                                         *om_u(i,j)*on_u(i,j)
#  endif
#  ifdef TS_HADV_RSUP3
     &                     +diff3d_u(i,j,k)
#  endif
     &                           )
# endif
#endif
          enddo
        enddo
        do j=jmin,jmax+1
          do i=imin,imax
            diff3v(i,j)= 
#ifdef TS_DIF2
     &                      0.5*(diff2(i,j,itrc)+diff2(i,j-1,itrc))
# if defined DIF_COEF_3D && defined TS_DIF_SMAGO
     &                     +0.5*(diff3d_r(i,j,k)+diff3d_r(i,j-1,k))
# endif
#elif defined TS_DIF4
     &                     +sqrt(
     &                      0.5*(diff4(i,j,itrc)+diff4(i,j-1,itrc))
# ifdef DIF_COEF_3D
#  ifdef TS_DIF_SMAGO
     &                     +0.5*(diff3d_r(i,j,k)+diff3d_r(i,j-1,k))
     &                                         *om_v(i,j)*on_v(i,j)
#  endif
#  ifdef TS_HADV_RSUP3
     &                     +diff3d_v(i,j,k)
#  endif
     &                           )
# endif
#endif
          enddo
        enddo
 
#ifdef TS_DIF2
!
!--------------------------------------------------------------------
!  Add in horizontal Laplacian diffusion along constant S-surfaces.
!--------------------------------------------------------------------
!
!  Compute XI- and ETA-components of diffusive tracer flux.
!
        do j=jstr,jend
          do i=istr,iend+1
              FX(i,j)=0.5*diff3u(i,j)
     &                   *pmon_u(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))*(
     &                     t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc)
# if defined CLIMAT_TS_MIXH || defined CLIMAT_TS_MIXH_FINE
     &                        -TCLM(i,j,k,itrc)+TCLM(i-1,j,k,itrc)
# endif
     &                                       ) SWITCH umask(i,j)
# if defined WET_DRY
     &                   *umask_wet(i,j)
# endif
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend
              FE(i,j)=0.5*diff3v(i,j)
     &                   *pnom_v(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))*(
     &                    t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc)
# if defined CLIMAT_TS_MIXH || defined CLIMAT_TS_MIXH_FINE
     &                       -TCLM(i,j,k,itrc)+TCLM(i,j-1,k,itrc)
# endif
     &                                        ) SWITCH vmask(i,j)
# if defined WET_DRY
     &                   *vmask_wet(i,j)
# endif
          enddo
        enddo
!
!  Add in horizontal diffusion of tracer [T m^3/s].
!
        do j=jstr,jend
          do i=istr,iend
            cff1=pm(i,j)*pn(i,j)
            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+dt*cff1
     &                 *(FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))
     &                                             /Hz(i,j,k)
          enddo
        enddo

# if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_PV
!
! Tracer diagnostics
!
        do j=jstr,jend
          do i=istr,iend
            cff1=pm(i,j)*pn(i,j)
            THmix(i,j,k,itrc)=FX(i+1,j)-FX(i,j)
     &                       +FE(i,j+1)-FE(i,j)
#  ifdef MASKING
     &                                        * rmask(i,j)
#  endif

            Trate(i,j,k,itrc)=(Hz(i,j,k)*t(i,j,k,nnew,itrc)
     &                    -Hz_bak(i,j,k)*t(i,j,k,nstp,itrc))
     &                                            /(dt*cff1)
#  ifdef MASKING
     &                                        * rmask(i,j)
#  endif
!
! Divide all diagnostic terms by the cell volume
! (Hz(i,j,k,itrc)/(pm(i,j).*pn(i,j)). There after the unit
! of diag terms are: (unit of tracers)* s-1.
!
            THmix(i,j,k,itrc)=THmix(i,j,k,itrc)*cff1/Hz(i,j,k)
            Trate(i,j,k,itrc)=Trate(i,j,k,itrc)*cff1/Hz(i,j,k)
          enddo
        enddo
# endif /* DIAGNOSTICS */

#endif /* TS_DIF2 */
!
!
#ifdef TS_DIF4
!
!--------------------------------------------------------------------
!  Compute/Add in horizontal biharmonic diffusion along constant
!  S-surfaces.
!--------------------------------------------------------------------
!  The biharmonic operator is computed by applying the Laplacian
!  operator twice.
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
!  Compute horizontal tracer flux in the XI-direction at U-points.
!
          do j=jmin,jmax
            do i=imin,imax+1
            FX(i,j)=0.5*diff3u(i,j)
     &                 *pmon_u(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))*(
     &                   t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc)
# if defined CLIMAT_TS_MIXH || defined CLIMAT_TS_MIXH_FINE
     &                      -TCLM(i,j,k,itrc)+TCLM(i-1,j,k,itrc)
# endif
     &                                        ) SWITCH umask(i,j)
          enddo
        enddo
!
!  Compute horizontal tracer flux in the ETA-direction at V-points.
!
        do j=jmin,jmax+1
          do i=imin,imax
            FE(i,j)=0.5*diff3v(i,j)
     &                 *pnom_v(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))*(
     &                     t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc)
# if defined CLIMAT_TS_MIXH || defined CLIMAT_TS_MIXH_FINE
     &                        -TCLM(i,j,k,itrc)+TCLM(i,j-1,k,itrc)
# endif
     &                                        ) SWITCH vmask(i,j)
          enddo
        enddo
!
!  Compute first Laplacian, without mixing coefficient.
!  Multiply by the metrics of the second Laplacian.
!  Save into work array "LapT".
!
         do j=jmin,jmax
            do i=imin,imax
            LapT(i,j)=(FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))
     &                          *pm(i,j)*pn(i,j)/Hz(i,j,k)
          enddo
        enddo
!
!  Apply boundary conditions (except periodic; closed or gradient)
!  to the first Laplacian.
!
# ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jmin,jmax
#  ifndef OBC_WEST
            LapT(istr-1,j)=0.
#  else
            LapT(istr-1,j)=LapT(istr,j)
#  endif
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jmin,jmax
#  ifndef OBC_EAST
            LapT(iend+1,j)=0.
#  else
            LapT(iend+1,j)=LapT(iend,j)
#  endif
          enddo
        endif
# endif /* !EW_PERIODIC */
# ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=imin,imax
#  ifndef OBC_SOUTH
            LapT(i,jstr-1)=0.
#  else
            LapT(i,jstr-1)=LapT(i,jstr)
#  endif
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=imin,imax
#  ifndef OBC_NORTH
            LapT(i,jend+1)=0.
#  else
            LapT(i,jend+1)=LapT(i,jend)
#  endif
          enddo
        endif
# endif /* !NS_PERIODIC */
!
!  Compute FX=d(LapT)/d(xi) and FE=d(LapT)/d(eta) terms
!  After that cmpute and add in biharmonic mixing [T m^3/s].
!  Multiply by mixing coefficient.
!
        do j=jstr,jend
          do i=istr,iend+1
              FX(i,j)=-0.5*diff3u(i,j)
     &                    *pmon_u(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))
     &                                *(LapT(i,j)-LapT(i-1,j)) 
# ifdef MASKING
     &                                        * umask(i,j)
# endif
# ifdef WET_DRY
     &                                        * umask_wet(i,j)
# endif
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend
              FE(i,j)=-0.5*diff3v(i,j)
     &                    *pnom_v(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))
     &                                *(LapT(i,j)-LapT(i,j-1))
# ifdef MASKING
     &                                         * vmask(i,j)
# endif
# ifdef WET_DRY
     &                                         * vmask_wet(i,j)
# endif
          enddo
        enddo
!
!  Add in horizontal diffusion of tracer [T m^3/s].
!
        do j=jstr,jend
          do i=istr,iend
            cff1=pm(i,j)*pn(i,j)
            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+dt*cff1
     &                 *(FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))
     &                                             /Hz(i,j,k)
          enddo
        enddo

# if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_PV
!
! Tracer diagnostics
!
        do j=jstr,jend
          do i=istr,iend
            cff1=pm(i,j)*pn(i,j)
            cff=Hz(i,j,k)/(pm(i,j)*pn(i,j))
            THmix(i,j,k,itrc)=THmix(i,j,k,itrc)*cff
     &                       +FX(i+1,j)-FX(i,j)
     &                       +FE(i,j+1)-FE(i,j)
#  ifdef MASKING
     &                                        * rmask(i,j)
#  endif
            Trate(i,j,k,itrc)=(Hz(i,j,k)*t(i,j,k,nnew,itrc)
     &                    -Hz_bak(i,j,k)*t(i,j,k,nstp,itrc))
     &                                            /(dt*cff1)
#  ifdef MASKING
     &                                        * rmask(i,j)
#  endif
!
! Divide all diagnostic terms by the cell volume
! (Hz(i,j,k,itrc)/(pm(i,j).*pn(i,j)). There after the unit
! of diag terms are: (unit of tracers)* s-1.
!

            cff=cff1/Hz(i,j,k)

            THmix(i,j,k,itrc)=THmix(i,j,k,itrc)*cff
            Trate(i,j,k,itrc)=Trate(i,j,k,itrc)*cff
            
#  if defined DIAGNOSTICS_TSVAR && !defined SPONGE_DIF2

#   include "finalize_diagnostics_tsadv.h"

#  endif
          enddo
        enddo
# endif /* DIAGNOSTICS_TS */

#endif /* TS_DIF4 */
      enddo   ! --> k
!---------------------------------------------------------------------
!
#ifdef DIAGNOSTICS_TS
# if defined DIAGNOSTICS_TS_MLD && !defined SPONGE_DIF2
!==========================================================
! Tracer diagnostics averaged over the MLD
!==========================================================
!
#  define T_mld_nnew FX
#  define T_mld_nstp FE
#   ifdef DIAGNOSTICS_TS_MLD_CRIT


        if (itrc.eq.iCRT4) then
          rho0i = 1.0/rho0
          zN2_c = g * mld_crit_D * rho0i
          do j=Jstr,Jend
            do i=Istr,Iend
              kbl(i,j)  = N                  ! Initialization to the number of w ocean point
              hmlp(i,j)  = 0.0                  ! here hmlp used as a dummy variable, integrating vertically N^2
            enddo
          enddo
        endif

        do j=Jstr,Jend
          do i=Istr,Iend
            indx_depth = N
            nearest_diff = abs(z_r(i,j,N)+mld_depth_ref)
! Find the index corresponding to the reference depth mld_depth_ref
            do k=N-1,1,-1
              if (abs(z_r(i,j,k) + mld_depth_ref) < nearest_diff) then
                nearest_diff = abs(z_r(i,j,k) + mld_depth_ref)
                indx_depth = k
              endif
            enddo

            if (itrc.eq.iCRT4) then
! for criterion based on brunt-vaisala frequency !
              do k=indx_depth,1,-1
                if (k.eq.N) then
                  hmlp(i,j) = hmlp(i,j) +  max(bvf(i,j,k),0.0) *
     &                        (2.0* (z_w(i,j,k+1) - z_r(i,j,k)))
                else
                  hmlp(i,j) = hmlp(i,j) + max(bvf(i,j,k),0.0) *
     &                        (z_r(i,j,k+1) - z_r(i,j,k))
                endif
                if (hmlp(i,j) < zN2_c) then
                  kbl(i,j) = k
                endif
              enddo


            else
! for criteria based on density and temperature !
              kbl(i,j)=N
              cff=0
              do k=N,1,-1
                if (itrc.eq.iCRT2) then
                  if (cff.eq.0.and.t(i,j,k,nnew,itemp).gt.
     &            t(i,j,indx_depth,nnew,itemp)-mld_crit_T) then
                    kbl(i,j)=k
                  endif
                  if (t(i,j,k,nnew,itemp).le.
     &            t(i,j,indx_depth,nnew,itemp)-mld_crit_T)  cff=1

                elseif (itrc.eq.iCRT3) then
                  if (cff.eq.0.and.t(i,j,k,nnew,itemp).gt.
     &            t(i,j,indx_depth,nnew,itemp)-mld_crit_T2) then
                    kbl(i,j)=k
                  endif
                  if (t(i,j,k,nnew,itemp).le.
     &            t(i,j,indx_depth,nnew,itemp)-mld_crit_T2)  cff=1
            
                else
                  if (cff.eq.0.and.rho(i,j,k).lt.
     &            rho(i,j,indx_depth)+mld_crit_D) then
                    kbl(i,j)=k
                  endif
                  if (rho(i,j,k).ge.
     &            rho(i,j,indx_depth)+mld_crit_D)  cff=1
                endif
              enddo
            endif
           
! CALCULATE MIXED LAYER DEPTH            
            Tcrit_mld(i,j,itrc)=z_w(i,j,N)-z_r(i,j,kbl(i,j))
#    ifdef MASKING
     &                        * rmask(i,j)
#    endif
          enddo
        enddo

#    endif /*  DIAGNOSTICS_TS_MLD_CRIT    */

       do j=Jstr,Jend
         do i=Istr,Iend
           THmix_mld(i,j,itrc)=0.
           T_mld_nnew(i,j)=0.
           T_mld_nstp(i,j)=0.
         enddo
       enddo
       do j=Jstr,Jend
         do i=Istr,Iend
#  ifdef LMD_SKPP
           kmld=kbl(i,j)
#  else
           kmld=N-5
#  endif
           do k=N,kmld,-1
             cff=Hz(i,j,k)/(z_w(i,j,N)-z_w(i,j,kmld-1))
#   ifdef DIAGNOSTICS_TS_MLD_CRIT
             if (itrc.eq.iCRT2 .or. itrc.eq.iCRT3 
     &       .or. itrc.eq.iCRT4) then 
               THmix_mld(i,j,itrc)=THmix_mld(i,j,itrc)+
     &                           THmix(i,j,k,itemp)*cff
               T_mld_nnew(i,j)=T_mld_nnew(i,j)+
     &                       t(i,j,k,nnew,itemp)*cff

             else
# endif
               THmix_mld(i,j,itrc)=THmix_mld(i,j,itrc)+
     &                           THmix(i,j,k,itrc)*cff
               T_mld_nnew(i,j)=T_mld_nnew(i,j)+
     &                       t(i,j,k,nnew,itrc)*cff
#   ifdef DIAGNOSTICS_TS_MLD_CRIT
             endif
# endif
           enddo
         enddo
       enddo

       do j=Jstr,Jend
         do i=Istr,Iend
#  if defined LMD_SKPP || defined GLS_MIXING
           if (kbl_nstp(i,j).eq.0) kbl_nstp(i,j)=kbl(i,j)
           kmld=kbl_nstp(i,j)
#  else
           kmld=N-5
#  endif
           do k=N,kmld,-1
             cff=Hz_bak(i,j,k)/(z_w(i,j,N)-z_w(i,j,kmld-1))
#   ifdef DIAGNOSTICS_TS_MLD_CRIT
             if (itrc.eq.iCRT2 .or. itrc.eq.iCRT3 
     &       .or. itrc.eq.iCRT4) then 
               T_mld_nstp(i,j)=T_mld_nstp(i,j)+
     &                       t(i,j,k,nstp,itemp)*cff
             else
# endif
               T_mld_nstp(i,j)=T_mld_nstp(i,j)+
     &                       t(i,j,k,nstp,itrc)*cff
#   ifdef DIAGNOSTICS_TS_MLD_CRIT
             endif
# endif
           enddo
           if (itrc .eq. NT) kbl_nstp(i,j)=kbl(i,j)
         enddo
       enddo
       do j=Jstr,Jend
         do i=Istr,Iend
            Trate_mld(i,j,itrc)=(T_mld_nnew(i,j)-T_mld_nstp(i,j))/dt
            Tentr_mld(i,j,itrc)=Trate_mld(i,j,itrc)-
     &                          TXadv_mld(i,j,itrc)-
     &                          TYadv_mld(i,j,itrc)-
     &                          TVadv_mld(i,j,itrc)-
     &                          TVmix_mld(i,j,itrc)-
     &                          THmix_mld(i,j,itrc)-
     &                          TForc_mld(i,j,itrc)
         enddo
       enddo
#  undef T_mld_nnew 
#  undef T_mld_nstp
# endif /* DIAGNOSTICS_TS_MLD && !SPONGE_DIF2 */

#endif /* DIAGNOSTICS_TS */

!$acc end kernels

!------------------------------------------------------------------
#if defined BD_EXCHANGE
# ifdef THREE_GHOST_POINTS_TS
        call exchange_r3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                               t(START_2D_ARRAY,1,nnew,itrc))
# else
!#   ifndef BASIN        
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          t(START_2D_ARRAY,1,nnew,itrc))
!#    endif        
# endif
#endif
!
      return
      end

#ifndef CHILD_SPG
# define CHILD_SPG
# ifdef AGRIF
#  include "t3dmix_S.F"
# endif
# undef CHILD_SPG
#endif  /* !CHILD_SPG */
