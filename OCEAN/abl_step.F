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
      subroutine abl_step (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call abl_step_tile (Istr,Iend,Jstr,Jend,
     &                               T2d(1,1,trd),T2d(1,2,trd),
     &                               T2d(1,3,trd),T2d(1,4,trd) )
      return
      end

      subroutine abl_step_tile (Istr,Iend,Jstr,Jend,
     &                                  elem_a,elem_b,elem_c,CF)
      !
      !======================================================================
      !                 ***  SUBROUTINE  bulk_flux  ***
      !  This routine computes the turbulent and radiative components of
      !  air-sea fluxes for the specification of surface boundary conditions
      !=======================================================================
      !
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "forces.h"
# include "scalars.h"
# include "params_bulk.h"
# ifdef ABL1D
#  include "abl1d.h"
# endif
      REAL           :: elem_a(PRIVATE_1D_SCRATCH_ARRAY,N_ABL)
      REAL           :: elem_b(PRIVATE_1D_SCRATCH_ARRAY,N_ABL)
      REAL           :: elem_c(PRIVATE_1D_SCRATCH_ARRAY,N_ABL)
      REAL           :: CF    (PRIVATE_1D_SCRATCH_ARRAY,N_ABL)
!++ Local integers
      INTEGER        :: i, j, k
      INTEGER        :: Istr, Iend, Jstr, Jend
      INTEGER        :: imin, imax, jmin, jmax
      INTEGER        :: kt
!++ Local scalars
      REAL           :: cff  , cff1 , cff2, sig
      REAL           :: rho0i, cpi  , Hlv , msk
      REAL           :: TseaC, TseaK, Qsea
      REAL           :: TairC, TairK, Qair , rhoAir
      REAL           :: WstarTstar, WstarQstar, evap
      REAL           :: hfsen, hflat, hflw , upvel
      REAL           :: uwnd_r,vwnd_r
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
! 1 - Advance tracers to time n+1
!======================================================================
       do j=jmin,jmax

         do k=3,N_abl-1
           do i=imin,imax
             elem_a(i,k) = - dt * avt_abl(i,j,k-1,nt_a) / Hzw_abl(k-1)
             elem_c(i,k) = - dt * avt_abl(i,j,k  ,nt_a) / Hzw_abl(k  )
             elem_b(i,k) = Hzr_abl(k) - elem_a(i,k) - elem_c(i,k)
           enddo
         enddo
         ! boundary conditions
         do i=imin,imax
           elem_a(i,2    ) = 0.
           elem_c(i,2    ) = -dt*avt_abl(i,j,    2,nt_a)/Hzw_abl(2    )
           elem_a(i,N_abl) = -dt*avt_abl(i,j,N_abl,nt_a)/Hzw_abl(N_abl)
           elem_b(i,N_abl) = Hzr_abl(N_abl) - elem_a(i,N_abl)
           elem_c(i,N_abl) = 0.
         enddo

         do kt=1,jptq
           do k=3,N_abl-1
             do i=imin,imax
               tq_abl(i,j,k,nt_a,kt)=Hzr_abl(k)*tq_abl(i,j,k,nt_n,kt)
             enddo
           enddo

           if(kt==jp_ta) then
             do i=imin,imax
               cff1 = Ch_du(i,j)
               cff2 = Ch_du(i,j)*sst_abl(i,j)
               elem_b(i,2) = Hzr_abl(2) - elem_c(i,2) + dt * cff1
               tq_abl(i,j,2    ,nt_a,kt)=
     &         Hzr_abl(2    )*tq_abl(i,j,    2,nt_n,kt) + dt * cff2
               tq_abl(i,j,N_abl,nt_a,kt)=
     &            Hzr_abl(N_abl)*tq_abl(i,j,N_abl,nt_n,kt)
             enddo
           else
             do i=imin,imax
               cff1 = Ce_du(i,j)
               cff2 = Ce_du(i,j)*ssq_abl(i,j)
               elem_b(i,2) = Hzr_abl(2) - elem_c(i,2) + dt * cff1
               tq_abl(i,j,2    ,nt_a,kt)=
     &            Hzr_abl(2    ) * tq_abl(i,j,    2,nt_n,kt) + dt * cff2
               tq_abl(i,j,N_abl,nt_a,kt)=
     &            Hzr_abl(N_abl)*tq_abl(i,j,N_abl,nt_n,kt)
             enddo
           endif
           !! Matrix inversion
           !! ----------------------------------------------------------
           do i=imin,imax
             cff                       =  1. / elem_b(i,2)
             CF    (i,   2           ) = - cff * elem_c(i,2)
             tq_abl(i,j,2,nt_a,kt) =    cff * tq_abl(i,j,2,nt_a,kt)
           enddo

           do k=3,N_abl
             do i=imin,imax
               cff = 1. / (elem_b(i,k) + elem_a(i,k)*CF(i,k-1))
               CF(i,k) = - cff * elem_c(i,k)
               tq_abl(i,j,k,nt_a,kt)=  cff * ( tq_abl(i,j,k  ,nt_a,kt)
     &                         - elem_a(i,k) * tq_abl(i,j,k-1,nt_a,kt) )
             enddo
           enddo

           do k=N_abl-1,2,-1
             do i=imin,imax
               tq_abl(i,j,k,nt_a,kt) = tq_abl(i,j,k  ,nt_a,kt)
     &                     + CF(i,k) * tq_abl(i,j,k+1,nt_a,kt)
             enddo
           enddo
         enddo   !<-- loop on tracers
       enddo   ! end outer loop
!======================================================================
! 2 - Compute Coriolis term with geostrophic guide
!======================================================================
       if( MOD( iic, 2 ) == 0 ) then
         ! Advance u_abl & v_abl to time n+1
         do k=2,N_abl    ! outer loop
           do j=jmin,jmax
             do i=imin,imax
               cff = f(i,j) * ( v_abl (i,j,k,nt_n) - pgv_dta(i,j,k) )
               u_abl(i,j,k,nt_a) = u_abl(i,j,k,nt_n) + dt * cff
               cff = f(i,j) * ( u_abl (i,j,k,nt_a) - pgu_dta(i,j,k)  )
               v_abl(i,j,k,nt_a) = Hzr_abl(k)*( v_abl(i,j,k,nt_n)
     &                                             - dt * cff )
               u_abl(i,j,k,nt_a) = Hzr_abl(k)*  u_abl(i,j,k,nt_a)
             enddo
           enddo
         enddo ! end outer loop
       else
           ! Advance v_abl & u_abl to time n+1
         do k=2,N_abl    ! outer loop
           do j=jmin,jmax
             do i=imin,imax
               cff = f(i,j) * ( u_abl(i,j,k,nt_n) - pgu_dta(i,j,k) )
               v_abl(i,j,k,nt_a) = v_abl(i,j,k,nt_n) - dt * cff
               cff = f(i,j) * ( v_abl(i,j,k,nt_a) - pgv_dta(i,j,k)  )
               u_abl(i,j,k,nt_a) = Hzr_abl(k) *( u_abl(i,j,k,nt_n)
     &                                                + dt * cff )
               v_abl(i,j,k,nt_a) = Hzr_abl(k) *  v_abl(i,j,k,nt_a)
             enddo
           enddo
         enddo ! end outer loop
       endif
!======================================================================
! 3 - Advance u to time n+1
!======================================================================
      do j=jmin,jmax
        !
        do k=3,N_abl
          do i=imin,imax
            elem_a(i,k) = -dt * Avm_abl(i,j,k-1,nt_a) / Hzw_abl(k-1)  ! lower-diagonal
            elem_c(i,k) = -dt * Avm_abl(i,j,k  ,nt_a) / Hzw_abl(k  )  ! upper-diagonal
            elem_b(i,k) = Hzr_abl(k) - elem_a(i,k) - elem_c(i,k) !       diagonal
          enddo
        enddo
        !
        do i=imin,imax    ! boundary conditions
          elem_a(i,2) = 0.
          elem_c(i,2) = - dt * Avm_abl(i,j,2,nt_a) / Hzw_abl(2)
           !
          cff1  = cd_du(i,j)
# ifdef KILPATRICK
          cff2  = 0.
# else
          cff2  = Cd_du(i,j) * 0.5*( u(i+1,j,N,nrhs)+u(i,j,N,nrhs) )
# endif
          elem_b(i,2) = Hzr_abl( 2 ) - elem_c(i,2) + dt * cff1
          u_abl(i,j,2,nt_a) = u_abl(i,j,2,nt_a)    + dt * cff2
          !++ Top Dirichlet B.C.
          elem_a(i,     N_abl )       = 0.
          elem_c(i,     N_abl )       = 0.
          elem_b(i,     N_abl )       = Hzr_abl( N_abl )
          u_abl(i,j,N_abl,nt_a)= Hzr_abl( N_abl )*pu_dta(i,j,N_abl)
        enddo
        !! Matrix inversion
        !! ----------------------------------------------------------
        do i=imin,imax
          cff                =   1. / elem_b(i,2)
          CF    (i,   2    ) =  - cff * elem_c(i,2)
          u_abl (i,j,2,nt_a) =    cff * u_abl(i,j,2,nt_a)
        enddo

        do k=3,N_abl
          do i=imin,imax
             cff = 1. / ( elem_b(i,k) + elem_a(i,k)*CF(i,k-1) )
             CF(i,k) = - cff * elem_c(i,k)
             u_abl(i,j,k,nt_a) = cff * ( u_abl(i,j,k,nt_a)
     &            - elem_a(i,k) * u_abl(i,j,k-1,nt_a) )
          enddo
        enddo

        do k=N_abl-1,2,-1
          do i=imin,imax
             u_abl(i,j,k,nt_a) = u_abl(i,j,k  ,nt_a)
     &                 + CF(i,k)*u_abl(i,j,k+1,nt_a)
          enddo
        enddo
         !-------------
      enddo             ! end outer loop
!======================================================================
! 4 - Advance v to time n+1
!======================================================================
      do j=jmin,jmax   ! outer loop
        !
        do k=3,N_abl-1
          do i=imin,imax
            elem_a(i,k) = -dt * Avm_abl(i,j,k-1,nt_a) / Hzw_abl(k-1)   ! lower-diagonal
            elem_c(i,k) = -dt * Avm_abl(i,j,k  ,nt_a) / Hzw_abl(k  )   ! upper-diagonal
            elem_b(i,k) = Hzr_abl(k) - elem_a(i,k) - elem_c(i,k)  !       diagonal
          enddo
        enddo

        do i=imin,imax
          !++ Surface boundary condition
          elem_a(i,2) = 0.
          elem_c(i,2) = - dt * Avm_abl(i,j,2,nt_a) / Hzw_abl(2)
          cff1        = cd_du(i,j)
# ifdef KILPATRICK
          cff2        = 0.
# else
          cff2        = Cd_du(i, j) * 0.5*( v(i,j+1,N,nrhs)+v(i,j,N,nrhs) )
# endif
          elem_b(i,2) = Hzr_abl(2) - elem_c(i,2) + dt * cff1
          v_abl(i,j,2,nt_a) =  v_abl(i,j,2,nt_a) + dt * cff2
          !++ Top Dirichlet B.C.
          elem_a(i,N_abl) = 0.
          elem_c(i,N_abl) = 0.
          elem_b(i,N_abl) = Hzr_abl( N_abl )
          v_abl(i,j,N_abl,nt_a)=Hzr_abl(N_abl)*pv_dta(i,j,N_abl)
        enddo
                !! Matrix inversion
                !! ----------------------------------------------------------
        do i=imin,imax
          cff               =  1. / elem_b(i,2)
          CF   (i,   2    ) =   - cff * elem_c(i,  2     )
          v_abl(i,j,2,nt_a) =     cff * v_abl (i,j,2,nt_a)
        enddo

        do k=3,N_abl
          do i=imin,imax
            cff = 1. / ( elem_b(i,k) + elem_a(i,k) * CF(i,k-1) )
            CF(i,k) = - cff * elem_c(i,k)
            v_abl(i,j,k,nt_a) = cff * ( v_abl(i,j,k,nt_a)
     &                - elem_a(i,k) * v_abl(i,j,k-1,nt_a) )
          enddo
        enddo

        do k=N_abl-1,2,-1
          do i=imin,imax
            v_abl(i,j,k,nt_a) = v_abl(i,j,k  ,nt_a)
     &               + CF(i,k) * v_abl(i,j,k+1,nt_a)
          enddo
        enddo
         !-------------
      enddo             ! end outer loop
!======================================================================
! 5 - Apply nudging on the dynamics and the tracers
!======================================================================
# ifdef ABL_NUDGING
!-------------------------
        do j=jmin,jmax   ! outer loop
          do i=imin,imax
!-------------------------
            cff = 1./MAX( pblh_min,  MIN(  pblh_max, ablh(i,j,nt_a)  ) )
#  ifdef MASKING
            msk = rmask(i,j)
#  else
            msk = 1.
#  endif
!-------------------------
            do k=2,N_abl
!-------------------------
              sig  = MIN(bmax,MAX(cff*zr_abl(k),bmin))
#  ifdef ABL_NUDGING_DYN
              cff1 = alp3_dyn*sig**3 + alp2_dyn*sig**2
     &             + alp1_dyn*sig    + alp0_dyn
              cff2 = (1.-msk) + msk*dt*cff1   ! cff2 = 1 for masked points
#   ifdef ABL_DYN_RESTORE_EQ
              cff2 = cff2 * rest_eq(i,j)
#   endif
              u_abl(i,j,k,nt_a)=(1.-cff2)* u_abl(i,j,k,nt_a)
     &                              + cff2 *pu_dta(i,j,k     )
              v_abl(i,j,k,nt_a)=(1.-cff2)* v_abl(i,j,k,nt_a)
     &                               +cff2 *pv_dta(i,j,k     )
#  endif
#  ifdef ABL_NUDGING_TRA
              cff1 = alp3_tra*sig**3 + alp2_tra*sig**2
     &             + alp1_tra*sig    + alp0_tra
              cff2 = (1.-msk) + msk*dt*cff1   ! cff2 = 1 for masked points
              tq_abl(i,j,k,nt_a,jp_ta)=(1.-cff2)*tq_abl(i,j,k,nt_a,jp_ta)
     &                                    +cff2 *(pt_dta(i,j,k)+273.16)
              tq_abl(i,j,k,nt_a,jp_qa)=(1.-cff2)*tq_abl(i,j,k,nt_a,jp_qa)
     &                                    +cff2 *pq_dta(i,j,k           )
#  endif
            enddo
          enddo
!-------------
        enddo             ! end outer loop
!-------------
# endif



!===============================================================================
!   @@@@@@@@   @@@@@@@@@@      @@@@@@@@@@ @@          @@        @@  @@     @@
!  @@      @@  @@      @@      @@         @@          @@        @@    @@  @@
!  @@      @@  @@@@@@@@@@  &&  @@@@@@@    @@          @@        @@      @@
!  @@      @@  @@      @@      @@         @@          @@        @@    @@  @@
!   @@@@@@@@   @@      @@      @@         @@@@@@@@@@   @@@@@@@@@@   @@      @@
!===============================================================================



!======================================================================
! Initialization of various constants
!======================================================================
!++  Inverse seawater density
      rho0i=1.0/rho0
!++  Inverse of specific heat for seawater (kg-degC/Joule) [cp is defined in scalars.h]
      cpi=1.0/cp
!----------------------
      DO j=jmin,jmax
        DO i=imin,imax
!----------------------
!++  Compute turbulent sensible heat flux (W/m2).
          TairK      = tq_abl(i,j,2,nt_a,jp_ta)
          TseaK      = sst_abl(i,j); TseaC = TseaK - CtoK
          Qair       = tq_abl(i,j,2,nt_a,jp_qa)
          Qsea       = ssq_abl(i,j)
          WstarTstar = Ch_du(i,j)*(TairK-TseaK)
          WstarQstar = Ce_du(i,j)*(Qair -Qsea )
          rhoAir     = rho_abl(i,j)
!
          hfsen = - blk_Cpa*rhoAir*WstarTstar
!++ Latent heat of vaporisation (J/kg)
          Hlv   = (2.5008 - 0.0023719*TseaC)*1.0E+6
!++  Compute turbulent latent heat flux (W/m2).
          hflat = - Hlv*rhoAir*WstarQstar
!++ Longwave Radiation.
# ifndef BULK_LW
          hflw   = -radlw(i,j) ! positive downward, i.e., warming the ocean
# else
!++  Treat input longwave data as downwelling radiation only and add
!++  outgoing IR from model sea surface temperature.
          hflw   = radlw(i,j)  ! positive downward
     &           - emiss_lw*rho0i*cpi*SigmaSB*TseaK*TseaK*TseaK*TseaK
# endif
!++  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
          upvel=-1.61*WstarQstar-(1.0+1.61*Qair)*WstarTstar/TairK
!++  Compute turbulent latent heat flux (W/m2).
          hflat=hflat+rhoAir*Hlv*upvel*Qair
!======================================================================
!  Latent and sensible heat fluxes
!======================================================================
          hflat=-hflat*rho0i*cpi
          hfsen=-hfsen*rho0i*cpi
!======================================================================
!  Compute kinematic, surface wind stress (m2/s2).
!======================================================================
          cff        = rhoAir*rho0i
          uwnd_r     = u_abl(i,j,2,nt_a)
# ifndef KILPATRICK
     &               - 0.5*( u(i+1,j,N,nrhs)+u(i,j,N,nrhs) )
# endif
          vwnd_r     = v_abl(i,j,2,nt_a)
# ifndef KILPATRICK
     &               - 0.5*( v(i,j+1,N,nrhs)+v(i,j,N,nrhs) )
# endif
          sustr(i,j) = cff*Cd_du(i,j) * uwnd_r
          svstr(i,j) = cff*Cd_du(i,j) * vwnd_r
# ifdef MASKING
          sustr(i,j)=sustr(i,j)*rmask(i,j)
          svstr(i,j)=svstr(i,j)*rmask(i,j)
# endif
!++  Total surface heat flux
          stflx(i,j,itemp)=srflx(i,j)+hflw+hflat+hfsen
# ifdef SALINITY
!++  Total surface salt flux
          evap=-cp*hflat/Hlv
          stflx(i,j,isalt)=(evap-prate(i,j))*t(i,j,N,nrhs,isalt)
# endif
!++  Store fluxes for diagnostics, when storing in netCDF file, multiplied by rho0.Cp to get W/m2
          shflx_rsw(i,j)=srflx(i,j)
          shflx_lat(i,j)=hflat
          shflx_sen(i,j)=hfsen
          shflx_rlw(i,j)=hflw
        ENDDO !<-- terminate i-loop
      ENDDO   !<-- terminate j-loop
!----------------------
      nt_n = 1 + MOD( iic  , 2)
      nt_a = 1 + MOD( iic+1, 2)
!----------------------
      return
      end
#else
      subroutine abl_step_empty
      return
      end
#endif /* !ABL1D */
