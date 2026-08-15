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
#ifdef BULK_FLUX
!======================================================================
      subroutine bulk_flux (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call bulk_flux_tile (Istr,Iend,Jstr,Jend,
     &                               A2d(1,1,trd),A2d(1,2,trd))
      return
      end

      subroutine bulk_flux_tile (Istr,Iend,Jstr,Jend, aer,cer)
# ifdef BULK_FLUX_MLP
      use bulk_flux_mlp, only: bulk_flux_mlp_tile, mlp_init, mlp_tick
# endif
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
!++ Local integers
      INTEGER              :: i, j, m, iter
      INTEGER              :: Istr, Iend, Jstr, Jend
      INTEGER              :: imin, imax, jmin, jmax
      INTEGER, PARAMETER   :: IterFl = 3
# ifdef CFB_WIND_TRA
      INTEGER, PARAMETER :: mb = 2
# else
      INTEGER, PARAMETER :: mb = 1
# endif

!++ Local scalars
      REAL           :: cff  , cff2
# ifdef BULK_FLUX_MLP
!  Tile-sized arrays for MLP: allocated after tile bounds are known
      DOUBLE PRECISION, ALLOCATABLE :: mlp_Cd(:,:), mlp_Ch(:,:), mlp_Ce(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: mlp_wspd(:,:), mlp_delT(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: mlp_delQ(:,:), mlp_TairC(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: mlp_RH(:,:), mlp_patm(:,:)
# endif
      REAL           :: rho0i, cpi  , RH   , Hlv     , wspd0
      REAL           :: TseaC, TseaK, Qsea , spec_hum, Tvstar
      REAL           :: TairC, TairK, Q    , rhoAir  , Bf
      REAL           :: delW(2) , delT , delQ , ddelW(2) , ddelT , ddelQ
      REAL           :: Cd   , Ch   , Ce   , qsat    , iexns , iexna
      REAL           :: Wstar(2), Tstar, Qstar, Wstar0  , Tstar0, Qstar0
      REAL           :: psi_u(2), ZoLu(2) , psi_t, ZoLt    , patm
      REAL           :: hfsen, hflat, hflw , upvel   , evap
!$acc declare create( cff  , cff2
!$acc&               ,rho0i, cpi  , RH   , Hlv     , wspd0
!$acc&               ,TseaC, TseaK, Qsea , Tvstar
!$acc&               ,TairC, TairK, Q    , rhoAir  , Bf
!$acc&               ,delW , delT , delQ , ddelW , ddelT , ddelQ
!$acc&               ,Cd   , Ch   , Ce   , iexns , iexna
!$acc&               ,Wstar, Tstar, Qstar, Wstar0  , Tstar0, Qstar0
!$acc&               ,psi_u, ZoLu, psi_t, ZoLt    , patm
!$acc&               ,hfsen, hflat, hflw , upvel   , evap
!$acc&               )
# ifdef CFB_WIND_TRA
      REAL           :: wspd0_cfb
# endif
# ifdef BULK_ECUMEV0
      REAL            :: Cdn(2), Chn, Cen
# endif
# ifdef BULK_ECUMEV6
      REAL            :: Cdn(2), partn, parqn
# endif
# ifdef BULK_WASP
      REAL            :: chi,aa,bb
      REAL            :: charn1, charn2, charn3, charn4, charn5
# endif
# ifdef ABL1D
      REAL            :: chosqrcd
      REAL            :: uwnd_r,vwnd_r
# endif
!++ Parameter values (blk_ZW and blk_ZT should not be parameters with ABL)
      REAL            :: blk_ZW, blk_ZT, blk_ZToZW
# ifndef ABL1D
      PARAMETER ( blk_ZW    =  10.0 )
      PARAMETER ( blk_ZT    =   2.0 )
      PARAMETER ( blk_ZToZW = blk_ZT/blk_ZW )
# endif
      REAL            ::  logus10,logts10
# if defined BULK_ECUMEV0 || defined BULK_ECUMEV6
#  ifndef ABL1D
      PARAMETER ( logus10   = LOG(blk_ZW/10.0) )
      PARAMETER ( logts10   = LOG(blk_ZT/10.0) )
#  endif
      REAL           :: bulk_psiu_ecume, bulk_psit_ecume
      REAL           :: deltau10n(2), deltat10n, deltaq10n
# else
      REAL           :: Log10oLogZw
#  ifndef ABL1D
      PARAMETER ( Log10oLogZw = LOG(  10.0*10000.)/
     &                          LOG(blk_ZW*10000.) )
#  endif
      REAL           :: bulk_psiu_coare, bulk_psit_coare, air_visc
      REAL           :: iZoW , iZoT  , CC, Ch10, Ribcu , charn
      REAL           :: iZo10, iZoT10, Ri, Rr  , VisAir
# endif
!++ Local arrays
      REAL           ::  aer(PRIVATE_2D_SCRATCH_ARRAY)
      REAL           ::  cer(PRIVATE_2D_SCRATCH_ARRAY)
# ifdef CFB_STRESS
      REAL           ::  stau(PRIVATE_2D_SCRATCH_ARRAY)
# endif
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
! Initialization of various constants
!======================================================================
# ifdef ABL1D
      blk_ZW    = zr_abl(2)
      blk_ZT    = zr_abl(2)
      blk_ZToZW = blk_ZT/blk_ZW
#  if defined BULK_ECUMEV0 || defined BULK_ECUMEV6
      logus10   = LOG(blk_ZW/10.0)
      logts10   = LOG(blk_ZT/10.0)
#  else
      Log10oLogZw = LOG(10.0*10000.)/LOG(blk_ZW*10000.)
#  endif
# endif
!++  Inverse seawater density
      rho0i=1.0/rho0
!++  Inverse of specific heat for seawater (kg-degC/Joule) [cp is defined in scalars.h]
      cpi=1.0/cp
!++
# ifndef READ_PATM
!++  Inverse of exner function at the surface
      iexns  = (psurf*ip00)**(-rdocpd)
!$acc update device( iexns, iexna, rho0i, hfsen, cpi )     
# endif

# ifdef BULK_FLUX_MLP
!======================================================================
!  BULK_FLUX_MLP — Pass 1: collect inputs across the full tile,
!  call the MLP ONCE for the whole tile, then use results in Pass 2.
!  This fixes the bug where bulk_flux_mlp_tile was called point-by-point
!  inside the i,j loop, causing step_counter to increment NX*NY times
!  per model timestep instead of once, which prevented mlp_active from
!  ever being set correctly and made both runs identical.
!======================================================================
      ALLOCATE( mlp_wspd (imin:imax,jmin:jmax),
     &          mlp_delT (imin:imax,jmin:jmax),
     &          mlp_delQ (imin:imax,jmin:jmax),
     &          mlp_TairC(imin:imax,jmin:jmax),
     &          mlp_RH   (imin:imax,jmin:jmax),
     &          mlp_patm (imin:imax,jmin:jmax),
     &          mlp_Cd   (imin:imax,jmin:jmax),
     &          mlp_Ch   (imin:imax,jmin:jmax),
     &          mlp_Ce   (imin:imax,jmin:jmax) )

!  --- Sub-loop: extract state variables for the full tile ---
      DO j=jmin,jmax
        DO i=imin,imax
#  ifdef READ_PATM
          psurf  = patm2d(i,j)
          iexns  = (psurf*ip00)**(-rdocpd)
#  endif
#  ifndef ABL1D
          wspd0  = wspd(i,j)
          wspd0  = MAX(wspd0, 0.1*MIN(10., blk_ZW))
          TairC  = tair(i,j)
          TairK  = TairC + CtoK
          RH     = rhum(i,j)
          Q      = spec_hum(RH,psurf,TairC)
#  else
          uwnd_r = u_abl(i,j,2,nt_n)
#   ifndef KILPATRICK
     &           - 0.5*(u(i+1,j,N,nrhs)+u(i,j,N,nrhs))
#   endif
          vwnd_r = v_abl(i,j,2,nt_n)
#   ifndef KILPATRICK
     &           - 0.5*(v(i,j+1,N,nrhs)+v(i,j,N,nrhs))
#   endif
          wspd0  = SQRT(uwnd_r*uwnd_r+vwnd_r*vwnd_r)
          TairK  = tq_abl(i,j,2,nt_n,jp_ta)
          TairC  = TairK - CtoK
          Q      = tq_abl(i,j,2,nt_n,jp_qa)
          wspd0  = MAX(wspd0, 0.1*MIN(10., blk_ZW))
#  endif
#  ifdef SST_SKIN
          TseaC = sst_skin(i,j)
#  else
          TseaC = t(i,j,N,nrhs,itemp)
#  endif
          TseaK = TseaC + CtoK
          CALL exner_patm_from_tairabs(iexna,patm,Q,TairK,blk_ZT,psurf)
          Qsea  = qsat(TseaK,psurf,0.98)
#  ifdef BULK_GUSTINESS
          delW(1) = SQRT(wspd0*wspd0+0.25)
#  else
          delW(1) = wspd0
#  endif
#  ifdef KILPATRICK
          delQ = 0.0
#  else
          delQ = Q - Qsea
#  endif
          cff  = CtoK*(iexna-iexns)
          delT = TairC*iexna - TseaC*iexns + cff
          ddelW(1) = delW(1)

          mlp_wspd (i,j) = DBLE(ddelW(1))
          mlp_delT (i,j) = DBLE(delT)
          mlp_delQ (i,j) = DBLE(delQ)
          mlp_TairC(i,j) = DBLE(TairC)
          mlp_RH   (i,j) = DBLE(RH)
          mlp_patm (i,j) = DBLE(patm)
        END DO
      END DO

!  --- Call MLP ONCE for the entire tile (fixes step_counter bug) ---
      CALL bulk_flux_mlp_tile(imin, imax, jmin, jmax,
     &     mlp_wspd,  mlp_delT,  mlp_delQ,
     &     mlp_TairC, mlp_RH,    mlp_patm,
     &     mlp_Cd,    mlp_Ch,    mlp_Ce)

!  --- DIAG: print centre-tile MLP output (once per tile call on rank 0) ---
#ifdef MPI
      if (mynode.eq.0) then
#else
      if (.true.) then
#endif
        write(*,'(A,3(2X,A,E12.5),2X,A,E12.5)')
     &    ' [BLK-DIAG] MLP output at tile centre:',
     &    'Cd=', mlp_Cd((imin+imax)/2,(jmin+jmax)/2),
     &    'Ch=', mlp_Ch((imin+imax)/2,(jmin+jmax)/2),
     &    'Ce=', mlp_Ce((imin+imax)/2,(jmin+jmax)/2),
     &    'wspd=', mlp_wspd((imin+imax)/2,(jmin+jmax)/2)
      endif

# endif /* BULK_FLUX_MLP */

!$acc kernels if(compute_on_device) default(present)
!----------------------
!$acc loop seq
      DO j=jmin,jmax
!$acc loop seq
        DO i=imin,imax
!----------------------
!======================================================================
! Initialization of bulk quantities (common to all bulk algorithms)
!======================================================================

# ifdef READ_PATM
!++  Get pressure from ATM file in Pa
          psurf  = patm2d(i,j)
!++  Inverse of exner function at the surface
          iexns  = (psurf*ip00)**(-rdocpd)
# endif

# ifndef ABL1D
          wspd0     = wspd(i,j)
#  ifdef CFB_WIND_TRA
          wspd0_cfb = wspd_cfb(i,j)
          wspd0_cfb = MAX ( wspd0_cfb , 0.1 * MIN(10., blk_ZW) )
#  endif
          TairC    = tair(i,j)
          TairK    = TairC + CtoK
          RH       = rhum(i,j)
          Q        = spec_hum (RH,psurf,TairC)
# else
          uwnd_r   = u_abl(i,j,2,nt_n)
#  ifndef KILPATRICK
     &             - 0.5*( u(i+1,j,N,nrhs)+u(i,j,N,nrhs) )
#  endif
          vwnd_r   = v_abl(i,j,2,nt_n)
#  ifndef KILPATRICK
     &             - 0.5*( v(i,j+1,N,nrhs)+v(i,j,N,nrhs) )
#  endif
          wspd0    = SQRT ( uwnd_r*uwnd_r+vwnd_r*vwnd_r )
          TairK    = tq_abl ( i,j,2,nt_n,jp_ta )
          TairC    = TairK - CtoK
          Q        = tq_abl ( i,j,2,nt_n,jp_qa )
# endif
          wspd0     = MAX ( wspd0 , 0.1 * MIN(10., blk_ZW) )

# ifdef SST_SKIN
          TseaC = sst_skin(i,j)
# else
          TseaC = t(i,j,N,nrhs,itemp)
# endif
          TseaK = TseaC + CtoK
!
          CALL exner_patm_from_tairabs(iexna,patm,Q,TairK,blk_ZT,psurf)
!++  Air density
          rhoAir = patm*(1.+Q) / ( blk_Rgas*TairK*(1.+MvoMa*Q) )
!++  Specific humidity at saturation
          Qsea   = qsat(TseaK,psurf,0.98)
# ifdef ABL1D
          sst_abl(i,j) = TseaK*iexns
#  ifdef KILPATRICK
          ssq_abl(i,j) = 0.0
#  else
          ssq_abl(i,j) = Qsea
#  endif
          rho_abl(i,j) = rhoAir
# endif
!++  Air-sea gradients delW, delT, delQ through the interface
# ifdef BULK_GUSTINESS
          delW(1) = SQRT(wspd0*wspd0+0.25)
#  ifdef CFB_WIND_TRA
          delW(2) = SQRT(wspd0_cfb*wspd0_cfb+0.25)
#  endif
# else
          delW(1) = wspd0
#  ifdef CFB_WIND_TRA
          delW(2) = wspd0_cfb
#  endif
# endif
# ifdef KILPATRICK
          delQ = 0.0
# else
          delQ = Q-Qsea
# endif
          cff  = CtoK*(iexna-iexns)
          delT = TairC*iexna - TseaC*iexns + cff

# if defined BULK_ECUMEV0 || defined BULK_ECUMEV6
!
!======================================================================
!   @@@@@@@@@@  @@@@@@@@@@  @@      @@  @@@    @@@  @@@@@@@@@@
!   @@          @@          @@      @@  @@ @  @ @@  @@
!   @@@@@@@     @@          @@      @@  @@  @@  @@  @@@@@@@
!   @@          @@          @@      @@  @@      @@  @@
!   @@@@@@@@@@  @@@@@@@@@@   @@@@@@@@   @@      @@  @@@@@@@@@@
!======================================================================
!
!++ Initial guess for ECUME algorithms
          ddelW(1:mb)     = SIGN(MAX(ABS(delW(1:mb)),10.0*dWstar0),
     &                           delW(1:mb))
          ddelT           = SIGN(MAX(ABS(delT),10.0*dTstar0),delT)
          ddelQ           = SIGN(MAX(ABS(delQ),10.0*dQstar0),delQ)
          deltau10n(1:mb) = ddelW(1:mb)
          deltat10n       = ddelT
          deltaq10n       = ddelQ
!======================================================================
!  ITERATIVE LOOP TO COMPUTE U*, T*, Q*.
!======================================================================
!----------------------
!$acc loop seq
          DO Iter=1,IterFl
!----------------------
#  ifdef BULK_ECUMEV0
!$acc loop seq
            DO m=1,mb
              cff    = deltau10n(m)*deltau10n(m)
!++  Neutral coefficient for wind speed cdn (ECUME_v0 formulation)
              IF ( deltau10n(m) <= utu1 ) THEN
                Cdn(m) = coefu10+coefu11*deltau10n(m) +
     &                   (coefu12+coefu13*deltau10n(m))*cff
              ELSEIF ( deltau10n(m) <= utu2 ) THEN
                Cdn(m) = coefu20+coefu21*deltau10n(m)
     &                 + (coefu22+coefu23*deltau10n(m))*cff
     &                 + coefu24 * cff * cff
              ELSE
                Cdn(m) = Cdn0
              ENDIF
            ENDDO
!++ Neutral coefficient for temperature chn (ECUME_v0 formulation)
            IF ( deltau10n(mb) <= utt ) THEN
              Chn =  coeft0+coeft1*deltau10n(mb)
     &            +  (coeft2+coeft3*deltau10n(mb))*cff
     &            +  (coeft4+coeft5*deltau10n(mb)) * cff * cff
            ELSE
              Chn = Chn0
            ENDIF
!++ Neutral coefficient for humidity cen (ECUME_v0 formulation)
            IF ( deltau10n(mb) <= utq1 ) THEN
              Cen = coefq10+coefq11*deltau10n(mb)
     &            + (coefq12+coefq13*deltau10n(mb))*cff
     &            + (coefq14                      )* cff* cff
            ELSEIF ( deltau10n(mb) <= utq2 ) then
              Cen = coefq20 + coefq21 * deltau10n(mb) +  ( coefq22 )*cff
            ELSE
              Cen = Cen0
            ENDIF
!++ Estimate Monin-Obukhov similarity scales (ECUME_v0 formulation)
            Wstar(1:mb) = SQRT(Cdn(1:mb))*deltau10n(1:mb)
            Tstar       = Chn/SQRT(Cdn(mb))*deltat10n
            Qstar       = Cen/SQRT(Cdn(mb))*deltaq10n
#  else
!$acc loop seq
            DO m=1,mb
              cff    = deltau10n(m)*deltau10n(m)
              cff2   = 1./deltau10n(m)
!++  Neutral parameter for wind speed (ECUME_V6 formulation)
              IF (deltau10n(m) <= utu) THEN
                Wstar(m) = (coefu0+coefu1*deltau10n(m))
     &                   + (coefu2+coefu3*deltau10n(m))*cff
     &                   + (coefu4+coefu5*deltau10n(m))*cff*cff
              ELSE
                Wstar(m) = cdiru*(deltau10n(m)-utu) + ordou
              ENDIF
            ENDDO
            Cdn(1:mb) = (Wstar(1:mb)/deltau10n(1:mb))**2.
!++  Neutral parameter for temperature (ECUME_V6 formulation)
            IF (deltau10n(mb) <= utt) then
              partn = (coeft0+coeft1*deltau10n(mb))
     &              + (coeft2+coeft3*deltau10n(mb))*cff
     &              +  coeft4*cff*cff
            ELSE
              partn = cdirt*(deltau10n(mb)-utt) + ordot
            ENDIF
            Tstar = partn*deltat10n*cff2
!++  Neutral parameter for humidity (ECUME_V6 formulation)
            IF (deltau10n(mb) <= utq) then
              parqn =   coefq0+coefq1*deltau10n(mb)+coefq2*cff
            ELSE
              parqn = cdirq*(deltau10n(mb)-utq) + ordoq
            ENDIF
            Qstar = parqn*deltaq10n*cff2
#  endif
! !++ Obukhovs stability param. z/l
!$acc loop seq
            DO m=1,mb
              ZoLu(m)  = blk_ZW*g*vonKar*(Tstar*(1.+cpvir*Q)
     &                   +cpvir*TairK*Qstar)/
     &                   (TairK*Wstar(m)*Wstar(m)*(1.+cpvir*Q)+eps)
              ZoLu(m)  = MAX(MIN(ZoLu(m),LMOmax),LMOmin)
              psi_u(m) = bulk_psiu_ecume ( ZoLu(m) )
            ENDDO
            ZoLt = ZoLu(mb)*blk_ZToZW
            ZoLt = MAX(MIN(ZoLt,LMOmax),LMOmin)
!++ Stability function psi (see Liu et al, 1979 ; Dyer and Hicks, 1970)
!++  Modified to include convective form following Fairall (unpublished)
            psi_t = bulk_psit_ecume ( ZoLt )
!
#  ifdef BULK_GUSTINESS
!++ Add gustiness contribution to wind speed in unstable conditions
            Tvstar   = Tstar*(1.0+cpvir*Q)+cpvir*TairK*Qstar
            Bf       = MAX(0.0,-g/TairK*Wstar(1)*Tvstar)
            cff2     = blk_beta*blk_beta*(Bf*blk_Zabl)**(2./3.)
            delW (1) = SQRT(wspd0*wspd0+cff2)
            ddelW(1) = SIGN(MAX(ABS(delW(1)),10.0*dWstar0),delW(1))
#   ifdef CFB_WIND_TRA
            Bf       = MAX(0.0,-g/TairK*Wstar(2)*Tvstar)
            cff2     = blk_beta*blk_beta*(Bf*blk_Zabl)**(2./3.)
            delW (2) = SQRT(wspd0_cfb*wspd0_cfb+cff2)
            ddelW(2) = SIGN(MAX(ABS(delW(2)),10.0*dWstar0),delW(2))
#   endif
#  endif
!++ Update 10 meter Neutral quantities
            cff          = 1./vonKar
            cff2         = ddelW(1)-Wstar(1)*(logus10-psi_u(1))*cff
            deltau10n(1) = SIGN(MAX(ABS(cff2),10.0*dWstar0),cff2)
#  ifdef CFB_WIND_TRA
            cff2         = ddelW(2)-Wstar(2)*(logus10-psi_u(2))*cff
            deltau10n(2) = SIGN(MAX(ABS(cff2),10.0*dWstar0),cff2)
#  endif
            cff2         = ddelT-Tstar*(logts10-psi_t)*cff
            deltat10n    = SIGN(MAX(ABS(cff2),10.0*dTstar0),cff2)
            cff2         = ddelQ-Qstar*(logts10-psi_t)*cff
            deltaq10n    = SIGN(MAX(ABS(cff2),10.0*dQstar0),cff2)
!----------------------
      ENDDO  !<-- terminate iterations-loop
!----------------------
#  ifdef ABL1D
            chosqrcd  = vonKar/(logts10-psi_t)
#  endif
# else
!
!======================================================================
!   @@@@@@@@@@   @@@@@@@@@@   @@@@@@@@@@  @@@@@@@@@@  @@@@@@@@@@
!   @@          @@        @@  @@      @@  @@      @@  @@
!   @@          @@        @@  @@@@@@@@@@  @@@@@@@@@@  @@@@@@@
!   @@          @@        @@  @@      @@  @@    @@    @@
!   @@@@@@@@@@   @@@@@@@@@@   @@      @@  @@      @@  @@@@@@@@@@
!
!                                     @@
!                                   @@
!                                 @@
!                               @@
!                             @@
!
!           @@      @@  @@@@@@@@@@  @@@@@@@@@@  @@@@@@@@@@
!           @@      @@  @@      @@  @@          @@      @@
!           @@  @@  @@  @@@@@@@@@@  @@@@@@@@@@  @@@@@@@@@@
!            @@ @@ @@   @@      @@          @@  @@
!             @@  @@    @@      @@  @@@@@@@@@@  @@
!======================================================================
!
# ifndef BULK_FLUX_MLP
!++ Initial guess for COARE algorithm
          Wstar(1:mb)  = 0.035*delW(1:mb)*Log10oLogZw   ! Friction velocity
          VisAir = air_visc ( TairC )       ! Molecular viscosity
          charn  = 0.011                    ! Charnock parameter
          Ch10   = 0.00115                                                ! 10m exchange coefficient for heat
          Ribcu  = -    blk_ZW / ( blk_Zabl * 0.004 * blk_beta**3 )
          
!$acc loop seq
          DO m=1,mb
            iZo10  = g*Wstar(m) /
     &               (charn*Wstar(m)*Wstar(m)*Wstar(m)+0.11*g*VisAir)    ! inverse of roughness length
            iZoT10 = 0.1 * exp(vonKar*vonKar /
     &                     ( Ch10*LOG( 10.0*iZo10 ) ) ) ! inverse of thermal roughness length
!++ Obukhovs stability param. z/l
            CC     = LOG( blk_ZW*iZo10 )*LOG( blk_ZW*iZo10 )/
     &               LOG( blk_ZT*iZoT10 )
#ifdef OPENACCNUMERIC
            CC = real(CC,kind=4)
#endif
     
            Ri     =  g * blk_ZW * ( delT+cpvir*TairK*delQ )/
     &                             ( TairK*delW(m)*delW(m) )
            IF ( Ri < 0.0 ) THEN
              ZoLu(m)=CC*Ri/(1.0+Ri/Ribcu)       ! Unstable
            ELSE
              ZoLu(m)=CC*Ri/(1.0+3.0*Ri/CC)      ! Stable
            ENDIF
            psi_u(m)     = bulk_psiu_coare ( ZoLu(m) )
            logus10   = LOG(blk_ZW* iZo10)
#ifdef OPENACCNUMERIC
            logus10 = real(logus10,kind=4)
#endif
            Wstar(m)  = delW(m)*vonKar/(logus10-psi_u(m))
          ENDDO
          ZoLt      = ZoLu(mb)*blk_ZToZW        ! rescale in case blk_ZW != blk_ZT
!++ Stability functions
          psi_t     = bulk_psit_coare ( ZoLt    )
!++ Initial guess for Monin-Obukhov similarity scales
          logts10   = LOG(blk_ZT*iZoT10)
#ifdef OPENACCNUMERIC
          logts10 = real( logts10, kind=4)  
#endif
          cff       =      vonKar/(logts10-psi_t)
          Tstar     = delT*cff
          Qstar     = delQ*cff
!++  Compute Charnock coefficient
#  ifdef BULK_WASP
          !
          !++ Compute wave age
          chi = cwage_wasp*(delW(1)/Wstar(1))
          !
          !++ Compute coefficient 1
          charn1 = coefu_wasp(1)*(delW(1)**coefu_wasp(2))
          !
          !++ Compute coefficient 2
          aa = coefa_wasp(1)
     &       + coefa_wasp(2)*delW(1)
     &       + coefa_wasp(3)*delW(1)**2
     &       + coefa_wasp(4)*delW(1)**3
          bb = coefb_wasp(1)
     &       + coefb_wasp(2)*delW(1)
     &       + coefb_wasp(3)*delW(1)**2
     &       + coefb_wasp(4)*delW(1)**3
          charn2 = aa*(chi**bb)
          !
          !++ Compute coefficient 3
          aa = coefa2_wasp(1)+coefa2_wasp(2)*delW(1)
          bb = coefb2_wasp(1)+coefb2_wasp(2)*delW(1)
          charn3 = aa*(chi**bb)
          IF (charn3 < limcharn_wasp) charn3 = limcharn_wasp
          !
          !++ Compute coefficient 4
          charn4 = polyu_wasp(1)
     &           + polyu_wasp(2)*delW(1)
     &           + polyu_wasp(3)*delW(1)**2
     &           + polyu_wasp(4)*delW(1)**3
          !
          !++ Compute coefficient 5
          charn5 = polyu2_wasp(1)
     &           + polyu2_wasp(2)*delW(1)
     &           + polyu2_wasp(3)*delW(1)**2
          !
          !++ Compute final charnock coefficient
          charn = (
     &          charn1*(0.5*(1.0+TANH((thr1_wasp-delW(1))/lon1_wasp)))
     &        + charn2*(0.5*(1.0+TANH((delW(1)-thr1_wasp)/lon1_wasp)))
     &            )
     &        * (       0.5*(1.0+TANH((thr2_wasp-delW(1))/lon2_wasp)))
     &        + charn3*(0.5*(1.0+TANH((delW(1)-thr2_wasp)/lon2_wasp)))
     &        * (       0.5*(1.0+TANH((thr3_wasp-delW(1))/lon3_wasp)))
     &        + charn4*(0.5*(1.0+TANH((delW(1)-thr3_wasp)/lon3_wasp)))
     &        * (       0.5*(1.0+TANH((thr4_wasp-delW(1))/lon4_wasp)))
     &        + charn5*(0.5*(1.0+TANH((delW(1)-thr4_wasp)/lon4_wasp)))
          !
!++  Limit Charnock parameter
          IF(charn > limcharn1_wasp) charn = limcharn1_wasp
#  else
          charn = 0.011
          IF     ( delW(1) > 18.0 ) then
            charn = 0.018
          ELSEIF ( delW(1) > 10.0 ) then
            charn = 0.011+0.125*(0.018-0.011)*(delW(1)-10.)
          ENDIF
#  endif
!
!======================================================================
!  ITERATIVE LOOP TO COMPUTE U*, T*, Q*.
!======================================================================
!----------------------
!$acc loop seq
          DO Iter=1,IterFl
!----------------------
!$acc loop seq
            DO m=1,mb
!++  Inverse of roughness length
              iZoW    = g*Wstar(m) / ( charn*Wstar(m)*Wstar(m)*Wstar(m)
     &                                 +0.11*g*VisAir )
!++  Inverse of thermal roughness length
#  ifdef BULK_WASP
              iZoT    = EXP(vonKar*vonKar/(Ch10*LOG(blk_ZW*iZoW)) )/
     &                      blk_ZW
#  else
              Rr      = Wstar(m)/(iZow*VisAir)
              iZoT    = MAX(8695.65,18181.8*(Rr**0.6))
#  endif
!++  Obukhovs stability param. z/l
              ZoLu(m)    = vonKar*g*blk_ZW*
     &             (Tstar*(1.0+cpvir*Q)+cpvir*TairK*Qstar)/
     &             (TairK*Wstar(m)*Wstar(m)*(1.0+cpvir*Q)+eps)
!++ Stability functions
              psi_u(m)     = bulk_psiu_coare ( ZoLu(m) )
              logus10 = LOG(blk_ZW*iZoW)
#ifdef OPENACCNUMERIC
              logus10 = real(logus10, kind=4)
#endif
              Wstar(m)   = delW(m)*vonKar/(logus10-psi_u(m))
            ENDDO

            ZoLt      = ZoLu(mb)*blk_ZToZW
            psi_t     = bulk_psit_coare ( ZoLt )

!++  Compute Monin-Obukhov similarity scales

            logts10   = LOG(blk_ZT*iZoT)
            cff       = vonKar/(logts10-psi_t)
            Tstar     = delT*cff
            Qstar     = delQ*cff
!++  Add BULK_GUSTINESS
#  ifdef BULK_GUSTINESS
              Bf=-g/TairK*Wstar(1)*(Tstar+cpvir*TairK*Qstar)
              if (Bf.gt.0.0) then
                cff=blk_beta*(Bf*blk_Zabl)**r3
              else
                cff=0.2
              endif
              delW(1)  = SQRT(wspd0*wspd0+cff*cff)
#ifdef OPENACCNUMERIC
              delW(1)  = real(delW(1), kind=4)  !y revenir
#endif

#   ifdef CFB_WIND_TRA
              Bf=-g/TairK*Wstar(2)*(Tstar+cpvir*TairK*Qstar)
              if (Bf.gt.0.0) then
                cff=blk_beta*(Bf*blk_Zabl)**r3
              else
                cff=0.2
              endif
              delW(2)  = SQRT(wspd0_cfb*wspd0_cfb+cff*cff)
#   endif

#  else
            delW(1)  = wspd0
#   ifdef CFB_WIND_TRA
            delW(2)  = wspd0_cfb
#   endif
#  endif
            ddelW(1:mb) = delW(1:mb)
!----------------------
          ENDDO !<-- terminate iterations-loop
# else
!======================================================================
!  BULK_FLUX_MLP: use pre-computed tile results (Cd/Ch/Ce already set
!  by the single tile-level call above, before the i,j loop).
!  Back-calculate Wstar, Tstar, Qstar so flux assembly below is unchanged.
!    Cd = (Wstar/delW)^2  =>  Wstar = sqrt(Cd)*delW
!    Ch = Wstar*Tstar/(delW*delT) => Tstar = Ch*delT/sqrt(Cd)
!    Ce = Wstar*Qstar/(delW*delQ) => Qstar = Ce*delQ/sqrt(Cd)
!======================================================================
          ddelW(1)  = delW(1)
          Wstar(1)  = SQRT(REAL(mlp_Cd(i,j))) * ddelW(1)
          Tstar     = REAL(mlp_Ch(i,j)) * delT /
     &                MAX(SQRT(REAL(mlp_Cd(i,j))), 1.E-10)
          Qstar     = REAL(mlp_Ce(i,j)) * delQ /
     &                MAX(SQRT(REAL(mlp_Cd(i,j))), 1.E-10)
          ddelW(1)  = MAX(ABS(ddelW(1)), 1.E-4)
# endif
! ---------------------
#  ifdef ABL1D
          chosqrcd  = vonKar/(logts10-psi_t)
#  endif
# endif

!===============================================================================
!   @@@@@@@@   @@@@@@@@@@      @@@@@@@@@@ @@          @@        @@  @@     @@
!  @@      @@  @@      @@      @@         @@          @@        @@    @@  @@
!  @@      @@  @@@@@@@@@@  &&  @@@@@@@    @@          @@        @@      @@
!  @@      @@  @@      @@      @@         @@          @@        @@    @@  @@
!   @@@@@@@@   @@      @@      @@         @@@@@@@@@@   @@@@@@@@@@   @@      @@
!===============================================================================
# ifdef ABL1D
          Cd          = (Wstar(1)/ddelW(1))**2
          Cd_du(i,j)  =                  Cd  * delW(1)
          Ch_du(i,j)  =  chosqrcd * SQRT(Cd) * delW(1)
          Ce_du(i,j)  =  chosqrcd * SQRT(Cd) * delW(1)
          ustar2(i,j) =  Wstar(1)*Wstar(1)
#  if defined BULK_ECUMEV0 || defined BULK_ECUMEV6
          z0_abl(i,j) = zr_abl(2) * EXP(-vonKar/SQRT(Cdn(1)))
#  else
          z0_abl(i,j) =  1./iZoW
#  endif
# else
!++  Compute turbulent sensible heat flux (W/m2).
          hfsen = - blk_Cpa*rhoAir*Wstar(mb)*Tstar
!++ Latent heat of vaporisation (J/kg)
          Hlv   = (2.5008 - 0.0023719*TseaC)*1.0E+6
!++  Compute turbulent latent heat flux (W/m2).
          hflat = - Hlv*rhoAir*Wstar(mb)*Qstar
!++ Longwave Radiation.
#  ifndef BULK_LW
          hflw   = -radlw(i,j) ! positive downward, i.e., warming the ocean
#  else
!++  Treat input longwave data as downwelling radiation only and add
!++  outgoing IR from model sea surface temperature.
          hflw   = radlw(i,j)  ! positive downward
     &           - emiss_lw*rho0i*cpi*SigmaSB*TseaK*TseaK*TseaK*TseaK
#  endif
!++  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
          upvel=-1.61*Wstar(mb)*Qstar-(1.0+1.61*Q)*Wstar(mb)
     &                                            *Tstar/TairK
!++  Compute turbulent latent heat flux (W/m2).
          hflat=hflat+rhoAir*Hlv*upvel*Q
!
!======================================================================
!  Latent and sensible heat fluxes
!======================================================================
!
          hflat=-hflat*rho0i*cpi
          hfsen=-hfsen*rho0i*cpi
!++  Total surface heat flux
          stflx(i,j,itemp)=srflx(i,j)+hflw+hflat+hfsen
#  ifdef SALINITY
          ! evap from latent heat flux
          evap=-cp*hflat/Hlv
#   ifdef RAIN_FLUX
          ! E-P bilan in m/s
          EmP(i,j)=evap-prate(i,j) 
          ! add temperature flux from rain
          stflx(i,j,itemp)=stflx(i,j,itemp)+prate(i,j)*t(i,j,N,nrhs,itemp)
#   else
          !++  Total surface salt flux
          stflx(i,j,isalt)=(evap-prate(i,j))*t(i,j,N,nrhs,isalt)
#   endif
#  endif
!
#  ifdef MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
#  endif
#  ifdef WET_DRY
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask_wet(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask_wet(i,j)
#  endif
!++  Save for the calculation of wind stress
          aer(i,j)  = rhoAir*delW(1)
          Cd        = (Wstar(1)/ddelW(1))**2
          cer(i,j)  = Cd
#  ifdef CFB_STRESS
          stau(i,j) = cfb_slope * delW(1) + cfb_offset
#  endif
!++  Store fluxes for diagnostics, when storing in netCDF file, multiplied by rho0.Cp to get W/m2
          shflx_rsw(i,j)=srflx(i,j)
          shflx_lat(i,j)=hflat
          shflx_sen(i,j)=hfsen
          shflx_rlw(i,j)=hflw
!----------------------
# endif    /* ABL1D */
!----------------------
        ENDDO !<-- terminate i-loop
      ENDDO   !<-- terminate j-loop
!----------------------

# ifdef BULK_FLUX_MLP
      DEALLOCATE( mlp_wspd, mlp_delT, mlp_delQ, mlp_TairC,
     &            mlp_RH, mlp_patm, mlp_Cd, mlp_Ch, mlp_Ce )
!  Advance the MLP time-step counter once per tile call.
!  mlp_tick() checks MOD(step_counter, MLP_SYNC_STEPS) and fires
!  mlp_sync() (MPI_Allreduce) when due.  Because bulk_flux_tile is
!  called once per OMP tile, and CROCO serialises tiles with a single
!  master thread when MPI is used without OpenMP (or with OMP MASTER),
!  this is safe.  With OpenMP tiling, only thread 0 should tick:
!$OMP MASTER
      CALL mlp_tick()
!$OMP END MASTER
# endif

# ifndef ABL1D
!
!======================================================================
!  SWITCH OFF shortwave radiation over "sea-ice"
!======================================================================
#  ifdef SEA_ICE_NOFLUX
!$acc loop seq
      do j=jmin,jmax
!$acc loop seq
        do i=imin,imax
              if( t(i,j,N,nrhs,itemp) .le. -1.8 ) then
                   stflx(i,j,itemp)=0.
#   if defined LMD_SKPP || defined LMD_BKPP || defined GLS_MIXING
                   srflx(i,j)=0.
#   endif
              endif
        enddo
      enddo
#  endif
!$acc loop seq
#  ifndef WAVE_SMFLUX
!
!======================================================================
!  Compute kinematic, surface wind stress (m2/s2).
!======================================================================
!
!$acc loop seq
      do j=jmin,jmax
!$acc loop seq
        do i=imin+1,imax
          cff =0.5*(cer(i-1,j)+cer(i,j))
          cff2=0.5*(aer(i-1,j)+aer(i,j))
          sustr(i,j)= ( cff*cff2*uwnd(i,j)
#   ifdef CFB_STRESS
     &     + 0.5*(stau(i,j)+stau(i-1,j))*u(i,j,N,nrhs)
#   endif
     &                             )*rho0i
#   ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
#   endif
#   ifdef WET_DRY
          sustr(i,j)=sustr(i,j)*umask_wet(i,j)
#   endif
        enddo
      enddo
!$acc end kernels      
!
!$acc kernels if(compute_on_device) default(present)
!$acc loop seq
      do j=jmin+1,jmax
!$acc loop seq
        do i=imin,imax
          cff =0.5*(cer(i,j-1)+cer(i,j))
          cff2=0.5*(aer(i,j-1)+aer(i,j))
          svstr(i,j)=( cff*cff2*vwnd(i,j)
#   ifdef CFB_STRESS
     &   + 0.5*(stau(i,j)+stau(i,j-1))*v(i,j,N,nrhs)
#   endif
     &                                        )*rho0i
#   ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
#   endif
#   ifdef WET_DRY
          svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
#   endif
        enddo
      enddo
!$acc end kernels      
#  endif /* !WAVE_SMFLUX */
# endif /* !ABL1D */
      return
      end


!======================================================================
!                     EXTERNAL FUNCTIONS
!======================================================================


!=======================================================================
!   Compute specific humidity from relative humidity
!=======================================================================
      FUNCTION spec_hum (RH,psfc,TairC)
!$acc routine(spec_hum) seq
      IMPLICIT NONE
# include "params_bulk.h"
      REAL       ::  RH   , spec_hum , cff
      REAL       ::  psfc, TairC
!$acc declare create(  RH,psfc,TairC )     
! !++  Compute air saturation vapor pressure (mb), using Teten formula.
      cff=(1.0007+3.46e-6*0.01*psfc)*6.1121*
     &        exp(17.502*TairC/(240.97+TairC))
!++  Compute specific humidity, Q (kg/kg).
      IF (RH.lt.2.0) then                             ! RH fraction
         cff=cff*RH                                   ! Vapor pres (mb)
         spec_hum=MvoMa*(cff/(psfc*0.01-0.378*cff))  ! Spec hum (kg/kg)
      ELSE             !RH input was actually specific humidity in g/kg
         spec_hum=0.001*RH                            ! Spec Hum (kg/kg)
      ENDIF
      END FUNCTION spec_hum
!=======================================================================
!    Compute Exner function from absolute air temperature
!=======================================================================
      SUBROUTINE exner_patm_from_tairabs (iexn,pair,q,tairabs,z,psfc)
!$acc routine(exner_patm_from_tairabs) seq
      IMPLICIT NONE
# include "params_bulk.h"
      REAL,INTENT(  out)   :: iexn,pair
      REAL,INTENT(in   )   :: q, tairabs, z, psfc
!++  Local variable declarations.
      REAL                 :: xm,q_sat,qsat
      REAL, PARAMETER      ::  g    = 9.80665
      INTEGER              :: iter
      INTEGER, PARAMETER   :: Niter = 3
!
      pair = psfc
!$acc loop seq
      DO Iter = 1, Niter
        q_sat = qsat(tairabs, pair, 1.)
        xm    =  mm_dryair + (q/q_sat) * ( mm_water - mm_dryair )
        pair  = psfc * EXP( -g * xm * z / ( r_gas * tairabs ) )
#ifdef OPENACCNUMERIC
        pair=real(pair,kind=4)
#endif
      ENDDO
      iexn =  (pair*ip00)**(-rdocpd)
!
      return
      END SUBROUTINE exner_patm_from_tairabs
!=======================================================================
!    Compute Saturation humidity
!=======================================================================
      FUNCTION qsat (TairK, patm, coeff)
!$acc routine(qsat) seq
      IMPLICIT NONE
# include "params_bulk.h"
      REAL                 ::  qsat
      REAL                 ::  TairK, patm, coeff
!++  Local variables declarations.
      REAL                 ::  psat
      REAL, PARAMETER      ::  alpw    = 60.2227554
      REAL, PARAMETER      ::  betaw   = 6822.40088
      REAL, PARAMETER      ::  gamw    = 5.13926744
      REAL, PARAMETER      ::  alpi    = 32.62117980819471
      REAL, PARAMETER      ::  betai   = 6295.421338904806
      REAL, PARAMETER      ::  gami    = 0.5631331575423155
!++  Compute Saturation Vapor Pressure
      IF (TairK .LE. CtoK) then
        psat = EXP( alpi - betai/TairK - gami*LOG(TairK) )
      ELSE
        psat = EXP( alpw - betaw/TairK - gamw*LOG(TairK) )
      ENDIF
      psat = coeff * psat
!++ Compute Saturation Humidity
      qsat = (MvoMa*psat)/(patm+(MvoMa-1.0)*psat)
#ifdef OPENACCNUMERIC
      qsat = real(qsat,kind=4)
#endif
      return
      END FUNCTION qsat

# if defined BULK_ECUMEV0 || defined BULK_ECUMEV6
!=======================================================================
!    Ecume stability function for velocity
!=======================================================================
      FUNCTION bulk_psiu_ecume (ZoL)
!$acc routine(bulk_psiu_ecume) seq
      IMPLICIT NONE
#  include "params_bulk.h"
      REAL                 ::  bulk_psiu_ecume
      REAL                 ::  ZoL
!++  Local variables declarations.
      REAL                 ::  chik, psik
      REAL                 ::  chic, psic
!
      IF (ZoL >= 0.0) THEN
        bulk_psiu_ecume = -7.0*ZoL
      ELSE
        chik  = (1.0-16.0*ZoL)**0.25
        psik  = 2.0*LOG(0.5*(1.0+chik))+LOG(0.5*(1.0+chik**2))
     &              -2.0*ATAN(chik)+pis2
        chic  = (1.0-12.87*ZoL)**r3    !for very unstable conditions
        psic  = 1.5*LOG(r3*(chic**2+chic+1.0))
     &               - sqr3*ATAN((2.0*chic+1.0)/sqr3)+2.*pis2osqr3
        bulk_psiu_ecume = psic+(psik-psic)/(1.0+ZoL**2) !match Kansas & free-conv. forms
#ifdef OPENACCNUMERIC
        bulk_psiu_ecume = real(bulk_psiu_ecume, kind=4)
#endif
      ENDIF
      return
      END FUNCTION bulk_psiu_ecume
!=======================================================================
!    Ecume stability function for tracers
!=======================================================================
      FUNCTION bulk_psit_ecume (ZoL)
!$acc routine(bulk_psit_ecume) seq
      IMPLICIT NONE
#  include "params_bulk.h"
      REAL                 ::  bulk_psit_ecume
      REAL                 ::  ZoL
      !  Local variables declarations.
      REAL                 ::  chik, psik
      REAL                 ::  chic, psic
      !
      IF (ZoL >= 0.0) THEN
        bulk_psit_ecume = -7.0*ZoL
      ELSE
        chik  = (1.0-16.0*ZoL)**0.25
        psik  = 2.0*LOG(0.5*(1.0+chik**2))
        chic  = (1.0-12.87*ZoL)**r3     !for very unstable conditions
        psic  = 1.5*LOG((chic**2+chic+1.0)*r3)
     &               -sqr3*ATAN((2.0*chic+1.0)/sqr3)
     &               +2.*pis2osqr3
        bulk_psit_ecume = psic+(psik-psic)/(1.0+ZoL**2)    !match Kansas & free-conv. forms
      ENDIF
      return
      END FUNCTION bulk_psit_ecume
!
# else
!=======================================================================
!    Compute molecular viscosity as a function of air temperature
!=======================================================================
      FUNCTION air_visc(TairC)
!$acc routine(air_visc) seq
      REAL                 :: air_visc,cff
      REAL, PARAMETER      :: c0 = 1.326E-5
      REAL, PARAMETER      :: c1 = 6.542E-3
      REAL, PARAMETER      :: c2 = 8.301E-6
      REAL, PARAMETER      :: c3 = 4.84E-9
      cff      = TairC*TairC
      air_visc = c0*(1.+c1*TairC+c2*cff-c3*cff*TairC)
      return
      END FUNCTION air_visc
!=======================================================================
!    Coare stability function for velocity
!=======================================================================
      FUNCTION bulk_psiu_coare (ZoL)
!$acc routine(bulk_psiu_coare) seq
      IMPLICIT NONE
#  include "params_bulk.h"
      REAL                 ::  bulk_psiu_coare
      REAL                 ::  ZoL
!  Local variables declarations.
      REAL                 ::  chik, psik
      REAL                 ::  chic, psic
!
      IF (ZoL <= 0.0) then  !  Unstable conditions.
        chik = (1.0-15.0*ZoL)**0.25
        psik = 2.0*LOG(0.5*(1.0+chik))+LOG(0.5*(1.0+chik**2))
     &                                  -2.0*ATAN(chik)+pis2
        chic = (1.0-10.15*ZoL)**r3
        psic  = 1.5*LOG(r3*(chic**2+chic+1.0))
     &          - sqr3*ATAN((2.0*chic+1.0)/sqr3)+2.*pis2osqr3
        bulk_psiu_coare=psic+(psik-psic)/(1.0+ZoL**2)
      ELSE                    !  Stable conditions
        chic=-MIN(50.0,0.35*ZoL)
        bulk_psiu_coare=-((1.0+ZoL)+0.6667*(ZoL-14.28)*EXP(chic)+8.525)
      ENDIF
#ifdef OPENACCNUMERIC
       bulk_psiu_coare = real( bulk_psiu_coare, kind=4 )
#endif
      return
      END FUNCTION bulk_psiu_coare
!=======================================================================
!    Coare stability function for tracers
!=======================================================================
      FUNCTION bulk_psit_coare (ZoL)
!$acc routine(bulk_psit_coare) seq
      IMPLICIT NONE
#  include "params_bulk.h"
      REAL                 ::  bulk_psit_coare
      REAL                 ::  ZoL
!  Local variables declarations.
      REAL                 ::  chik, psik
      REAL                 ::  chic, psic
!
      IF (ZoL < 0.0) THEN
        chik  = (1.0-15.0*ZoL)**0.25
        psik  = 2.0*LOG(0.5*(1.0+chik**2))
        chic  = (1.0-34.15*ZoL)**r3     !for very unstable conditions
        psic  = 1.5*LOG((chic**2+chic+1.0)*r3)
     &               -sqr3*ATAN((2.0*chic+1.0)/sqr3)
     &               +2.*pis2osqr3
        bulk_psit_coare = psic+(psik-psic)/(1.0+ZoL**2) !match Kansas & free-conv. forms
#ifdef OPENACCNUMERIC        
        bulk_psit_coare = real(bulk_psit_coare,kind=4)
#endif        
      ELSE
        chic=-MIN(50.0,0.35*ZoL)
        bulk_psit_coare = -((1.0+2.0*ZoL/3.0)**1.5+
     &            0.6667*(ZoL-14.28)*EXP(chic)+8.525)
#ifdef OPENACCNUMERIC        
        bulk_psit_coare = real(bulk_psit_coare,kind=4)
#endif        
      ENDIF
      return
      END FUNCTION bulk_psit_coare
!=======================================================================
# endif

#else
      subroutine bulk_flux_empty
      return
      end
#endif
