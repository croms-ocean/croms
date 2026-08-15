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
!
!  This is the "online_interpolate_3d_bulk.F" script
!------------------------------------------------------------------------------
!  This file contains the subfunctions enabling the ONLINE extraction of the
!  forcing from an atmospheric dataset. A spatial and time interpolation are 
!  applied to the extracted data in order to adapt these to the considered 
!  simulation domain as well as the associated MPI/OPENMP discretisation (if 
!  defined MPI/OPENMP).
!------------------------------------------------------------------------------
#include "cppdefs.h"

#if defined BULK_FLUX && defined ONLINE && defined ABL1D


!*******************************************************************************
      subroutine interpolate_3D_bulk_online(NXref,NYref,NZref,
     &                                      lonref,latref,varref,
     &                                      blkvar_id,iblkrec)
!------------------------------------------------------------------------------
!  This subfunction complete the call to the interpolation and the recalculation
!  of the forcing fields to generate the roms bulk forcing (tairg, rhumg, prateg,
!  radlwg, radswg, uwndg, vwndg) on the current tile.
!------------------------------------------------------------------------------
!  The main steps of this interpolation and recalculation are:
!  - Interpolation from the varref matrix of the dum_array matrix on the
!  simulation domain for the blkvar_id variable.
!  - Recalculation of the roms bulk forcings (tairg, rhumg, prateg, radlwg,
!  radswg, uwndg, vwndg) from the interpolated dum_array fields.

!  INPUTS
!  NXref (Vertical size of the input matrix used for the interpolation)
!  NYref (Horizontal size of the input matrix used for the interpolation)
!  lonref (longitudes of the references data used for the interpolation)
!  latref (latitudes of the references data used for the interpolation)
!  varref (reference data used for the interpolation)
!  Istr, Iend, Jstr, Jend (limits index of the interpolation)
!  blkvar_id (index of the extrapolated bulk variable)
!  iblkrec (index of the leap alternation)
!
!  OUTPUTS
!  dum_array (interpolated forcing fields)
!  tairg, rhumg, prateg, radlwg, radswg, uwndg, vwndg (recalculated roms
!                                                             bulk forcings)
      implicit none
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "online.h"
# include "abl1d.h"
      integer, intent(in) :: NXref, NYref, NZref, blkvar_id, iblkrec
      real(kind=8), intent(in) :: lonref(1:NXref), latref(1:NYref)
      real(kind=8), intent(in) :: varref(1:NXref,1:NYref,1:NZref)
      real blk_dum1(GLOBAL_2D_ARRAY)
      real blk_dum2(GLOBAL_2D_ARRAY)
      real Pref, ew, Qsat
      integer i,j,k


      Pref=1020.         ! default air pressure [mbars]

      !--------
      ! Limits:
      !--------
      ! Local loop ranges
      ! This is not parallel with openMP

# ifdef MPI
#  define I_RANGE 0,Lmmpi+1
#  define J_RANGE 0,Mmmpi+1
# else
#  define I_RANGE 0,LLm+1
#  define J_RANGE 0,MMm+1
# endif

      !===============
      ! Interpolation:
      !===============

# ifdef CUBIC_INTERP
#  define MYINTERP cinterp2d
# else
#  define MYINTERP linterp2d
# endif

      if (blkvar_id.eq.1) then
!
! 1 process Tair
! Temperature: Convert from Kelvin to Celsius
!
        do k=1, N_abl
        call MYINTERP(1, NXref, 1, NYref,
     &                lonref, latref, varref(:,:,k),
     &                I_RANGE,J_RANGE,
     &                lonr(GLOBAL_2D_ARRAY),
     &                latr(GLOBAL_2D_ARRAY),
     &                pt_dtag(GLOBAL_2D_ARRAY,k,iblkrec))
        enddo

        do j=J_RANGE
          do i=I_RANGE
            pt_dtag(i,j,:,iblkrec)=pt_dtag(i,j,:,iblkrec)-273.15
            tairg(i,j,iblkrec)=pt_dtag(i,j,1,iblkrec)-273.15
          enddo
        enddo

# ifdef MPI
        call exchange_r2d_tile(1,Lm,1,Mm,
     &                         tairg(START_2D_ARRAY,iblkrec))
        do k=1, N_abl
          call exchange_r2d_tile (1,Lm,1,Mm,
     &                          pt_dtag(START_2D_ARRAY,k,iblkrec))
        enddo
# endif

      elseif  (blkvar_id.eq.2) then
!
! 2 process Rhum
! Relative humidity: Convert specific humidity to relative
! (except from METEO-FRANCE fluxes where it's already RH)
!
        do k=1, N_abl
        call MYINTERP(1, NXref, 1, NYref,
     &                lonref, latref, varref(:,:,k),
     &                I_RANGE,J_RANGE,
     &                lonr(GLOBAL_2D_ARRAY),
     &                latr(GLOBAL_2D_ARRAY),
     &                pq_dtag(GLOBAL_2D_ARRAY,k,iblkrec))
        enddo

        do j=J_RANGE
          do i=I_RANGE
            rhumg(i,j,iblkrec)=pq_dtag(i,j,1,iblkrec)
          enddo
        enddo

# ifndef AROME
        do j=J_RANGE
          do i=I_RANGE
            ew=6.1121*(1.0007+3.46e-6*Pref)*
     &         exp((17.502*tairg(i,j,iblkrec))/
     &             (240.97+tairg(i,j,iblkrec)))
            Qsat=0.62197*(ew/(Pref-0.378*ew))
            rhumg(i,j,iblkrec)=rhumg(i,j,iblkrec)/Qsat
          enddo
        enddo
# endif

# ifdef MPI
        call exchange_r2d_tile(1,Lm,1,Mm,
     &                         rhumg(START_2D_ARRAY,iblkrec))
        do k=1, N_abl
          call exchange_r2d_tile (1,Lm,1,Mm,
     &                          pq_dtag(START_2D_ARRAY,k,iblkrec))
        enddo
# endif

      elseif (blkvar_id.eq.8) then
!
! 8 process uwnd
! U-component_of_wind
!
        do k=1, N_abl
        call MYINTERP(1, NXref, 1, NYref,
     &                lonref, latref, varref(:,:,k),
     &                I_RANGE,J_RANGE,
     &                lonr(GLOBAL_2D_ARRAY),
     &                latr(GLOBAL_2D_ARRAY),
     &                pu_dtag(GLOBAL_2D_ARRAY,k,iblkrec))
        enddo

        do j=J_RANGE
          do i=I_RANGE
            uwndg(i,j,iblkrec)=pu_dtag(i,j,1,iblkrec)
          enddo
        enddo

      elseif (blkvar_id.eq.9) then
!
! 9 process vwnd
! V-component_of_wind
!
        do k=1, N_abl
        call MYINTERP(1, NXref, 1, NYref,
     &                lonref, latref, varref(:,:,k),
     &                I_RANGE,J_RANGE,
     &                lonr(GLOBAL_2D_ARRAY),
     &                latr(GLOBAL_2D_ARRAY),
     &                pv_dtag(GLOBAL_2D_ARRAY,k,iblkrec))
        enddo

        do j=J_RANGE
          do i=I_RANGE
            vwndg(i,j,iblkrec)=pv_dtag(i,j,1,iblkrec)
          enddo
        enddo

        do j=J_RANGE
          do i=I_RANGE
            wspdg(i,j,iblkrec)= sqrt(uwndg(i,j,iblkrec)
     &                              *uwndg(i,j,iblkrec)
     &                              +vwndg(i,j,iblkrec)
     &                              *vwndg(i,j,iblkrec))
           enddo
        enddo

# ifdef MPI
        call exchange_r2d_tile(1,Lm,1,Mm,
     &                         uwndg(START_2D_ARRAY,iblkrec))
        call exchange_r2d_tile(1,Lm,1,Mm,
     &                         vwndg(START_2D_ARRAY,iblkrec))
        call exchange_r2d_tile(1,Lm,1,Mm,
     &                         wspdg(START_2D_ARRAY,iblkrec))

        do k=1, N_abl
          call exchange_r2d_tile(1,Lm,1,Mm,
     &                         pu_dtag(START_2D_ARRAY,k,iblkrec))
          call exchange_r2d_tile(1,Lm,1,Mm,
     &                         pv_dtag(START_2D_ARRAY,k,iblkrec))
        enddo
# endif

      elseif (blkvar_id.eq.11) then

!
! 11 process uhpg
! U geostrophic wind
!
        do k=1, N_abl
        call MYINTERP(1, NXref, 1, NYref,
     &                lonref, latref, varref(:,:,k),
     &                I_RANGE,J_RANGE,
     &                lonr(GLOBAL_2D_ARRAY),
     &                latr(GLOBAL_2D_ARRAY),
     &                pgu_dtag(GLOBAL_2D_ARRAY,k,iblkrec))
        enddo

# ifdef MPI
        do k=1, N_abl
          call exchange_r2d_tile(1,Lm,1,Mm,
     &                         pgu_dtag(START_2D_ARRAY,k,iblkrec))
        enddo
# endif

      elseif (blkvar_id.eq.12) then

!
! 12 process vhpg
! V geostrophic wind
!
        do k=1, N_abl
        call MYINTERP(1, NXref, 1, NYref,
     &                lonref, latref, varref(:,:,k),
     &                I_RANGE,J_RANGE,
     &                lonr(GLOBAL_2D_ARRAY),
     &                latr(GLOBAL_2D_ARRAY),
     &                pgv_dtag(GLOBAL_2D_ARRAY,k,iblkrec))
        enddo

# ifdef MPI
        do k=1, N_abl
          call exchange_r2d_tile(1,Lm,1,Mm,
     &                         pgv_dtag(START_2D_ARRAY,k,iblkrec))
        enddo
# endif


      endif

      return
      end


#endif /* BULK_FLUX && ONLINE && ABL1D */
