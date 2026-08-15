! MUSTANG - This software is governed by the CeCILL-C license
! see LICENSE_MUSTANG.txt

MODULE lateral_erosion
   !============================================================================
   !!                   ***  MODULE  lateral_erosion  ***
   !!
   !! This module manage lateral erosion effect
   !============================================================================

#include "cppdefs.h"

#if defined MUSTANG

   USE comMUSTANG, ONLY: coef_erolat, htncrit_eros, coef_tauskin_lat, &
                         l_erolat_wet_cell, flx_s2w, &
                         flx_s2w_corip1, flx_s2w_corim1, &
                         flx_s2w_corjp1, flx_s2w_corjm1, &
                         phieau_s2w_corip1, phieau_s2w_corim1, &
                         phieau_s2w_corjp1, phieau_s2w_corjm1, &
                         h0fond
   USE comsubstance, ONLY: nv_adv, rsh

   IMPLICIT NONE
   PRIVATE

   PUBLIC lateral_erosion_reset
   PUBLIC lateral_erosion_get
   PUBLIC lateral_erosion_compute
   PUBLIC lateral_erosion_saveflx
   PUBLIC lateral_erosion_apply
   ! TYPE
   PUBLIC lateral_erosion_type

   TYPE lateral_erosion_type
      LOGICAL        :: l_drycell ! True if computation was in a drycell
      REAL(KIND=rsh) :: hgt_e ! height use in computation at east side
      REAL(KIND=rsh) :: hgt_w ! height use in computation at west side
      REAL(KIND=rsh) :: hgt_n ! height use in computation at north side
      REAL(KIND=rsh) :: hgt_s ! height use in computation at south side
      REAL(KIND=rsh) :: factor_surf_e ! surfcell divide by surfcell at east side
      REAL(KIND=rsh) :: factor_surf_w ! surfcell divide by surfcell at west side
      REAL(KIND=rsh) :: factor_surf_n ! surfcell divide by surfcell at north side
      REAL(KIND=rsh) :: factor_surf_s ! surfcell divide by surfcell at south side
      REAL(KIND=rsh) :: tau_e ! equivalent stress used in computation at east side
      REAL(KIND=rsh) :: tau_w ! equivalent stress used in computation at west side
      REAL(KIND=rsh) :: tau_n ! equivalent stress used in computation at north side
      REAL(KIND=rsh) :: tau_s ! equivalent stress used in computation at south side
      REAL(KIND=rsh) :: ero ! sum of lateral erosion fluxes (all 4 sides)
      REAL(KIND=rsh) :: flx_e ! fluxes correction factor at east side
      REAL(KIND=rsh) :: flx_w ! fluxes correction factor at west side
      REAL(KIND=rsh) :: flx_n ! fluxes correction factor at north side
      REAL(KIND=rsh) :: flx_s ! fluxes correction factor atsouth side
      REAL(KIND=rsh) :: phieau_e ! water fluxes correction factor at east side
      REAL(KIND=rsh) :: phieau_w ! water fluxes correction factor at west side
      REAL(KIND=rsh) :: phieau_n ! water fluxes correction factor at north side
      REAL(KIND=rsh) :: phieau_s ! water fluxes correction factor at south side
   END TYPE

CONTAINS

   !============================================================================
   SUBROUTINE lateral_erosion_reset()
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE lateral_erosion_reset ***
      !!
      !! ** Purpose : Reset LATERAL EROSION arrays to zero at begining of
      !!              time step
      !! ** Called by : MUSTANG_update
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE

      IF (coef_erolat .NE. 0.0_rsh) THEN
         flx_s2w_corip1(:, :, :) = 0.0_rsh
         flx_s2w_corim1(:, :, :) = 0.0_rsh
         flx_s2w_corjp1(:, :, :) = 0.0_rsh
         flx_s2w_corjm1(:, :, :) = 0.0_rsh
#if ! defined key_nofluxwat_IWS
         phieau_s2w_corip1(:, :) = 0.0_rsh
         phieau_s2w_corim1(:, :) = 0.0_rsh
         phieau_s2w_corjp1(:, :) = 0.0_rsh
         phieau_s2w_corjm1(:, :) = 0.0_rsh
#endif
      END IF  ! end if coef_erolat
   END SUBROUTINE lateral_erosion_reset

   !============================================================================
   SUBROUTINE lateral_erosion_get(htotij, htote, htotw, htots, htotn, surf_ij, &
                                  surf_e, surf_w, surf_s, surf_n, &
                                  V_NEAR_E, V_NEAR_W, U_NEAR_S, U_NEAR_N, &
                                  erolat_res)
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE lateral_erosion_get ***
      !!
      !! ** Purpose : prepare variable from htot and current values at
      !!              the begining of erosion computation to avoid to redo it
      !!              at each sub time step
      !! ** Called by : MUSTANG_erosion
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=rsh), INTENT(IN) :: htotij
      REAL(KIND=rsh), INTENT(IN) :: htote, htotw, htots, htotn
      REAL(KIND=rsh), INTENT(IN) :: V_NEAR_E, V_NEAR_W, U_NEAR_S, U_NEAR_N
      REAL(KIND=rsh), INTENT(IN) :: surf_ij
      REAL(KIND=rsh), INTENT(IN) :: surf_e, surf_w, surf_s, surf_n
      TYPE(lateral_erosion_type), INTENT(INOUT) :: erolat_res

      REAL(KIND=rsh) :: hwat

      if (htotij > h0fond) then
         erolat_res%l_drycell = .false.
         hwat = htotij - h0fond
         erolat_res%hgt_e = max(0.0_rsh, (htote - htncrit_eros) - hwat)
         erolat_res%hgt_w = max(0.0_rsh, (htotw - htncrit_eros) - hwat)
         erolat_res%hgt_n = max(0.0_rsh, (htotn - htncrit_eros) - hwat)
         erolat_res%hgt_s = max(0.0_rsh, (htots - htncrit_eros) - hwat)
      else
         erolat_res%l_drycell = .true.
         erolat_res%hgt_e = max(0.0_rsh, htote - htncrit_eros)
         erolat_res%hgt_w = max(0.0_rsh, htotw - htncrit_eros)
         erolat_res%hgt_n = max(0.0_rsh, htotn - htncrit_eros)
         erolat_res%hgt_s = max(0.0_rsh, htots - htncrit_eros)
      end if

      erolat_res%tau_e = coef_tauskin_lat*V_NEAR_E**2
      erolat_res%tau_w = coef_tauskin_lat*V_NEAR_W**2
      erolat_res%tau_n = coef_tauskin_lat*U_NEAR_N**2
      erolat_res%tau_s = coef_tauskin_lat*U_NEAR_S**2

      erolat_res%factor_surf_e = surf_ij/surf_e
      erolat_res%factor_surf_w = surf_ij/surf_w
      erolat_res%factor_surf_n = surf_ij/surf_n
      erolat_res%factor_surf_s = surf_ij/surf_s

      erolat_res%ero = 0.0_rsh
      erolat_res%flx_e = 0.0_rsh
      erolat_res%flx_w = 0.0_rsh
      erolat_res%flx_n = 0.0_rsh
      erolat_res%flx_s = 0.0_rsh
      erolat_res%phieau_e = 0.0_rsh
      erolat_res%phieau_w = 0.0_rsh
      erolat_res%phieau_n = 0.0_rsh
      erolat_res%phieau_s = 0.0_rsh

   END SUBROUTINE lateral_erosion_get

   !============================================================================
   SUBROUTINE lateral_erosion_compute(toce, erolat_res)
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE lateral_erosion_compute ***
      !!
      !! ** Purpose : compute erolat fluxes
      !! ** Called by : MUSTANG_erosion
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=rsh), INTENT(IN)    :: toce
      TYPE(lateral_erosion_type), INTENT(INOUT) :: erolat_res

      REAL(KIND=rsh) :: ero, eroe, erow, eron, eros

      eroe = 0.0_rsh
      erow = 0.0_rsh
      eros = 0.0_rsh
      eron = 0.0_rsh
      IF (coef_erolat .NE. 0.0_rsh .AND. &
          (l_erolat_wet_cell .OR. erolat_res%l_drycell)) THEN
         eroe = lateral_erosion_value(erolat_res%hgt_e, erolat_res%tau_e, toce)
         erow = lateral_erosion_value(erolat_res%hgt_w, erolat_res%tau_w, toce)
         eros = lateral_erosion_value(erolat_res%hgt_s, erolat_res%tau_s, toce)
         eron = lateral_erosion_value(erolat_res%hgt_n, erolat_res%tau_n, toce)
      END IF
      ero = eroe + erow + eros + eron
      erolat_res%ero = ero

      erolat_res%flx_e = eroe/ero*erolat_res%factor_surf_e
      erolat_res%flx_w = erow/ero*erolat_res%factor_surf_w
      erolat_res%flx_s = eros/ero*erolat_res%factor_surf_s
      erolat_res%flx_n = eron/ero*erolat_res%factor_surf_n

      erolat_res%phieau_e = eroe/ero
      erolat_res%phieau_w = erow/ero
      erolat_res%phieau_s = eros/ero
      erolat_res%phieau_n = eron/ero

   END SUBROUTINE lateral_erosion_compute

   !============================================================================
   REAL FUNCTION lateral_erosion_value(height, tau, toce)
      !!------------------------------------------------------------------------
      !!                ***  FUNCTION lateral_erosion_value ***
      !!
      !! ** Purpose : Compute lateral fluxes from height and equivalent stress
      !!------------------------------------------------------------------------
      IMPLICIT NONE
          !! * Arguments
      REAL(KIND=rsh), INTENT(IN) :: height
      REAL(KIND=rsh), INTENT(IN) :: tau
      REAL(KIND=rsh), INTENT(IN) :: toce

      lateral_erosion_value = 0.0_rsh

      IF (coef_erolat .NE. 0.0_rsh) THEN
         lateral_erosion_value = coef_erolat &
                                 *max(0.0_rsh, tau - toce)*height
      END IF

   END FUNCTION lateral_erosion_value

   !============================================================================
   SUBROUTINE lateral_erosion_saveflx(erolat_res, flx_s2w_eroij &
                                      , flx_s2w_corip1_ij, flx_s2w_corim1_ij &
                                      , flx_s2w_corjm1_ij, flx_s2w_corjp1_ij &
#if ! defined key_nofluxwat_IWS
                                      , phieau_ero_ij &
                                      , phieau_s2w_corip1_ij &
                                      , phieau_s2w_corim1_ij &
                                      , phieau_s2w_corjm1_ij &
                                      , phieau_s2w_corjp1_ij &
#endif
                                      )
      !!------------------------------------------------------------------------
      !!                ***  ROUTINE lateral_erosion_saveflx ***
      !!
      !! ** Purpose : memorisation of lateral erosion fluxes from dry cell
      !! ** Called by : MUSTANG_erosion
      !!------------------------------------------------------------------------

      implicit none

      TYPE(lateral_erosion_type), INTENT(IN) :: erolat_res
      REAL(KIND=rsh), DIMENSION(-1:nv_adv), INTENT(INOUT)  ::  flx_s2w_eroij
      REAL(KIND=rsh), DIMENSION(-1:nv_adv), INTENT(INOUT)  ::  flx_s2w_corip1_ij
      REAL(KIND=rsh), DIMENSION(-1:nv_adv), INTENT(INOUT)  ::  flx_s2w_corim1_ij
      REAL(KIND=rsh), DIMENSION(-1:nv_adv), INTENT(INOUT)  ::  flx_s2w_corjm1_ij
      REAL(KIND=rsh), DIMENSION(-1:nv_adv), INTENT(INOUT)  ::  flx_s2w_corjp1_ij
#if ! defined key_nofluxwat_IWS
      REAL(KIND=rsh), INTENT(INOUT)  ::  phieau_ero_ij
      REAL(KIND=rsh), INTENT(INOUT)  ::  phieau_s2w_corip1_ij
      REAL(KIND=rsh), INTENT(INOUT)  ::  phieau_s2w_corim1_ij
      REAL(KIND=rsh), INTENT(INOUT)  ::  phieau_s2w_corjm1_ij
      REAL(KIND=rsh), INTENT(INOUT)  ::  phieau_s2w_corjp1_ij
#endif

      INTEGER :: iv

      ! memorisation of lateral erosion for dry cell
      IF (erolat_res%l_drycell) THEN
         DO iv = -1, nv_adv
            flx_s2w_corip1_ij(iv) = flx_s2w_corip1_ij(iv) + &
                                    flx_s2w_eroij(iv)*erolat_res%flx_e
            flx_s2w_corim1_ij(iv) = flx_s2w_corim1_ij(iv) + &
                                    flx_s2w_eroij(iv)*erolat_res%flx_w
            flx_s2w_corjm1_ij(iv) = flx_s2w_corjm1_ij(iv) + &
                                    flx_s2w_eroij(iv)*erolat_res%flx_s
            flx_s2w_corjp1_ij(iv) = flx_s2w_corjp1_ij(iv) + &
                                    flx_s2w_eroij(iv)*erolat_res%flx_n
            flx_s2w_eroij(iv) = 0.0_rsh
         END DO
#if ! defined key_nofluxwat_IWS
         phieau_s2w_corip1_ij = phieau_s2w_corip1_ij + &
                                phieau_ero_ij*erolat_res%phieau_e
         phieau_s2w_corim1_ij = phieau_s2w_corim1_ij + &
                                phieau_ero_ij*erolat_res%phieau_w
         phieau_s2w_corjm1_ij = phieau_s2w_corjm1_ij + &
                                phieau_ero_ij*erolat_res%phieau_s
         phieau_s2w_corjp1_ij = phieau_s2w_corjp1_ij + &
                                phieau_ero_ij*erolat_res%phieau_n
         phieau_ero_ij = 0.0_rsh
#endif
      END IF

   END SUBROUTINE lateral_erosion_saveflx

   !============================================================================
   SUBROUTINE lateral_erosion_apply(ifirst, ilast, jfirst, jlast, dtinv)
      !!------------------------------------------------------------------------
      !!                ***  ROUTINE lateral_erosion_apply ***
      !!
      !! ** Purpose : Apply fluxes correction if LATERAL EROSION
      !! ** Called by : MUSTANG_update
      !!------------------------------------------------------------------------

      implicit none
          !! * Arguments
      INTEGER, INTENT(IN) :: ifirst, ilast, jfirst, jlast
      REAL(KIND=rsh) :: dtinv

      INTEGER :: iv, i, j

      IF (coef_erolat .NE. 0.0_rsh) THEN
#if defined MPI
         call lateral_erosion_exchange(ifirst, ilast, jfirst, jlast)
#endif
         ! correction : neighboring cells of eroded laterally dry cell receive
         ! one fraction of eroded sediment
         DO j = jfirst, jlast
            DO i = ifirst, ilast
               DO iv = -1, nv_adv
                  flx_s2w(iv, i, j) = flx_s2w(iv, i, j) + &
                                      dtinv*(flx_s2w_corip1(iv, i - 1, j) + &
                                             flx_s2w_corim1(iv, i + 1, j) + &
                                             flx_s2w_corjp1(iv, i, j - 1) + &
                                             flx_s2w_corjm1(iv, i, j + 1))
               END DO
#if ! defined key_nofluxwat_IWS
               phieau_s2w(i, j) = phieau_s2w(i, j) + &
                                  phieau_s2w_corip1(i - 1, j) + &
                                  phieau_s2w_corim1(i + 1, j) + &
                                  phieau_s2w_corjp1(i, j - 1) + &
                                  phieau_s2w_corjm1(i, j + 1)
#endif
            END DO
         END DO

      END IF  ! end if coef_erolat
   END SUBROUTINE lateral_erosion_apply

   !============================================================================
#if defined MPI
   SUBROUTINE lateral_erosion_exchange(ifirst, ilast, jfirst, jlast)
      !-------------------------------------------------------------------------
      !!                 ***  ROUTINE lateral_erosion_exchange ***
      !!
      !! ** Purpose : MPI exchange of lateral erosion flux between processors
      !! ** Called by : lateral_erosion_apply
      !-------------------------------------------------------------------------

      USE module_substance ! for dimensions

      !! * Arguments
      INTEGER, INTENT(IN) :: ifirst, ilast, jfirst, jlast

      REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch
      INTEGER :: iv

      do iv = -1, nv_adv
         workexch(:, :) = flx_s2w_corip1(iv, :, :)
         call exchange_r2d_tile(ifirst, ilast, jfirst, jlast,  &
                 &          workexch(START_2D_ARRAY))
         flx_s2w_corip1(iv, :, :) = workexch(:, :)

         workexch(:, :) = flx_s2w_corim1(iv, :, :)
         call exchange_r2d_tile(ifirst, ilast, jfirst, jlast,  &
                 &          workexch(START_2D_ARRAY))
         flx_s2w_corim1(iv, :, :) = workexch(:, :)

         workexch(:, :) = flx_s2w_corjp1(iv, :, :)
         call exchange_r2d_tile(ifirst, ilast, jfirst, jlast,  &
                 &          workexch(START_2D_ARRAY))
         flx_s2w_corjp1(iv, :, :) = workexch(:, :)

         workexch(:, :) = flx_s2w_corjm1(iv, :, :)
         call exchange_r2d_tile(ifirst, ilast, jfirst, jlast,  &
                 &          workexch(START_2D_ARRAY))
         flx_s2w_corjm1(iv, :, :) = workexch(:, :)
      end do

      !! TODO missing exchange for phieau* ??

   END SUBROUTINE lateral_erosion_exchange
#endif /* defined MPI  */

#endif /* ifdef MUSTANG */
   !============================================================================
END MODULE lateral_erosion
