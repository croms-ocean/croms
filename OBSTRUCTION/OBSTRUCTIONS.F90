MODULE OBSTRUCTIONS

#include "cppdefs.h"

#ifdef OBSTRUCTION
   !!==========================================================================
   !!                   ***  MODULE  OBSTRUCTIONS  ***
   !!
   !! ** Purpose : obstructions interactions with hydrodynamics,
   !!              link from hydro model to OBSTRUCTIONS1DV module
   !!
   !!              Strong dependancy with hydro model for mesh size,
   !!              position of u,v, name of variable, order of indexes (i,j,k),
   !!              MPI exchange...
   !!
   !! ** Description :
   !!
   !!     subroutine obst_update  ! Main subroutine controlling updates
   !!                               of obstructions characteristics
   !!     subroutine obst_rd_timeserie ! Read file for time-varying
   !!                                  ! obstructions characteristics
   !!==========================================================================
   USE com_OBSTRUCTIONS
   USE OBSTRUCTIONS1DV

   IMPLICIT NONE

   PUBLIC obst_update
   PRIVATE

CONTAINS
   !!==========================================================================
   SUBROUTINE obst_update(limin, limax, ljmin, ljmax)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_update  ***
      !!
      !! ** Purpose : Update of obstruction parameters at each time step
      !!---------------------------------------------------------------------
      USE module_OBSTRUCTIONS  ! needed for zob, h, u,v, cm0, Zt_avg1, time
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: limin, limax, ljmin, ljmax

      TYPE(o1dv_out_type) :: obst_output_ij
      REAL(KIND=rsh)      :: hwat
      REAL(KIND=rsh), DIMENSION(obst_kmax)   :: obst_uz
      REAL(KIND=rsh), DIMENSION(obst_kmax)   :: obst_vz
      REAL(KIND=rsh), DIMENSION(obst_kmax)   :: obst_zc
      REAL(KIND=rsh), DIMENSION(obst_kmax)   :: obst_dz
      INTEGER :: i, j, k

      ! ******************************
      ! * READING CHARACTERISTICS FILE
      ! ******************************
      ! update obst_height_inst, obst_width_inst,
      ! obst_thick_inst, obst_dens_inst
      CALL obst_rd_timeserie(time)

      DO j = ljmin, ljmax
         DO i = limin, limax
            ! **************************************
            ! * COMPUTES 3D VELOCITIES AT RHO POINT
            ! **************************************
            obst_uz(:) = (u(i, j, :, nstp) + u(i + 1, j, :, nstp))/2.
            obst_vz(:) = (v(i, j, :, nstp) + v(i, j + 1, :, nstp))/2.

            ! ***********************************************
            ! * COMPUTES CELL THICKNESS AND HEIGHT AT CENTER
            ! ***********************************************
            obst_dz(:) = Hz(i, j, :)
            obst_zc(1) = Hz(i, j, 1)/2.
            DO k = 2, obst_kmax
               obst_zc(k) = obst_zc(k - 1) + Hz(i, j, k - 1)/2.+Hz(i, j, k)/2.
            END DO
            ! at this phase of time stepping, total water height is obtain with Zt_avg1
            hwat = h(i, j) + Zt_avg1(i, j)

            ! ************************
            ! * UPDATING OBSTRUCTION
            ! ************************
            CALL o1dv_main(hwat, cm0, obst_z0bed(i, j), obst_uz, obst_vz, obst_dz, obst_zc, &
                           obst_height_inst(:, i, j), obst_dens_inst(:, i, j), &
                           obst_thick_inst(:, i, j), obst_width_inst(:, i, j), &
                           obst_position(:, i, j), obst_height(:, i, j), &
                           obst_output_ij)

            ! ************************
            ! * SAVE FOR OUTPUT
            ! ************************
            obst_height(:, i, j) = obst_output_ij%height(:)
            obst_fuz(i, j, :) = obst_output_ij%fuz(:)
            obst_fvz(i, j, :) = obst_output_ij%fvz(:)
            obst_a3d(:, :, i, j) = obst_output_ij%a3d(:, :)
            obst_tau(i, j, :) = obst_output_ij%tau(:)
            obst_t(i, j, :) = obst_output_ij%t(:)
            IF (obst_nv_noturb > 0) THEN
               zob(i, j) = obst_output_ij%z0obst(obst_nbvar + 1)
            END IF
            ! Non mandatory variables, allocated only if wanted in output
            IF (l_obstout_dens_e) then
               obst_dens3d(:, :, i, j) = obst_output_ij%dens3d(:, :)
            END IF
            IF (l_obstout_width_e) then
               obst_width3d(:, :, i, j) = obst_output_ij%width3d(:, :)
            END IF
            IF (l_obstout_thick_e) then
               obst_thick3d(:, :, i, j) = obst_output_ij%thick3d(:, :)
            END IF
            IF (l_obstout_theta) then
               obst_theta3d(:, :, i, j) = obst_output_ij%theta3d(:, :)
            END IF
            IF (l_obstout_frac_xy) then
               obst_fracxy(:, i, j) = obst_output_ij%fracxy(:)
            END IF
            IF (l_obstout_frac_z) then
               obst_fracz3d(:, :, i, j) = obst_output_ij%fracz3d(:, :)
            END IF
            IF (l_obstout_a2d) then
               obst_a2d(:, i, j) = obst_output_ij%a2d(:)
            END IF
            IF (l_obstout_s2d) then
               obst_s2d(:, i, j) = obst_output_ij%s2d(:)
            END IF
            IF (l_obstout_s3d) then
               obst_s3d(:, :, i, j) = obst_output_ij%s3d(:, :)
            END IF
            IF (l_obstout_drag) then
               obst_drag3d(:, :, i, j) = obst_output_ij%drag3d(:, :)
            END IF
         END DO
      END DO

#ifdef MPI /* exchange needed to retrieve fuz/fvz at cells edges */
      CALL exchange_r3d_tile(limin, limax, ljmin, ljmax, obst_fuz(START_2D_ARRAY, 1))
      CALL exchange_r3d_tile(limin, limax, ljmin, ljmax, obst_fvz(START_2D_ARRAY, 1))
#endif /* MPI */

   END SUBROUTINE obst_update

   !!==========================================================================================================

   SUBROUTINE obst_rd_timeserie(loc_time)
   !!---------------------------------------------------------------------
   !!                 *** SUBROUTINE OBSTRUCTIONS_readfile_timeserie  ***
   !!
   !! ** Purpose : Read file for time-varying obstructions characteristics
   !!
   !!---------------------------------------------------------------------

      IMPLICIT NONE
      REAL(KIND=rlg), INTENT(IN) :: loc_time

      REAL a, b, c
      INTEGER iv

      DO iv = 1, obst_nbvar
         IF (obst_l_filetimeserie(iv)) THEN
            ! update t_before and t_after given loc_time
            IF (loc_time > obst_ts_time(iv, obst_ts_tafter(iv))) THEN
               obst_ts_tbefore(iv) = obst_ts_tafter(iv)
               obst_ts_tafter(iv) = min(obst_ts_tafter(iv) + 1, obst_ts_tmax(iv))
            END IF
            IF (obst_ts_tbefore(iv) == obst_ts_tafter(iv)) THEN ! particular case
               obst_height_inst(iv, :, :) = obst_ts_height(iv, obst_ts_tbefore(iv))
               obst_dens_inst(iv, :, :) = obst_ts_dens(iv, obst_ts_tbefore(iv))
               obst_width_inst(iv, :, :) = obst_ts_width(iv, obst_ts_tbefore(iv))
               obst_thick_inst(iv, :, :) = obst_ts_thick(iv, obst_ts_tbefore(iv))
            else
               ! time linear interpolation
               a = real(loc_time - obst_ts_time(iv, obst_ts_tbefore(iv)))
               b = real(obst_ts_time(iv, obst_ts_tafter(iv)) - loc_time)
               c = a + b
               ! compute interpolation for each var height, dens, width, thick
               obst_height_inst(iv, :, :) = obst_ts_height(iv, obst_ts_tbefore(iv))*b/c &
                                            + obst_ts_height(iv, obst_ts_tafter(iv))*a/c
               obst_dens_inst(iv, :, :) = obst_ts_dens(iv, obst_ts_tbefore(iv))*b/c &
                                          + obst_ts_dens(iv, obst_ts_tafter(iv))*a/c
               obst_width_inst(iv, :, :) = obst_ts_width(iv, obst_ts_tbefore(iv))*b/c &
                                           + obst_ts_width(iv, obst_ts_tafter(iv))*a/c
               obst_thick_inst(iv, :, :) = obst_ts_thick(iv, obst_ts_tbefore(iv))*b/c &
                                           + obst_ts_thick(iv, obst_ts_tafter(iv))*a/c
            END IF

            ! application of position mask
            WHERE (obst_position(iv, :, :) == 0.)
               obst_height_inst(iv, :, :) = 0.
               obst_dens_inst(iv, :, :) = 0.
               obst_width_inst(iv, :, :) = 0.
               obst_thick_inst(iv, :, :) = 0.
            END WHERE

         END IF
      END DO

      ! *********************************
   END SUBROUTINE obst_rd_timeserie

#endif
END MODULE OBSTRUCTIONS
