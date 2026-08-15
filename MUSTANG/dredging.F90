! MUSTANG - This software is governed by the CeCILL-C license
! see LICENSE_MUSTANG.txt

MODULE dredging
   !============================================================================
   !!                   ***  MODULE  dredging  ***
   !!
   !! This module manage dredging effect on sediment.
   !! This is done by removing sediments in areas where a given sediment height
   !! is exceeded and optionnaly dumping them in other areas of the domain.
   !!
   !! This module is linked to modules : module_substance, 
   !! comsubstance, comMUSTANG
   !! This module also use netcdf
   !!
   !! This module has been adapted and generalized in 2024 by F. Grasso and
   !! S. Le Gac from the work of J.P. Lemoine in a MARS version for ARES
   !! hindcast in the SEINE estuary.
   !============================================================================

#include "cppdefs.h"

#if defined MUSTANG
   USE module_substance ! for dimension
   USE comsubstance, ONLY: lchain, rsh, rlg, riosh, nvp, isand1, isand2, &
                           igrav1, igrav2, surf_cell, l_subs2D, name_var, nv_tot
   USE comMUSTANG, ONLY: dredging_location_file, dredging_settings_file, &
                         dredging_out_file, dredging_dt, dredging_dt_out, &
                         dredging_dumping_layer

   IMPLICIT NONE
   PRIVATE

   PUBLIC dredging_init_param ! called by initMUSTANG
   PUBLIC dredging_main ! called by MUSTANG_update
   PUBLIC l_dredging ! boolean , TRUE to activate dredging
   PUBLIC dump_gravel_flx ! to put gravel directly in deposit (not in water)
   PUBLIC dredg_hsed_init ! to keep initial state info

   ! Shared module variables
   LOGICAL :: l_dredging ! boolean , TRUE to activate dredging
   REAL(KIND=rsh), DIMENSION(:, :, :), ALLOCATABLE :: dump_gravel_flx
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: dredg_hsed_init

   ! Local variables
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: dredg_mass_byclass
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: dredg_mass_byclass_byloc
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: dump_mass_byclass_byloc
   REAL(KIND=rsh), DIMENSION(:, :), ALLOCATABLE :: dredg_mass_byclass_byloc_cum
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: dump_surface
   INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: dredg_loc
   INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: dump_loc
   CHARACTER(lchain), DIMENSION(:), ALLOCATABLE :: dredg_name
   CHARACTER(lchain), DIMENSION(:), ALLOCATABLE :: dump_name
   INTEGER, DIMENSION(:, :), ALLOCATABLE :: dredg_flag
   INTEGER, DIMENSION(:), ALLOCATABLE :: dredg_dump_flag
   REAL(KIND=rsh), DIMENSION(:), ALLOCATABLE :: dredg_depth
   INTEGER :: n_dredg_loc
   INTEGER :: n_dump_loc
   REAL(KIND=rsh) :: dredging_time
   REAL(KIND=rsh) :: dredging_time_out
   INTEGER :: dredging_ncid
   INTEGER :: dredging_t_dimid, dredging_t_varid
   INTEGER :: dredging_next_record
   INTEGER :: dredging_area_dimid, dredging_area_varid
   INTEGER :: dredging_class_dimid, dredging_class_varid
   INTEGER :: dredging_mass_varid
   INTEGER, DIMENSION(:), ALLOCATABLE :: dump_layer

CONTAINS

   !============================================================================
   SUBROUTINE dredging_init_param(imin, imax, jmin, jmax)
      !!------------------------------------------------------------------------
      !!                 *** SUBROUTINE dredging_init_param  ***
      !!
      !! ** Purpose : Initialization of dredging parameters and variables
      !!
      !!------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax

      l_dredging = .TRUE.
      IF (TRIM(dredging_location_file) == "" .OR. &
          TRIM(dredging_settings_file) == "") THEN
         l_dredging = .FALSE.
      ELSE
         !! reading association dredging/dumping areas and dredging depths
         CALL dredging_get_dimension_from_settings_file
         IF (n_dredg_loc == 0) THEN
            l_dredging = .FALSE.
         ELSE
            !! check_param
            CALL dredging_check_param
            !! allocate all arrays
            CALL dredging_alloc
            !! reading association dredging/dumping areas and dredging depths
            CALL dredging_read_settings_file
            !! reading dredging_location and dumping_location
            CALL dredging_read_location_file
            !! initialize dump_surface and dredg_flag
            CALL dredging_init_var(imin, imax, jmin, jmax)
            !! initialize output file
            CALL dredging_def_output
            !! mpi case for dumping_surface variable
            CALL dumping_mpi_surface
         END IF ! test if there is at least one dredging area

      END IF ! test if there is dreddging
      RETURN

   END SUBROUTINE dredging_init_param
   !============================================================================

   SUBROUTINE dredging_get_dimension_from_settings_file
      !!------------------------------------------------------------------------
      !!    *** SUBROUTINE dredging_get_dimension_from_settings_file ***
      !!
      !! ** Purpose : Get n_dredg_loc and n_dump_loc values from file
      !!
      !!------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: i, ios, idump
      CHARACTER(LEN=1000) :: line
      CHARACTER(LEN=lchain) :: temp_dredging, temp_dumping
      CHARACTER(LEN=lchain), DIMENSION(:), ALLOCATABLE :: temp_dumping_name
      REAL(KIND=rsh) :: temp_depth
      LOGICAL :: new_dump_zone, line_ok

      OPEN (unit=51, file=dredging_settings_file, &
            status='old', action='read', iostat=ios)
      IF (ios /= 0) THEN
         WRITE (stdout, *) "Error opening dredging_settings_file : ", &
            dredging_settings_file
         STOP 1
      END IF

      READ (51, '(A)') ! reading header line
      ! retrieve n_dredg_loc
      n_dredg_loc = 0
      ios = 0
      line_ok = .TRUE.
      DO WHILE (line_ok)
         READ (51, '(A)', iostat=ios) line
         IF (ios == 0) THEN
            READ (line, *, iostat=ios) temp_dredging, temp_depth, temp_dumping
            IF (ios == 0) THEN
               ! If reading is successful (iostat is 0), count this line
               n_dredg_loc = n_dredg_loc + 1
            ELSE
               line_ok = .FALSE.
            END IF
         ELSE
            line_ok = .FALSE.
         END IF
      END DO

      REWIND (51)
      READ (51, '(A)') ! reading header line
      ! retrieve n_dump_loc
      n_dump_loc = 0
      ALLOCATE (temp_dumping_name(n_dredg_loc))
      temp_dumping_name(:) = "none"
      DO i = 1, n_dredg_loc
         READ (51, '(A)', iostat=ios) line
         READ (line, *, iostat=ios) temp_dredging, temp_depth, temp_dumping

         new_dump_zone = .true.
         IF (temp_dumping == "none") THEN
            new_dump_zone = .false.
         ELSE
            DO idump = 1, i - 1
               IF (TRIM(temp_dumping) == TRIM(temp_dumping_name(idump))) THEN
                  new_dump_zone = .false.
               END IF
            END DO
         END IF
         IF (new_dump_zone) THEN
            n_dump_loc = n_dump_loc + 1
         END IF
         temp_dumping_name(i) = temp_dumping
      END DO

      CLOSE (51)

      DEALLOCATE (temp_dumping_name)

   END SUBROUTINE dredging_get_dimension_from_settings_file
   !============================================================================

   SUBROUTINE dredging_read_settings_file
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dredging_read_settings_file  ***
      !!
      !! ** Purpose : Read text file containing dredged depths and
      !!              correspondance between dredging zones and dumping zones
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: i, ios, iz, idump
      CHARACTER(LEN=1000) :: line
      CHARACTER(LEN=lchain) :: temp_dredging, temp_dumping
      REAL(KIND=rsh) :: temp_depth
      LOGICAL :: new_dump_zone

      OPEN (unit=51, file=dredging_settings_file, &
            status='old', action='read', iostat=ios)
      IF (ios /= 0) THEN
         WRITE (stdout, *) "Error opening dredging_settings_file : ", &
            dredging_settings_file
         STOP 1
      END IF

      READ (51, '(A)') ! reading header line

      idump = 0
      dump_name(:) = ""
      DO i = 1, n_dredg_loc
         READ (51, '(A)', iostat=ios) line
         READ (line, *, iostat=ios) temp_dredging, temp_depth, temp_dumping

         dredg_name(i) = TRIM(temp_dredging)
         dredg_depth(i) = temp_depth

         ! finding index of corresponding dumping zone
         IF (temp_dumping == "none") THEN
            dredg_dump_flag(i) = 0
         ELSE
            new_dump_zone = .true.
            DO iz = 1, idump
               IF (TRIM(temp_dumping) == TRIM(dump_name(iz))) THEN
                  new_dump_zone = .false.
                  dredg_dump_flag(i) = iz
               END IF
            END DO
            IF (new_dump_zone) THEN
               idump = idump + 1
               dump_name(idump) = TRIM(temp_dumping)
               dredg_dump_flag(i) = idump
            END IF

         END IF
      END DO
      CLOSE (51)

   END SUBROUTINE dredging_read_settings_file
   !============================================================================

   SUBROUTINE dredging_read_location_file
      !!------------------------------------------------------------------------
      !!                 *** SUBROUTINE dredging_read_location_file ***
      !!
      !! ** Purpose : Reading netcdf file containing dradging and dumping
      !!              locations
      !!
      !!------------------------------------------------------------------------
      USE netcdf
      IMPLICIT NONE

      INTEGER :: iz, ncid, varid, status
      CHARACTER(lchain) :: name
      REAL(KIND=rsh) :: tmp(GLOBAL_2D_ARRAY)
      tmp(:, :) = 0.

      ! Open input NetCDF file.
      status = NF90_OPEN(dredging_location_file, NF90_NOWRITE, ncid)
      IF (status /= NF90_NOERR) THEN
         WRITE (stdout, *) "Error opening dredging location file"
         WRITE (stdout, *) TRIM(NF90_STRERROR(status))
         STOP 1
      END IF

      DO iz = 1, n_dredg_loc
         name = TRIM(dredg_name(iz))
         status = NF90_INQ_VARID(ncid, name, varid)
         IF (status /= NF90_NOERR) THEN
            WRITE (stdout, *) "Unable to retrieve dredg variable id: ", &
               iz, " ", name
            WRITE (stdout, *) TRIM(NF90_STRERROR(status))
            STOP 1
         END IF
         CALL dredging_ncget2D(ncid, varid, tmp)
         dredg_loc(iz, GLOBAL_2D_ARRAY) = INT(tmp(GLOBAL_2D_ARRAY))
      END DO

      DO iz = 1, n_dump_loc
         name = TRIM(dump_name(iz))
         IF (name == "none") THEN
            WRITE (stdout, *) "Invalid dump variable name: ", name
            WRITE (stdout, *) "none can not be used as a dump zone name"
            STOP 1
         END IF
         status = NF90_INQ_VARID(ncid, name, varid)
         IF (status /= NF90_NOERR) THEN
            WRITE (stdout, *) "Unable to retrieve dump variable id: ", &
               iz, " ", name
            WRITE (stdout, *) TRIM(NF90_STRERROR(status))
            STOP 1
         END IF
         CALL dredging_ncget2D(ncid, varid, tmp)
         dump_loc(iz, GLOBAL_2D_ARRAY) = INT(tmp(GLOBAL_2D_ARRAY))
      END DO

      ! Close input NetCDF file.
      status = NF90_CLOSE(ncid)
      IF (status /= NF90_NOERR) THEN
         WRITE (stdout, *) "Error closing dredging location file"
         WRITE (stdout, *) TRIM(NF90_STRERROR(status))
         STOP 1
      END IF

      RETURN

   END SUBROUTINE dredging_read_location_file
   !============================================================================

   SUBROUTINE dredging_init_var(imin, imax, jmin, jmax)
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dredging_init_var  ***
      !!
      !! ** Purpose :
      !! Initialize dredg_flag from all dredg_loc.
      !! In case of overlapping dredging zone, the deepest dredging is keeped
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax

      INTEGER :: i, j, iz, iv, tmp_flag
      LOGICAL :: warning_out

      dredging_time = time
      dredging_time_out = time
      warning_out = .FALSE.

      DO j = jmin, jmax
         DO i = imin, imax
            tmp_flag = 0
            DO iz = 1, n_dredg_loc
               IF (dredg_loc(iz, i, j) > 0) THEN
                  IF (tmp_flag == 0) THEN
                     tmp_flag = iz
                  ELSE IF (dredg_depth(iz) > dredg_depth(tmp_flag)) THEN
                     ! keeping the deepest
                     tmp_flag = iz
                  END IF
               END IF
            END DO
            dredg_flag(i, j) = tmp_flag

            DO iz = 1, n_dump_loc
               IF (dump_loc(iz, i, j) > 0) THEN
                  dump_surface(iz) = dump_surface(iz) + surf_cell(i, j)
               END IF
            END DO
         END DO !i
      END DO !j

      IF (warning_out) THEN
         WRITE (stdout, *) "WARNING - dredging"
         WRITE (stdout, *) "Overlapping of dredging zone"
         WRITE (stdout, *) "The deepest one is keeped"
      END IF

      dump_layer(:) = dredging_dumping_layer
#ifdef key_sand2D
      DO iv = isand1, isand2
         IF (l_subs2D(iv)) THEN
            dump_layer(iv) = 1
         END IF
      END DO
#endif

   END SUBROUTINE dredging_init_var
   !============================================================================

   SUBROUTINE dredging_check_param
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dredging_check_param  ***
      !!
      !! ** Purpose : Check values of param in namelist
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE
      IF (dredging_dumping_layer > N) THEN
         WRITE (stdout, *) "Error dredging_dumping_layer should be <= N"
         WRITE (stdout, *) "N = ", N
         WRITE (stdout, *) "dredging_dumping_layer = ", dredging_dumping_layer
         STOP 1
      END IF
      IF (dredging_dumping_layer < 1) THEN
         WRITE (stdout, *) "Error dredging_dumping_layer should be >= 1"
         WRITE (stdout, *) "dredging_dumping_layer = ", dredging_dumping_layer
         STOP 1
      END IF
      IF (dredging_dt <= 0.) THEN
         WRITE (stdout, *) "Error dredging_dt should be > 0."
         WRITE (stdout, *) "dredging_dt = ", dredging_dt
         STOP 1
      END IF
      IF (dredging_dt_out < dredging_dt) THEN
         WRITE (stdout, *) "Error dredging_dt_out should be >= dredging_dt"
         WRITE (stdout, *) "dredging_dt = ", dredging_dt
         WRITE (stdout, *) "dredging_dt_out = ", dredging_dt_out
         STOP 1
      END IF
      IF (dredging_dt_out <= 0.) THEN
         WRITE (stdout, *) "Error dredging_dt_out should be > 0."
         WRITE (stdout, *) "dredging_dt_out = ", dredging_dt_out
         STOP 1
      END IF

   END SUBROUTINE dredging_check_param
   !============================================================================

   SUBROUTINE dredging_alloc
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dredging_alloc ***
      !!
      !! ** Purpose : Allocate arrays for the dredging module
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE
      ALLOCATE (dredg_name(n_dredg_loc))
      ALLOCATE (dump_name(n_dump_loc))

      ALLOCATE (dredg_mass_byclass(nvp))
      ALLOCATE (dredg_mass_byclass_byloc(n_dredg_loc, nvp))
      ALLOCATE (dump_mass_byclass_byloc(n_dump_loc, nvp))
      ALLOCATE (dredg_mass_byclass_byloc_cum(n_dredg_loc, nvp))
      dredg_mass_byclass(:) = 0.0_rsh
      dredg_mass_byclass_byloc(:, :) = 0.0_rsh
      dump_mass_byclass_byloc(:, :) = 0.0_rsh
      dredg_mass_byclass_byloc_cum(:, :) = 0.0_rsh

      ALLOCATE (dredg_loc(n_dredg_loc, GLOBAL_2D_ARRAY))
      ALLOCATE (dump_loc(n_dump_loc, GLOBAL_2D_ARRAY))
      dredg_loc(:, :, :) = 0
      dump_loc(:, :, :) = 0

      ALLOCATE (dump_gravel_flx(igrav1:igrav2, GLOBAL_2D_ARRAY))
      dump_gravel_flx(:, :, :) = 0.

      ALLOCATE (dump_layer(nvp))
      dump_layer(:) = 1

      ALLOCATE (dredg_depth(n_dredg_loc))
      ALLOCATE (dredg_dump_flag(n_dredg_loc))
      ALLOCATE (dredg_flag(GLOBAL_2D_ARRAY))
      ALLOCATE (dredg_hsed_init(GLOBAL_2D_ARRAY))
      ALLOCATE (dump_surface(n_dump_loc))
      dredg_depth(:) = 0.0_rsh
      dredg_dump_flag(:) = 0
      dredg_flag(:, :) = 0
      dredg_hsed_init(:, :) = 0.0_rsh
      dump_surface(:) = 0.0_rsh

   END SUBROUTINE dredging_alloc

   !============================================================================

   SUBROUTINE dredging_main(imin, imax, jmin, jmax, watconc, z_w, hmod, hsed, &
                            dzs, ksmi, ksma, cv_sed, c_sedtot)
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dredging_main ***
      !!
      !! ** Purpose : Main subroutine called by MUSTANG_update, main steps:
      !!                    - Check time
      !!                    - Compute dredged mass
      !!                    - Compute dumping mass
      !!                    - Transfer dredged mass to dumping area
      !!                      (water and deposit)
      !!                    - Output
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax
      REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY, N, 3, NT), &
         INTENT(INOUT) :: watconc
      REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY, 0:N), INTENT(IN) ::  z_w
      REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY), INTENT(IN) :: hmod
      REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY), INTENT(INOUT) :: hsed
      REAL(KIND=rsh), DIMENSION(ksdmin:ksdmax, GLOBAL_2D_ARRAY), &
         INTENT(INOUT) :: dzs
      INTEGER, DIMENSION(GLOBAL_2D_ARRAY), INTENT(IN) :: ksmi
      INTEGER, DIMENSION(GLOBAL_2D_ARRAY), INTENT(INOUT) :: ksma
      REAL(KIND=rsh), DIMENSION(-1:nv_tot, ksdmin:ksdmax, GLOBAL_2D_ARRAY), &
         INTENT(INOUT) :: cv_sed
      REAL(KIND=rsh), DIMENSION(ksdmin:ksdmax, GLOBAL_2D_ARRAY), &
         INTENT(INOUT) :: c_sedtot

      IF (time .GE. dredging_time) THEN

         CALL dredging_compute_mass(imin, imax, jmin, jmax, hmod, hsed, &
                                    dzs, ksmi, ksma, cv_sed, c_sedtot)

         CALL dredging_mpi_mass

         CALL dumping_compute_mass

         dredg_mass_byclass_byloc_cum(:, :) = &
            dredg_mass_byclass_byloc_cum(:, :) + &
            dredg_mass_byclass_byloc(:, :)

         CALL dredging_dump_mass(watconc, z_w)

         CALL dredging_mpi_waterconcentration(imin, imax, jmin, jmax, watconc)

         IF (time .GE. dredging_time_out) THEN
            MPI_master_only CALL dredging_write_output
            ! update next time for output
            dredging_time_out = time + dredging_dt_out
         END IF
         ! update next time
         dredging_time = time + dredging_dt
      END IF

   END SUBROUTINE dredging_main
   !============================================================================

   SUBROUTINE dredging_compute_mass(imin, imax, jmin, jmax, hmod, hsed, &
                                    dzs, ksmi, ksma, cv_sed, c_sedtot)
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dredging_compute_mass ***
      !!
      !! ** Purpose : computation of dredged masses
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax
      REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY), INTENT(IN) :: hmod
      REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY), INTENT(INOUT) :: hsed
      REAL(KIND=rsh), DIMENSION(ksdmin:ksdmax, GLOBAL_2D_ARRAY), &
         INTENT(INOUT) :: dzs
      INTEGER, DIMENSION(GLOBAL_2D_ARRAY), INTENT(IN) :: ksmi
      INTEGER, DIMENSION(GLOBAL_2D_ARRAY), INTENT(INOUT) :: ksma
      REAL(KIND=rsh), DIMENSION(-1:nv_tot, ksdmin:ksdmax, GLOBAL_2D_ARRAY), &
         INTENT(INOUT) :: cv_sed
      REAL(KIND=rsh), DIMENSION(ksdmin:ksdmax, GLOBAL_2D_ARRAY), &
         INTENT(INOUT) :: c_sedtot

      INTEGER :: i, j, iv, iz, k
      REAL(KIND=rsh) :: dzs_dredged, tmp_depth

      dredg_mass_byclass_byloc(:, :) = 0.0_rsh
      DO j = jmin, jmax
         DO i = imin, imax
            hsed(i, j) = 0.0_rsh
            DO k = ksmi(i, j), ksma(i, j)
               hsed(i, j) = hsed(i, j) + dzs(k, i, j)
            END DO
            iz = dredg_flag(i, j)
            IF (iz > 0) THEN ! in dredging area
#if defined MORPHODYN
               ! retrieve initial depth, as hmod is updated during simulation
               tmp_depth = hmod(i, j) + hsed(i, j) - dredg_hsed_init(i, j)
#else
               ! initial depth
               tmp_depth = hmod(i, j)
#endif /* MORPHODYN */
               DO WHILE (((tmp_depth - (hsed(i, j) - dredg_hsed_init(i, j))) &
                          .LT. dredg_depth(iz)) .AND. ksma(i, j) > 0)

                  k = ksma(i, j)
                  IF ((tmp_depth - &
                       (hsed(i, j) - dzs(k, i, j) - dredg_hsed_init(i, j))) &
                      .LE. dredg_depth(iz)) THEN
                     ! whole layer can be dredged
                     dzs_dredged = dzs(k, i, j)
                     DO iv = 1, nvp
                        dredg_mass_byclass_byloc(iz, iv) = &
                           dredg_mass_byclass_byloc(iz, iv) + &
                           dzs_dredged*cv_sed(iv, k, i, j)*surf_cell(i, j)
                        cv_sed(iv, k, i, j) = 0.0_rsh
                     END DO
                     hsed(i, j) = hsed(i, j) - dzs_dredged
                     dzs(k, i, j) = 0.0_rsh
                     c_sedtot(k, i, j) = 0.0_rsh
                     ksma(i, j) = ksma(i, j) - 1
                  ELSE
                     ! dredged only partial layer
                     dzs_dredged = min(dzs(k, i, j), &
                                       dredg_depth(iz) - &
                                       (tmp_depth - &
                                        (hsed(i, j) - dredg_hsed_init(i, j))))
                     DO iv = 1, nvp
                        dredg_mass_byclass_byloc(iz, iv) = &
                           dredg_mass_byclass_byloc(iz, iv) + &
                           dzs_dredged*cv_sed(iv, k, i, j)*surf_cell(i, j)
                     END DO
                     hsed(i, j) = hsed(i, j) - dzs_dredged
                     dzs(k, i, j) = dzs(k, i, j) - dzs_dredged
                  END IF
               END DO
            END IF
         END DO
      END DO

   END SUBROUTINE dredging_compute_mass
   !============================================================================

   SUBROUTINE dumping_compute_mass
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dumping_compute_mass ***
      !!
      !! ** Purpose : Computation of dumping masses
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: iz, idump

      dump_mass_byclass_byloc(:, :) = 0.0_rsh
      DO iz = 1, n_dredg_loc
         DO idump = 1, n_dump_loc
            IF (dredg_dump_flag(iz) == idump) THEN
               dump_mass_byclass_byloc(idump, :) = &
                  dredg_mass_byclass_byloc(iz, :)
            END IF
         END DO
      END DO

   END SUBROUTINE dumping_compute_mass

   !============================================================================

   SUBROUTINE dredging_dump_mass(watconc, z_w)
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dredging_dump_mass ***
      !!
      !! ** Purpose : Transfer of dredged masses to dumping areas
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE
      REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY, N, 3, NT), &
         INTENT(INOUT) :: watconc
      REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY, 0:N), INTENT(IN) :: z_w
      INTEGER :: i, j, idump, iv

      dump_gravel_flx(:, :, :) = 0.0_rsh

      DO idump = 1, n_dump_loc
         IF (sum(dump_mass_byclass_byloc(idump, :)) > 0.0_rsh) THEN
            ! if there is mass to dump
            DO iv = isand1, nvp ! gravels are treated separately
               watconc(:, :, dump_layer(iv), nstp, itsubs1 - 1 + iv) = &
                  watconc(:, :, dump_layer(iv), nstp, itsubs1 - 1 + iv) + &
                  REAL(dump_loc(idump, :, :), rsh) &
                  *dump_mass_byclass_byloc(idump, iv) &
                  /max(1e-30_rsh, dump_surface(idump)) &
                  /max(1e-30_rsh, z_w(:, :, dredging_dumping_layer) &
                       - z_w(:, :, dredging_dumping_layer - 1))
            END DO

            DO iv = igrav1, igrav2
               dump_gravel_flx(iv, :, :) = dump_gravel_flx(iv, :, :) &
                                           + REAL(dump_loc(idump, :, :), rsh) &
                                           *dump_mass_byclass_byloc(idump, iv) &
                                           /dump_surface(idump)
            END DO

         END IF
      END DO

   END SUBROUTINE dredging_dump_mass
   !============================================================================

   SUBROUTINE dumping_mpi_surface
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dumping_mpi_surface ***
      !!
      !! ** Purpose : MPI treatment of dredged and dumped mass
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE

# ifdef MPI
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE), blank, ierr
#  ifdef XIOS
#include "mpi_cpl.h"
#  endif /* XIOS */

      INTEGER :: iz, ierror
      REAL(KIND=rsh) :: tmp

      DO iz = 1, n_dump_loc
         CALL MPI_ALLREDUCE(dump_surface(iz), &
                            tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                            MPI_COMM_WORLD, ierror)
         dump_surface(iz) = tmp
      END DO

#endif /*MPI*/

   END SUBROUTINE dumping_mpi_surface
   !============================================================================

   SUBROUTINE dredging_mpi_mass
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dredging_mpi_mass ***
      !!
      !! ** Purpose : MPI treatment of dredged and dumped mass
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE

# ifdef MPI
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE), blank, ierr
#  ifdef XIOS
#include "mpi_cpl.h"
#  endif /* XIOS */

      INTEGER :: iv, iz, ierror
      REAL(KIND=rsh) :: tmp

      DO iv = 1, nvp
         DO iz = 1, n_dredg_loc
            CALL MPI_ALLREDUCE(dredg_mass_byclass_byloc(iz, iv), &
                               tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               MPI_COMM_WORLD, ierror)
            dredg_mass_byclass_byloc(iz, iv) = tmp
         END DO
      END DO

#endif /*MPI*/

   END SUBROUTINE dredging_mpi_mass
   !============================================================================

   SUBROUTINE dredging_mpi_waterconcentration(Istr, Iend, Jstr, Jend, watconc)
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dredging_mpi_waterconcentration ***
      !!
      !! ** Purpose : MPI exchange of water concentration
      !!
      !!------------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY, N, 3, NT), &
         INTENT(INOUT) :: watconc
      INTEGER, INTENT(IN) :: Istr, Iend, Jstr, Jend

# ifdef MPI
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE), blank, ierr
#  ifdef XIOS
#include "mpi_cpl.h"
#  endif

      INTEGER :: itrc

      do itrc = 1, NT
         call exchange_r3d_3pts_tile(Istr, Iend, Jstr, Jend, &
                                     watconc(START_2D_ARRAY, 1, nstp, itrc))
      end do
#endif /*MPI*/

   END SUBROUTINE dredging_mpi_waterconcentration
   !============================================================================

   SUBROUTINE dredging_def_output
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dredging_def_output ***
      !!
      !! ** Purpose : Define netcdf file to output cummulative dredging mass
      !!
      !!------------------------------------------------------------------------

      USE netcdf

      IMPLICIT NONE
      INTEGER :: j, lchain_dimid

      INTEGER :: name_area_varid
      REAL(KIND=rlg), DIMENSION(:), ALLOCATABLE :: area_var
      CHARACTER(LEN=lchain), DIMENSION(:), ALLOCATABLE :: area_name
      INTEGER, DIMENSION(2) :: dimids2_area

      INTEGER :: name_class_varid
      REAL(KIND=rlg), DIMENSION(:), ALLOCATABLE :: class_var
      CHARACTER(LEN=lchain), DIMENSION(:), ALLOCATABLE :: class_name
      INTEGER, DIMENSION(2) :: dimids2_class

      INTEGER, DIMENSION(3) :: dimids3_mass_t

#ifdef MPI
      IF (mynode .eq. 0) THEN
#endif
         CALL dredging_check( &
            nf90_create(dredging_out_file, NF90_SHARE, dredging_ncid))

         CALL dredging_check( &
            nf90_def_dim(dredging_ncid, 'lchain', lchain, lchain_dimid))
         CALL dredging_check( &
            nf90_def_dim(dredging_ncid, 'time', NF90_UNLIMITED, &
                         dredging_t_dimid))
         CALL dredging_check( &
            nf90_def_var(dredging_ncid, "time", NF90_DOUBLE, &
                         dredging_t_dimid, dredging_t_varid))
         CALL dredging_check( &
            nf90_put_att(dredging_ncid, dredging_t_varid, &
                         "units", "seconds since 1900-01-01"))

         CALL dredging_check( &
            nf90_def_dim(dredging_ncid, 'area', &
                         n_dredg_loc, dredging_area_dimid))
         CALL dredging_check( &
            nf90_def_var(dredging_ncid, "area", NF90_DOUBLE, &
                         dredging_area_dimid, dredging_area_varid))
         dimids2_area = (/lchain_dimid, dredging_area_dimid/)
         CALL dredging_check( &
            nf90_def_var(dredging_ncid, 'area_name', NF90_CHAR, &
                         dimids2_area, name_area_varid))

         CALL dredging_check( &
            nf90_def_dim(dredging_ncid, 'class', &
                         nvp, dredging_class_dimid))
         CALL dredging_check( &
            nf90_def_var(dredging_ncid, "class", NF90_DOUBLE, &
                         dredging_class_dimid, dredging_class_varid))
         dimids2_class = (/lchain_dimid, dredging_class_dimid/)
         CALL dredging_check( &
            nf90_def_var(dredging_ncid, 'class_name', NF90_CHAR, &
                         dimids2_class, name_class_varid))

         dimids3_mass_t = (/dredging_area_dimid, dredging_class_dimid, &
                            dredging_t_dimid/)
         CALL dredging_check( &
            nf90_def_var(dredging_ncid, 'cummulative_mass', NF90_DOUBLE, &
                         dimids3_mass_t, dredging_mass_varid))
         CALL dredging_check( &
            nf90_put_att(dredging_ncid, dredging_mass_varid, &
                         "description", "Cummulative mass"))

         CALL dredging_check(nf90_enddef(dredging_ncid))
         CALL dredging_check(nf90_sync(dredging_ncid))

         ALLOCATE (area_var(n_dredg_loc))
         ALLOCATE (area_name(n_dredg_loc))
         DO j = 1, n_dredg_loc
            area_var(j) = j
            area_name(j) = dredg_name(j)
         END DO
         CALL dredging_check( &
            nf90_put_var(dredging_ncid, dredging_area_varid, area_var(:), &
                         start=(/1/), count=(/n_dredg_loc/)))
         CALL dredging_check( &
            nf90_put_var(dredging_ncid, name_area_varid, area_name, &
                         start=(/1, 1/), count=(/lchain, n_dredg_loc/)))

         CALL dredging_check(nf90_sync(dredging_ncid))

         IF (nvp .gt. 0) THEN
            ALLOCATE (class_var(nvp))
            ALLOCATE (class_name(nvp))
            DO j = 1, nvp
               class_var(j) = j
               class_name(j) = name_var(j)
            end do
            CALL dredging_check( &
               nf90_put_var(dredging_ncid, dredging_class_varid, class_var(:), &
                            start=(/1/), count=(/nvp/)))
            CALL dredging_check( &
               nf90_put_var(dredging_ncid, name_class_varid, &
                            class_name, &
                            start=(/1, 1/), count=(/lchain, nvp/)))
         END IF

         CALL dredging_check(nf90_sync(dredging_ncid))

         dredging_next_record = 1

#ifdef MPI
      END IF
#endif

   END SUBROUTINE dredging_def_output
   !============================================================================

   SUBROUTINE dredging_write_output
      !!------------------------------------------------------------------------
      !!       *** SUBROUTINE dredging_write_output ***
      !!
      !! ** Purpose : Output cummulative dredging mass
      !!
      !!------------------------------------------------------------------------

      USE netcdf

      IMPLICIT NONE
      integer :: iv, iz
      integer, dimension(3) :: start

      ! write time
      call dredging_check( &
         nf90_put_var(dredging_ncid, dredging_t_varid, &
                      time, (/dredging_next_record/)))

      ! write open line
      do iz = 1, n_dredg_loc
         do iv = 1, nvp
            start = (/iz, iv, dredging_next_record/)
            call dredging_check( &
               nf90_put_var(dredging_ncid, &
                            dredging_mass_varid, &
                            dredg_mass_byclass_byloc_cum(iz, iv), start))
         end do
      end do

      call dredging_check(nf90_sync(dredging_ncid))

      dredging_next_record = dredging_next_record + 1

   END SUBROUTINE dredging_write_output
   !============================================================================

   SUBROUTINE dredging_check(status)
      !!------------------------------------------------------------------------
      !!                 *** SUBROUTINE dredging_check  ***
      !!
      !! ** Purpose : Check netcdf function
      !!------------------------------------------------------------------------

      USE netcdf
      INTEGER, INTENT(IN) :: status

      IF (status /= nf90_noerr) THEN
         WRITE (stdout, *) 'dredging_check(): '
         WRITE (stdout, *) TRIM(nf90_strerror(status))
         STOP 1
      END IF

   END SUBROUTINE dredging_check

   SUBROUTINE dredging_ncget2D(ncid, varid, tmp)
      !!------------------------------------------------------------------------
      !!                 *** SUBROUTINE dredging_ncget2D  ***
      !!
      !! ** Purpose : Get var from netcdf
      !!              (from nf_fread.F translate in f90)
      !!------------------------------------------------------------------------
      USE netcdf
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ncid, varid
      REAL(KIND=rsh), INTENT(INOUT) ::  tmp(GLOBAL_2D_ARRAY)

      INTEGER :: imin, imax, jmin, jmax, start(2), count(2), status

      jmin = 0
      imin = 0
      start(1) = 1
      start(2) = 1

#ifdef MPI
      IF (ii .gt. 0) THEN
         start(1) = 1 - imin + iminmpi
         imin = 1
      END IF
      IF (ii .eq. NP_XI - 1) THEN
         imax = Lmmpi + 1
      ELSE
         imax = Lmmpi
      END IF
      IF (jj .gt. 0) THEN
         start(2) = 1 - jmin + jminmpi
         jmin = 1
      END IF
      IF (jj .eq. NP_ETA - 1) THEN
         jmax = Mmmpi + 1
      ELSE
         jmax = Mmmpi
      END IF
#else
      imax = Lm + 1
      jmax = Mm + 1
#endif

      count(1) = imax - imin + 1
      count(2) = jmax - jmin + 1

      status = NF90_GET_VAR(ncid, varid, tmp(imin:imax, jmin:jmax), &
                            start, count)
      IF (status /= NF90_NOERR) THEN
         WRITE (stdout, *) "Error NF90_GET_VAR in dredging"
         WRITE (stdout, *) TRIM(NF90_STRERROR(status))
         STOP 1
      END IF

#ifdef MPI
# define LOCALLM Lmmpi
# define LOCALMM Mmmpi
#else
# define LOCALLM Lm
# define LOCALMM Mm
#endif /* MPI */
#if defined EW_PERIODIC || defined NS_PERIODIC  || defined MPI
   !!/* exchange needed */
      CALL exchange_r2d_tile(1, LOCALLM, 1, LOCALMM, tmp)
#endif /* MPI */

   END SUBROUTINE dredging_ncget2D

#endif /* ifdef MUSTANG */
   !============================================================================
END MODULE dredging
