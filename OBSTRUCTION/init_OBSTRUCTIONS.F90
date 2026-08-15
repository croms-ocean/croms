MODULE init_OBSTRUCTIONS

#include "cppdefs.h"

#ifdef OBSTRUCTION
   !!==========================================================================
   !!                   ***  MODULE  init_OBSTRUCTIONS  ***
   !!
   !!
   !! ** Purpose : concerns all subroutines related to obstruction module
   !!              initialization
   !!
   !! ** Description :
   !!
   !!     subroutine obst_init   ! Defines obstructions outputs parameters
   !!                            ! and initialize some variables
   !!
   !!     subroutine obst_vname  ! Defines output vname vector for def_his and
   !!                            ! wrt_his
   !!
   !!     subroutine obst_readvar ! Read the file (variables_obstructions.dat)
   !!
   !!     subroutine obst_write_summary ! Write a summary of options
   !!
   !!     subroutine obst_compatibility ! Check compatibility of options/params
   !!
   !!     subroutine obst_nccheck ! For netcdf error logging
   !!
   !!     subroutine obst_ncget2D ! Read 2D variables in netcdf
   !!
   !!     subroutine obst_init_pos ! Reads input file for obstruction position
   !!
   !!     subroutine obst_init_spatial ! Read the spatial initialization file
   !!                                  ! or apply constant values
   !!
   !!     subroutine obst_init_timeseries ! Read timeseries files and
   !!                                     ! initialize obst_ts_* variables
   !!
   !!     subroutine obst_readfile_distri    ! Read input file for vertical
   !!                                        ! distribution
   !!
   !!     subroutine obst_alloc_nbvar  ! Allocates variables depending on
   !!                                  ! number of obstructions variables
   !!
   !!     subroutine obst_alloc_xyz ! Allocates spatial tables for obstructions
   !!
   !!==========================================================================
   !! * Modules used
   USE com_OBSTRUCTIONS
   USE OBSTRUCTIONS1DV

   IMPLICIT NONE

   ! function & routines of this module, called outside :
   PUBLIC obst_init

   PRIVATE

CONTAINS

   !!==========================================================================
   SUBROUTINE obst_init(h0fond_in)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_init  ***
      !!
      !! ** Purpose : Initialization of obstruction parameters
      !!
      !!---------------------------------------------------------------------

      USE module_OBSTRUCTIONS ! for GLOBAL_2D_ARRAY, N, stdout, mynode
      IMPLICIT NONE

      REAL(KIND=rsh), INTENT(IN) :: h0fond_in

      INTEGER                    :: iv, k, lstr, lenstr

      NAMELIST /obst_main/ obst_nbvar, obst_c_paramhuv
      NAMELIST /obst_input/ obst_fn_position, obst_fn_var
      NAMELIST /obst_output/ l_obstout_pos, l_obstout_height_f, l_obstout_height_e, &
         l_obstout_dens_f, l_obstout_dens_e, l_obstout_width_f, l_obstout_width_e, &
         l_obstout_thick_f, l_obstout_thick_e, &
         l_obstout_theta, l_obstout_frac_xy, l_obstout_frac_z, &
         l_obstout_fuzvz, l_obstout_a2d, l_obstout_a3d, l_obstout_s2d, &
         l_obstout_s3d, l_obstout_drag, l_obstout_tau

      obst_kmax = N
      iscreenlog = stdout
      ierrorlog = stdout
      iwarnlog = stdout

      ! ************************
      ! * READING NAMELIST
      ! ************************
      lstr = lenstr(obstname)
      OPEN (50, file=obstname(1:lstr), status='old', form='formatted', access='sequential')
      READ (50, obst_main)
      ALLOCATE (obst_fn_var(1:obst_nbvar))
      READ (50, obst_input)
      READ (50, obst_output)
      CLOSE (50)
      MPI_master_only WRITE (iscreenlog, *) &
         ' '
      MPI_master_only WRITE (iscreenlog, *) &
         '*****************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         '***** module OBSTRUCTIONS, subroutine OBST_INIT *****'
      MPI_master_only WRITE (iscreenlog, *) &
         '*****************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         ' Reading file ', TRIM(obstname)

      ! ************************************
      ! * ALLOCATION OF VARIABLES PARAMETERS
      ! ************************************
      CALL obst_alloc_nbvar

      ! ***********************
      ! * READING VARIABLES
      ! ***********************
      DO iv = 1, obst_nbvar
         CALL obst_readvar(iv)
      END DO
      ! **********************
      ! * TABLES ALLOCATION
      ! ***********************
      CALL obst_alloc_xyz

      ! **********************
      ! * CHECKS AND SUMMARY
      ! ***********************
      CALL obst_compatibility
      CALL obst_write_summary

      ! *******************************
      ! * INITIALIZATION OF OUTPUT
      ! *******************************
      CALL obst_vname

      ! ***********************
      ! * READING POSITION FILE
      ! ***********************
      CALL obst_init_pos

      ! ********************************************************
      ! * READING SPATIAL INIT FILE or init from constant values
      ! ********************************************************
      CALL obst_init_spatial

      ! ***************************
      ! * READING DISTRIBUTION FILE
      ! ***************************
      CALL obst_readfile_distri

      ! *******************************
      ! * INITIALIZATION OF TIME SERIES
      ! *******************************
      CALL obst_init_timeseries

      ! ***********************************
      ! * INITIALIZATION OF ROUGHNESS SEDIM
      ! ***********************************
      obst_l_z0bstress_tot = .FALSE.
      DO iv = 1, obst_nbvar
         IF (obst_l_z0bstress(iv)) THEN
            obst_l_z0bstress_tot = .TRUE. ! Only one variable used z0sed
         END IF
      END DO

      ! ***********************
      ! * OTHER INITIALIZATIONS
      ! ***********************
      obst_z0bed(:, :) = zob(:, :) ! Saving the surface z0 tables from hydraulical initialization (without obstructions)

      CALL o1dv_alloc(obst_nbvar, obst_kmax, obst_nb_max_hnorm)
      CALL o1dv_init(h0fond_in, obst_c_paramhuv, &
                     obst_varname, obst_type, obst_l_flexible, obst_l_abdelposture, &
                     obst_l_param_height, obst_l_cylinder, obst_l_noturb, obst_l_drag_cste, &
                     obst_l_abdelrough_cste, obst_l_fracxy, obst_l_z0bstress, &
                     obst_fracxy_type, obst_c_abdel_nmax, obst_z0bstress_option, &
                     obst_c_rho, obst_c_height_x0, obst_c_height_x1, obst_c_shelter, &
                     obst_c_lift, obst_c_drag, obst_c_lz, obst_c_crough_x0, obst_c_crough_x1, &
                     obst_c_fracxy_k0, obst_c_fracxy_k1, obst_c_fracxy_l, obst_c_z0bstress, &
                     obst_c_z0bstress_x0, obst_c_z0bstress_x1, &
                     obst_l_filedistri, obst_nbhnorm, obst_height_norm, obst_dens_norm, &
                     stdout)

   END SUBROUTINE obst_init

   !!==========================================================================================================

   SUBROUTINE obst_vname
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_vname  ***
      !!
      !! ** Purpose : initialize output vname and variables
      !!
      !!---------------------------------------------------------------------
      USE module_OBSTRUCTIONS
      IMPLICIT NONE

      INTEGER :: iv, indvar

      obst_nout_pos = 'pos'       ! Name obstruction position
      obst_nout_height_f = 'height_f'  ! Name 2D obstruction height (forcing) (iv,i,j)
      obst_nout_height_e = 'height_e'  ! Name 2D obstruction height (effective) (iv,i,j)
      obst_nout_dens_f = 'dens_f'    ! Name 2D obstruction density (forcing) (iv,i,j)
      obst_nout_dens_e = 'dens_e'    ! Name 3D obstruction density (3D effective) (iv,k,i,j)
      obst_nout_width_f = 'width_f'   ! Name 2D obstruction width (forcing) (iv,i,j)
      obst_nout_width_e = 'width_e'   ! Name 3D obstruction width (3D effective) (iv,k,i,j)
      obst_nout_thick_f = 'thick_f'   ! Name 2D obstruction thick (forcing) (iv,i,j)
      obst_nout_thick_e = 'thick_e'   ! Name 3D obstruction thick (3D effective) (iv,k,i,j)
      obst_nout_theta = 'theta'     ! Name 3D obstruction bending angle (iv,k,i,j)
      obst_nout_frac_xy = 'frac_xy'   ! Name 2D obstruction fragmentation correction (iv,i,j)
      obst_nout_frac_z = 'frac_z'    ! Name 2D obstruction fraction of sigma layer (iv,i,j)
      obst_nout_fuzvz = 'fuzvz'     ! Name 3D obstruction friction force (k,i,j)
      obst_nout_a2d = 'a2d'       ! Name 2D obstruction horizontal area (iv+3,i,j)
      obst_nout_a3d = 'a3d'       ! Name 3D obstruction horizontal area (iv+3,k,i,j)
      obst_nout_s2d = 's2d'       ! Name 2D obstruction vertical area (iv+3,i,j)
      obst_nout_s3d = 's3d'       ! Name 3D obstruction vertical area (iv+3,k,i,j)
      obst_nout_drag = 'cd3d'      ! Name 3D obstruction drag coefficient (iv,k,i,j)
      obst_nout_tau = 'tau3d'     ! Name 3D obstruction turbulence dissipation scale (k,i,j)

      ALLOCATE (hisObst(1:17*obst_nbvar + 15))
      ALLOCATE (avgObst(1:17*obst_nbvar + 15))
      ALLOCATE (outObst(1:17*obst_nbvar + 15))
      ALLOCATE (out2DObst(1:17*obst_nbvar + 15))
      ALLOCATE (vname_obst(20, 1:17*obst_nbvar + 15))

      vname_obst(:,:) = ""
      outObst(:) = .FALSE.
      out2DObst(:) = .FALSE.

      DO iv = 1, obst_nbvar
         indvar = iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_pos)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction occupation rate for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = '-                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_pos) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .TRUE.
         END IF

         indvar = obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_height_f)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction forcing height for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = 'm                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_height_f) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .TRUE.
         END IF

         indvar = 2*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_height_e)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction effective height for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = 'm                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_height_e) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .TRUE.
         END IF

         indvar = 3*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_dens_f)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction forcing density for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = 'm-2                                  '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_dens_f) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .TRUE.
         END IF

         indvar = 4*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_dens_e)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction effective density for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = 'm-2                                  '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_dens_e) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .FALSE.
         END IF

         indvar = 5*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_width_f)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction forcing width for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = 'm                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_width_f) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .TRUE.
         END IF

         indvar = 6*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_width_e)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction effective width for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = 'm                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_width_e) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .FALSE.
         END IF

         indvar = 7*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_thick_f)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction forcing thickness for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = 'm                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_thick_f) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .TRUE.
         END IF

         indvar = 8*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_thick_e)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction effective thickness for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = 'm                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_thick_e) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .FALSE.
         END IF

         indvar = 9*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_theta)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction bending angle for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = 'deg                                  '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_theta) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .FALSE.
         END IF

         indvar = 10*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_frac_xy)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction fragmentation correction factor for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = '-                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_frac_xy) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .TRUE.
         END IF

         indvar = 11*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_frac_z)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction sigma fraction for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = '-                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_frac_z) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .FALSE.
         END IF

         indvar = 12*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_drag)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            'Obstruction drag coefficient for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = '-                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_drag) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .FALSE.
         END IF

         indvar = 13*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_a2d)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            '2D Obstruction horizontal area for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = '-                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_a2d) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .TRUE.
         END IF

         indvar = 14*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_a3d)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            '3D Obstruction horizontal area for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = '-                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_a3d) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .FALSE.
         END IF

         indvar = 15*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_s2d)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            '2D Obstruction vertical area for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = '-                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_s2d) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .TRUE.
         END IF

         indvar = 16*obst_nbvar + iv
         vname_obst(1, indvar) = &
            TRIM(obst_nout_s3d)//'_'//TRIM(obst_varname(iv))
         vname_obst(2, indvar) = &
            '3D Obstruction vertical area for '//TRIM(obst_varname(iv))
         vname_obst(3, indvar) = '-                                    '
         vname_obst(4, indvar) = '                                     '
         vname_obst(5, indvar) = '                                     '
         vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
         vname_obst(7, indvar) = '                                     '
         IF (l_obstout_s3d) THEN
            outObst(indvar) = .TRUE.
            out2DObst(indvar) = .FALSE.
         END IF
      END DO

      indvar = 17*obst_nbvar + 1
      vname_obst(1, indvar) = &
         TRIM(obst_nout_a2d)//'_NoTurb'
      vname_obst(2, indvar) = &
         '2D Obstruction horizontal area for NoTurb variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_a2d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .TRUE.
      END IF

      indvar = 17*obst_nbvar + 2
      vname_obst(1, indvar) = &
         TRIM(obst_nout_a2d)//'_Turb'
      vname_obst(2, indvar) = &
         '2D Obstruction horizontal area for Turb variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_a2d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .TRUE.
      END IF

      indvar = 17*obst_nbvar + 3
      vname_obst(1, indvar) = &
         TRIM(obst_nout_a2d)//'_All'
      vname_obst(2, indvar) = &
         '2D Obstruction horizontal area for All variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_a2d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .TRUE.
      END IF

      indvar = 17*obst_nbvar + 4
      vname_obst(1, indvar) = &
         TRIM(obst_nout_a3d)//'_NoTurb'
      vname_obst(2, indvar) = &
         '3D Obstruction horizontal area for NoTurb variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_a3d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .FALSE.
      END IF

      indvar = 17*obst_nbvar + 5
      vname_obst(1, indvar) = &
         TRIM(obst_nout_a3d)//'_Turb'
      vname_obst(2, indvar) = &
         '3D Obstruction horizontal area for Turb variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_a3d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .FALSE.
      END IF

      indvar = 17*obst_nbvar + 6
      vname_obst(1, indvar) = &
         TRIM(obst_nout_a3d)//'_All'
      vname_obst(2, indvar) = &
         '3D Obstruction horizontal area for All variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_a3d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .FALSE.
      END IF

      indvar = 17*obst_nbvar + 7
      vname_obst(1, indvar) = &
         TRIM(obst_nout_s2d)//'_NoTurb'
      vname_obst(2, indvar) = &
         '2D Obstruction vertical area for NoTurb variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_s2d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .TRUE.
      END IF

      indvar = 17*obst_nbvar + 8
      vname_obst(1, indvar) = &
         TRIM(obst_nout_s2d)//'_Turb'
      vname_obst(2, indvar) = &
         '2D Obstruction vertical area for Turb variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_s2d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .TRUE.
      END IF

      indvar = 17*obst_nbvar + 9
      vname_obst(1, indvar) = &
         TRIM(obst_nout_s2d)//'_All'
      vname_obst(2, indvar) = &
         '2D Obstruction vertical area for All variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time lat_rho lon_rho                 '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_s2d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .TRUE.
      END IF

      indvar = 17*obst_nbvar + 10
      vname_obst(1, indvar) = &
         TRIM(obst_nout_s3d)//'_NoTurb'
      vname_obst(2, indvar) = &
         '3D Obstruction vertical area for NoTurb variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_s3d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .FALSE.
      END IF

      indvar = 17*obst_nbvar + 11
      vname_obst(1, indvar) = &
         TRIM(obst_nout_s3d)//'_Turb'
      vname_obst(2, indvar) = &
         '3D Obstruction vertical area for Turb variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_s3d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .FALSE.
      END IF

      indvar = 17*obst_nbvar + 12
      vname_obst(1, indvar) = &
         TRIM(obst_nout_s3d)//'_All'
      vname_obst(2, indvar) = &
         '3D Obstruction vertical area for All variables'
      vname_obst(3, indvar) = '-                                    '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_s3d) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .FALSE.
      END IF

      indvar = 17*obst_nbvar + 13
      vname_obst(1, indvar) = &
         TRIM(obst_nout_tau)
      vname_obst(2, indvar) = &
         'Obstruction 3D turbulent dissipation'
      vname_obst(3, indvar) = 'N.m-2                                '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_tau) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .FALSE.
      END IF

      indvar = 17*obst_nbvar + 14
      vname_obst(1, indvar) = &
         TRIM(obst_nout_fuzvz)//'_uz'
      vname_obst(2, indvar) = &
         'Obstruction 3D friction force FUZ'
      vname_obst(3, indvar) = 'N.m-2                                '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_fuzvz) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .FALSE.
      END IF

      indvar = 17*obst_nbvar + 15
      vname_obst(1, indvar) = &
         TRIM(obst_nout_fuzvz)//'_vz'
      vname_obst(2, indvar) = &
         'Obstruction 3D friction force FVZ'
      vname_obst(3, indvar) = 'N.m-2                                '
      vname_obst(4, indvar) = '                                     '
      vname_obst(5, indvar) = '                                     '
      vname_obst(6, indvar) = 'time N lat_rho lon_rho               '
      vname_obst(7, indvar) = '                                     '
      IF (l_obstout_fuzvz) THEN
         outObst(indvar) = .TRUE.
         out2DObst(indvar) = .FALSE.
      END IF

   END SUBROUTINE obst_vname

   !!==========================================================================================================

   SUBROUTINE obst_readvar(iv)
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_readvar  ***
      !!
      !! ** Purpose : read "obstruction_variables.dat" file describing obstructions
      !!
      !!---------------------------------------------------------------------
      USE module_OBSTRUCTIONS ! for mynode
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iv

      CHARACTER(len=lchain) :: filepc
      ! For obst_var_main
      CHARACTER(len=lchain) :: r_obst_varname
      CHARACTER(LEN=2)      :: r_obst_type
      LOGICAL               :: r_l_obst_cylinder
      ! For obst_var_option
      LOGICAL               :: r_l_obst_flexible, r_l_obst_noturb
      LOGICAL               :: r_l_obst_filetimeserie
      CHARACTER(len=lchain) :: r_obst_fn_timeserie
      ! For obst_var_init
      LOGICAL               :: r_l_obst_init_spatial, r_l_obst_filedistri
      CHARACTER(len=lchain) :: r_obst_fn_initspatial, r_obst_fn_distrib
      REAL(KIND=rsh)        :: r_obst_i_height, r_obst_i_width, r_obst_i_thick, r_obst_i_dens
      ! For obst_var_flexibility
      LOGICAL               :: r_l_obst_abdelposture, r_l_obst_param_height
      INTEGER               :: r_obst_c_abdel_nmax
      REAL(KIND=rsh)        :: r_obst_c_rho, r_obst_c_lift, r_obst_c_shelter, r_obst_c_height_x0, r_obst_c_height_x1
      ! For obst_var_roughdrag
      LOGICAL               :: r_l_obst_drag_cste, r_l_obst_abdelrough_cste
      REAL(KIND=rsh)        :: r_obst_c_crough_x0, r_obst_c_crough_x1, r_obst_c_drag, r_obst_c_lz
      ! For obst_var_fracxy
      LOGICAL               :: r_l_obst_fracxy
      INTEGER               :: r_obst_fracxy_type
      REAL(KIND=rsh)        :: r_obst_c_fracxy_k0, r_obst_c_fracxy_k1, r_obst_c_fracxy_l
      ! For obst_var_bstress
      LOGICAL               :: r_l_obst_z0bstress
      INTEGER               :: r_obst_z0bstress_option
      REAL(KIND=rsh)        :: r_obst_c_z0bstress, r_obst_c_z0bstress_x0, r_obst_c_z0bstress_x1
      !! * Namelists
      NAMELIST /obst_var_main/ r_obst_varname, r_obst_type, r_l_obst_cylinder
      NAMELIST /obst_var_option/ r_l_obst_flexible, r_l_obst_noturb, &
         r_l_obst_filetimeserie, r_obst_fn_timeserie, r_l_obst_filedistri, r_obst_fn_distrib
      NAMELIST /obst_var_init/ r_l_obst_init_spatial, &
         r_obst_fn_initspatial, &
         r_obst_i_height, r_obst_i_width, r_obst_i_thick, r_obst_i_dens
      NAMELIST /obst_var_flexibility/ r_l_obst_abdelposture, r_l_obst_param_height, &
         r_obst_c_abdel_nmax, r_obst_c_rho, r_obst_c_lift, &
         r_obst_c_shelter, r_obst_c_height_x0, r_obst_c_height_x1
      NAMELIST /obst_var_roughdrag/ r_l_obst_drag_cste, r_l_obst_abdelrough_cste, &
         r_obst_c_crough_x0, r_obst_c_crough_x1, r_obst_c_drag, r_obst_c_lz
      NAMELIST /obst_var_fracxy/ r_l_obst_fracxy, r_obst_fracxy_type, r_obst_c_fracxy_k0, &
         r_obst_c_fracxy_k1, r_obst_c_fracxy_l
      NAMELIST /obst_var_bstress/ r_l_obst_z0bstress, r_obst_z0bstress_option, r_obst_c_z0bstress, &
         r_obst_c_z0bstress_x0, r_obst_c_z0bstress_x1
      !!----------------------------------------------------------------------
      !! * Executable part

      ! save into simu.log
      !-------------------
      MPI_master_only WRITE (iscreenlog, *) &
         ' '
      MPI_master_only WRITE (iscreenlog, *) &
         '****************************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         '***** module OBSTRUCTIONS, subroutine obst_READVAR *****'
      MPI_master_only WRITE (iscreenlog, *) &
         '****************************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         ' Reading file ', TRIM(obst_fn_var(iv))
      MPI_master_only WRITE (iscreenlog, *) &
         ' defining obstructions parameters'
      ! ******************************
      ! * Start reading variables file
      ! ******************************
      filepc = obst_fn_var(iv)
      OPEN (55, file=filepc, status='old', form='formatted', access='sequential')
      READ (55, obst_var_main)
      READ (55, obst_var_option)
      READ (55, obst_var_init)
      READ (55, obst_var_flexibility)
      READ (55, obst_var_roughdrag)
      READ (55, obst_var_fracxy)
      READ (55, obst_var_bstress)
      CLOSE (55)
      ! ************************************
      ! * Alocate to corresponding parameter
      ! ************************************
      ! * For namelist obst_var_main
      obst_varname(iv) = r_obst_varname
      obst_l_cylinder(iv) = r_l_obst_cylinder
      obst_type(iv) = r_obst_type
      ! * For namelist obst_var_option
      obst_l_flexible(iv) = r_l_obst_flexible
      obst_l_noturb(iv) = r_l_obst_noturb
      obst_l_filetimeserie(iv) = r_l_obst_filetimeserie
      obst_fn_timeserie(iv) = r_obst_fn_timeserie
      ! * For namelist obst_var_init
      obst_l_init_spatial(iv) = r_l_obst_init_spatial
      obst_l_filedistri(iv) = r_l_obst_filedistri
      obst_fn_initspatial(iv) = r_obst_fn_initspatial
      obst_fn_distrib(iv) = r_obst_fn_distrib
      obst_i_height(iv) = r_obst_i_height
      obst_i_width(iv) = r_obst_i_width
      obst_i_thick(iv) = r_obst_i_thick
      obst_i_dens(iv) = r_obst_i_dens
      ! * For namelist obst_var_flexibility
      obst_l_abdelposture(iv) = r_l_obst_abdelposture
      obst_l_param_height(iv) = r_l_obst_param_height
      obst_c_abdel_nmax(iv) = r_obst_c_abdel_nmax
      obst_c_rho(iv) = r_obst_c_rho
      obst_c_lift(iv) = r_obst_c_lift
      obst_c_shelter(iv) = r_obst_c_shelter
      obst_c_height_x0(iv) = r_obst_c_height_x0
      obst_c_height_x1(iv) = r_obst_c_height_x1
      ! * For namelist obst_var_roughdrag
      obst_l_drag_cste(iv) = r_l_obst_drag_cste
      obst_l_abdelrough_cste(iv) = r_l_obst_abdelrough_cste
      obst_c_crough_x0(iv) = r_obst_c_crough_x0
      obst_c_crough_x1(iv) = r_obst_c_crough_x1
      obst_c_drag(iv) = r_obst_c_drag
      obst_c_lz(iv) = r_obst_c_lz
      ! * For namelist obst_var_fracxy
      obst_l_fracxy(iv) = r_l_obst_fracxy
      obst_fracxy_type(iv) = r_obst_fracxy_type
      obst_c_fracxy_k0(iv) = r_obst_c_fracxy_k0
      obst_c_fracxy_k1(iv) = r_obst_c_fracxy_k1
      obst_c_fracxy_l(iv) = r_obst_c_fracxy_l
      ! * For namelist obst_var_bstress
      obst_l_z0bstress(iv) = r_l_obst_z0bstress
      obst_z0bstress_option(iv) = r_obst_z0bstress_option
      obst_c_z0bstress(iv) = r_obst_c_z0bstress
      obst_c_z0bstress_x0(iv) = r_obst_c_z0bstress_x0
      obst_c_z0bstress_x1(iv) = r_obst_c_z0bstress_x1
      !!**********************************
   END SUBROUTINE obst_readvar

   !!==========================================================================================================

   SUBROUTINE obst_write_summary
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_write_summary  ***
      !!
      !! ** Purpose : Write a summary of obstructions parameters
      !!
      !!---------------------------------------------------------------------
      USE module_OBSTRUCTIONS ! for mynode
      IMPLICIT NONE

      INTEGER          :: iv

      MPI_master_only WRITE (iscreenlog, *) &
         ' '
      MPI_master_only WRITE (iscreenlog, *) &
         '***********************************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         '********************** module OBSTRUCTIONS ****************************'
      MPI_master_only WRITE (iscreenlog, *) &
         '************** subroutine obst_WRITE_SUMMARY ******************'
      MPI_master_only WRITE (iscreenlog, *) &
         '***********************************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         'LISTING OF OBSTRUCTION VARIABLES :'
      MPI_master_only WRITE (iscreenlog, *) &
         '------------------------------------------------'
      MPI_master_only WRITE (iscreenlog, *) &
         'TOTAL NUMBER OF OBSTRUCTION VARIABLES : ', obst_nbvar
      MPI_master_only WRITE (iscreenlog, *) &
         'Number of RIGID_UP                    : ', obst_nv_rigid_up
      MPI_master_only WRITE (iscreenlog, *) &
         'Number of RIGID_DO                    : ', obst_nv_rigid_do
      MPI_master_only WRITE (iscreenlog, *) &
         'Number of FLEXI_UP                    : ', obst_nv_flexi_up
      MPI_master_only WRITE (iscreenlog, *) &
         'Number of FLEXI_DO                    : ', obst_nv_flexi_do
      MPI_master_only WRITE (iscreenlog, *) &
         'Number of 3DVARS                      : ', obst_nv_3d
      MPI_master_only WRITE (iscreenlog, *) &
         '------------------------------------------------'
      MPI_master_only WRITE (iscreenlog, *) &
         'File for obstruction position is : ', TRIM(obst_fn_position)
      DO iv = 1, obst_nbvar
         MPI_master_only WRITE (iscreenlog, *) &
            '***********************************************************************************'
         MPI_master_only WRITE (iscreenlog, *) &
            '!===========================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '! NAMELIST : obst_var_main  !'
         MPI_master_only WRITE (iscreenlog, *) &
            '!===========================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '  - Name (identifier) of the variable         : ', obst_varname(iv)
         MPI_master_only WRITE (iscreenlog, *) &
            '  - If the current variable is cylinder-like  : ', obst_l_cylinder(iv)
         MPI_master_only WRITE (iscreenlog, *) &
            '  - Current variable itype                    : ', obst_type(iv)
         MPI_master_only WRITE (iscreenlog, *) &
            '!=============================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '! NAMELIST : obst_var_option  !'
         MPI_master_only WRITE (iscreenlog, *) &
            '!=============================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '  - If the current variable is flexible                     : ', obst_l_flexible(iv)
         MPI_master_only WRITE (iscreenlog, *) &
            '  - If the current variable is considered a macro-roughness : ', obst_l_noturb(iv)
         MPI_master_only WRITE (iscreenlog, *) &
            '  - Use a time-series file for obstructions characteristics : ', obst_l_filetimeserie(iv)
         IF (obst_l_filetimeserie(iv)) THEN
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Name of temporal file for obstructions characteristics      : ', obst_fn_timeserie(iv)
         ELSE
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Name of temporal file for obstructions characteristics      : NOT USED'
         END IF
         MPI_master_only WRITE (iscreenlog, *) &
            '!===========================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '! NAMELIST : obst_var_init  !'
         MPI_master_only WRITE (iscreenlog, *) &
            '!===========================!'

         MPI_master_only WRITE (iscreenlog, *) &
            '  - Use a spatially variable file obstructions characteristics             : ', obst_l_init_spatial(iv)
         MPI_master_only WRITE (iscreenlog, *) &
            '  - Use a file describing the vertical distribution of obstruction density : ', obst_l_filedistri(iv)
         IF (obst_l_init_spatial(iv)) THEN
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Name of spatial file for obstructions characteristics      : ', obst_fn_initspatial(iv)
         ELSE
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Name of spatial file for obstructions characteristics      : NOT USED'
         END IF
         IF (obst_l_filedistri(iv)) THEN
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Name of file for the vertical distribution of obstruction density      : ', obst_fn_distrib(iv)
         ELSE
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Name of file for the vertical distribution of obstruction density      : NOT USED'
         END IF
         IF (obst_l_init_spatial(iv)) THEN
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Initial height (unbent, eg. leaf-length for segrasses) of obstructions : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Initial width (or diameter for cylindric obstructions) of obstructions : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Initial thick of obstructions (along the flow)                         : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Initial density of obstructions (or maximum density)                   : NOT USED'
         ELSE
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Initial height (unbent, eg. leaf-length for segrasses) of obstructions : ', obst_i_height(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Initial width (or diameter for cylindric obstructions) of obstructions : ', obst_i_width(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Initial thick of obstructions (along the flow)                         : ', obst_i_thick(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Initial density of obstructions (or maximum density)                   : ', obst_i_dens(iv)
         END IF
         MPI_master_only WRITE (iscreenlog, *) &
            '!==================================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '! NAMELIST : obst_var_flexibility  !'
         MPI_master_only WRITE (iscreenlog, *) &
            '!==================================!'
         IF (obst_l_flexible(iv)) THEN
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Use Abdhelhrmans (2007) procedure to compute bending                   : ', obst_l_abdelposture(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Use empirical (exponential desrease) formulation to compute bending    : ', obst_l_param_height(iv)
            IF (obst_l_abdelposture(iv)) THEN
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Number of segments for Abdhelhrmans (2007) procedure                   : ', obst_c_abdel_nmax(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : ', obst_c_rho(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : ', obst_c_lift(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : ', obst_c_shelter(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - First parameter for empirical formulation                              : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Second parameter for empirical formulation                             : NOT USED'
            ELSEIF (obst_l_param_height(iv)) THEN
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Number of segments for Abdhelhrmans (2007) procedure                   : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - First parameter for empirical formulation                              : ', obst_c_height_x0(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Second parameter for empirical formulation                             : ', obst_c_height_x1(iv)
            ELSE
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Number of segments for Abdhelhrmans (2007) procedure                   : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - First parameter for empirical formulation                              : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Second parameter for empirical formulation                             : NOT USED'
            END IF
         ELSE
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Use Abdhelhrmans (2007) procedure to compute bending                   : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Use empirical (exponential desrease) formulation to compute bending    : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Number of segments for Abdhelhrmans (2007) procedure                   : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Volumic mass of obstructions for Abdhelhrmans (2007) procedure         : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Lift coefficient for Abdhelhrmans (2007) procedure                     : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Sheltering coefficient Afor bdhelhrmans (2007) procedure               : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - First parameter for empirical formulation                              : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Second parameter for empirical formulation                             : NOT USED'
         END IF
         MPI_master_only WRITE (iscreenlog, *) &
            '!================================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '! NAMELIST : obst_var_roughdrag  !'
         MPI_master_only WRITE (iscreenlog, *) &
            '!================================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '  - Use a constant drag coefficient for obstructions in hydrodynamics      : ', obst_l_drag_cste(iv)
         MPI_master_only WRITE (iscreenlog, *) &
            '  - Use a constant drag during reconstruction of velocity profiles         : ', obst_l_abdelrough_cste(iv)
         IF (obst_l_abdelrough_cste(iv)) THEN
            MPI_master_only WRITE (iscreenlog, *) &
               '  - First coefficient for CD during velocity profiles reconstruction       : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Second coefficient for CD during velocity profiles reconstruction      : NOT USED'
         ELSE
            MPI_master_only WRITE (iscreenlog, *) &
               '  - First coefficient for CD during velocity profiles reconstruction       : ', obst_c_crough_x0(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Second coefficient for CD during velocity profiles reconstruction      : ', obst_c_crough_x1(iv)
         END IF
         MPI_master_only WRITE (iscreenlog, *) &
            '  - Drag coefficient (max value if not constant) for obstructions elements : ', obst_c_drag(iv)
         IF (obst_l_noturb(iv)) THEN
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Coef. turbulent dissipation time-scale between obstructions elements   : NOT USED'
         ELSE
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Coef. turbulent dissipation time-scale between obstructions elements   : ', obst_c_lz(iv)
         END IF
         MPI_master_only WRITE (iscreenlog, *) &
            '!=============================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '! NAMELIST : obst_var_fracxy  !'
         MPI_master_only WRITE (iscreenlog, *) &
            '!=============================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '  - Take account for non-linear patchiness correction                      : ', obst_l_fracxy(iv)
         IF (.NOT. obst_l_fracxy(iv)) THEN
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Kind of non-linear correction method                                   : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Coefficient for the exponential correction                             : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - First parameter for correction of the exponential coefficient          : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Second parameter for correction of the exponential coefficient         : NOT USED'
         ELSE
            IF (obst_fracxy_type(iv) .EQ. 1) THEN
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Kind of non-linear correction method                                   : Exponential'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Coefficient for the exponential correction                             : ', obst_c_fracxy_k0(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - First parameter for correction of the exponential coefficient          : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Second parameter for correction of the exponential coefficient         : NOT USED'
            ELSEIF (obst_fracxy_type(iv) .EQ. 2) THEN
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Kind of non-linear correction method                                   : Exponential K0-variable'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Coefficient for the exponential correction                             : ', obst_c_fracxy_k0(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - First parameter for correction of the exponential coefficient          : ', obst_c_fracxy_k1(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Second parameter for correction of the exponential coefficient         : ', obst_c_fracxy_l(iv)
            END IF
         END IF
         MPI_master_only WRITE (iscreenlog, *) &
            '!==============================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '! NAMELIST : obst_var_bstress  !'
         MPI_master_only WRITE (iscreenlog, *) &
            '!==============================!'
         MPI_master_only WRITE (iscreenlog, *) &
            '  - To activate the impact of obstruction on Z0 (BSS computation)          : ', obst_l_z0bstress(iv)
         IF (obst_l_z0bstress(iv)) THEN
            IF (obst_z0bstress_option(iv) .EQ. 0) THEN
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Option to compute the obstruction induced roughness length             : Constant'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Constant (corrected value of roughness length)                         : ', obst_c_z0bstress(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - First parameter for rouhgness length computation (in 3D)               : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Second parameter for rouhgness length computation (in 3D)              : NOT USED'
            ELSEIF (obst_z0bstress_option(iv) .EQ. 1) THEN
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Option to compute the obstruction induced roughness length             : Parameterized'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Constant (corrected value of roughness length)                         : NOT USED'
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - First parameter for rouhgness length computation (in 3D)               : ', obst_c_z0bstress_x0(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  '  - Second parameter for rouhgness length computation (in 3D)              : ', obst_c_z0bstress_x1(iv)
            END IF
         ELSE
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Option to compute the obstruction induced roughness length             : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Constant (corrected value of roughness length)                         : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - First parameter for rouhgness length computation (in 3D)               : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Second parameter for rouhgness length computation (in 3D)              : NOT USED'
            MPI_master_only WRITE (iscreenlog, *) &
               '  - Coefficient to correct 3D roughness length into 2D roughness length    : NOT USED'
         END IF
      END DO
      !-------------------------------------------
   END SUBROUTINE obst_write_summary

   !!==========================================================================================================

   SUBROUTINE obst_compatibility
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_compatibility ***
      !!
      !! ** Purpose : Check compatibility of parameterization used
      !!
      !!---------------------------------------------------------------------
      USE module_OBSTRUCTIONS ! for mynode
      IMPLICIT NONE

      INTEGER          :: iv, nv_tot
      LOGICAL          :: l_turb

      ! ********************************************
      ! COUNT NUMBER OF VARIABLES OF DIFFERENT TYPES
      ! ********************************************
      obst_nv_up = 0
      obst_nv_DO = 0
      obst_nv_3d = 0
      obst_nv_rigid_up = 0
      obst_nv_rigid_do = 0
      obst_nv_flexi_up = 0
      obst_nv_flexi_do = 0
      obst_nv_turb = 0
      obst_nv_noturb = 0
      DO iv = 1, obst_nbvar
         IF (obst_l_noturb(iv)) THEN
            obst_nv_noturb = obst_nv_noturb + 1
         ELSE
            obst_nv_turb = obst_nv_turb + 1
         END IF
         IF (obst_type(iv) == "3D") THEN
            obst_nv_3d = obst_nv_3d + 1
         ELSE
            IF (obst_type(iv) == "DO") THEN
               obst_nv_DO = obst_nv_DO + 1
               IF (obst_l_flexible(iv)) THEN
                  obst_nv_flexi_do = obst_nv_flexi_do + 1
               ELSE
                  obst_nv_rigid_do = obst_nv_rigid_do + 1
               END IF
            ELSE
               obst_nv_up = obst_nv_up + 1
               IF (obst_l_flexible(iv)) THEN
                  obst_nv_flexi_up = obst_nv_flexi_up + 1
               ELSE
                  obst_nv_rigid_up = obst_nv_rigid_up + 1
               END IF
            END IF
         END IF
      END DO
      ! *************************************
      ! TESTS AND WRITE ON SIMULOG WARNLOG...
      ! *************************************
      !--------------------------------
      ! TESTING THE NUMBER OF VARIABLES
      !--------------------------------
      nv_tot = obst_nv_rigid_up + obst_nv_rigid_do + obst_nv_flexi_up + obst_nv_flexi_do + obst_nv_3d
      IF (nv_tot /= obst_nbvar) THEN
         MPI_master_only WRITE (ierrorlog, *) &
            ' '
         MPI_master_only WRITE (ierrorlog, *) &
            '**********************************************************************'
         MPI_master_only WRITE (ierrorlog, *) &
            '********************** module OBSTRUCTIONS ***************************'
         MPI_master_only WRITE (ierrorlog, *) &
            '************* subroutine obst_compatibility ******************'
         MPI_master_only WRITE (ierrorlog, *) &
            ' ERROR : The total number of obstruction variables read from'
         MPI_master_only WRITE (ierrorlog, *) &
            '         file', TRIM(obst_fn_var(iv))
         MPI_master_only WRITE (ierrorlog, *) &
            '         is DIFFERENT from the number of variables defined'
         MPI_master_only WRITE (ierrorlog, *) &
            '         nbvar', obst_nbvar
         MPI_master_only WRITE (ierrorlog, *) &
            '         nv_tot', nv_tot
         MPI_master_only WRITE (ierrorlog, *) &
            ' --> THE SIMULATION IS STOPPED !!! '
         STOP
      END IF

      !---------------------------
      ! TEST FOR 2D and GLS_KEPSILON
      !---------------------------
#ifndef SOLVE3D
      MPI_master_only WRITE (ierrorlog, *) &
         ' '
      MPI_master_only WRITE (ierrorlog, *) &
         '**********************************************************************'
      MPI_master_only WRITE (ierrorlog, *) &
         '********************** module OBSTRUCTIONS ***************************'
      MPI_master_only WRITE (ierrorlog, *) &
         '************* subroutine obst_compatibility ******************'
      MPI_master_only WRITE (ierrorlog, *) &
         ' ERROR : OBSTRUCTION only available with SOLVE3D'
      MPI_master_only WRITE (ierrorlog, *) &
         ' --> THE SIMULATION IS STOPPED !!!'
      STOP
#endif
      !---------------------------------------------------
      ! TEST FOR VARIABLES USING FULL TURBULENCE PROCEDURE
      !---------------------------------------------------
#ifndef GLS_KEPSILON
      MPI_master_only WRITE (ierrorlog, *) &
         ' '
      MPI_master_only WRITE (ierrorlog, *) &
         '**********************************************************************'
      MPI_master_only WRITE (ierrorlog, *) &
         '************************* module OBSTRUCTIONS ************************'
      MPI_master_only WRITE (ierrorlog, *) &
         '************* subroutine obst_compatibility ******************'
      MPI_master_only WRITE (ierrorlog, *) &
         ' ERROR : OBSTRUCTION only available with GLS_KEPSILON'
      MPI_master_only WRITE (ierrorlog, *) &
         ' --> THE SIMULATION IS STOPPED !!!'
      STOP
#endif

      !--------------------------------
      ! TESTING THE NUMBER OF VARIABLES
      !--------------------------------
      IF (obst_nbvar > 6) THEN
         MPI_master_only WRITE (ierrorlog, *) &
            ' '
         MPI_master_only WRITE (ierrorlog, *) &
            '**********************************************************************'
         MPI_master_only WRITE (ierrorlog, *) &
            '********************** module OBSTRUCTIONS ***************************'
         MPI_master_only WRITE (ierrorlog, *) &
            '************* subroutine obst_compatibility ******************'
         MPI_master_only WRITE (ierrorlog, *) &
            ' ERROR : The total number of obstruction variables '
         MPI_master_only WRITE (ierrorlog, *) &
            '         is > 6 which is the maximum allowed for output reasons (see hisObst) '
         MPI_master_only WRITE (ierrorlog, *) &
            '         nbvar = ', obst_nbvar
         MPI_master_only WRITE (ierrorlog, *) &
            ' --> THE SIMULATION IS STOPPED !!! '
         STOP
      END IF

      !--------------------------------------------------------------
      ! TEST CONSISTENCY FOR SCHEME AND DESCRIPTION FOR EACH VARIALBE
      !--------------------------------------------------------------
      DO iv = 1, obst_nbvar
         IF ((obst_type(iv) == "DO") .AND. (obst_l_noturb(iv))) THEN
            MPI_master_only WRITE (ierrorlog, *) &
               ' '
            MPI_master_only WRITE (ierrorlog, *) &
               '**********************************************************************'
            MPI_master_only WRITE (ierrorlog, *) &
               '********************** module OBSTRUCTIONS ***************************'
            MPI_master_only WRITE (ierrorlog, *) &
               '************* subroutine obst_compatibility ******************'
            MPI_master_only WRITE (ierrorlog, *) &
               ' ERROR : Inconsistency for variable :', obst_varname(iv)
            MPI_master_only WRITE (ierrorlog, *) &
               '         This variable is defined as a downward one, while'
            MPI_master_only WRITE (ierrorlog, *) &
               '         it uses simplified (obst_l_noturb) procedure'
            MPI_master_only WRITE (ierrorlog, *) &
               ' --> THE SIMULATION IS STOPPED !!!'
            STOP
         END IF
         IF ((obst_type(iv) == "3D") .AND. (obst_l_noturb(iv))) THEN
            MPI_master_only WRITE (ierrorlog, *) &
               ' '
            MPI_master_only WRITE (ierrorlog, *) &
               '**********************************************************************'
            MPI_master_only WRITE (ierrorlog, *) &
               '************************* module OBSTRUCTIONS ************************'
            MPI_master_only WRITE (ierrorlog, *) &
               '************* subroutine obst_compatibility ******************'
            MPI_master_only WRITE (ierrorlog, *) &
               ' ERROR : Inconsistency for variable :', obst_varname(iv)
            MPI_master_only WRITE (ierrorlog, *) &
               '         This variable is defined as a full 3d one, while'
            MPI_master_only WRITE (ierrorlog, *) &
               '         it uses simplified (obst_l_noturb) procedure'
            MPI_master_only WRITE (ierrorlog, *) &
               ' --> THE SIMULATION IS STOPPED !!!'
            STOP
         END IF
         IF ((obst_type(iv) == "3D") .AND. (obst_l_flexible(iv))) THEN
            MPI_master_only WRITE (ierrorlog, *) &
               ' '
            MPI_master_only WRITE (ierrorlog, *) &
               '**********************************************************************'
            MPI_master_only WRITE (ierrorlog, *) &
               '********************** module OBSTRUCTIONS ***************************'
            MPI_master_only WRITE (ierrorlog, *) &
               '************* subroutine obst_compatibility ******************'
            MPI_master_only WRITE (ierrorlog, *) &
               ' ERROR : Inconsistency for variable :', obst_varname(iv)
            MPI_master_only WRITE (ierrorlog, *) &
               '         This variable is defined as a full 3d one, while'
            MPI_master_only WRITE (ierrorlog, *) &
               '         it is also defined as a flexible one'
            MPI_master_only WRITE (ierrorlog, *) &
               '         This kind of variable is not available yet...'
            MPI_master_only WRITE (ierrorlog, *) &
               ' --> THE SIMULATION IS STOPPED !!!'
            STOP
         END IF

         IF ((obst_type(iv) == "DO") .AND. (obst_l_flexible(iv)) .AND. (obst_l_abdelposture(iv))) THEN
            MPI_master_only WRITE (ierrorlog, *) &
               ' '
            MPI_master_only WRITE (ierrorlog, *) &
               '**********************************************************************'
            MPI_master_only WRITE (ierrorlog, *) &
               '********************** module OBSTRUCTIONS ***************************'
            MPI_master_only WRITE (ierrorlog, *) &
               '************* subroutine obst_compatibility ******************'
            MPI_master_only WRITE (ierrorlog, *) &
               ' ERROR : Inconsistency for variable :', obst_varname(iv)
            MPI_master_only WRITE (ierrorlog, *) &
               '         This variable is defined as a flexible and downward one'
            MPI_master_only WRITE (ierrorlog, *) &
               '         and uses Abdelhrman procedure'
            MPI_master_only WRITE (ierrorlog, *) &
               '         This procedure is not available yet for downward variables...'
            MPI_master_only WRITE (ierrorlog, *) &
               ' --> THE SIMULATION IS STOPPED !!!'
            STOP
         END IF
         IF (obst_l_init_spatial(iv) .AND. obst_l_filetimeserie(iv)) THEN
            MPI_master_only WRITE (ierrorlog, *) &
               ' '
            MPI_master_only WRITE (ierrorlog, *) &
               '**********************************************************************'
            MPI_master_only WRITE (ierrorlog, *) &
               '********************** module OBSTRUCTIONS ***************************'
            MPI_master_only WRITE (ierrorlog, *) &
               '************* subroutine obst_compatibility ******************'
            MPI_master_only WRITE (ierrorlog, *) &
               ' ERROR : Inconsistency for variable :', obst_varname(iv)
            MPI_master_only WRITE (ierrorlog, *) &
               '         Spatial initialization is defined with also temporal forcing'
            MPI_master_only WRITE (ierrorlog, *) &
               '         for obstruction characteristics : This is not ready yet...'
            MPI_master_only WRITE (ierrorlog, *) &
               '         Choose between : spatial initialization (but no time-varying characteristics)'
            MPI_master_only WRITE (ierrorlog, *) &
               '         and homogeneous initialization (but with time-varying characteristics)'
            STOP
         END IF
      END DO
      !-------------------------------------------
   END SUBROUTINE obst_compatibility

   !!==========================================================================

   SUBROUTINE obst_nccheck(status, msg)
      !!----------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_nccheck  ***
      !!
      !! ** Purpose : Purpose : check netcdf function
      !!----------------------------------------------------------------------
      USE netcdf
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: status
      CHARACTER*(*), INTENT(IN) :: msg

      IF (status /= NF90_NOERR) THEN
         WRITE (ierrorlog, *) 'obstruction_nc_check(): '
         WRITE (ierrorlog, *) msg
         WRITE (ierrorlog, *) trim(NF90_STRERROR(status))
         STOP
      END IF

   END SUBROUTINE obst_nccheck

   !!==========================================================================

   SUBROUTINE obst_ncget2D(ncid, varid, tmp, msg)
      !!----------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_nccheck  ***
      !!
      !! ** Purpose : Purpose : get 2D from netcdf (from nf_fread.F translate in f90)
      !!----------------------------------------------------------------------
      USE module_OBSTRUCTIONS ! for GLOBAL_2D_ARRAY
      USE netcdf
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ncid, varid
      REAL, INTENT(INOUT) ::  tmp(GLOBAL_2D_ARRAY)
      CHARACTER*(*), INTENT(IN) :: msg

      INTEGER :: imin, imax, jmin, jmax, start(3), count(3)

      jmin = 0
      imin = 0
      start(1) = 1
      start(2) = 1
      start(3) = 1

#ifdef MPI
      if (ii .gt. 0) then
         start(1) = 1 - imin + iminmpi
         imin = 1
      end if
      if (ii .eq. NP_XI - 1) then
         imax = Lmmpi + 1
      else
         imax = Lmmpi
      end if
      if (jj .gt. 0) then
         start(2) = 1 - jmin + jminmpi
         jmin = 1
      end if
      if (jj .eq. NP_ETA - 1) then
         jmax = Mmmpi + 1
      else
         jmax = Mmmpi
      end if
#else
      imax = Lm + 1
      jmax = Mm + 1
#endif

      count(1) = imax - imin + 1
      count(2) = jmax - jmin + 1
      count(3) = 1

      CALL obst_nccheck(NF90_GET_VAR(ncid, varid, tmp(imin:imax, jmin:jmax), start, count), msg)
#ifdef MPI
# define LOCALLM Lmmpi
# define LOCALMM Mmmpi
#else
# define LOCALLM Lm
# define LOCALMM Mm
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC  || defined MPI /* exchange needed */
      CALL exchange_r2d_tile(1, LOCALLM, 1, LOCALMM, tmp)
#endif /* MPI */

   END SUBROUTINE obst_ncget2D
   !!==========================================================================================================

   SUBROUTINE obst_init_pos
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_init_pos  ***
      !!
      !! ** Purpose : read the obstruction position file
      !!
      !!---------------------------------------------------------------------
      USE module_OBSTRUCTIONS ! for GLOBAL_2D_ARRAY
      USE netcdf
      IMPLICIT NONE

      INTEGER iv, ncid, varid
      CHARACTER(30) name
      REAL tmp(GLOBAL_2D_ARRAY)

      CALL obst_nccheck(NF90_OPEN(obst_fn_position, NF90_NOWRITE, ncid), &
                        "obst_init_pos nf90_open "//obst_fn_position)

      DO iv = 1, obst_nbvar
         name = TRIM(obst_nout_pos)//'_'//TRIM(obst_varname(iv))
         CALL obst_nccheck(NF90_INQ_VARID(ncid, name, varid), &
                           "obst_init_pos nf90_inq_varid "//name)
         CALL obst_ncget2D(ncid, varid, tmp, &
                           "obst_init_pos getvar "//name)
         obst_position(iv, GLOBAL_2D_ARRAY) = tmp(GLOBAL_2D_ARRAY)
      END DO

      ! Close input NetCDF file.
      CALL obst_nccheck(NF90_CLOSE(ncid), &
                        "obst_init_pos nf90_close "//obst_fn_position)

      RETURN

   END SUBROUTINE obst_init_pos

   !!==========================================================================================================

   SUBROUTINE obst_init_spatial
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_init_spatial  ***
      !!
      !! ** Purpose : read the obstruction spatial initialization files
      !!              or apply constant values
      !!
      !!---------------------------------------------------------------------
      USE module_OBSTRUCTIONS ! for GLOBAL_2D_ARRAY
      USE netcdf
      IMPLICIT NONE

      INTEGER iv, ncid, varid
      CHARACTER(30) name
      REAL tmp(GLOBAL_2D_ARRAY)

      DO iv = 1, obst_nbvar
         IF (obst_l_init_spatial(iv)) THEN
            CALL obst_nccheck(NF90_OPEN(obst_fn_initspatial(iv), NF90_NOWRITE, ncid), &
                              "obst_init_spatial nf90_open "//obst_fn_initspatial(iv))

            name = TRIM(obst_nout_height_f)//'_'//TRIM(obst_varname(iv))
            CALL obst_nccheck(NF90_INQ_VARID(ncid, name, varid), &
                              "obst_init_spatial nf90_inq_varid "//obst_fn_initspatial(iv)//" "//name)
            CALL obst_ncget2D(ncid, varid, tmp, &
                              "obst_init_spatial getvar "//obst_fn_initspatial(iv)//" "//name)
            obst_height_inst(iv, GLOBAL_2D_ARRAY) = tmp(GLOBAL_2D_ARRAY)

            name = TRIM(obst_nout_dens_f)//'_'//TRIM(obst_varname(iv))
            CALL obst_nccheck(NF90_INQ_VARID(ncid, name, varid), &
                              "obst_init_spatial nf90_inq_varid "//obst_fn_initspatial(iv)//" "//name)
            CALL obst_ncget2D(ncid, varid, tmp, &
                              "obst_init_spatial getvar "//obst_fn_initspatial(iv)//" "//name)
            obst_dens_inst(iv, GLOBAL_2D_ARRAY) = tmp(GLOBAL_2D_ARRAY)

            name = TRIM(obst_nout_width_f)//'_'//TRIM(obst_varname(iv))
            CALL obst_nccheck(NF90_INQ_VARID(ncid, name, varid), &
                              "obst_init_spatial nf90_inq_varid "//obst_fn_initspatial(iv)//" "//name)
            CALL obst_ncget2D(ncid, varid, tmp, &
                              "obst_init_spatial getvar "//obst_fn_initspatial(iv)//" "//name)
            obst_width_inst(iv, GLOBAL_2D_ARRAY) = tmp(GLOBAL_2D_ARRAY)

            name = TRIM(obst_nout_thick_f)//'_'//TRIM(obst_varname(iv))
            CALL obst_nccheck(NF90_INQ_VARID(ncid, name, varid), &
                              "obst_init_spatial nf90_inq_varid "//obst_fn_initspatial(iv)//" "//name)
            CALL obst_ncget2D(ncid, varid, tmp, &
                              "obst_init_spatial getvar "//obst_fn_initspatial(iv)//" "//name)
            obst_thick_inst(iv, GLOBAL_2D_ARRAY) = tmp(GLOBAL_2D_ARRAY)

            ! Close input NetCDF file.
            CALL obst_nccheck(NF90_CLOSE(ncid), &
                              "obst_init_spatial nf90_close "//obst_fn_initspatial(iv))
         ELSE
            obst_height_inst(iv, :, :) = obst_i_height(iv)
            obst_dens_inst(iv, :, :) = obst_i_dens(iv)
            obst_width_inst(iv, :, :) = obst_i_width(iv)
            obst_thick_inst(iv, :, :) = obst_i_thick(iv)
         END IF
      END DO

      WHERE (obst_position == 0.)
         obst_height_inst = 0.
         obst_dens_inst = 0.
         obst_width_inst = 0.
         obst_thick_inst = 0.
      END WHERE

      RETURN

   END SUBROUTINE obst_init_spatial

   !!==========================================================================================================

   SUBROUTINE obst_init_timeseries
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_init_timeseries  ***
      !!
      !! ** Purpose : Read timeseries files and initialize obst_ts_* variables
      !!
      !!---------------------------------------------------------------------
      USE netcdf
      IMPLICIT NONE

      INTEGER iv, ncid, varid, ts_tmax, dimid, lendim
      CHARACTER(30) name
      REAL(KIND=rlg) origin_date_in_sec
      REAL, DIMENSION(:), ALLOCATABLE :: tmp

      ! allocate vars depending on obst_nbvar
      ALLOCATE (obst_ts_tmax(1:obst_nbvar))
      obst_ts_tmax(:) = 0
      ALLOCATE (obst_ts_tbefore(1:obst_nbvar))
      obst_ts_tbefore(:) = 1 ! initialize tbefore to the first index
      ALLOCATE (obst_ts_tafter(1:obst_nbvar))
      obst_ts_tafter(:) = 1 ! initialize tafter to the first index

      ! get max number of time dim in all files
      DO iv = 1, obst_nbvar
         IF (obst_l_filetimeserie(iv)) THEN
            CALL obst_nccheck(NF90_OPEN(obst_fn_timeserie(iv), NF90_NOWRITE, ncid), &
                              "obst_init_timeseries nf90_open "//obst_fn_timeserie(iv))
            CALL obst_nccheck(NF90_INQ_DIMID(ncid, 'time', dimid), &
                              "obst_init_timeseries nf90_inq_dimid "//obst_fn_timeserie(iv)//' time')
            CALL obst_nccheck(NF90_INQUIRE_DIMENSION(ncid, dimid, len=lendim), &
                              "obst_init_timeseries nf90_inquire_dimension "//obst_fn_timeserie(iv)//' time')
            CALL obst_nccheck(NF90_CLOSE(ncid), &
                              "obst_init_timeseries nf90_close "//obst_fn_timeserie(iv))
            obst_ts_tmax(iv) = lendim
         END IF
      END DO

      ts_tmax = maxval(obst_ts_tmax(:))
      IF (ts_tmax > 0) THEN ! there is at least one time series to allocate
         ! allocate vars depending on obst_nbvar and ts_tmax
         ALLOCATE (obst_ts_time(1:obst_nbvar, 1:ts_tmax))
         obst_ts_time(:, :) = 0.
         ALLOCATE (obst_ts_height(1:obst_nbvar, 1:ts_tmax))
         obst_ts_height(:, :) = 0.
         ALLOCATE (obst_ts_dens(1:obst_nbvar, 1:ts_tmax))
         obst_ts_dens(:, :) = 0.
         ALLOCATE (obst_ts_width(1:obst_nbvar, 1:ts_tmax))
         obst_ts_width(:, :) = 0.
         ALLOCATE (obst_ts_thick(1:obst_nbvar, 1:ts_tmax))
         obst_ts_thick(:, :) = 0.
         ALLOCATE (tmp(1:ts_tmax))
         tmp(:) = 0.

         ! read height, dens, width, thick vars
         DO iv = 1, obst_nbvar
            IF (obst_l_filetimeserie(iv)) THEN
               CALL obst_nccheck(NF90_OPEN(obst_fn_timeserie(iv), NF90_NOWRITE, ncid), &
                                 "obst_init_timeseries nf90_open "//obst_fn_timeserie(iv))

               CALL obst_nccheck(NF90_INQ_DIMID(ncid, 'time', dimid), &
                                 "obst_init_timeseries nf90_inq_dimid "//obst_fn_timeserie(iv)//" time")

               CALL obst_nccheck(NF90_INQ_VARID(ncid, 'time', varid), &
                                 "obst_init_timeseries nf90_inq_varid "//obst_fn_timeserie(iv)//" time")

               CALL tool_origindate(ncid, varid, origin_date_in_sec)

               CALL obst_nccheck(NF90_GET_VAR(ncid, varid, tmp(1:obst_ts_tmax(iv))), &
                                 "obst_init_timeseries nf90_get_var "//obst_fn_timeserie(iv)//" time")
               obst_ts_time(iv, 1:obst_ts_tmax(iv)) = tmp(1:obst_ts_tmax(iv)) + origin_date_in_sec

               name = TRIM(obst_nout_height_f)//'_'//TRIM(obst_varname(iv))
               CALL obst_nccheck(NF90_INQ_VARID(ncid, name, varid), &
                                 "obst_init_timeseries nf90_inq_varid "//obst_fn_timeserie(iv)//" "//name)
               CALL obst_nccheck(NF90_GET_VAR(ncid, varid, tmp(1:obst_ts_tmax(iv))), &
                                 "obst_init_timeseries nf90_getvar "//obst_fn_timeserie(iv)//" "//name)
               obst_ts_height(iv, 1:obst_ts_tmax(iv)) = tmp(1:obst_ts_tmax(iv))

               name = TRIM(obst_nout_dens_f)//'_'//TRIM(obst_varname(iv))
               CALL obst_nccheck(NF90_INQ_VARID(ncid, name, varid), &
                                 "obst_init_timeseries nf90_inq_varid "//obst_fn_timeserie(iv)//" "//name)
               CALL obst_nccheck(NF90_GET_VAR(ncid, varid, tmp(1:obst_ts_tmax(iv))), &
                                 "obst_init_timeseries nf90_getvar "//obst_fn_timeserie(iv)//" "//name)
               obst_ts_dens(iv, 1:obst_ts_tmax(iv)) = tmp(1:obst_ts_tmax(iv))

               name = TRIM(obst_nout_width_f)//'_'//TRIM(obst_varname(iv))
               CALL obst_nccheck(NF90_INQ_VARID(ncid, name, varid), &
                                 "obst_init_timeseries nf90_inq_varid "//obst_fn_timeserie(iv)//" "//name)
               CALL obst_nccheck(NF90_GET_VAR(ncid, varid, tmp(1:obst_ts_tmax(iv))), &
                                 "obst_init_timeseries nf90_getvar "//obst_fn_timeserie(iv)//" "//name)
               obst_ts_width(iv, 1:obst_ts_tmax(iv)) = tmp(1:obst_ts_tmax(iv))

               name = TRIM(obst_nout_thick_f)//'_'//TRIM(obst_varname(iv))
               CALL obst_nccheck(NF90_INQ_VARID(ncid, name, varid), &
                                 "obst_init_timeseries nf90_inq_varid "//obst_fn_timeserie(iv)//" "//name)
               CALL obst_nccheck(NF90_GET_VAR(ncid, varid, tmp(1:obst_ts_tmax(iv))), &
                                 "obst_init_timeseries nf90_getvar "//obst_fn_timeserie(iv)//" "//name)
               obst_ts_thick(iv, 1:obst_ts_tmax(iv)) = tmp(1:obst_ts_tmax(iv))

               CALL obst_nccheck(NF90_CLOSE(ncid), &
                                 "obst_init_timeseries nf90_close "//obst_fn_timeserie(iv))
            END IF
         END DO
      END IF

      RETURN

   END SUBROUTINE obst_init_timeseries

   !!==========================================================================================================

   SUBROUTINE obst_readfile_distri
      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_readfile_distri  ***
      !!
      !! ** Purpose : Initialization of obstruction parameters
      !!
      !!---------------------------------------------------------------------
      USE module_OBSTRUCTIONS ! for mynode
      IMPLICIT NONE

      INTEGER :: kk, iv
      LOGICAL :: ex

      MPI_master_only WRITE (iscreenlog, *) &
         '************************************************************************'
      MPI_master_only WRITE (iscreenlog, *) &
         '***** module OBSTRUCTIONS, subroutine obst_READFILE_DISTRI *****'
      MPI_master_only WRITE (iscreenlog, *) &
         '************************************************************************'
      DO iv = 1, obst_nbvar
         IF (obst_l_filedistri(iv)) THEN
            MPI_master_only WRITE (iscreenlog, *) &
               ' '
            MPI_master_only WRITE (iscreenlog, *) &
               'Obstruction variable : ', obst_varname(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               'file defining vertical distribution of obstuctions :'
            MPI_master_only WRITE (iscreenlog, *) TRIM(obst_fn_distrib(iv))
            INQUIRE (file=obst_fn_distrib(iv), exist=ex)
            IF (ex) THEN
               OPEN (53, file=obst_fn_distrib(iv), form='formatted')
               READ (53, *) ! Filename
               READ (53, *) ! Title
               READ (53, *) obst_nbhnorm(iv)
               MPI_master_only WRITE (iscreenlog, *) &
                  ' '
               MPI_master_only WRITE (iscreenlog, *) &
                  'Number of vertical discretization : ', obst_nbhnorm(iv)
               CLOSE (53)
            ELSE
               MPI_master_only WRITE (ierrorlog, *) &
                  ' '
               MPI_master_only WRITE (ierrorlog, *) &
                  ' '
               MPI_master_only WRITE (ierrorlog, *) &
                  '************************************************************************'
               MPI_master_only WRITE (ierrorlog, *) &
                  '***** module OBSTRUCTIONS, subroutine obst_READFILE_DISTRI *****'
               MPI_master_only WRITE (ierrorlog, *) &
                  '************************************************************************'
               MPI_master_only WRITE (ierrorlog, *) &
                  ' ERROR : File ', TRIM(obst_fn_distrib(iv)), 'DOes not exist'
               MPI_master_only WRITE (ierrorlog, *) &
                  ' --> THE SIMULATION IS STOPPED !!! '
               MPI_master_only WRITE (ierrorlog, *) &
                  '************************************************************************'
               STOP
            END IF
         ELSE
            obst_nbhnorm(iv) = 11
            MPI_master_only WRITE (iscreenlog, *) &
               ' '
            MPI_master_only WRITE (iscreenlog, *) &
               'Obstruction variable : ', obst_varname(iv)
            MPI_master_only WRITE (iscreenlog, *) &
               'An homogene distribution is applied'
            MPI_master_only WRITE (iscreenlog, *) &
               ' '
            MPI_master_only WRITE (iscreenlog, *) &
               'Number of vertical discretization : ', obst_nbhnorm(iv)
         END IF
      END DO

      obst_nb_max_hnorm = 0
      IF (obst_nbvar .EQ. 1) THEN
         obst_nb_max_hnorm = obst_nbhnorm(1)
      ELSE
         DO iv = 1, obst_nbvar
            obst_nb_max_hnorm = MAX(obst_nb_max_hnorm, obst_nbhnorm(iv))
         END DO
      END IF
      ALLOCATE (obst_dens_norm(1:obst_nbvar, 1:obst_nb_max_hnorm))
      ALLOCATE (obst_height_norm(1:obst_nbvar, 1:obst_nb_max_hnorm))

      DO iv = 1, obst_nbvar
         IF (obst_l_filedistri(iv)) THEN
            OPEN (53, file=obst_fn_distrib(iv), form='formatted')
            READ (53, *) ! Filename
            READ (53, *) ! Title
            READ (53, *) ! NBhnorm
            READ (53, *) ! Column titles
            kk = 1
            DO WHILE (kk .LE. obst_nbhnorm(iv))
               READ (53, *) obst_height_norm(iv, kk), obst_dens_norm(iv, kk)
               MPI_master_only WRITE (iscreenlog, *) &
                  'kk', kk
               MPI_master_only WRITE (iscreenlog, *) &
                  'height : ', obst_height_norm(iv, kk), ' dens : ', obst_dens_norm(iv, kk)
               kk = kk + 1
            END DO
            CLOSE (53)
         ELSE
            ! Homogene distribution of obstructions distribution
            obst_dens_norm(iv, :) = 100.0_rsh
            obst_height_norm(iv, :) = 0.001_rsh
            DO kk = 1, obst_nbhnorm(iv)
               obst_height_norm(iv, kk) = (REAL(kk - 1))*10.0_rsh
               obst_dens_norm(iv, kk) = 100.0_rsh
            END DO
            obst_height_norm(iv, obst_nbhnorm(iv)) = 100.01_rsh
         END IF
      END DO

   END SUBROUTINE obst_readfile_distri

   !!==========================================================================================================

   SUBROUTINE obst_alloc_nbvar

      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_alloc_nbvar  ***
      !!
      !! ** Purpose : Allocation of tables depending on number of obstructions variables
      !!
      !!---------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------
      ! Variables on (iv)
      !--------------------
      ALLOCATE (obst_fracxy_type(1:obst_nbvar))
      obst_fracxy_type(:) = 0
      ALLOCATE (obst_nbhnorm(1:obst_nbvar))
      obst_nbhnorm(:) = 0
      ALLOCATE (obst_c_abdel_nmax(1:obst_nbvar))
      obst_c_abdel_nmax(:) = 1

      ALLOCATE (obst_l_filetimeserie(1:obst_nbvar))
      obst_l_filetimeserie(:) = .FALSE.
      ALLOCATE (obst_l_filedistri(1:obst_nbvar))
      obst_l_filedistri(:) = .FALSE.
      ALLOCATE (obst_l_init_spatial(1:obst_nbvar))
      obst_l_init_spatial(:) = .FALSE.
      ALLOCATE (obst_l_flexible(1:obst_nbvar))
      obst_l_flexible(:) = .FALSE.
      ALLOCATE (obst_l_cylinder(1:obst_nbvar))
      obst_l_cylinder(:) = .FALSE.
      ALLOCATE (obst_l_noturb(1:obst_nbvar))
      obst_l_noturb(:) = .FALSE.
      ALLOCATE (obst_l_abdelrough_cste(1:obst_nbvar))
      obst_l_abdelrough_cste = .FALSE.
      ALLOCATE (obst_l_fracxy(1:obst_nbvar))
      obst_l_fracxy(:) = .FALSE.
      ALLOCATE (obst_l_abdelposture(1:obst_nbvar))
      obst_l_abdelposture(:) = .FALSE.
      ALLOCATE (obst_l_param_height(1:obst_nbvar))
      obst_l_param_height(:) = .FALSE.
      ALLOCATE (obst_l_drag_cste(1:obst_nbvar))
      obst_l_drag_cste(:) = .TRUE.
      ALLOCATE (obst_l_z0bstress(1:obst_nbvar))
      obst_l_z0bstress(:) = .FALSE.

      ALLOCATE (obst_z0bstress_option(1:obst_nbvar))
      obst_z0bstress_option(:) = 0

      ALLOCATE (obst_varname(1:obst_nbvar))
      obst_varname(:) = '.'
      ALLOCATE (obst_type(1:obst_nbvar))
      obst_type(:) = '..'
      ALLOCATE (obst_fn_timeserie(1:obst_nbvar))
      obst_fn_timeserie(:) = '.'
      ALLOCATE (obst_fn_initspatial(1:obst_nbvar))
      obst_fn_initspatial(:) = '.'
      ALLOCATE (obst_fn_distrib(1:obst_nbvar))
      obst_fn_distrib(:) = '.'

      ALLOCATE (obst_i_height(1:obst_nbvar))
      obst_i_height(:) = 0.0_rsh
      ALLOCATE (obst_i_width(1:obst_nbvar))
      obst_i_width(:) = 0.0_rsh
      ALLOCATE (obst_i_thick(1:obst_nbvar))
      obst_i_thick(:) = 0.0_rsh
      ALLOCATE (obst_i_dens(1:obst_nbvar))
      obst_i_dens(:) = 0.0_rsh

      ALLOCATE (obst_c_rho(1:obst_nbvar))
      obst_c_rho(:) = 0.0_rsh
      ALLOCATE (obst_c_drag(1:obst_nbvar))
      obst_c_drag(:) = 0.0_rsh
      ALLOCATE (obst_c_lift(1:obst_nbvar))
      obst_c_lift(:) = 0.0_rsh
      ALLOCATE (obst_c_z0bstress(1:obst_nbvar))
      obst_c_z0bstress(:) = 0.0_rsh
      ALLOCATE (obst_c_fracxy_k0(1:obst_nbvar))
      obst_c_fracxy_k0(:) = 0.0_rsh
      ALLOCATE (obst_c_fracxy_k1(1:obst_nbvar))
      obst_c_fracxy_k1(:) = 0.0_rsh
      ALLOCATE (obst_c_fracxy_l(1:obst_nbvar))
      obst_c_fracxy_l(:) = 0.0_rsh
      ALLOCATE (obst_c_crough_x0(1:obst_nbvar))
      obst_c_crough_x0(:) = 0.0_rsh
      ALLOCATE (obst_c_crough_x1(1:obst_nbvar))
      obst_c_crough_x1(:) = 0.0_rsh
      ALLOCATE (obst_c_lz(1:obst_nbvar))
      obst_c_lz(:) = 0.0_rsh
      ALLOCATE (obst_c_shelter(1:obst_nbvar))
      obst_c_shelter(:) = 1.0_rsh
      ALLOCATE (obst_c_height_x0(1:obst_nbvar))
      obst_c_height_x0 = 0.0_rsh
      ALLOCATE (obst_c_height_x1(1:obst_nbvar))
      obst_c_height_x1(:) = 0.0_rsh
      ALLOCATE (obst_c_z0bstress_x0(1:obst_nbvar))
      obst_c_z0bstress_x0(:) = 0.0_rsh
      ALLOCATE (obst_c_z0bstress_x1(1:obst_nbvar))
      obst_c_z0bstress_x1(:) = 0.0_rsh
      !-------------------------
   END SUBROUTINE obst_alloc_nbvar

   !!==========================================================================================================

   SUBROUTINE obst_alloc_xyz

      !!---------------------------------------------------------------------
      !!                 *** SUBROUTINE obst_alloc_xyz  ***
      !!
      !! ** Purpose : Allocation of spatial obstruction tables
      !!
      !!---------------------------------------------------------------------
      USE module_OBSTRUCTIONS ! for GLOBAL_2D_ARRAY
      IMPLICIT NONE

      !----------------------
      ! Variables on (iv,i,j)
      !----------------------
      ALLOCATE (obst_dens_inst(1:obst_nbvar, GLOBAL_2D_ARRAY)) ! always needed
      obst_dens_inst(:, :, :) = 0.0_rsh
      ALLOCATE (obst_width_inst(1:obst_nbvar, GLOBAL_2D_ARRAY)) ! always needed
      obst_width_inst(:, :, :) = 0.0_rsh
      ALLOCATE (obst_thick_inst(1:obst_nbvar, GLOBAL_2D_ARRAY)) ! always needed
      obst_thick_inst(:, :, :) = 0.0_rsh
      ALLOCATE (obst_height_inst(1:obst_nbvar, GLOBAL_2D_ARRAY)) ! always needed
      obst_height_inst(:, :, :) = 0.0_rsh

      ALLOCATE (obst_position(1:obst_nbvar, GLOBAL_2D_ARRAY)) ! always needed
      obst_position(:, :, :) = 0.0_rsh
      ALLOCATE (obst_height(1:obst_nbvar, GLOBAL_2D_ARRAY)) ! always needed for previous value
      obst_height(:, :, :) = 0.0_rsh

      IF (l_obstout_frac_xy) THEN
         ALLOCATE (obst_fracxy(1:obst_nbvar, GLOBAL_2D_ARRAY)) ! optionnal
         obst_fracxy(:, :, :) = 0.0_rsh
      END IF
      IF (l_obstout_a2d) THEN
         ALLOCATE (obst_a2d(1:obst_nbvar + 3, GLOBAL_2D_ARRAY)) ! optionnal
         obst_a2d(:, :, :) = 0.0_rsh
      END IF
      IF (l_obstout_s2d) THEN
         ALLOCATE (obst_s2d(1:obst_nbvar + 3, GLOBAL_2D_ARRAY)) ! optionnal
         obst_s2d(:, :, :) = 0.0_rsh
      END IF

      ALLOCATE (obst_z0obst(1:obst_nbvar + 3, GLOBAL_2D_ARRAY))
      obst_z0obst(:, :, :) = 0.0_rsh

      !-------------------
      ! Variables on (i,j)
      !-------------------
      ALLOCATE (obst_z0bed(GLOBAL_2D_ARRAY))
      obst_z0bed(:, :) = 0.0_rsh

      ALLOCATE (obst_z0bstress(GLOBAL_2D_ARRAY))
      obst_z0bstress(:, :) = 0.0_rsh

      !------------------------
      ! Variables on (iv,k,i,j)
      !------------------------
      IF (l_obstout_dens_e) THEN
         ALLOCATE (obst_dens3d(1:obst_nbvar, 1:obst_kmax, GLOBAL_2D_ARRAY)) ! optionnal
         obst_dens3d(:, :, :, :) = 0.0_rsh
      END IF
      IF (l_obstout_width_e) THEN
         ALLOCATE (obst_width3d(1:obst_nbvar, 1:obst_kmax, GLOBAL_2D_ARRAY)) ! optionnal
         obst_width3d(:, :, :, :) = 0.0_rsh
      END IF
      IF (l_obstout_thick_e) THEN
         ALLOCATE (obst_thick3d(1:obst_nbvar, 1:obst_kmax, GLOBAL_2D_ARRAY)) ! optionnal
         obst_thick3d(:, :, :, :) = 0.0_rsh
      END IF
      IF (l_obstout_theta) THEN
         ALLOCATE (obst_theta3d(1:obst_nbvar, 1:obst_kmax, GLOBAL_2D_ARRAY)) ! optionnal
         obst_theta3d(:, :, :, :) = 0.0_rsh
      END IF
      IF (l_obstout_frac_z) THEN
         ALLOCATE (obst_fracz3d(1:obst_nbvar, 1:obst_kmax, GLOBAL_2D_ARRAY)) ! optionnal
         obst_fracz3d(:, :, :, :) = 0.0_rsh
      END IF
      IF (l_obstout_drag) THEN
         ALLOCATE (obst_drag3d(1:obst_nbvar, 1:obst_kmax, GLOBAL_2D_ARRAY)) ! optionnal
         obst_drag3d(:, :, :, :) = 0.0_rsh
      END IF
      ALLOCATE (obst_a3d(1:obst_nbvar + 3, 1:obst_kmax, GLOBAL_2D_ARRAY)) ! always needed
      obst_a3d(:, :, :, :) = 0.0_rsh
      IF (l_obstout_s3d) THEN
         ALLOCATE (obst_s3d(1:obst_nbvar + 3, 1:obst_kmax, GLOBAL_2D_ARRAY)) ! optionnal
         obst_s3d(:, :, :, :) = 0.0_rsh
      END IF
      !---------------------
      ! Variables on (i,j,k)
      !---------------------
      ALLOCATE (obst_fuz(GLOBAL_2D_ARRAY, 1:obst_kmax)) ! always needed
      obst_fuz(:, :, :) = 0.0_rsh
      ALLOCATE (obst_fvz(GLOBAL_2D_ARRAY, 1:obst_kmax)) ! always needed
      obst_fvz(:, :, :) = 0.0_rsh
      ALLOCATE (obst_t(GLOBAL_2D_ARRAY, 1:obst_kmax)) ! always needed
      obst_t(:, :, :) = 0.0_rsh
      ALLOCATE (obst_tau(GLOBAL_2D_ARRAY, 1:obst_kmax)) ! always needed
      obst_tau(:, :, :) = 0.0_rsh
      !------------------------------
   END SUBROUTINE obst_alloc_xyz

      !!==========================================================================================================

#endif
END MODULE init_OBSTRUCTIONS
