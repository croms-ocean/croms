#include "cppdefs.h"

MODULE iom

#if defined key_pisces

#ifdef XIOS
   USE xios
#endif

   IMPLICIT NONE
   PUBLIC

  PUBLIC iom_open, iom_close, iom_get, iom_put
  PUBLIC iom_use, iom_getszuld, iom_varid, iom_rstput

  INTERFACE iom_put
     MODULE PROCEDURE iom_p0d, iom_p1d, iom_p2d, iom_p3d
  END INTERFACE
  INTERFACE iom_get
     MODULE PROCEDURE iom_g0d, iom_g1d, iom_g2d, iom_g3d
  END INTERFACE
  INTERFACE iom_rstput
     MODULE PROCEDURE iom_rst_p0d, iom_rst_p1d, iom_rst_p2d, iom_rst_p3d
  END INTERFACE
  INTEGER :: jpdom_data = 1
  INTEGER :: jpdom_auto = 1
  INTEGER :: jptra_sms  = 1

#  include "ocean2pisces.h90"
#  include "do_loop_substitute.h90"

CONTAINS

   !!----------------------------------------------------------------------
   !!                   INTERFACE iom
   !!----------------------------------------------------------------------

   SUBROUTINE iom_open( cdname, kiomid )
      CHARACTER(len=*), INTENT(in   )    ::   cdname   ! File name
      INTEGER         , INTENT(  out)    ::   kiomid   ! iom identifier of the opened file
      kiomid = -1
   END SUBROUTINE iom_open

   SUBROUTINE iom_close( kiomid )
      INTEGER         , INTENT(in)       ::   kiomid   ! iom identifier of the opened file
   END SUBROUTINE iom_close

   SUBROUTINE iom_p0d( cdname, pfield0d )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdname
      REAL,  INTENT(in) ::   pfield0d
#ifdef XIOS
      CALL xios_send_field(cdname, (/pfield0d/))
#endif
   END SUBROUTINE iom_p0d

   SUBROUTINE iom_p1d( cdname, pfield1d )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdname
      REAL,     DIMENSION(:), INTENT(in) ::   pfield1d
#ifdef XIOS
      CALL xios_send_field( cdname, RESHAPE( (/pfield1d/), (/1,1,SIZE(pfield1d)/) ) )
#endif
   END SUBROUTINE iom_p1d

   SUBROUTINE iom_p2d( cdname, pfield2d )
      CHARACTER(LEN=*)            , INTENT(in) ::   cdname
      REAL,     DIMENSION(:,:), INTENT(in) ::   pfield2d
#ifdef XIOS
      CALL xios_send_field(cdname, pfield2d)
#endif
   END SUBROUTINE iom_p2d

   SUBROUTINE iom_p3d( cdname, pfield3d )
      CHARACTER(LEN=*)                , INTENT(in) ::   cdname
      REAL,       DIMENSION(:,:,:), INTENT(in) ::   pfield3d
#ifdef XIOS
      CALL xios_send_field(cdname, pfield3d)
#endif
   END SUBROUTINE iom_p3d

   LOGICAL FUNCTION iom_use( cdname )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      !!----------------------------------------------------------------------
#ifdef XIOS
      iom_use = xios_field_is_active( cdname )
#else
      iom_use = .FALSE.
#endif
   END FUNCTION iom_use

   FUNCTION iom_getszuld ( kiomid )
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  iom_getszuld  ***
      !!
      !! ** Purpose : get the size of the unlimited dimension in a file
      !!              (return -1 if not found)
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kiomid   ! file Identifier
      !
      INTEGER                ::   iom_getszuld
      !!-----------------------------------------------------------------------
      iom_getszuld = -1
   END FUNCTION iom_getszuld

  FUNCTION iom_varid ( kiomid, cdvar, kdimsz, kndims, lduld, ldstop )
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  iom_varid  ***
      !!
      !! ** Purpose : get the id of a variable in a file (return 0 if not found)
      !!-----------------------------------------------------------------------
      INTEGER              , INTENT(in   )           ::   kiomid   ! file Identifier
      CHARACTER(len=*)     , INTENT(in   )           ::   cdvar    ! name of the variable
      INTEGER, DIMENSION(:), INTENT(  out), OPTIONAL ::   kdimsz   ! size of each dimension
      INTEGER              , INTENT(  out), OPTIONAL ::   kndims   ! number of dimensions
      LOGICAL              , INTENT(  out), OPTIONAL ::   lduld    ! true if the last dimension is unlimited (time)
      LOGICAL              , INTENT(in   ), OPTIONAL ::   ldstop   ! stop if looking for non-existing variable (default = .TRUE.)
      !
      INTEGER                        ::   iom_varid, iiv, i_nvd
      LOGICAL                        ::   ll_fnd
      CHARACTER(LEN=100)             ::   clinfo                   ! info character
      LOGICAL                        ::   llstop                   ! local definition of ldstop
      !!-----------------------------------------------------------------------
      iom_varid = 0                         ! default definition
      !
   END FUNCTION iom_varid

   SUBROUTINE iom_rst_p0d( kt, kwrite, kiomid, cdvar,  pvdp0d )
      INTEGER                    , INTENT(in) ::   kt       ! ocean time-step
      INTEGER                    , INTENT(in) ::   kwrite   ! writing time-step
      INTEGER                    , INTENT(in) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*)           , INTENT(in) ::   cdvar    ! time axis name
      REAL                       , INTENT(in) ::   pvdp0d    ! read field (0D case), double precision
   END SUBROUTINE iom_rst_p0d

   SUBROUTINE iom_rst_p1d( kt, kwrite, kiomid, cdvar, pvdp1d )
      INTEGER           , INTENT(in) ::   kt       ! ocean time-step
      INTEGER           , INTENT(in) ::   kwrite   ! writing time-step
      INTEGER           , INTENT(in) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*)  , INTENT(in) ::   cdvar    ! time axis name
      REAL, DIMENSION(:), INTENT(in) ::   pvdp1d    ! read field (1D case), double precision
   END SUBROUTINE iom_rst_p1d

   SUBROUTINE iom_rst_p2d( kt, kwrite, kiomid, cdvar, pvdp2d )
      INTEGER             , INTENT(in) ::   kt       ! ocean time-step
      INTEGER             , INTENT(in) ::   kwrite   ! writing time-step
      INTEGER             , INTENT(in) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*)    , INTENT(in) ::   cdvar    ! time axis name
      REAL, DIMENSION(:,:), INTENT(in) ::   pvdp2d    ! read field (1D case), double precision
   END SUBROUTINE iom_rst_p2d

   SUBROUTINE iom_rst_p3d( kt, kwrite, kiomid, cdvar, pvdp3d )
      INTEGER               , INTENT(in) ::   kt       ! ocean time-step
      INTEGER               , INTENT(in) ::   kwrite   ! writing time-step
      INTEGER               , INTENT(in) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in) ::   cdvar    ! time axis name
      REAL, DIMENSION(:,:,:), INTENT(in) ::   pvdp3d    ! read field (1D case), double precision
   END SUBROUTINE iom_rst_p3d

   SUBROUTINE iom_g0d( kiomid, cdvar,  pvdp0d )
      INTEGER                    , INTENT(in) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*)           , INTENT(in) ::   cdvar    ! time axis name
      REAL                       , INTENT(in) ::   pvdp0d    ! read field (0D case), double precision
   END SUBROUTINE iom_g0d

   SUBROUTINE iom_g1d( kiomid, kdom, cdvar, pvdp1d )
      INTEGER           , INTENT(in) ::   kiomid   ! Identifier of the file
      INTEGER           , INTENT(in) ::   kdom   ! writing time-step
      CHARACTER(len=*)  , INTENT(in) ::   cdvar    ! time axis name
      REAL, DIMENSION(:), INTENT(in) ::   pvdp1d    ! read field (1D case), double precision
   END SUBROUTINE iom_g1d

   SUBROUTINE iom_g2d( kiomid, kdom, cdvar, pvdp2d )
      INTEGER             , INTENT(in) ::   kiomid   ! Identifier of the file
      INTEGER             , INTENT(in) ::   kdom   ! writing time-step
      CHARACTER(len=*)    , INTENT(in) ::   cdvar    ! time axis name
      REAL, DIMENSION(:,:), INTENT(in) ::   pvdp2d    ! read field (1D case), double precision
   END SUBROUTINE iom_g2d

   SUBROUTINE iom_g3d( kiomid, kdom, cdvar, pvdp3d )
      INTEGER               , INTENT(in) ::   kiomid   ! Identifier of the file
      INTEGER               , INTENT(in) ::   kdom   ! writing time-step
      CHARACTER(len=*)      , INTENT(in) ::   cdvar    ! time axis name
      REAL, DIMENSION(:,:,:), INTENT(in) ::   pvdp3d    ! read field (1D case), double precision
   END SUBROUTINE iom_g3d

#endif

END MODULE iom
