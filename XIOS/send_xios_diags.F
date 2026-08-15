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

#ifdef XIOS

      MODULE xios_module
      USE xios
#ifdef key_pisces
      USE trc, ONLY : profsed
#endif
      IMPLICIT NONE

      PUBLIC set_grid
      PUBLIC iom_set_domain_attr
      PUBLIC iom_set_axis_attr
      PUBLIC iom_set_field_attr
      PUBLIC iom_set_file_attr
      PUBLIC iom_set_file_global_attr_char
      PUBLIC iom_set_file_global_attr_float
      PUBLIC iom_get_file_attr
      PUBLIC set_scalar
#if defined ONLINE_ANALYSIS
      PUBLIC set_scalogram_dom
      PUBLIC set_scalogram_axis
#endif    

      CONTAINS


      SUBROUTINE set_grid( cdgrd, plon, plat)
# include "param.h"
# include "ncscrum.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_grid  ***
      !!
      !! ** Purpose :   set domain definition to XIOS library
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                A. Ponte    , adaptation for CROCO
      !!                S. Masson   , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)            , INTENT(in) ::   cdgrd
      REAL plon(GLOBAL_2D_ARRAY)
      REAL plat(GLOBAL_2D_ARRAY)
      !
      INTEGER vtype, horiz_type, vert_type
      INTEGER imin, imax, jmin, jmax, ni, nj
      INTEGER nx_glo, ny_glo 
      INTEGER istart, jstart
      INTEGER data_ni,data_nj 
      INTEGER :: ji
      !!----------------------------------------------------------------------
      !
      ! See ncscrum.h for r2dvar ...
      SELECT CASE (TRIM(cdgrd))
      CASE ("rho")
         vtype=r2dvar
         nx_glo=LLm+2
         ny_glo=MMm+2
      CASE ("b")
         vtype=r2dvar
         nx_glo=LLm+2
         ny_glo=MMm+2
      CASE ("u")
         vtype=u2dvar
         nx_glo=LLm+1
         ny_glo=MMm+2
      CASE ("v")
         vtype=v2dvar
         nx_glo=LLm+2
         ny_glo=MMm+1
      CASE ("w")
         vtype=r2dvar
         nx_glo=LLm+2
         ny_glo=MMm+2
      CASE ("psi")
         vtype=p2dvar
         nx_glo=LLm+1
         ny_glo=MMm+1
      END SELECT

      ! From nf_fread.F
      vert_type=vtype/4         ! vertical and horizontal
      horiz_type=vtype-4*vert_type ! grid types, then calculate
      jmin=horiz_type/2         ! starting indices indices
      imin=horiz_type-2*jmin    ! in horizontal directions.
      !     
      istart=1
      jstart=1

# ifdef MPI
      IF (ii.gt.0) then
         istart=1-imin+iminmpi
         imin=1
      ENDIF
      IF (ii.eq.NP_XI-1) then
         imax=Lmmpi+1
      ELSE
         imax=Lmmpi
      ENDIF
      IF (jj.gt.0) then
         jstart=1-jmin+jminmpi
         jmin=1
      ENDIF
      IF (jj.eq.NP_ETA-1) then
         jmax=Mmmpi+1
      ELSE
         jmax=Mmmpi
      ENDIF
# else
      imin=1
      imax=Lm+1
      jmin=1
      jmax=Mm+1
# endif

      ! there is probably already a variable for ni and nj 
      ni=imax-imin+1
      nj=jmax-jmin+1

      ! from set_global_definitions.h
# ifdef THREE_GHOST_POINTS
#  ifdef MPI
      data_ni = Lm+3+padd_X+2+1
      data_nj = Mm+3+padd_E+2+1
#  endif
# else
#  ifdef MPI
      data_ni = Lm+2+padd_X+1+1
      data_nj = Mm+2+padd_E+1+1
#  endif
# endif

      ! notes iminmpi, jminmpi are computed in MPI_Setup.F
      CALL iom_set_domain_attr(TRIM(cdgrd),ni_glo=nx_glo, nj_glo=ny_glo,
     &     ibegin=istart-1, ni=ni,jbegin=jstart-1,nj=nj)
      ! to do: below should have be rewritten to match all configs (periodic, number of ghost points, ...)
      ! see set_global_definitions.h
      CALL iom_set_domain_attr(TRIM(cdgrd),data_dim=2,
     &     data_ibegin=START_1D_ARRAYXI-imin,
     &     data_ni=data_ni,
     &     data_jbegin=START_1D_ARRAYETA-jmin, 
     &     data_nj=data_nj)

      CALL iom_set_domain_attr(TRIM(cdgrd),
     &     lonvalue=RESHAPE(plon(imin:imax,jmin:jmax), (/ni*nj/)),
     &     latvalue=RESHAPE(plat(imin:imax,jmin:jmax), (/ni*nj/)))
      !      
      END SUBROUTINE set_grid


#ifdef MASKING
      SUBROUTINE set_mask(cdgrd)
# include "param.h"
# include "ncscrum.h"
# include "scalars.h"
# include "grid.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_mask  ***
      !!
      !! ** Purpose :   set mask in domain definition to XIOS library
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Petton   , adapt from set_grid
      !!                S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                A. Ponte    , adaptation for CROCO
      !!                S. Masson   , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)            , INTENT(in) ::   cdgrd
      LOGICAL zmask(GLOBAL_2D_ARRAY)
      !
      INTEGER vtype, horiz_type, vert_type
      INTEGER imin, imax, jmin, jmax, ni, nj
      !!----------------------------------------------------------------------
      !
      zmask(:,:) = .false.
      SELECT CASE (TRIM(cdgrd))
      CASE ("rho")
         vtype=r2dvar
         where(rmask(:,:) > 0) zmask(:,:) = .true.
      CASE ("b")
         vtype=r2dvar
         where(rmask(:,:) > 0) zmask(:,:) = .true.
      CASE ("u")
         vtype=u2dvar
         where(umask(:,:) > 0) zmask(:,:) = .true.
      CASE ("v")
         vtype=v2dvar
         where(vmask(:,:) > 0) zmask(:,:) = .true.
      CASE ("w")
         vtype=r2dvar
         where(rmask(:,:) > 0) zmask(:,:) = .true.
      CASE ("psi")
         vtype=p2dvar
         where(pmask(:,:) > 0) zmask(:,:) = .true.
      END SELECT

      ! From nf_fread.F
      vert_type=vtype/4            ! vertical and horizontal
      horiz_type=vtype-4*vert_type ! grid types, then calculate
      jmin=horiz_type/2            ! starting indices indices
      imin=horiz_type-2*jmin       ! in horizontal directions.

# ifdef MPI
      IF (ii.gt.0) then
         imin=1
      ENDIF
      IF (ii.eq.NP_XI-1) then
         imax=Lmmpi+1
      ELSE
         imax=Lmmpi
      ENDIF
      IF (jj.gt.0) then
         jmin=1
      ENDIF
      IF (jj.eq.NP_ETA-1) then
         jmax=Mmmpi+1
      ELSE
         jmax=Mmmpi
      ENDIF
# else
      imin=1
      imax=Lm+1
      jmin=1
      jmax=Mm+1
# endif

      ! there is probably already a variable for ni and nj 
      ni=imax-imin+1
      nj=jmax-jmin+1

      CALL iom_set_domain_attr(TRIM(cdgrd),
     &     mask=RESHAPE(zmask(imin:imax,jmin:jmax), (/ni*nj/)))
      
      END SUBROUTINE set_mask
#endif

      SUBROUTINE iom_set_domain_attr( cdid, ni_glo, nj_glo, ibegin,
     & jbegin, ni, nj, zoom_ibegin, zoom_jbegin, zoom_ni, zoom_nj,
     & data_dim, data_ibegin, data_ni, data_jbegin, data_nj, 
     & lonvalue, latvalue, mask,nvertex, bounds_lon, bounds_lat, area )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_domain_attr  ***
      !!
      !! ** Purpose :   set domain definition to XIOS library
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                R. Benshila , adaptation for CROCO
      !!                S. Masson   , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*) , INTENT(in) :: cdid
      INTEGER , OPTIONAL, INTENT(in) :: ni_glo, nj_glo, ibegin, jbegin 
      INTEGER , OPTIONAL, INTENT(in) :: ni, nj
      INTEGER , OPTIONAL, INTENT(in) :: data_dim, data_ibegin, data_ni
      INTEGER , OPTIONAL, INTENT(in) :: data_jbegin, data_nj
      INTEGER , OPTIONAL, INTENT(in) :: zoom_ibegin, zoom_jbegin
      INTEGER , OPTIONAL, INTENT(in) :: zoom_ni, zoom_nj, nvertex
      REAL, DIMENSION(:)   , OPTIONAL, INTENT(in) :: lonvalue, latvalue
      REAL, DIMENSION(:,:) , OPTIONAL, INTENT(in) :: bounds_lon, 
     &                                               bounds_lat, area
      LOGICAL,  DIMENSION(:) , OPTIONAL, INTENT(in) ::   mask
      !!----------------------------------------------------------------------
      !
      IF ( xios_is_valid_domain(cdid) ) THEN
         CALL xios_set_domain_attr( cdid, ni_glo=ni_glo, nj_glo=nj_glo, 
     &    ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,
     &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni,
     &    data_jbegin=data_jbegin, data_nj=data_nj,   
     &    lonvalue_1D=lonvalue, latvalue_1D=latvalue, mask_1D=mask,
     &    nvertex=nvertex, bounds_lon_1D=bounds_lon,       
     &    bounds_lat_1D=bounds_lat, area=area, type='curvilinear')
      ENDIF
      !
      IF ( xios_is_valid_domaingroup(cdid) ) THEN
         CALL xios_set_domaingroup_attr( cdid, ni_glo=ni_glo,
     &    nj_glo=nj_glo, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,
     &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, 
     &    data_jbegin=data_jbegin, data_nj=data_nj ,   
     &    lonvalue_1D=lonvalue, latvalue_1D=latvalue, mask_1D=mask,
     &    nvertex=nvertex, bounds_lon_1D=bounds_lon,                  
     &    bounds_lat_1D=bounds_lat, area=area, type='curvilinear')
      ENDIF
      !
      CALL xios_solve_inheritance()
      ! 
      END SUBROUTINE iom_set_domain_attr

      SUBROUTINE iom_set_zoom_domain_attr( cdid, ibegin, jbegin, ni, nj)
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      INTEGER         , OPTIONAL, INTENT(in) ::   ibegin, jbegin, ni, nj

      IF ( xios_is_valid_zoom_domain(cdid) ) THEN
         CALL xios_set_zoom_domain_attr(cdid, ibegin=ibegin-1, 
     &                                  jbegin=jbegin-1, ni=ni, nj=nj)
      ENDIF
      END SUBROUTINE iom_set_zoom_domain_attr
      
      SUBROUTINE iom_set_axis_attr( cdid, paxis,bounds )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_axis_attr  ***
      !!
      !! ** Purpose :   set axis attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                A. Ponte    , adaptation for CROCO
      !!                S. Masson   , original version (NEMO)
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)      , INTENT(in) ::   cdid
      REAL, DIMENSION(:)  , OPTIONAL, INTENT(in) ::   paxis
      REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) ::   bounds
      !!----------------------------------------------------------------------
      !
      IF ( PRESENT(paxis) ) THEN
      IF ( xios_is_valid_axis (cdid) ) then
       CALL xios_set_axis_attr(cdid,n_glo=SIZE(paxis), value=paxis )
      ENDIF
      !
      IF ( xios_is_valid_axisgroup(cdid) ) then
       CALL xios_set_axisgroup_attr(cdid, n_glo=SIZE(paxis),value=paxis)
      ENDIF
      !
      ENDIF
      IF ( xios_is_valid_axis     (cdid) )  then 
         CALL xios_set_axis_attr     ( cdid, bounds=bounds )
      ENDIF
      IF ( xios_is_valid_axisgroup(cdid) )  then
         CALL xios_set_axisgroup_attr( cdid, bounds=bounds )
      ENDIF
      
      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_axis_attr
 
      SUBROUTINE iom_set_field_attr( cdid, freq_op, freq_offset )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_field_attr  ***
      !!
      !! ** Purpose :   set output field attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                A. Ponte    , adaptation for CROCO 
      !!                S. Masson   , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      TYPE(xios_file) :: file_hdl
      TYPE(xios_variable) :: var_hdl
      LOGICAL ok
      TYPE(xios_duration),OPTIONAL , INTENT(in) ::   freq_op
      TYPE(xios_duration),OPTIONAL , INTENT(in) ::   freq_offset
      !!----------------------------------------------------------------------
      !
      IF ( xios_is_valid_field     (cdid) ) 
     &   CALL xios_set_field_attr   
     &  ( cdid, freq_op=freq_op, freq_offset=freq_offset )
      !
      IF ( xios_is_valid_fieldgroup(cdid) )   
     & CALL xios_set_fieldgroup_attr( 
     &     cdid, freq_op=freq_op, freq_offset=freq_offset )
      !
      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_field_attr

      SUBROUTINE iom_set_file_global_attr_char(cdid, type, name, value)
# include "param.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_file_global_attr_char  ***
      !!
      !! ** Purpose :   set file global attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten : add global attributes , May, 2019

      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*)          , INTENT(in) ::   type
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name, value
      LOGICAL :: ok
      TYPE(xios_file) :: file_hdl
      TYPE(xios_variable) :: var_hdl
      !!----------------------------------------------------------------------

      CALL xios_solve_inheritance()

      IF ( xios_is_valid_file     (cdid) )   THEN
          CALL xios_get_handle(cdid,file_hdl)

          CALL xios_add_child(file_hdl, var_hdl, trim(name))
          CALL xios_set_variable_attr(trim(name), name=name, type=type)
          ok = xios_setvar(trim(name),value)
      ENDIF

      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_file_global_attr_char


      SUBROUTINE iom_set_file_global_attr_float(cdid, type, name, value)
# include "param.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_file_global_attr_char  ***
      !!
      !! ** Purpose :   set file global attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten : add global attributes , May, 2019

      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*)          , INTENT(in) ::   type
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name
      REAL                                   ::   value
      LOGICAL :: ok
      TYPE(xios_file) :: file_hdl
      TYPE(xios_variable) :: var_hdl
      !!----------------------------------------------------------------------

      CALL xios_solve_inheritance()

      IF ( xios_is_valid_file     (cdid) )   THEN
          CALL xios_get_handle(cdid,file_hdl)

          CALL xios_add_child(file_hdl, var_hdl, trim(name))
          CALL xios_set_variable_attr(trim(name), name=name, type=type)
          ok = xios_setvar(trim(name),value)
      ENDIF

      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_file_global_attr_float
  

      SUBROUTINE iom_set_file_global_attr_int(cdid, type, name, value)
# include "param.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_file_global_attr_int  ***
      !!
      !! ** Purpose :   set file global attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten : add global attributes , May, 2019

      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*)          , INTENT(in) ::   type
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name
      INTEGER                                ::   value
      LOGICAL :: ok
      TYPE(xios_file) :: file_hdl
      TYPE(xios_variable) :: var_hdl
      !!----------------------------------------------------------------------

      CALL xios_solve_inheritance()

      IF ( xios_is_valid_file     (cdid) )   THEN
          CALL xios_get_handle(cdid,file_hdl)

          CALL xios_add_child(file_hdl, var_hdl, trim(name))
          CALL xios_set_variable_attr(trim(name), name=name, type=type)
          ok = xios_setvar(trim(name),value)
      ENDIF

      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_file_global_attr_int
  
 
      SUBROUTINE iom_set_file_attr( cdid, name, name_suffix )
# include "param.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_file_attr  ***
      !!
      !! ** Purpose :   set file attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   A. Ponte   , adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name, name_suffix
      !!----------------------------------------------------------------------
      IF ( xios_is_valid_file     (cdid) )  CALL xios_set_file_attr  
     &    ( cdid, name=name, name_suffix=name_suffix )
      !
      IF ( xios_is_valid_filegroup(cdid) )  CALL xios_set_filegroup_attr
     & ( cdid, name=name, name_suffix=name_suffix )
      !
      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_file_attr

      
      SUBROUTINE iom_get_file_attr(cdid, name, name_suffix, output_freq)
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_get_file_attr  ***
      !!
      !! ** Purpose :   get attribute from xml file
      !!
      !! ** Method  :   
      !!
      !! ** History :   R. Benshila, adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in ) :: 
     &   cdid
      CHARACTER(LEN=*),OPTIONAL , INTENT(out) :: 
     &   name, name_suffix
      TYPE(xios_duration)   ,OPTIONAL , INTENT(out) :: output_freq
      LOGICAL                                 :: 
     &   llexist1,llexist2,llexist3
      !---------------------------------------------------------------------
      !
      IF(PRESENT(name       )) name = ''          ! default values
      IF(PRESENT(name_suffix)) name_suffix = ''
      IF(PRESENT(output_freq)) output_freq = xios_duration(0,0,0,0,0,0)
       IF ( xios_is_valid_file     (cdid) ) THEN
         CALL xios_solve_inheritance()
         CALL xios_is_defined_file_attr   
     &      ( cdid, name = llexist1, name_suffix = llexist2,
     &        output_freq = llexist3)
         IF(llexist1)   CALL xios_get_file_attr   
     &      ( cdid, name = name )
         IF(llexist2)   CALL xios_get_file_attr   
     &      ( cdid, name_suffix = name_suffix )
         IF(llexist3)   CALL xios_get_file_attr   
     &      ( cdid, output_freq = output_freq )
      ENDIF

      IF ( xios_is_valid_filegroup(cdid) ) THEN
         CALL xios_solve_inheritance()
         CALL xios_is_defined_filegroup_attr( cdid, name = llexist1, 
     &    name_suffix = llexist2, output_freq = llexist3)
         IF(llexist1)   CALL xios_get_filegroup_attr( cdid, 
     &    name = name )
         IF(llexist2)   CALL xios_get_filegroup_attr( cdid, 
     &    name_suffix = name_suffix )
         IF(llexist3)   CALL xios_get_filegroup_attr( cdid, 
     &    output_freq = output_freq )
      ENDIF
      !
      END SUBROUTINE iom_get_file_attr

      
      SUBROUTINE set_scalar
# include "param.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_scalar  ***
      !!
      !! ** Purpose :   fake axis for scalar output
      !!
      !! ** Method  :   
      !! ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                R. Benshila , adaptation for CROCO 
      !!                S. Masson   , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      REAL, DIMENSION(1) ::   zz = 1.
      !!----------------------------------------------------------------------
      ! 
      CALL iom_set_domain_attr('scalarpoint', ni_glo=NNODES, nj_glo=1,
     &                         ibegin=mynode, jbegin=0, ni=1, nj=1)
      CALL iom_set_domain_attr('scalarpoint', data_dim=2, data_ibegin=1,
     &                         data_ni=1, data_jbegin=1, data_nj=1)
      zz = REAL(mynode+1)
      CALL iom_set_domain_attr('scalarpoint', lonvalue=zz, latvalue=zz)
      !
      END SUBROUTINE set_scalar

#if defined ONLINE_ANALYSIS
      SUBROUTINE set_scalogram_axis( cdid
     &                         ,nsclg_glo
     &                         ,nsclg_loc
     &                         ,sclg_coord1D, isclg_beg
     &                         ,mynode )
      !
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_scalogram_axis  ***
      !!
      !! ** Purpose :   xios axis definition for ONLINE_ANALYSIS scalogram
      !!
      !! ** Method  :   XIOS version >= 2.5 
      !! ** History :   B. Lemieux-Dudon (2021) 
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)      , INTENT(in) ::   cdid
      !> Module OA variable names but passed as args IN 
      INTEGER, INTENT(IN) :: nsclg_glo
      INTEGER, INTENT(IN) :: nsclg_loc, isclg_beg
      INTEGER, INTENT(IN) :: mynode 
      !> Local MPI process/global scalogram indices index_s_oa/index_s_cr
      INTEGER, DIMENSION(0:nsclg_loc-1),   INTENT(IN) :: sclg_coord1D

      INTEGER            :: nj, jbeg !i, j, iper, jsclg, 
      INTEGER            :: dat_nj, dat_jbeg
      REAL, DIMENSION(:), POINTER :: coord_val
      REAL, DIMENSION(:), ALLOCATABLE, TARGET :: sclg_coord1D_r
      INTEGER            :: jsclg
      !!----------------------------------------------------------------------
      ! 

# ifndef XIOS2
      print*,'ERROR OA : in set_scalogram_dom XIOS version>=2.5'
      stop
# else
      if ( .not. allocated(sclg_coord1D_r) ) then 
        if (nsclg_loc>0) then
            allocate( sclg_coord1D_r(0:nsclg_loc-1) )
        else
            allocate( sclg_coord1D_r(0) )
        endif
      endif
      do jsclg  = 0,nsclg_loc-1
        if (nsclg_loc>0) then
         sclg_coord1D_r(jsclg) = REAL(sclg_coord1D(jsclg))
        endif
      enddo
      if ( nsclg_loc > 0 ) then
                jbeg=isclg_beg
                nj=size(sclg_coord1D(:))
                dat_nj=nj
                dat_jbeg=0
                coord_val => sclg_coord1D_r(:)
      else
                jbeg=0
      !initially set to nj=1 but nj=0 supported
                nj=0
                dat_nj=0
                dat_jbeg=0
                coord_val => sclg_coord1D_r(:)
      endif
      if ( xios_is_valid_axis     (cdid) ) then
         CALL xios_set_axis_attr( cdid
     &     ,n_glo=nsclg_glo
     &     ,n=nj
     &     ,begin=jbeg
     &     ,data_n=dat_nj
     &     ,data_begin=dat_jbeg
     &     ,standard_name="scalogram"
     &     ,long_name="scalogram index"
     &     ,value=coord_val)
   
      else
        print*,'ERROR OA : not a valid XIOS axis ',cdid
        stop
      endif
      if ( allocated(sclg_coord1D_r) ) deallocate( sclg_coord1D_r )
# endif
      !
      END SUBROUTINE set_scalogram_axis
#endif  /* ONLINE_ANALYSIS */

#if defined ONLINE_ANALYSIS
      SUBROUTINE set_scalogram_dom( cdid
     &                         ,nsclg_glo
     &                         ,nsclg_loc, nper_sclg_max
     &                         ,sclg_coord1D, isclg_beg
     &                         ,period_coord2D
     &                         ,mynode )
      !
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_scalogram_dom  ***
      !!
      !! ** Purpose :   XIOS domain definition for ONLINE_ANALYSIS scalogram
      !!
      !! ** Method  :   XIOS version >= 2.5 
      !! ** History :   B. Lemieux-Dudon (2021) 
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)      , INTENT(in) ::   cdid
      !> Module OA variable names but passed as args IN 
      INTEGER, INTENT(IN) :: nsclg_glo,  nper_sclg_max
      INTEGER, INTENT(IN) :: nsclg_loc, isclg_beg
      INTEGER, INTENT(IN) :: mynode 
      !> Local MPI process/global scalogram indices index_s_oa/index_s_cr
      INTEGER, DIMENSION(0:nsclg_loc-1),   INTENT(IN) :: sclg_coord1D
      REAL, DIMENSION(nper_sclg_max,nsclg_loc),
     &                        TARGET,INTENT(IN) :: period_coord2D

      INTEGER            :: nj, jbeg
      INTEGER            :: dat_nj, dat_jbeg
      REAL, DIMENSION(:,:), POINTER :: lonval, latval
      REAL, DIMENSION(:,:), ALLOCATABLE, TARGET :: sclg_coord2D
      INTEGER            :: iper, jsclg
      !!----------------------------------------------------------------------
      ! 
# ifndef XIOS2
      print*,'ERROR OA : in set_scalogram_dom XIOS version>=2.5'
      stop
# else
      if ( .not. allocated(sclg_coord2D) ) then 
        if (nsclg_loc>0) then
            allocate( sclg_coord2D(nper_sclg_max,0:nsclg_loc-1) )
        else
            allocate( sclg_coord2D(nper_sclg_max,0) )
        endif
      endif
      do jsclg  = 0,nsclg_loc-1
       do iper = 1,nper_sclg_max
        if (nsclg_loc>0) then
         sclg_coord2D(iper,jsclg) = REAL(sclg_coord1D(jsclg))
        endif
       enddo
      enddo
      !zero(nper_sclg_max,1)=0.
      if ( nsclg_loc > 0 ) then
                jbeg=isclg_beg
                nj=size(sclg_coord1D(:))
                dat_nj=nj
                dat_jbeg=0
                lonval => period_coord2D(:,:)
                latval => sclg_coord2D(:,:)
      else
                jbeg=0
      !initially set 2 nj=1 but nj=0 supported
                nj=0
                dat_nj=0
      !suppose to be nj-1
                dat_jbeg=0
                lonval => period_coord2D(:,:)
                latval => sclg_coord2D(:,:)
      endif
      if ( xios_is_valid_domain     (cdid) ) then
         CALL xios_set_domain_attr( cdid
     &     ,ni_glo=nper_sclg_max
     &     ,nj_glo=nsclg_glo
     &     ,ni=nper_sclg_max
     &     ,ibegin=0
     &     ,nj=nj
     &     ,jbegin=jbeg
     &     ,data_dim=2
     &     ,data_ni=nper_sclg_max
     &     ,data_ibegin=0
     &     ,data_nj=dat_nj
     &     ,data_jbegin=dat_jbeg
     &     ,lon_name="period"
     &     ,lat_name="scalogram"
     &     ,lonvalue_2d=lonval
     &     ,latvalue_2d=latval)
   
      else
        print*,'ERROR OA : not a valid XIOS domain ',cdid
        stop
      endif
      if ( allocated(sclg_coord2D) ) deallocate( sclg_coord2D )
# endif
      !
      END SUBROUTINE set_scalogram_dom
#endif  /* ONLINE_ANALYSIS */

      SUBROUTINE set_xml_global_att(cdid)
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "strings.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_xml_global_att  ***
      !!
      !! ** Purpose :   automatic definitions of global attributs
      !!
      !! ** Method  :
      !!
      !! ** History :   S. Theetten, (May, 2019)

      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=40)  ::   char_value

      char_value=trim(title)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='title', value=char_value)
      char_value=trim(rstname)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='rst_file', value=char_value)
#ifndef ANA_GRID
      char_value=trim(grdname)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='grd_file', value=char_value)
#endif
#ifndef ANA_INITIAL
      char_value=trim(ininame)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='ini_file', value=char_value)
#endif
#if !defined ANA_SMFLUX || !defined ANA_STFLUX \
      || !defined ANA_BTFLUX \
      || (defined BBL && !defined ANA_BSEDIM)  \
      || (defined BBL && !defined ANA_WWAVE)   \
      || (defined SALINITY && !defined ANA_SSFLUX)  \
      || ((defined LMD_SKPP || defined LMD_BKPP)  && !defined ANA_SRFLUX)
      char_value=trim(frcname)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='frc_file', value=char_value)
#endif
#ifdef PSOURCE_NCFILE
C qbar_file
#endif
#ifdef ASSIMILATION
C ass_file
C apar_file
#endif
#ifdef SOLVE3D
# ifdef  NEW_S_COORD
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='VertCoordType', value='NEW')
# endif
# ifdef  LMD_SKPP2005
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='skpp', value='2005')
# endif
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='theta_s', value=theta_s)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='theta_s_expl'
     &     ,value='S-coordinate surface control parameter')
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='theta_b', value=theta_b)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='theta_b_expl'
     &     , value='S-coordinate bottom control parameter')
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='Tcline', value=Tcline)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='Tcline_expl'
     &     , value='S-coordinate surface/bottom layer width')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='Tcline_units', value='meter')
      CALL iom_set_file_global_attr_float   (cdid
     &      ,type='float',name='hc', value=hc)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='hc_expl'
     &     , value='S-coordinate parameter, critical depth')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='hc_units', value='meter')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='Cs_w'
     &     , value='have a look at variable Cs_w in this file')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='Cs_rho'
     &     , value='have a look at variable Cs_rho in this file')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='s_w'
     &     , value='have a look at variable s_w in this file')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='s_rho'
     &     , value='have a look at variable s_rho in this file')     
#endif
      CALL iom_set_file_global_attr_int   (cdid
     &     ,type='int',name='ntimes', value=ntimes)
      CALL iom_set_file_global_attr_int   (cdid
     &     ,type='int',name='ndtfast', value=ndtfast)
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='dt', value=dt)
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='dtfast', value=dtfast)
      CALL iom_set_file_global_attr_int   (cdid
     &      ,type='int',name='nwrt', value=nwrt)
#ifdef UV_VIS2
#endif
#ifdef UV_VIS4
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
# endif
# ifdef TS_DIF4
C       CALL iom_set_file_global_attr_float   (cdid
C      &     ,type='float',name='tnu4', value=tnu4)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='tnu4_expl'
     &     , value='biharmonic mixing coefficient for tracers')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='units', value='meter4 second-1')
# endif
# if !defined LMD_MIXING
!
! Background vertical viscosity and tracer mixing coefficients.
!
# endif
#endif
!
! Bottom drag coefficients.
!
C      if (maxval(Zob).ne.0.) then
C 		  elseif (rdrg2.gt.0.) then
C 			  elseif (rdrg.ne.0) then
 
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='rdrg', value=rdrg)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='rdrg_expl'
     &     , value='linear drag coefficient')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='rdrg_units'
     &     , value='meter second-1')
C 				  endif
#ifdef SOLVE3D
!
! Equation of State parameters.
!
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='rho0', value=rho0)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='rho0_expl'
     &     ,value='Mean density used in Boussinesq approximation')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='rho0_units'
     &     ,value='kilogram meter-3')
# ifndef NONLIN_EOS
# endif

! Various parameters.
!
# ifdef BODYFORCE
# endif
#endif /* SOLVE3D */
!
! Slipperiness parameters.
!
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='gamma2', value=gamma2)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='gamma2_expl', value=
     &          'Slipperiness parameter')
# ifdef SPONGE
!
! Sponge parameters
!
      CALL iom_set_file_global_attr_float   (cdid
     &      ,type='float',name='x_sponge', value=x_sponge)
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='v_sponge', value=v_sponge)
      CALL iom_set_file_global_attr_char   (cdid,
     &     type='string',name='sponge_expl',value=
     &     'Sponge parameters : extent (m) & viscosity (m2.s-1)')
# endif
# ifdef SEDIMENT
!
! Sediment parameters
!
# endif
!
! List of Source Codes and Activated CPP-switches
!
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='SRCS', value=srcs)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='CPP-options', value=Coptions)

      END SUBROUTINE set_xml_global_att

      SUBROUTINE set_xmlatt
# include "param.h"
# include "scalars.h"

      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_xmlatt  ***
      !!
      !! ** Purpose :   automatic definitions of some of the xml attributs...
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten, adaptation for XIOS-2.5 (May, 2019)
      !!                R. Benshila, adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=1),DIMENSION( 3) ::   clgrd                    ! suffix name
      CHARACTER(len=256)             ::   clsuff                   ! suffix name
      CHARACTER(len=1)               ::   cl1                      ! 1 character
      CHARACTER(len=2)               ::   cl2                      ! 2 characters
      CHARACTER(len=3)               ::   cl3                      ! 3 characters
      CHARACTER(LEN=20)              ::   char_value
      INTEGER                        ::   ji, jg                   ! loop counters
      INTEGER                        ::   ix, iy                   ! i-,j- index
      REAL            ,DIMENSION(11) ::   zlontao                  ! longitudes of tao    moorings
      REAL            ,DIMENSION( 7) ::   zlattao                  ! latitudes  of tao    moorings
      REAL            ,DIMENSION( 4) ::   zlonrama                 ! longitudes of rama   moorings
      REAL            ,DIMENSION(11) ::   zlatrama                 ! latitudes  of rama   moorings
      REAL            ,DIMENSION( 3) ::   zlonpira                 ! longitudes of pirata moorings
      REAL            ,DIMENSION( 9) ::   zlatpira                 ! latitudes  of pirata moorings
      !!----------------------------------------------------------------------
      !
      ! frequency of the call of iom_put (attribut: freq_op)
      WRITE(cl1,'(i1)') 1
      ! output file names (attribut: name)
      DO ji = 1, 9
         WRITE(cl1,'(i1)') ji
         CALL iom_update_file_name('file'//cl1)
      END DO
      DO ji = 1, 99
         WRITE(cl2,'(i2.2)') ji
         CALL iom_update_file_name('file'//cl2)
      END DO
      DO ji = 1, 999
         WRITE(cl3,'(i3.3)') ji
         CALL iom_update_file_name('file'//cl3)
      END DO

      ! Zooms...
      clgrd = (/ 'T', 'U', 'W' /)
      DO jg = 1, SIZE(clgrd)     ! grid type
         cl1 = clgrd(jg)
         ! Equatorial section (attributs: jbegin, ni, name_suffix)
         CALL dom_ngb( 0., 0., ix, iy, cl1 )
         CALL iom_set_domain_attr ('Eq'//cl1, 
     &    zoom_jbegin=iy, zoom_ni=LLm0)
         CALL iom_get_file_attr   ('Eq'//cl1, 
     &    name_suffix = clsuff             )
         CALL iom_set_file_attr   ('Eq'//cl1, 
     &    name_suffix = TRIM(clsuff)//'_Eq')
         CALL iom_update_file_name('Eq'//cl1)
      END DO
      ! TAO moorings (attributs: ibegin, jbegin, name_suffix)
      zlontao = (/ 137.0, 147.0, 156.0, 165.0, -180.0, -170.0,
     &  -155.0, -140.0, -125.0, -110.0, -95.0 /)
      zlattao = (/  -8.0,  -5.0,  -2.0,   0.0,    2.0,    5.0,
     &     8.0 /)
      CALL set_mooring( zlontao, zlattao )
      ! RAMA moorings (attributs: ibegin, jbegin, name_suffix)
      zlonrama = (/  55.0,  67.0, 80.5, 90.0 /)
      zlatrama = (/ -16.0, -12.0, -8.0, -4.0, -1.5, 0.0, 1.5,
     &  4.0, 8.0, 12.0, 15.0 /)
      CALL set_mooring( zlonrama, zlatrama )
      ! PIRATA moorings (attributs: ibegin, jbegin, name_suffix)
      zlonpira = (/ -38.0, -23.0, -10.0 /)
      zlatpira = (/ -19.0, -14.0,  -8.0, 0.0, 4.0, 8.0, 12.0,15.0,20.0/)
      CALL set_mooring( zlonpira, zlatpira )
      !
      END SUBROUTINE set_xmlatt


      SUBROUTINE set_mooring( plon, plat)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_mooring  ***
      !!
      !! ** Purpose :   automatic definitions of moorings xml attributs...
      !!
      !! ** Method  :   provide a list of coordinates
      !!
      !! ** History :   R. Benshila, adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      REAL, DIMENSION(:), INTENT(in) ::  plon, plat           ! longitudes/latitudes oft the mooring
      !
      CHARACTER(len=1),DIMENSION(1) ::   clgrd = (/ 'T' /)        ! suffix name
      CHARACTER(len=256)            ::   clname                   ! file name
      CHARACTER(len=256)            ::   clsuff                   ! suffix name
      CHARACTER(len=1)              ::   cl1                      ! 1 character
      CHARACTER(len=6)              ::   clon,clat                ! name of longitude, latitude
      INTEGER                       ::   ji, jj, jg               ! loop counters
      INTEGER                       ::   ix, iy                   ! i-,j- index
      REAL                          ::   zlon, zlat
      !!----------------------------------------------------------------------
      !
      DO jg = 1, SIZE(clgrd)
         cl1 = clgrd(jg)
         DO ji = 1, SIZE(plon)
            DO jj = 1, SIZE(plat)
               zlon = plon(ji)
               zlat = plat(jj)
               ! modifications for RAMA moorings
               IF( zlon ==  67. .AND. zlat ==  15. )   zlon =  65.
               IF( zlon ==  90. .AND. zlat <=  -4. )   zlon =  95.
               IF( zlon ==  95. .AND. zlat ==  -4. )   zlat =  -5.
               ! modifications for PIRATA moorings
               IF( zlon == -38. .AND. zlat == -19. )   zlon = -34.
               IF( zlon == -38. .AND. zlat == -14. )   zlon = -32.
               IF( zlon == -38. .AND. zlat ==  -8. )   zlon = -30.
               IF( zlon == -38. .AND. zlat ==   0. )   zlon = -35.
               IF( zlon == -23. .AND. zlat ==  20. )   zlat =  21.
               IF( zlon == -10. .AND. zlat == -14. )   zlat = -10.
               IF( zlon == -10. .AND. zlat ==  -8. )   zlat =  -6.
               IF( zlon == -10. .AND. zlat ==   4. ) THEN  
                   zlon = 0.   ;   zlat = 0.   
               ENDIF
               CALL dom_ngb( zlon, zlat, ix, iy, cl1 )
               IF( zlon >= 0. ) THEN  
                  IF( zlon == REAL(NINT(zlon)) ) THEN   
                  WRITE(clon, '(i3,  a)') NINT( zlon), 'e'
                  ELSE                                      
                  WRITE(clon, '(f5.1,a)')       zlon , 'e'
                  ENDIF
               ELSE             
                  IF( zlon == REAL(NINT(zlon)) ) THEN  
                  WRITE(clon, '(i3,  a)') NINT(-zlon), 'w'
                  ELSE                                     
                  WRITE(clon, '(f5.1,a)')      -zlon , 'w'
                  ENDIF
               ENDIF
               IF( zlat >= 0. ) THEN  
                  IF( zlat == REAL(NINT(zlat)) ) THEN  
                  WRITE(clat, '(i2,  a)') NINT( zlat), 'n'
                  ELSE                                      
                  WRITE(clat, '(f4.1,a)')       zlat , 'n'
                  ENDIF
               ELSE             
                  IF( zlat == REAL(NINT(zlat)) ) THEN   
                  WRITE(clat, '(i2,  a)') NINT(-zlat), 's'
                  ELSE                                      
                  WRITE(clat, '(f4.1,a)')      -zlat , 's'
                  ENDIF
               ENDIF
               clname = TRIM(ADJUSTL(clat))//TRIM(ADJUSTL(clon))
               CALL iom_set_domain_attr 
     &           (TRIM(clname)//cl1, zoom_ibegin= ix, zoom_jbegin= iy)
               CALL iom_get_file_attr   
     &           (TRIM(clname)//cl1, name_suffix = clsuff)
               CALL iom_set_file_attr   
     &           (TRIM(clname)//cl1, name_suffix = TRIM(clsuff)
     &            //'_'//TRIM(clname))
               CALL iom_update_file_name
     &           (TRIM(clname)//cl1)
            END DO
         END DO
      END DO
      ! 
      END SUBROUTINE set_mooring


      SUBROUTINE dom_ngb( plon, plat, kii, kjj, cdgrid )
# include "param.h"
# include "grid.h"
# ifdef MPI
#      include "mpi_cpl.h"
       include 'mpif.h'
# endif
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dom_ngb  ***
      !!
      !! ** Purpose :   find the closest grid point from a given lon/lat position
      !!
      !! ** Method  :   look for minimum distance in cylindrical projection 
      !!                -> not good if located at too high latitude...
      !!
      !!** History  : R. Benshila, adaptation for CROCO 
      !!              S. Masson  , original version (NEMO)  
      !!
      !----------------------------------------------------------------------
      REAL            , INTENT(in   ) ::   plon, plat   ! longitude,latitude of the point
      INTEGER         , INTENT(  out) ::   kii, kjj     ! i-,j-index of the closes grid point
      CHARACTER(len=1), INTENT(in   ) ::   cdgrid       ! grid name 'T', 'U', 'V', 'W'
      !
      INTEGER , DIMENSION(2) ::   iloc
      REAL                   ::   zlon, zmini
      REAL,DIMENSION(GLOBAL_2D_ARRAY) ::  zglam, zgphi, 
     &                                    zmask, zdist
      REAL, DIMENSION(2,1) ::   zain,zaout
      INTEGER :: ierror
      !!--------------------------------------------------------------------
      !
      zmask(:,:) = 0.
# ifdef SPHERICAL
      SELECT CASE( cdgrid )
      CASE( 'U' ) 
      zglam(:,:) = lonu(:,:) 
      zgphi(:,:) = latu(:,:)
#  ifdef MASKING
      zmask(:,:) = umask(:,:)
#  endif
      CASE( 'V' ) 
      zglam(:,:) = lonv(:,:) 
      zgphi(:,:) = latv(:,:) 
#  ifdef MASKING
      zmask(:,:) = vmask(:,:)
#  endif
      CASE DEFAULT 
      zglam(:,:) = lonr(:,:) 
      zgphi(:,:) = latr(:,:) ; 
#  ifdef MASKING
      zmask(:,:) = rmask(:,:)
#  endif
      END SELECT
# else
      SELECT CASE( cdgrid )
      CASE( 'U' ) 
      zglam(:,:) = yr(:,:) 
      zgphi(:,:) = xp(:,:)
#  ifdef MASKING
      zmask(:,:) = umask(:,:)
#  endif
      CASE( 'V' ) 
      zglam(:,:) = yp(:,:) 
      zgphi(:,:) = xr(:,:) 
# ifdef MASKING
      zmask(:,:) = vmask(:,:)
# endif
      CASE DEFAULT 
      zglam(:,:) = yr(:,:) 
      zgphi(:,:) = xr(:,:) ; 
#  ifdef MASKING
      zmask(:,:) = rmask(:,:)
#  endif
      END SELECT
# endif

# ifdef SPHERICAL
      zlon       = MOD( plon       + 720., 360. )                       ! plon between    0 and 360
      zglam(:,:) = MOD( zglam(:,:) + 720., 360. )                       ! glam between    0 and 360
      IF( zlon > 270. ) zlon = zlon - 360.                              ! zlon between  -90 and 270
      IF( zlon <  90. ) 
     &      WHERE( zglam(:,:) > 180. ) zglam(:,:) = zglam(:,:) - 360.   ! glam between -180 and 180
# endif
      zglam(:,:) = zglam(:,:) - zlon
      zgphi(:,:) = zgphi(:,:) - plat
      zdist(:,:) = zglam(:,:) * zglam(:,:) + zgphi(:,:) * zgphi(:,:)
      
# ifdef MPI
      zmini = MINVAL( zdist(:,:) , mask= zmask == 1.e0 )
      iloc = MINLOC( zdist(:,:) , mask= zmask == 1.e0 )
      !
      kii = iloc(1) + iminmpi - 1
      kjj = iloc(2) + jminmpi - 1
      !
      zain(1,:)=zmini
      zain(2,:)=kii+10000.*kjj
      !
      CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,
     &                    MPI_MINLOC,MPI_COMM_WORLD,ierror)
      !
      zmini = zaout(1,1)
      kii = INT(zaout(2,1)/10000.)
      kjj = INT(zaout(2,1) - 10000.*kjj )
# else
         iloc(:) = MINLOC( zdist(:,:), mask = zmask(:,:) == 1.e0 )
         kii = iloc(1) - 1
         kjj = iloc(2) - 1
# endif
      !
      END SUBROUTINE dom_ngb


      SUBROUTINE iom_update_file_name( cdid )
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "strings.h"
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE iom_update_file_name  ***
      !!
      !! ** Purpose :   update output file name from model information
      !!
      !! ** Method  :   request xml attribute and modify it
      !!
      !! ** History :   S. Theetten, adaptation for XIOS-2.5 (May, 2019)
      !!                R. Benshila, adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!---------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      !
      CHARACTER(LEN=256) ::   clname
      CHARACTER(LEN=20)  ::   clfreq
      CHARACTER(LEN=20)  ::   cldate
      CHARACTER(LEN=20)  ::   type, name, value, nom, valeur
      CHARACTER(LEN=20)  ::   char_value
      INTEGER            ::   idx
      INTEGER            ::   jn
      INTEGER            ::   itrlen
      INTEGER            ::   lvar,lenstr
      INTEGER            ::   iyear, imonth, iday, isec
      REAL               ::   zsec, jday0,rvalue
      LOGICAL            ::   llexist
      TYPE(xios_duration)   ::   output_freq
      !!----------------------------------------------------------------------
      !
# if !defined USE_CALENDAR
      ! start_date with format yyyy-mm-dd
       READ(xios_origin_date(1:4),fmt='(i4)') iyear
       READ(xios_origin_date(6:7),fmt='(i2)') imonth
       READ(xios_origin_date(9:10),fmt='(i2)') iday
       CALL ymds2ju( iyear, imonth,iday, 0., jday0 )
# else
       CALL ymds2ju(1900, 1, 1, 0., jday0 ) 
# endif
      DO jn = 1,2

         IF( jn == 1 )   CALL iom_get_file_attr( cdid, name 
     &           = clname, output_freq = output_freq )
         IF( jn == 2 )   CALL iom_get_file_attr( cdid, 
     &                                           name_suffix = clname )

         IF ( TRIM(clname) /= '' ) THEN 

            idx = INDEX(clname,'@expname@') + INDEX(clname,'@EXPNAME@')
            DO WHILE ( idx /= 0 ) 
               clname = clname(1:idx-1)//
     &                  TRIM(replace_blank(TRIM(title),"_"))//
     &                  clname(idx+9:LEN_TRIM(clname))
               idx = INDEX(clname,'@expname@') +  
     &               INDEX(clname,'@EXPNAME@')
            END DO

            idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            DO WHILE ( idx /= 0 )
              IF ( output_freq%timestep /= 0) THEN
                  WRITE(clfreq,'(I18,A2)')INT(output_freq%timestep),'ts'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%second /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%second),'s'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%minute /= 0 ) THEN
                  WRITE(clfreq,'(I18,A2)')INT(output_freq%minute),'mi'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%hour /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%hour),'h'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%day /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%day),'d'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%month /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%month),'m'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%year /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%year),'y'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE
                  if(mynode.eq.0) print*,'error in the name of file id '
     &           //TRIM(cdid),
     & ' attribute output_freq is undefined -> cannot replace @freq@'//
     & ' in '//TRIM(clname) 
                  STOP
              ENDIF
              clname = clname(1:idx-1)//TRIM(ADJUSTL(clfreq))//
     &        clname(idx+6:LEN_TRIM(clname))
              idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            END DO

            idx = INDEX(clname,'@startdate@') + 
     &            INDEX(clname,'@STARTDATE@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( tdays+jday0  )
               clname = clname(1:idx-1)//TRIM(cldate)//
     &          clname(idx+11:LEN_TRIM(clname))
               idx = INDEX(clname,'@startdate@') + 
     &               INDEX(clname,'@STARTDATE@')
            END DO

            idx = INDEX(clname,'@startdatefull@') + 
     &            INDEX(clname,'@STARTDATEFULL@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate(tdays+jday0 , ldfull = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)//
     &          clname(idx+15:LEN_TRIM(clname))
               idx = INDEX(clname,'@startdatefull@') +
     &               INDEX(clname,'@STARTDATEFULL@')
            END DO

            idx = INDEX(clname,'@enddate@') +
     &            INDEX(clname,'@ENDDATE@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( tdays+jday0 + dt
     &          / day2sec * REAL( ntimes - ntstart ), ld24 = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)
     &          //clname(idx+9:LEN_TRIM(clname))
               idx = INDEX(clname,'@enddate@') + 
     &               INDEX(clname,'@ENDDATE@')
            END DO

            idx = INDEX(clname,'@enddatefull@') + 
     &            INDEX(clname,'@ENDDATEFULL@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( tdays+jday0 + dt  
     &          /day2sec * REAL( ntimes - ntstart ), ld24 = .TRUE., 
     &           ldfull = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)
     &          //clname(idx+13:LEN_TRIM(clname))
               idx = INDEX(clname,'@enddatefull@') +
     &               INDEX(clname,'@ENDDATEFULL@')
            END DO

            IF( jn == 1 ) CALL iom_set_file_attr( cdid,
     &                                            name  = clname )
            IF( jn == 2 ) CALL iom_set_file_attr( cdid,
     &                                            name_suffix = clname )
 
            CALL set_xml_global_att(cdid)

         ENDIF
         !
      END DO
      !
      END SUBROUTINE iom_update_file_name


      FUNCTION iom_sdate( pjday, ld24, ldfull )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_sdate  ***
      !!
      !! ** Purpose :   send back the date corresponding to the given julian day
      !!
      !! ** Method  :   
      !!
      !! ** History :   R. Benshila, adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!---------------------------------------------------------------------
      REAL    , INTENT(in   )           ::   pjday         ! julian day
      LOGICAL , INTENT(in   ), OPTIONAL ::   ld24          ! true to force 24:00 instead of 00:00
      LOGICAL , INTENT(in   ), OPTIONAL ::   ldfull        ! true to get the compleate date: yyyymmdd_hh:mm:ss
      !
      CHARACTER(LEN=20) ::   iom_sdate
      CHARACTER(LEN=50) ::   clfmt                         !  format used to write the date
      INTEGER           ::   iyear, imonth, iday, ihour, iminute, isec
      REAL              ::   zsec
      LOGICAL           ::   ll24, llfull
      !!---------------------------------------------------------------------
      !
      IF( PRESENT(ld24) ) THEN   ;   ll24 = ld24
      ELSE                       ;   ll24 = .FALSE.
      ENDIF

      IF( PRESENT(ldfull) ) THEN   ;   llfull = ldfull
      ELSE                         ;   llfull = .FALSE.
      ENDIF

      CALL ju2ymds( pjday, iyear, imonth, iday, zsec )
      isec = NINT(zsec)

      IF( ll24 .AND. isec == 0 ) THEN   ! 00:00 of the next day -> move to 24:00 of the current day
         CALL ju2ymds( pjday - 1., iyear, imonth, iday, zsec )
         isec = 86400
      ENDIF

      IF( iyear < 10000 ) THEN  
      clfmt = "i4.4,2i2.2"                ! format used to write the date
      ELSE                   
      WRITE(clfmt, "('i',i1,',2i2.2')") INT(LOG10(REAL(iyear))) + 1
      ENDIF

      IF( llfull ) THEN
         clfmt = TRIM(clfmt)//",'_',3i2.2"
         ihour   = isec / 3600
         isec    = MOD(isec, 3600)
         iminute = isec / 60
         isec    = MOD(isec, 60)
         WRITE(iom_sdate, '('//TRIM(clfmt)//')') 
     &    iyear, imonth, iday, ihour, iminute, isec    ! date of the end of run
      ELSE
         WRITE(iom_sdate, '('//TRIM(clfmt)//')')
     &     iyear, imonth, iday                          ! date of the end of run
      ENDIF

      END FUNCTION iom_sdate

 
      SUBROUTINE ymds2ju (year,month,day,sec,julian)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE ymds2ju  ***
      !!
      !! ** Purpose :   send back the julian day corresponding to the given date
      !!
      !! ** Method  :   
      !!
      !! ** History :   R. Benshila, adaptation for CROCO 
      !!                IPSL       , original version  
      !!
      !!---------------------------------------------------------------------
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: year,month,day
      REAL,INTENT(IN)    :: sec
      REAL,INTENT(OUT) :: julian
      !
      INTEGER :: julian_day
      REAL    :: julian_sec
      REAL,PARAMETER :: one_day = 86400.0    
      !---------------------------------------------------------------------
      !
      CALL ymds2ju_internal (year,month,day,sec,julian_day,julian_sec)
      julian = julian_day+julian_sec/one_day
      !
      END SUBROUTINE ymds2ju


      SUBROUTINE ymds2ju_internal (year,month,day,sec,julian_day,
     &                             julian_sec)
      !---------------------------------------------------------------------
      !- Converts year, month, day and seconds into a julian day
      !-
      !- In 1968 in a letter to the editor of Communications of the ACM
      !- (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
      !- and Thomas C. Van Flandern presented such an algorithm.
      !-
      !- See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm
      !-
      !- In the case of the Gregorian calendar we have chosen to use
      !- the Lilian day numbers. This is the day counter which starts
      !- on the 15th October 1582.
      !- This is the day at which Pope Gregory XIII introduced the
      !- Gregorian calendar.
      !- Compared to the true Julian calendar, which starts some
      !- 7980 years ago, the Lilian days are smaler and are dealt with
      !- easily on 32 bit machines. With the true Julian days you can only
      !- the fraction of the day in the real part to a precision of
      !- a 1/4 of a day with 32 bits.
      !---------------------------------------------------------------------
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: year,month,day
      REAL,INTENT(IN)    :: sec
      INTEGER,INTENT(OUT) :: julian_day
      REAL,INTENT(OUT)    :: julian_sec
      !
      REAL,PARAMETER :: one_day = 86400.0    
      REAL,PARAMETER :: one_year = 365.2425
      INTEGER :: mon_len(12)=(/31,28,31,30,3
     &     1,30,31,31,30,31,30,31/)
      INTEGER :: jd,m,y,d,ml
      !---------------------------------------------------------------------
      !
      m = month
      y = year
      d = day
      !
      ! We deduce the calendar from the length of the year as it
      ! is faster than an INDEX on the calendar variable.
      !
      IF ( (one_year > 365.0).AND.(one_year < 366.0) ) THEN
      !-- "Gregorian"
       jd = (1461*(y+4800+INT((m-14)/12)))/4 
     &      +(367*(m-2-12*(INT((m-14)/12))))/12 
     &      -(3*((y+4900+INT((m-14)/12))/100))/4 
     &      +d-32075
       jd = jd-2299160
      ELSE IF (    (ABS(one_year-365.0) <= EPSILON(one_year))  
     &  .OR.(ABS(one_year-366.0) <= EPSILON(one_year)) ) THEN
      !-- "No leap" or "All leap"
      ml = SUM(mon_len(1:m-1))
      jd = y*NINT(one_year)+ml+(d-1)
      ELSE
      !-- Calendar with regular month
       ml = NINT(one_year/12.)
       jd = y*NINT(one_year)+(m-1)*ml+(d-1)
      ENDIF
      !
      julian_day = jd
      julian_sec = sec
      !
      END SUBROUTINE ymds2ju_internal

     
      SUBROUTINE ju2ymds (julian,year,month,day,sec)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE ymds2ju  ***
      !!
      !! ** Purpose :   send back the date corresponding to the given julian day
      !!
      !! ** Method  :   
      !!
      !! ** History :   R. Benshila, adaptation for CROCO 
      !!                IPSL       , original version  
      !!
      !!---------------------------------------------------------------------
      IMPLICIT NONE
      !
      REAL,INTENT(IN) :: julian
      INTEGER,INTENT(OUT) :: year,month,day
      REAL,INTENT(OUT)    :: sec
      !
      INTEGER :: julian_day
      REAL    :: julian_sec
      REAL,PARAMETER :: one_day = 86400.0
      !---------------------------------------------------------------------
      !   
      julian_day = INT(julian)
      julian_sec = (julian-julian_day)*one_day
      !
      CALL ju2ymds_internal(julian_day,julian_sec,year,month,day,sec)
      !
      END SUBROUTINE ju2ymds


      SUBROUTINE ju2ymds_internal (julian_day,julian_sec,year,
     &  month,day,sec)
      !---------------------------------------------------------------------
      !- This subroutine computes from the julian day the year,
      !- month, day and seconds
      !-
      !- In 1968 in a letter to the editor of Communications of the ACM
      !- (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
      !- and Thomas C. Van Flandern presented such an algorithm.
      !-
      !- See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm
      !-
      !- In the case of the Gregorian calendar we have chosen to use
      !- the Lilian day numbers. This is the day counter which starts
      !- on the 15th October 1582. This is the day at which Pope
      !- Gregory XIII introduced the Gregorian calendar.
      !- Compared to the true Julian calendar, which starts some 7980
      !- years ago, the Lilian days are smaler and are dealt with easily
      !- on 32 bit machines. With the true Julian days you can only the
      !- fraction of the day in the real part to a precision of a 1/4 of
      !- a day with 32 bits.
      !---------------------------------------------------------------------
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: julian_day
      REAL,INTENT(IN)    :: julian_sec
      INTEGER,INTENT(OUT) :: year,month,day
      REAL,INTENT(OUT)    :: sec
      ! 
      INTEGER :: l,n,i,jd,j,d,m,y,ml
      INTEGER :: add_day
      REAL :: eps_day
      REAL,PARAMETER :: one_day = 86400.0
      REAL,PARAMETER :: one_year = 365.2425
      INTEGER :: mon_len(12)=(/31,28,31,30,3
     &     1,30,31,31,30,31,30,31/)
      !---------------------------------------------------------------------
      !	 
      eps_day = SPACING(one_day)
      !
      jd = julian_day
      sec = julian_sec
      IF (sec > (one_day-eps_day)) THEN
        add_day = INT(sec/one_day)
        sec = sec-add_day*one_day
        jd = jd+add_day
      ENDIF
      IF (sec < -eps_day) THEN
        sec = sec+one_day
        jd = jd-1
      ENDIF
      !
      IF ( (one_year > 365.0).AND.(one_year < 366.0) ) THEN
        !-- Gregorian
        jd = jd+2299160
        !
        l = jd+68569
        n = (4*l)/146097
        l = l-(146097*n+3)/4
        i = (4000*(l+1))/1461001
        l = l-(1461*i)/4+31
        j = (80*l)/2447
        d = l-(2447*j)/80
        l = j/11
        m = j+2-(12*l)
        y = 100*(n-49)+i+l
      ELSE IF ((ABS(one_year-365.0) <= EPSILON(one_year)) 
     &    .OR. (ABS(one_year-366.0) <= EPSILON(one_year))) THEN
        !-- No leap or All leap
        y = jd/NINT(one_year)
        l = jd-y*NINT(one_year)
        m = 1
        ml = 0
        DO WHILE (ml+mon_len(m) <= l)
          ml = ml+mon_len(m)
          m = m+1
        ENDDO
        d = l-ml+1
      ELSE
        !-- others
        ml = NINT(one_year/12.)
        y = jd/NINT(one_year)
        l = jd-y*NINT(one_year)
        m = (l/ml)+1
        d = l-(m-1)*ml+1
      ENDIF
      !
      day = d
      month = m
      year = y
      !
      END SUBROUTINE ju2ymds_internal

   
      FUNCTION replace_blank(string,rep) RESULT(outs)
      !!----------------------------------------------------------------------
      !!                     ***  FUNCTION replace_blank  ***
      !!
      !! ** Purpose :   replace blank strin by a given pattern
      !!
      !! ** Method  :   
      !!
      !! ** History :   R. Benshila
      !!
      !!---------------------------------------------------------------------
      !
      CHARACTER(len=*) :: string, rep
      CHARACTER(len=100) :: outs
      !!---------------------------------------------------------------------
      INTEGER :: stringLen 
      INTEGER :: last, actual
      !!---------------------------------------------------------------------
      ! 
      outs=string
      stringLen = LEN(string)
      last = 1
      actual = 1

      DO WHILE (actual < stringLen)
        IF (string(last:last) == ' ') THEN
            actual = actual + 1
            outs(last:last) = '_'
            last=actual
        ELSE
            last = last + 1
            IF (actual < last) actual = last
        ENDIF
      END DO
      !
      END FUNCTION replace_blank


      END MODULE xios_module

#include "cppdefs.h"
      SUBROUTINE init_xios(tile,cdnamin)

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::  tile
      CHARACTER(len=*),OPTIONAL,INTENT(IN) :: cdnamin

# include "param.h"

# include "compute_tile_bounds.h"
      
      CALL init_xios_tile (Istr,Iend,Jstr,Jend)
      !
      END SUBROUTINE init_xios

      SUBROUTINE init_xios_tile(Istr,Iend,Jstr,Jend)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_sdate  ***
      !!
      !! ** Purpose :   definition of grid for output
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten, adaptation for XIOS-2.5 (May, 2019)
      !!                R. Benshila, improvement for CROCO 
      !!                A. Ponte   , adaption for ROMS  
      !!                S. Masson  , original version (NEMO)  
      !!
      !!---------------------------------------------------------------------
      !
      USE xios     ! XIOS library
      USE xios_module

#if defined MUSTANG
      USE comMUSTANG, only: nk_nivsed_out
# endif
#if defined ONLINE_ANALYSIS
      use module_oa_variables, only : nzv_oa
      use module_interface_oa, only : nzupd0d_oa, nsclg_loc, isclg_beg
     &                               ,nper_sclg_max
     &                               ,model_mode_glo, mdg, nmod_oa
     &                               ,oa_sclg_grd,oa_sclg_dgrd
     &                               ,oa_sclg_axis, oa_sclg_dom

#ifdef MPI
     &                               ,model_mode_loc, mdl
#endif
     &                               ,if_mpi_oa
     &                               ,scalogram_analysis
     &                               ,if_record_sclg_ijpoints
     &                               ,if_xios_dom_grid

#endif
      !    
      IMPLICIT NONE
      !
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"
# include "ncscrum.h"
# include "mpi_cpl.h"
#ifdef ABL1D
# include "abl1d.h"
#endif

      !
      INTEGER :: Istr,Iend,Jstr,Jend
      !
      INTEGER tile,ji
      INTEGER ilocal_comm
      TYPE(xios_duration) :: dtime    = xios_duration(0,
     &  0, 0, 0, 0, 0)
      TYPE(xios_date)      :: time_origin
      TYPE(xios_context) :: ctx_hdl
      TYPE(xios_field) :: field_hdl
      TYPE(xios_fieldgroup) :: fieldgroup_hdl
      TYPE(xios_file) :: file_hdl
      TYPE(xios_variable) :: var_hdl
      LOGICAL ok
      !
      CHARACTER(len=19) :: cldate
      CHARACTER(len=10) :: clname
      CHARACTER(len=10) :: cdname
      INTEGER           :: imd
      INTEGER iday,imonth,iyear,ihour,iminute,isec
#ifdef USE_CALENDAR
      CHARACTER (len=19)    :: tool_sectodat
#endif

      !!---------------------------------------------------------------------
      !  
# include "compute_auxiliary_bounds.h"

      !
# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif

# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
      !
      cdname='crocox'
      clname = cdname
# ifdef AGRIF
      IF( TRIM(Agrif_CFixed()) /= '0' )
     & clname =TRIM(cdname)//"."//TRIM(Agrif_CFixed())
# endif
# if defined MPI
      CALL xios_context_initialize(TRIM(clname), MPI_COMM_WORLD)
# else
      CALL xios_context_initialize(TRIM(clname), 0)
# endif
      CALL iom_swap( cdname )

      ! calendar parameters

# if !defined USE_CALENDAR
! start_date with format yyyy-mm-dd
      READ(xios_origin_date(1:4),fmt='(i4)') iyear
      READ(xios_origin_date(6:7),fmt='(i2)') imonth
      READ(xios_origin_date(9:10),fmt='(i2)') iday
      ihour=0
      iminute=0
      isec=0
# else
      cldate = tool_sectodat(time)
      CALL tool_decompdate(cldate,iday,imonth,iyear,ihour,iminute,isec)
# endif
      CALL xios_define_calendar( TYPE = "Gregorian", 
     &      start_date  = xios_date(start_year,start_month,start_day,
     &                    start_hour,start_minute,start_second),
     &      time_origin = xios_date(origin_year,origin_month,origin_day,
     &                    origin_hour,origin_minute,origin_second))

      ! horizontal grid definition
      CALL set_scalar

# ifdef SPHERICAL
      CALL set_grid( "rho", lonr, latr )
      CALL set_grid( "b"  , lonr, latr )
      CALL set_grid( "u"  , lonu, latu )
      CALL set_grid( "v"  , lonv, latv )
      CALL set_grid( "w"  , lonr, latr )
      CALL set_grid( "psi", lonu, latv )
# else
      CALL set_grid( "rho", xr, yr )
      CALL set_grid( "b"  , xr, yr )
      CALL set_grid( "u"  , xp, yr )
      CALL set_grid( "v"  , xr, yp )
      CALL set_grid( "w"  , xr, yr )
      CALL set_grid( "psi", xp, yp )
# endif

# ifdef MASKING
      CALL set_mask("rho")
      CALL set_mask("b")
      CALL set_mask("u")
      CALL set_mask("v")
      CALL set_mask("w")
      CALL set_mask("psi")
# endif
 
      ! vertical grid definition
      CALL iom_set_axis_attr( "s_rho", sc_r(1:N) )
      CALL iom_set_axis_attr( "s_w"  , sc_w(0:N) )
# ifdef M3FAST_SEDLAYERS
      CALL iom_set_axis_attr( "s_rho_total", sc_r(-N_sl+1:N) )      
      CALL iom_set_axis_attr( "s_w_total"  , sc_w(-N_sl:N) )
# endif      
# ifdef SEDIMENT
      CALL iom_set_axis_attr( "s_b",(/ (REAL(ji), ji=1,NLAY) /))
# elif defined MUSTANG
      CALL iom_set_axis_attr( "s_b",(/ (REAL(ji), ji=1,nk_nivsed_out) /))
# else
      ! fake axis to avoid changing xml file 
      CALL iom_set_axis_attr( "s_b",(/ (REAL(ji), ji=1,1) /))
# endif

# ifdef ABL1D
      CALL iom_set_axis_attr( "N_abl", zr_abl)
# endif

#ifdef ONLINE_ANALYSIS
      oa_wth_sclg : if ( scalogram_analysis ) then

        oa_dom : if ( if_xios_dom_grid ) then

          oa_mmode_loop : do imd=0,nmod_oa

          if ( .not. xios_is_valid_grid(oa_sclg_dgrd(imd)) ) then
            print*,'ERROR OA : oa_sclg_dgrd XIOS INvalid grid'
            print*,oa_sclg_dgrd(imd)
            stop
          endif

          if  (if_record_sclg_ijpoints) then

           if ( .not. xios_is_valid_grid(oa_sclg_grd(imd)) ) then
             print*,'ERROR OA : oa_sclg_grd XIOS INvalid grid'
             print*,oa_sclg_grd(imd)
             stop
           endif

          endif

          enddo oa_mmode_loop


#ifdef MPI
          distributed_sclg_axis : if ( .not. if_mpi_oa ) then

          oa_mmode_loop0 : do imd=0,nmod_oa

            ! BLXD remove nzv_oa unused var
            CALL set_scalogram_dom( oa_sclg_dom(imd)
     &                           ,nzupd0d_oa(imd)
     &                           ,nsclg_loc(imd)
     &                           ,nper_sclg_max(imd)
     &                           ,mdl(imd)%index_s_oa, isclg_beg(imd)
     &                           ,mdl(imd)%per0d_oa
     &                           ,mynode  )

            CALL xios_solve_inheritance()

            ijkscal : if  (if_record_sclg_ijpoints) then

                CALL set_scalogram_axis( oa_sclg_axis(imd)
     &                           ,nzupd0d_oa(imd)
     &                           ,nsclg_loc(imd)
     &                           ,mdl(imd)%index_s_oa, isclg_beg(imd)
     &                           ,mynode  )
            else ijkscal

             if(xios_is_valid_axis(oa_sclg_axis(imd)))then
              call xios_set_axis_attr( oa_sclg_axis(imd)
     &                                 ,n_glo=1)
             else
!$OMP MASTER
              print*,'WARNING OA : INvalid XIOS axis'
              print*,oa_sclg_axis(imd)
!$OMP END MASTER
              !stop
             endif

            endif ijkscal

          enddo oa_mmode_loop0

          else distributed_sclg_axis  

          oa_mmode_loop1 : do imd=0,nmod_oa

            CALL set_scalogram_dom( oa_sclg_dom(imd)
     &                           ,nzupd0d_oa(imd)
     &                           ,nzupd0d_oa(imd)
     &                           ,nper_sclg_max(imd)
     &                           ,mdg(imd)%index_s_cr, 0
     &                           ,mdg(imd)%per0d_cr
     &                           ,mynode  )

            CALL xios_solve_inheritance()

            ijkscal2 : if  (if_record_sclg_ijpoints) then

                CALL set_scalogram_axis( oa_sclg_axis(imd)
     &                           ,nzupd0d_oa(imd)
     &                           ,nzupd0d_oa(imd)
     &                           ,mdg(imd)%index_s_cr, 0
     &                           ,mynode  )

            else ijkscal2

             if(xios_is_valid_axis(oa_sclg_axis(imd)))then
              call xios_set_axis_attr( oa_sclg_axis(imd)
     &                                 ,n_glo=1)
             else
!$OMP MASTER
              print*,'WARNING OA : INvalid XIOS axis'
              print*,oa_sclg_axis(imd)
!$OMP END MASTER
              !stop
             endif

            endif ijkscal2

          enddo oa_mmode_loop1

          endif distributed_sclg_axis  

          CALL xios_solve_inheritance()

#else

          print*,'ERROR OA : XIOS without MPI not tested'
          print*,'- test => comment stop/compile/restart'
          stop

          oa_mmode_loop2 : do imd=0,nmod_oa

          CALL set_scalogram_dom( oa_sclg_dom(imd)
     &                         ,nzupd0d_oa(imd)
     &                         ,nzupd0d_oa(imd)
     &                         ,nper_sclg_max(imd)
     &                         ,mdg(imd)%index_s_cr, 0
     &                         ,mdg(imd)%per0d_cr
     &                         ,mynode  )

          CALL xios_solve_inheritance()

          ijkscal3 : if  (if_record_sclg_ijpoints) then

            CALL set_scalogram_axis( oa_sclg_axis(imd)
     &                           ,nzupd0d_oa(imd)
     &                           ,nzupd0d_oa(imd)
     &                           ,mdg(imd)%index_s_cr, 0
     &                           ,mynode  )


          else ijkscal3

             if(xios_is_valid_axis(oa_sclg_axis(imd)))then
              call xios_set_axis_attr( oa_sclg_axis(imd)
     &                                 ,n_glo=1)
             else
!$OMP MASTER
              print*,'WARNING OA : INvalid XIOS axis'
              print*,oa_sclg_axis(imd)
!$OMP END MASTER
              !stop
             endif

          endif ijkscal3

          CALL xios_solve_inheritance()

          enddo oa_mmode_loop2

#endif /* MPI */

      !BLXD alternate option to oa_sclg_dom with axis grid 
        !else oa_dom
        endif oa_dom

      else oa_wth_sclg

      oa_mmode_loop4 : do imd=0,nmod_oa

      !BLXD Prevent failure with potential scalogram domain/axis grid
      !     incompletely def. in xml files
      null_sclg_dom : if ( xios_is_valid_domain(oa_sclg_dom(imd)) ) then
         CALL xios_set_domain_attr(oa_sclg_dom(imd)
     &     ,ni_glo=1
     &     ,nj_glo=1
     &     )
      endif null_sclg_dom 

      enddo oa_mmode_loop4


      endif oa_wth_sclg
#endif  /* ONLINE_ANALYSIS */

#ifdef key_pisces
      IF(ALLOCATED(profsed))  CALL iom_set_axis_attr("s_rhoS", profsed)
#endif
      
      ! time axis (for CROCO tools)

      ! automatic definitions of some of the xml attributs
      CALL set_xmlatt

      ! end file definition

      ! Academic test with dt > 1 second are non comptabile with
      ! XIOS, set xios time step unit to 1 model time_step
      if (dt < 1) then
        dtime%second =  1
      else
        dtime%second = dt
      endif
      CALL xios_set_timestep(dtime)

      ! end definition phase
      CALL xios_close_context_definition()
      !
      CALL xios_update_calendar(0)

      END SUBROUTINE init_xios_tile
      
      
      SUBROUTINE iom_swap( cdname )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_swap  ***
      !!
      !! ** Purpose :  swap context between different agrif grid for xmlio_server
      !!---------------------------------------------------------------------
      USE xios
      !
      CHARACTER(len=*), INTENT(in) :: cdname
      TYPE(xios_context) :: roms_hdl
# ifdef AGRIF
      IF( TRIM(Agrif_CFixed()) == '0' ) THEN
        CALL xios_get_handle(TRIM(cdname),roms_hdl)
      ELSE
        CALL xios_get_handle(TRIM(cdname)//"."//
     &                       TRIM(Agrif_CFixed()),roms_hdl)
      ENDIF
# else
      CALL xios_get_handle(TRIM(cdname),roms_hdl)

# endif
      !
      CALL xios_set_current_context(roms_hdl)
      !
      END SUBROUTINE iom_swap


      SUBROUTINE iom_context_finalize( cdname )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_context_finalize  ***
      !!
      !! ** Purpose :  finalization of output definition
      !!---------------------------------------------------------------------
      USE xios
      !
      CHARACTER(LEN=*), INTENT(in) :: cdname
      TYPE(xios_context) :: roms_hdl
      !
      IF( xios_is_valid_context(cdname) ) THEN
         CALL iom_swap( cdname )   ! swap to cdname context
         CALL xios_context_finalize() ! finalize the context
         IF( cdname /= "crocox" ) CALL iom_swap( "crocox" )   ! return back to croco context
      ENDIF
      !
      END SUBROUTINE iom_context_finalize


      subroutine send_xios_diags (tile)

      implicit none

     
# include "param.h"
# ifdef SOLVE3D
#  include "work.h"
#  include "ncscrum.h"
# endif
      integer tile
# include "compute_tile_bounds.h"

      call send_xios_diags_tile (tile,Istr,Iend,Jstr,Jend)

      return
      end


      subroutine send_xios_diags_tile (tile,Istr,Iend,Jstr,Jend)

!
!================================================== John M. Klinck ===
!  Copyright (c) 2000 Rutgers/UCLA                                   !
!================================================ Hernan G. Arango ===
!                                                                    !
!  This routine sends every diagnostics fields to Xios server         
!
!  ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
!                                                                    !
!=====================================================================
!

          USE xios     ! XIOS io
          USE xios_module
#if defined BIOLOGY && defined PISCES
          USE trcwri_pisces
#endif

#if defined MUSTANG
      USE comsubstance, ONLY : name_var, nvpc, ws_part, imud1, imud2
      USE sed_MUSTANG,  ONLY : sed_MUSTANG_outres
      USE comMUSTANG, only: nk_nivsed_out, l_outsed_hsed , 
     &    l_outsed_poro, l_outsed_dzs, l_outsed_temp_sed, l_outsed_salt_sed, 
     &    l_outsed_cv_sed, l_outsed_toce, l_outsed_flx_s2w_w2s,
     &    l_outsed_pephm_fcor, l_outsed_bedload, 
     &    l_outsed_fsusp, l_outsed_frmudsup, l_outsed_dzs_ksmax, 
     &    l_outsed_theoric_active_layer, l_outsed_ero_details,
     &    l_outsed_consolidation, l_outsed_ws,
     &    var3D_cvsed, var3D_dzs,
     &    var3D_TEMP, var3D_SAL, var3D_poro,
     &    var2D_hsed, ksma, tauskin, tauskin_c, tauskin_w,
     &    var2D_toce, var2D_flx_s2w, var2D_flx_w2s,
     &    var2D_pephm_fcor, var2D_flx_bx, var2D_flx_by, 
     &    var2D_bil_bedload, var2D_fsusp, 
     &    var2D_frmudsup, var2D_dzs_ksmax, 
     &    var2D_theoric_active_layer, 
     &    var2D_tero_noncoh, var2D_tero_coh, 
     &    var2D_pct_iter_noncoh, var2D_pct_iter_coh, 
     &    var2D_niter_ero, 
     &    var2D_flx_s2w_coh, var2D_flx_s2w_noncoh, 
     &    var2D_flx_w2s_coh, var2D_flx_w2s_noncoh, 
     &    var2D_flx_bx_int, var2D_flx_by_int, 
     &    var2D_bil_bedload_int, z0sed,
     &    l_dyn_insed, var3Dksed_dtsdzs, 
     &    var3Dksed_loadograv, var3Dksed_permeab, 
     &    var3Dksed_sed_rate, var3Dksed_sigmadjge, 
     &    var3Dksed_sigmapsg, var3Dksed_hinder, 
     &    var3Dksed_stateconsol

#ifdef key_sand2D
      USE comsubstance, ONLY : l_subs2D, nv_grav, nv_sand, l_outsandrouse, rsh
      USE comMUSTANG,   ONLY : sum_tmp, rouse2D
#endif
#endif

#if defined OBSTRUCTION
      USE com_OBSTRUCTIONS, ONLY : l_obstout_pos, l_obstout_height_f,
     &   l_obstout_height_e, l_obstout_dens_f, l_obstout_dens_e,
     &   l_obstout_dens_2d, l_obstout_width_f, l_obstout_width_e, 
     &   l_obstout_width_2d, l_obstout_thick_f, l_obstout_thick_e, 
     &   l_obstout_theta, l_obstout_frac_xy, l_obstout_frac_z, 
     &   l_obstout_fuzvz, l_obstout_a2d, l_obstout_a3d, l_obstout_s2d, 
     &   l_obstout_s3d, l_obstout_drag, l_obstout_tau, 
     &   obst_nout_pos, obst_nout_height_f, 
     &   obst_nout_height_e, obst_nout_dens_f , obst_nout_dens_e, 
     &   obst_nout_width_f, obst_nout_width_e , obst_nout_thick_f, 
     &   obst_nout_thick_e, obst_nout_theta, obst_nout_frac_xy, 
     &   obst_nout_frac_z, obst_nout_fuzvz, obst_nout_a2d, 
     &   obst_nout_a3d, obst_nout_s2d , obst_nout_s3d, 
     &   obst_nout_drag, obst_nout_tau,
     &   obst_varname, obst_position, obst_height_inst, obst_height, 
     &   obst_dens_inst, obst_dens3d, obst_width_inst, obst_width3d,
     &   obst_thick_inst, obst_thick3d, obst_theta3d, obst_fracxy, 
     &   obst_fracz3d, obst_drag3d, obst_a2d, obst_a3d, obst_s2d, obst_s3d, 
     &   obst_tau, obst_fuz, obst_fvz, obst_nbvar
#endif
          IMPLICIT NONE

# include "param.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "forces.h"
# include "mixing.h"
# include "coupling.h"
# ifdef SEDIMENT
#  include "sediment.h"
# endif
# ifdef BBL
#  include "bbl.h"
# endif
# ifdef WKB_WWAVE
#  include "wkb_wwave.h"
# endif
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_TS
#  include "diagnostics.h"
# endif
# if defined DIAGNOSTICS_VRT
#  include "diags_vrt.h"
# endif
# ifdef DIAGNOSTICS_EK
#  include "diags_ek.h"
# endif
# ifdef DIAGNOSTICS_PV
#  include "diags_pv.h"
# endif
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
# ifdef M3FAST
# include "nbq.h"
# endif
# include "scalars.h"
# include "work.h"
# include "private_scratch.h"
# ifdef ABL1D
# include "abl1d.h"
# endif

!
#ifdef DIAGNOSTICS_EDDY
      real workru(GLOBAL_2D_ARRAY,1:N)
      common /work3d_r/ workru
      real workrv(GLOBAL_2D_ARRAY,1:N)
      common /work3d_r/ workrv

      real workrsu(GLOBAL_2D_ARRAY)
      common /work2d_r/ workrsu
      real workrsv(GLOBAL_2D_ARRAY)
      common /work2d_r/ workrsv

      real workrbu(GLOBAL_2D_ARRAY)
      common /work2d_r/ workrbu
      real workrbv(GLOBAL_2D_ARRAY)
      common /work2d_r/ workrbv

      real workrug(GLOBAL_2D_ARRAY)
      common /work2d_r/ workrug
      real workrvg(GLOBAL_2D_ARRAY)
      common /work3d_r/ workrvg
      real GRho
#endif
      integer tile
      integer Istr,Iend,Jstr,Jend
      integer i,j,k
      integer ilc, itrc
      character*20 nametrc
      real eps
      parameter (eps=1.D-20)
      real stf_cff
      real zz(GLOBAL_2D_ARRAY,0:N)
      parameter(stf_cff=86400/0.01)
!
#if defined MUSTANG
      real workrsed(GLOBAL_2D_ARRAY,nk_nivsed_out),
     &     workrsed1(GLOBAL_2D_ARRAY,nk_nivsed_out)

#endif

#if defined OBSTRUCTION
      integer iv
#endif


# include "compute_auxiliary_bounds.h"
!
#  define I_EXT_RANGE Istr-1,Iend+1
#  define J_EXT_RANGE Jstr-1,Jend+1

#  ifdef EW_PERIODIC 
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
#  else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
#  endif

#  ifdef NS_PERIODIC 
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
#  else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
#  endif

#ifdef OPENACC
# include "copy_to_host_xios_diag.h"
#endif

       ilc=1+iic-ntstart   ! number of time step since restart
       ! To be conform with clasic (no xios) croco file first 
       ! time step must be writen
       if (ilc == 0) ilc=1
       call xios_update_calendar(ilc)

      if (xios_field_is_active("time")) then
       call xios_send_field("time",time)
      endif

! static attributes
# ifdef SOLVE3D
!
! Variable dimension on vertical : s_rho and s_w
!
! done automatically
!
! S-coordinate independent variables "sc_w", "sc_rho" and stretching
! curves "Cs_w", "Cs_rho" at W- and RHO-points.
!
      if (xios_field_is_active("Cs_rho")) then
      call xios_send_field("Cs_rho",Cs_r(1:N))
# ifdef M3FAST_SEDLAYERS   
      call xios_send_field("Cs_rho_total",Cs_r(-N_sl+1:N))
# endif      
      endif
      if (xios_field_is_active("Cs_w")) then
      call xios_send_field("Cs_w",Cs_w(0:N))
# ifdef M3FAST_SEDLAYERS   
      call xios_send_field("Cs_w_total",Cs_w(-N_sl:N))
# endif      
      endif
      if (xios_field_is_active("hc")) then
      call xios_send_field("hc",hc)
      endif
      if (xios_field_is_active("theta_s")) then
      call xios_send_field("theta_s",theta_s)
      endif
      if (xios_field_is_active("theta_b")) then
      call xios_send_field("theta_b",theta_b)
      endif
      if (xios_field_is_active("Tcline")) then
      call xios_send_field("Tcline",Tcline)
      endif
      if (xios_field_is_active("Vtransform")) then
      call xios_send_field("Vtransform",REAL(Vtransform))
      endif
c !$acc wait

      if (xios_field_is_active("levels_rho")) then
      call xios_send_field("levels_rho",z_r(:,:,:))
      endif
      if (xios_field_is_active("levels_w")) then
      call xios_send_field("levels_w"  ,z_w(:,:,:))
      endif
      if (xios_field_is_active("levels_u")) then
      do k=1,N
        do j=JstrR,JendR
          do i=Istr,IendR
            workr(i,j,k)=0.5*(z_r(i,j,k)+z_r(i-1,j,k))
          enddo
        enddo
      enddo
      call xios_send_field("levels_u"  ,workr(:,:,:))
      endif
      if (xios_field_is_active("levels_v")) then
      do k=1,N
       do j=JV_RANGE
          do i=IV_RANGE
            workr(i,j,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))
          enddo
       enddo
      enddo
      call xios_send_field("levels_v"  ,workr(:,:,:))
      endif
      
# endif
!
! Bathymetry.
!
      if (xios_field_is_active("h")) then
      call xios_send_field("h",h)
      endif
!
! Coriolis parameter.
!
      if (xios_field_is_active("f")) then
      call xios_send_field("f",f)
      endif
!
! Curvilinear transformation metrics.
!
      if (xios_field_is_active("pm")) then
      call xios_send_field("pm",pm)
      endif
      if (xios_field_is_active("pn")) then
      call xios_send_field("pn",pn)
      endif
!
! Longitude/latitude or cartezian coordinates of RHO-points
!
# ifdef SPHERICAL
      if (xios_field_is_active("lon_rho")) then
      call xios_send_field("lon_rho",lonr)
      endif
      if (xios_field_is_active("lat_rho")) then
      call xios_send_field("lat_rho",latr)
      endif
      if (xios_field_is_active("lon_u")) then
      call xios_send_field("lon_u",lonu)
      endif
      if (xios_field_is_active("lat_u")) then
      call xios_send_field("lat_u",latu)
      endif
      if (xios_field_is_active("lon_v")) then
      call xios_send_field("lon_v",lonv)
      endif
      if (xios_field_is_active("lat_v")) then
      call xios_send_field("lat_v",latv)
      endif
# else  
      if (xios_field_is_active("x_rho")) then
      call xios_send_field("x_rho",xr)
      endif
      if (xios_field_is_active("y_rho")) then
      call xios_send_field("y_rho",yr)
      endif
# endif
# ifdef CURVGRID
      if (xios_field_is_active("angle")) then
      call xios_send_field("angle",angler)
      endif
# endif
# ifdef MASKING
      if (xios_field_is_active("mask_rho")) then
      call xios_send_field("mask_rho",rmask)
      endif
# endif

!
! Critical Depth for Drying cells
!
# ifdef WET_DRY
      if (xios_field_is_active("Dcrit")) then
      call xios_send_field("Dcrit",Dcrit)
      endif
# endif
      if (xios_field_is_active("zeta")) then
!!$acc update host ( zeta )      
      call xios_send_field("zeta",zeta(:,:,fast_indx_out))
      endif
      if (xios_field_is_active("ubar")) then
!!$acc update host ( ubar )              
      call xios_send_field("ubar",ubar(:,:,fast_indx_out))
      endif
      if (xios_field_is_active("vbar")) then
!!$acc update host ( vbar )              
      call xios_send_field("vbar",vbar(:,:,fast_indx_out))
      endif

# ifdef MORPHODYN
      if (xios_field_is_active("hmorph")) then
      call xios_send_field("hmorph",h)
      endif
# endif
# ifdef WET_DRY
      if (xios_field_is_active("rmask_wet")) then
      call xios_send_field("rmask_wet",rmask_wet)
      endif
      if (xios_field_is_active("umask_wet")) then
      call xios_send_field("umask_wet",umask_wet)
      endif
      if (xios_field_is_active("vmask_wet")) then
      call xios_send_field("vmask_wet",vmask_wet)
      endif
#endif

      if (xios_field_is_active("bustr")) then
      call xios_send_field("bustr",bustr*rho0)
      endif
      if (xios_field_is_active("bvstr")) then
      call xios_send_field("bvstr",bvstr*rho0)
      endif
      if (xios_field_is_active("bostr")) then
      do j=0,Mm
         do i=0,Lm
            work2d(i,j)=0.5*sqrt((bustr(i,j)+bustr(i+1,j))**2
     &                           +(bvstr(i,j)+bvstr(i,j+1))**2)
     &                                                    *rho0
         enddo
      enddo
      call xios_send_field("bostr",work2d)
      endif

      if (xios_field_is_active("wstr")) then
      do j=1,Mm
         do i=1,Lm
            work2d(i,j)=0.5*sqrt((sustr(i,j)+sustr(i+1,j))**2
     &                           +(svstr(i,j)+svstr(i,j+1))**2)
     &                                                    *rho0
         enddo
      enddo
      call xios_send_field("wstr",work2d)
      endif

      if (xios_field_is_active("sustr")) then
      work2d=sustr*rho0
      call xios_send_field("sustr",work2d)
      endif
      if (xios_field_is_active("svstr")) then
      work2d=svstr*rho0 
      call xios_send_field("svstr",work2d)
      endif

# ifdef ABL1D
      if (xios_field_is_active("zr_abl")) then
      call xios_send_field("zr_abl",zr_abl)
      endif
      if (xios_field_is_active("zw_abl")) then
      call xios_send_field("zw_abl",zw_abl)
      endif
      if (xios_field_is_active("Hzr_abl")) then
      call xios_send_field("Hzr_abl",Hzr_abl)
      endif
      if (xios_field_is_active("Hzw_abl")) then
      call xios_send_field("Hzw_abl",Hzw_abl)
      endif
      if (xios_field_is_active("pu_dta")) then
      call xios_send_field("pu_dta",pu_dta)
      endif
      if (xios_field_is_active("pv_dta")) then
      call xios_send_field("pv_dta",pv_dta)
      endif
      if (xios_field_is_active("pt_dta")) then
      call xios_send_field("pt_dta",pt_dta)
      endif
      if (xios_field_is_active("pq_dta")) then
      call xios_send_field("pq_dta",pq_dta)
      endif
      if (xios_field_is_active("pgu_dta")) then
      call xios_send_field("pgu_dta",pgu_dta)
      endif
      if (xios_field_is_active("pgv_dta")) then
      call xios_send_field("pgv_dta",pgv_dta)
      endif
      if (xios_field_is_active("u_abl")) then
      call xios_send_field("u_abl",u_abl(:,:,:,nstp))
      endif
      if (xios_field_is_active("v_abl")) then
      call xios_send_field("v_abl",v_abl(:,:,:,nstp))
      endif
      if (xios_field_is_active("t_abl")) then
      call xios_send_field("t_abl",tq_abl(:,:,:,nstp,1))
      endif
      if (xios_field_is_active("q_abl")) then
      call xios_send_field("q_abl",tq_abl(:,:,:,nstp,2))
      endif
      if (xios_field_is_active("tke_abl")) then
      call xios_send_field("tke_abl",tke_abl(:,:,:,nstp))
      endif
      if (xios_field_is_active("mxlm_abl")) then
      call xios_send_field("mxlm_abl",mxlm_abl(:,:,:,nstp))
      endif
      if (xios_field_is_active("mxld_abl")) then
      call xios_send_field("mxld_abl",mxld_abl(:,:,:,nstp))
      endif
      if (xios_field_is_active("avm_abl")) then
      call xios_send_field("avm_abl",avm_abl(:,:,:,nstp))
      endif
      if (xios_field_is_active("avt_abl")) then
      call xios_send_field("avt_abl",avt_abl(:,:,:,nstp))
      endif
      if (xios_field_is_active("ablh")) then
      call xios_send_field("ablh",ablh(:,:,nstp))
      endif
# endif

# ifdef BULK_FLUX
      if (xios_field_is_active("tair")) then
#   ifdef ABL1D
      call xios_send_field("tair",tq_abl(:,:,2,nstp,1))
#   else
      call xios_send_field("tair",tair)
#   endif
      endif
      if (xios_field_is_active("rhum")) then
#   ifdef ABL1D
      call xios_send_field("rhum",tq_abl(:,:,2,nstp,2))
#   else
      call xios_send_field("rhum",rhum)
#   endif
      endif
      if (xios_field_is_active("prate")) then
      call xios_send_field("prate",prate)
      endif
#   ifdef READ_PATM
      if (xios_field_is_active("patm")) then
      call xios_send_field("patm",patm2d)
      endif
#   endif
# endif

# if defined BULK_SM_UPDATE || defined BULK_FLUX 
      if (xios_field_is_active("uwnd")) then
#   ifdef ABL1D
      call xios_send_field("uwnd",u_abl(:,:,2,nstp))
#   else
      call xios_send_field("uwnd",uwnd)
#   endif
      endif
      if (xios_field_is_active("vwnd")) then
#   ifdef ABL1D
      call xios_send_field("vwnd",v_abl(:,:,2,nstp))
#   else
      call xios_send_field("vwnd",vwnd)
#   endif
      endif
# endif
# ifdef WKB_WWAVE
      if (xios_field_is_active("hrm")) then
      call xios_send_field("hrm",hrm(:,:,wstp))
      endif
      if (xios_field_is_active("frq")) then
      call xios_send_field("frq",frq(:,:,wstp))
      endif
      if (xios_field_is_active("wac")) then
      call xios_send_field("wac",wac(:,:,wstp))
      endif
      if (xios_field_is_active("wkx")) then
      call xios_send_field("wkx",wkx(:,:,wstp))
      endif
      if (xios_field_is_active("wke")) then
      call xios_send_field("wke",wke(:,:,wstp))
      endif
      if (xios_field_is_active("epb")) then
      call xios_send_field("epb",wepb)
      endif
      if (xios_field_is_active("epd")) then
      call xios_send_field("epd",wepd)
      endif
#  ifdef WAVE_ROLLER
      if (xios_field_is_active("war")) then
      call xios_send_field("war",war(:,:,wstp))
      endif
      if (xios_field_is_active("epr")) then
      call xios_send_field("epr",wepr)
      endif

#  endif
# endif

# ifdef MRL_WCI
      if (xios_field_is_active("sup")) then
      call xios_send_field("sup",sup)
      endif
      if (xios_field_is_active("ust2d")) then
      call xios_send_field("ust2d",ust2d)
      endif
      if (xios_field_is_active("vst2d")) then
      call xios_send_field("vst2d",vst2d)
      endif
#  ifdef SOLVE3D
      if (xios_field_is_active("ust")) then
      call xios_send_field("ust",ust)
      endif
      if (xios_field_is_active("vst")) then
      call xios_send_field("vst",vst)
      endif
      if (xios_field_is_active("wst")) then
      call xios_send_field("wst",wst)
      endif
      if (xios_field_is_active("Akb")) then
      call xios_send_field("Akb",Akb)
      endif
      if (xios_field_is_active("Akw")) then
      call xios_send_field("Akw",Akw)
      endif
      if (xios_field_is_active("kvf")) then
      call xios_send_field("kvf",kvf)
      endif
      if (xios_field_is_active("calP")) then
      call xios_send_field("calP",calP)
      endif
      if (xios_field_is_active("Kapsrf")) then
      call xios_send_field("Kapsrf",Kapsrf)
      endif
#  endif
# endif

# ifdef OW_COUPLING
      if (xios_field_is_active("hrm")) then
      call xios_send_field("hrm",whrm)
      endif
      if (xios_field_is_active("frq")) then
      call xios_send_field("frq",wfrq)
      endif
      if (xios_field_is_active("wkx")) then
      call xios_send_field("wkx",wwkx)
      endif
      if (xios_field_is_active("wke")) then
      call xios_send_field("wke",wwke)
      endif
      if (xios_field_is_active("epb")) then
      call xios_send_field("epb",wepb)
      endif
      if (xios_field_is_active("epd")) then
      call xios_send_field("epd",wepd)
      endif
#  endif
! ------ 3D variables
! nstp vs nnew : present choice follows set_avg.F
#  ifdef SOLVE3D
      if (xios_field_is_active("u")) then
      call xios_send_field("u",u(:,:,:,nstp))
      endif
      if (xios_field_is_active("v")) then
      call xios_send_field("v",v(:,:,:,nstp)) 
      endif
      if (xios_field_is_active("u_surf")) then
      call xios_send_field("u_surf",u(:,:,N,nstp)) 
      endif
      if (xios_field_is_active("v_surf")) then
      call xios_send_field("v_surf",v(:,:,N,nstp)) 
      endif
# ifdef TEMPERATURE
      if (xios_field_is_active("temp")) then
      call xios_send_field("temp",t(:,:,:,nstp,itemp)) 
      endif
      if (xios_field_is_active("temp_surf")) then
      call xios_send_field("temp_surf",t(:,:,N,nstp,itemp)) 
      endif
# endif
# ifdef SALINITY
      if (xios_field_is_active("salt")) then
      call xios_send_field("salt",t(:,:,:,nstp,isalt))
      endif
      if (xios_field_is_active("salt_surf")) then
      call xios_send_field("salt_surf",t(:,:,N,nstp,isalt))
      endif
# endif /* SALINITY */
! WARNING here passive tracers are missing
! we could easily add them in using all the vname stuff
! but I found it quiet ugly

      if (xios_field_is_active("rho")) then
      call xios_send_field("rho",rho)
      endif
      if (xios_field_is_active("rho_surf")) then
      call xios_send_field("rho_surf",rho(:,:,N))
      endif

# if defined ANA_VMIX || defined LMD_MIXING \
|| defined LMD_SKPP   || defined LMD_BKPP \
|| defined GLS_MIXING || defined UV_VIS_SMAGO_3D
      if (xios_field_is_active("bvf")) then
      call xios_send_field("bvf",bvf)
      endif
# endif

      if (xios_field_is_active("omega")) then
      do k=0,N
         do j=0,Mm+1
            do i=0,Lm+1
               work(i,j,k)= ( We(i,j,k)
#  ifdef VADV_ADAPT_IMP
     &                    +   Wi(i,j,k)
#  endif                  
     &                           ) *pm(i,j)*pn(i,j) 
#  ifdef NBQ_MASS
     &                         /rho_nbq_avg1(i,j,k)
#  endif
            enddo
         enddo
      enddo
      call xios_send_field("omega",work)
      endif

#  ifdef M3FAST
      if (xios_field_is_active("u_nbq")) then
        worknbq=0.
          do k=-N_sl+1,N-1
          do j=0,Mm
            do i=0,Lm
            worknbq(i,j,k+N_sl)=2.*qdmu_nbq(i,j,k)
     &                /(Hz(min(Lm,i+1),j,k)+Hz(i,j,k))
            enddo
          enddo
          enddo
      call xios_send_field(" u_nbq",worknbq)
      endif
      if (xios_field_is_active("v_nbq")) then
        worknbq=0.
          do k=-N_sl+1,N-1
          do j=0,Mm
            do i=0,Lm
            worknbq(i,j,k+N_sl)=2.*qdmv_nbq(i,j,k)
     &                /(Hz(i,min(Mm,j+1),k)+Hz(i,j,k))
            enddo
          enddo
          enddo
      call xios_send_field(" v_nbq",worknbq)
      endif
      if (xios_field_is_active("w_nbq")) then
      !call xios_send_field("w_nbq",wz(:,:,:,nstp))
      
        worknbq=0.
          do k=-N_sl+1,N-1
          do j=0,Mm
            do i=0,Lm
            worknbq(i,j,k+N_sl)=2.*(
     &         qdmw_nbq(i,j,k)
     &        )/(Hz(i,j,k)+Hz(i,j,k+1))
#  ifdef NHINT_WH
!     &               +wzh_nbq(i,j,k)
#  endif
            enddo
          enddo
          enddo
          k=-N_sl
          do j=0,Mm
            do i=0,Lm
#   ifndef M3FAST_SEDLAYERS
            worknbq(i,j,k)=2.*(
     &         qdmw_nbq(i,j,k)
     &        )/Hz(i,j,k+1) 
#  ifdef NHINT_WH
!     &               +wzh_nbq(i,j,k)
#  endif
#   else
            worknbq(i,j,k)=2.*(
     &                qdmw_nbq(i,j,k)
     &                )/(Hz(i,j,k)+Hz(i,j,k+1))
#  ifdef NHINT_WH
!     &               +wzh_nbq(i,j,k)
#  endif
#   endif
            enddo
          enddo
          k=N
          do j=0,Mm
            do i=0,Lm 
            worknbq(i,j,k)=2.*(
     &                qdmw_nbq(i,j,k)
     &                )/Hz(i,j,k)
#  ifdef NHINT_WH
!     &               +wzh_nbq(i,j,k)
#  endif
            enddo
          enddo
     
      call xios_send_field(" w_nbq",worknbq)
      endif
      if (xios_field_is_active("rho_nbq")) then
!$acc update host ( rho_nbq
# ifdef M3FAST_DIAGACOUS              
!$acc&  , p_nbq,p_nbq_max 
# endif              
!$acc&                  )  
      call xios_send_field("rho_nbq",rho_nbq)
# ifdef M3FAST_DIAGACOUS              
      call xios_send_field("p_nbq",p_nbq)
      call xios_send_field("p_nbq_max",p_nbq_max)
# endif       
      endif
#  endif

      if (xios_field_is_active("w")) then
# ifndef M3FAST
      call Wvlcty (tile, workr) 
      call xios_send_field("w",workr)
# else        
      worknbq=0.
      do k=0,N
      do j=0,Mm
         do i=0,Lm
            worknbq(i,j,k)=
!          rw_int_nbq(i,j,k)
     &         wz     (i,j,k,nstp)
#  ifdef NHINT_WH
     &        +wzh_nbq(i,j,k)
#  endif
         enddo
      enddo
      enddo
      call xios_send_field("w",worknbq(:,:,0:N))
# endif
      endif

#  ifdef DIAGNOSTICS_EDDY

         GRho=-1*g/rho0

          do k=1,N
            do j=JstrR,JendR
              do i=IstrR,IendR

                  GRho=-1*g/rho0

                  workru(i,j,k) = 0.5*(u(i+1,j,k,nstp)+u(i,j,k,nstp))
                  workrv(i,j,k) = 0.5*(v(i,j+1,k,nstp)+v(i,j,k,nstp))
                if (k.eq.N) then
                  workrsu(i,j) = 0.5*(sustr(i+1,j)+sustr(i,j))
                  workrsv(i,j) = 0.5*(svstr(i,j+1)+svstr(i,j))
                  workrbu(i,j) = 0.5*(bustr(i+1,j)+bustr(i,j))
                  workrbv(i,j) = 0.5*(bvstr(i,j+1)+bvstr(i,j))

                  workrug(i,j) = 0.5*(
     &                    -g/f(i+1,j)*(vmask(i+1,j+1)
     &              *(zeta(i+1,j+1,nstp)-zeta(i+1,j,nstp))
     &                   *0.5*(pm(i+1,j+1)+pm(i+1,j))
     &                    ) + 
     &                    -g/f(i,j)*(vmask(i,j+1)
     &              *(zeta(i,j+1,nstp)-zeta(i,j,nstp))
     &                   *0.5*(pm(i,j+1)+pm(i,j))
     &                    ) 
     &                    )
                  workrvg(i,j) = 0.5*(
     &                     g/f(i,j+1)*(umask(i+1,j+1)
     &              *(zeta(i+1,j+1,nstp)-zeta(i,j+1,nstp))
     &                   *0.5*(pn(i+1,j+1)+pn(i,j+1)) ) + 
     &                     g/f(i,j)*(umask(i+1,j)
     &              *(zeta(i+1,j,nstp)-zeta(i,j,nstp))
     &                   *0.5*(pn(i+1,j)+pn(i,j)) )
     &                    )

                endif

                enddo
              enddo
           enddo

      if (xios_field_is_active("uu")) then
          call xios_send_field("uu", workru**2)
      endif
      if (xios_field_is_active("uv")) then
          call xios_send_field("uv", workru*workrv)
      endif
      if (xios_field_is_active("vv")) then
          call xios_send_field("vv", workrv**2)
      endif
      if (xios_field_is_active("ub")) then
          call xios_send_field("ub", workru*GRho*rho)
      endif
      if (xios_field_is_active("vb")) then
          call xios_send_field("vb", workrv*GRho*rho)
      endif
      if (xios_field_is_active("wb")) then
          call xios_send_field("wb", workr *GRho*rho)
      endif
      if (xios_field_is_active("uw")) then
          call xios_send_field("uw", workru*workr)
      endif
      if (xios_field_is_active("vw")) then
          call xios_send_field("vw", workrv*workr)
      endif
      if (xios_field_is_active("usustr")) then
          call xios_send_field("usustr", workru(:,:,N)*workrsu)
      endif
      if (xios_field_is_active("vsvstr")) then
          call xios_send_field("vsvstr", workrv(:,:,N)*workrsv)
      endif
      if (xios_field_is_active("ugsustr")) then
          call xios_send_field("ugsustr", workrug*workrsu)
      endif
      if (xios_field_is_active("vgsvstr")) then
          call xios_send_field("vgsvstr", workrvg*workrsv)
      endif
      if (xios_field_is_active("ubustr")) then
          call xios_send_field("ubustr", workru(:,:,1)*workrbu)
      endif
      if (xios_field_is_active("vbvstr")) then
          call xios_send_field("vbvstr", workrv(:,:,1)*workrbv)
      endif
      if (xios_field_is_active("uT")) then
          call xios_send_field("uT", workru*t(:,:,:,nstp,itemp))
      endif
      if (xios_field_is_active("vT")) then
          call xios_send_field("vT", workrv*t(:,:,:,nstp,itemp))
      endif
      if (xios_field_is_active("wT")) then
          call xios_send_field("wT", workr *t(:,:,:,nstp,itemp))
      endif
      if (xios_field_is_active("uS")) then
          call xios_send_field("uS", workru*t(:,:,:,nstp,isalt))
      endif
      if (xios_field_is_active("vS")) then
          call xios_send_field("vS", workrv*t(:,:,:,nstp,isalt))
      endif
      if (xios_field_is_active("wS")) then
          call xios_send_field("wS", workr *t(:,:,:,nstp,isalt))
      endif

#  endif


        
#  ifdef VIS_COEF_3D
      if (xios_field_is_active("visc3d")) then
      call xios_send_field("visc3d",visc3d_r)
      endif
#  endif

#  ifdef DIF_COEF_3D
      if (xios_field_is_active("diff3d")) then
      do k=1,N
         do j=1,Mm
            do i=1,Lm
               workr(i,j,k)=
#   ifdef TS_DIF2
     &              diff2(i,j,itemp)
#    ifdef TS_DIF_SMAGO
     &              +diff3d_r(i,j,k)
#    endif
#   elif defined TS_DIF4
     &              diff4(i,j,itemp)
#    ifdef TS_DIF_SMAGO
     &              +diff3d_r(i,j,k)*om_r(i,j)*on_r(i,j)
#    endif
     &              +0.25*(diff3d_u(i,j,k)+diff3d_u(i+1,j,k)
     &              +diff3d_v(i,j,k)+diff3d_v(i,j+1,k))
#   endif
            enddo
         enddo
      enddo
      call xios_send_field("diff3d",workr)
      endif
#  endif

      if (xios_field_is_active("AKv")) then
      call xios_send_field("AKv",Akv)
      endif
#  ifdef TEMPERATURE
      if (xios_field_is_active("AKt")) then
      work=Akt(:,:,:,itemp)
      call xios_send_field("AKt",work)
      endif
#  endif
#  ifdef SALINITY      
      if (xios_field_is_active("AKs")) then
      work=Akt(:,:,:,isalt)
      call xios_send_field("AKs",work)
      endif
#  endif
#  ifdef GLS_MIXING
      if (xios_field_is_active("AKk")) then
!      call xios_send_field("AKk",Akk)
      endif
      if (xios_field_is_active("AKp")) then
!      call xios_send_field("AKp",Akp)
      endif
#  endif
#  ifdef LMD_SKPP
      if (xios_field_is_active("hbl")) then
#   ifdef LMD_SKPP2005
      work2d=hbls(:,:,nstp)
#   else
      work2d=hbl
#   endif
      call xios_send_field("hbl",work2d)
      endif
#  endif
#  ifdef LMD_BKPP
      if (xios_field_is_active("hbbl")) then
      call xios_send_field("hbbl",work2d)
      endif
#  endif
#  ifdef GLS_MIXING
      if (xios_field_is_active("tke")) then
      call xios_send_field("tke",trb(:,:,:,nstp,itke))
      endif
      if (xios_field_is_active("gls")) then
      call xios_send_field("gls",trb(:,:,:,nstp,igls))
      endif
      if (xios_field_is_active("Lscale")) then
      call xios_send_field("Lscale",Lscale)
      endif
#  endif

#  ifdef TEMPERATURE
      if (xios_field_is_active("shflx")) then
      work2d=stflx(:,:,itemp)*rho0*Cp
     &       SWITCH rmask
      call xios_send_field("shflx",work2d)
      endif
#  endif      
#  ifdef SALINITY    
      if (xios_field_is_active("swflx")) then
         do j=0,Mm+1
            do i=0,Lm+1
               work2d(i,j)=stf_cff*stflx(i,j,isalt)/
     &              ( max(eps,t(i,j,N,nstp,isalt))) 
     &              SWITCH rmask(i,j)
            enddo
        enddo
      call xios_send_field("swflx",work2d)
      endif
#  endif

#  ifdef BULK_FLUX
      if (xios_field_is_active("radsw")) then
      work2d=shflx_rsw*rho0*Cp
      call xios_send_field("radsw",work2d)
      endif
#  else 
      if (xios_field_is_active("swrad")) then
      work2d=srflx*rho0*Cp
      call xios_send_field("swrad",work2d)
      endif
#  endif
      
#  ifdef BULK_FLUX
      if (xios_field_is_active("shflx_rlw")) then
      work2d=shflx_rlw*rho0*Cp
      call xios_send_field("shflx_rlw",work2d)
      endif

      if (xios_field_is_active("shflx_lat")) then
      work2d=shflx_lat*rho0*Cp
      call xios_send_field("shflx_lat",work2d)
      endif

      if (xios_field_is_active("shflx_sen")) then
      work2d=shflx_sen*rho0*Cp
      call xios_send_field("shflx_sen",work2d)
      endif
#  endif
#  ifdef SST_SKIN
      if (xios_field_is_active("sst_skin")) then
      call xios_send_field("sst_skin",sst_skin)
      endif
#  endif
     

#  ifdef DIAGNOSTICS_UV

      itrc = 1
      if (xios_field_is_active("u_rate")) then
          call xios_send_field("u_rate", Mrate(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_adv")) then
          call xios_send_field("u_adv", MXadv(:,:,:,itrc)
     &           + MYadv(:,:,:,itrc) + MVadv(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_Cor")) then
          call xios_send_field("u_Cor",MCor(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_Prsgrd")) then
          call xios_send_field("u_Prsgrd",MPrsgrd(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_Hmix")) then
          call xios_send_field("u_Hmix",MHmix(:,:,:,itrc,nstp))
      endif
      if (xios_field_is_active("u_Hdiff")) then
          call xios_send_field("u_Hdiff",MHdiff(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_Vmix")) then
          call xios_send_field("u_Vmix",MVmix(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_Vmix2")) then
          call xios_send_field("u_Vmix2",MVmix2(:,:,:,itrc))
      endif
#  ifdef DIAGNOSTICS_BARO
      if (xios_field_is_active("u_Baro")) then
          call xios_send_field("u_Baro",MBaro(:,:,:,itrc))
      endif
#  endif
#  ifdef M3FAST
      if (xios_field_is_active("u_fast")) then
          call xios_send_field("u_fast",Mfast(:,:,:,itrc))
      endif
#  endif
#  ifdef MRL_WCI
      if (xios_field_is_active("u_vf")) then
          call xios_send_field("u_vf",Mvf(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_brk")) then
          call xios_send_field("u_brk",Mbrk(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_StCo")) then
          call xios_send_field("u_StCo",MStCo(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_Vvf")) then
          call xios_send_field("u_Vvf",MVvf(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_Prscrt")) then
          call xios_send_field("u_Prscrt",MPrscrt(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_sbk")) then
          call xios_send_field("u_sbk",Msbk(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_bwf")) then
          call xios_send_field("u_bwf",Mbwf(:,:,:,itrc))
      endif
      if (xios_field_is_active("u_frc")) then
          call xios_send_field("u_frc",Mfrc(:,:,:,itrc))
      endif
#  endif

      itrc = 2
      if (xios_field_is_active("v_rate")) then
          call xios_send_field("v_rate", Mrate(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_adv")) then
          call xios_send_field("v_adv", MXadv(:,:,:,itrc)
     &            + MYadv(:,:,:,itrc) + MVadv(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_Cor")) then
          call xios_send_field("v_Cor",MCor(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_Prsgrd")) then
          call xios_send_field("v_Prsgrd",MPrsgrd(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_Hmix")) then
          call xios_send_field("v_Hmix",MHmix(:,:,:,itrc,nstp))
      endif
      if (xios_field_is_active("v_Hdiff")) then
          call xios_send_field("v_Hdiff",MHdiff(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_Vmix")) then
          call xios_send_field("v_Vmix",MVmix(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_Vmix2")) then
          call xios_send_field("v_Vmix2",MVmix2(:,:,:,itrc))
      endif
#  ifdef DIAGNOSTICS_BARO
      if (xios_field_is_active("v_Baro")) then
          call xios_send_field("v_Baro",MBaro(:,:,:,itrc))
      endif
#  endif
#  ifdef M3FAST
      if (xios_field_is_active("v_fast")) then
          call xios_send_field("v_fast",Mfast(:,:,:,itrc))
      endif
#  endif
#  ifdef MRL_WCI
      if (xios_field_is_active("v_vf")) then
          call xios_send_field("v_vf",Mvf(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_brk")) then
          call xios_send_field("v_brk",Mbrk(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_StCo")) then
          call xios_send_field("v_StCo",MStCo(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_Vvf")) then
          call xios_send_field("v_Vvf",MVvf(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_Prscrt")) then
          call xios_send_field("v_Prscrt",MPrscrt(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_sbk")) then
          call xios_send_field("v_sbk",Msbk(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_bwf")) then
          call xios_send_field("v_bwf",Mbwf(:,:,:,itrc))
      endif
      if (xios_field_is_active("v_frc")) then
          call xios_send_field("v_frc",Mfrc(:,:,:,itrc))
      endif
#  endif

#  endif



#  ifdef DIAGNOSTICS_TS
         itrc=1
#  ifdef TEMPERATURE
      if (xios_field_is_active("T_rate")) then
          call xios_send_field("T_rate", Trate(:,:,:,itrc))
      endif
      if (xios_field_is_active("T_adv")) then
          call xios_send_field("T_adv", TXadv(:,:,:,itrc)
     &           + TYadv(:,:,:,itrc) + TVadv(:,:,:,itrc))
      endif
      if (xios_field_is_active("T_Hmix")) then
          call xios_send_field("T_Hmix",THmix(:,:,:,itrc))
      endif
      if (xios_field_is_active("T_Vmix")) then
          call xios_send_field("T_Vmix",TVmix(:,:,:,itrc))
      endif
      if (xios_field_is_active("T_Forc")) then
          call xios_send_field("T_Forc",TForc(:,:,:,itrc))
      endif
#  endif          
#  ifdef SALINITY
        itrc = 2
      if (xios_field_is_active("S_rate")) then
          call xios_send_field("S_rate", Trate(:,:,:,itrc))
      endif
      if (xios_field_is_active("S_adv")) then
          call xios_send_field("S_adv", TXadv(:,:,:,itrc)
     &              + TYadv(:,:,:,itrc) + TVadv(:,:,:,itrc))
      endif
      if (xios_field_is_active("S_Hmix")) then
          call xios_send_field("S_Hmix",THmix(:,:,:,itrc))
      endif
      if (xios_field_is_active("S_Vmix")) then
          call xios_send_field("S_Vmix",TVmix(:,:,:,itrc))
      endif
      if (xios_field_is_active("S_Forc")) then
          call xios_send_field("S_Forc",TForc(:,:,:,itrc))
      endif
#  endif
#  endif




#  ifdef DIAGNOSTICS_VRT
      if (xios_field_is_active("vrtrate")) then
          call xios_send_field("vrtrate", vrtrate(:,:))
      endif
      if (xios_field_is_active("vrtadv")) then
          call xios_send_field("vrtadv", vrtXadv(:,:) + vrtYadv(:,:))
      endif
      if (xios_field_is_active("vrtCor")) then
          call xios_send_field("vrtCor",vrtCor(:,:))
      endif
      if (xios_field_is_active("vrtPrsgrd")) then
          call xios_send_field("vrtPrsgrd",vrtPrsgrd(:,:))
      endif
      if (xios_field_is_active("vrtHmix")) then
          call xios_send_field("vrtHmix",vrtHmix(:,:))
      endif
      if (xios_field_is_active("vrtHdiff")) then
          call xios_send_field("vrtHdiff",vrtHdiff(:,:))
      endif
      if (xios_field_is_active("vrtVmix")) then
!          call xios_send_field("vrtVmix",vrtVmix(:,:))
      endif
      if (xios_field_is_active("vrtVmix2")) then
          call xios_send_field("vrtVmix2",vrtVmix2(:,:))
      endif
#  ifdef DIAGNOSTICS_BARO
      if (xios_field_is_active("vrtBaro")) then
          call xios_send_field("vrtBaro",vrtBaro(:,:))
      endif
#  endif
      if (xios_field_is_active("vrtDrag")) then
          call xios_send_field("vrtDrag",vrtDrag(:,:))
      endif
      if (xios_field_is_active("vrtWind")) then
          call xios_send_field("vrtWind",vrtWind(:,:))
      endif
#  ifdef M3FAST
      if (xios_field_is_active("vrtfast")) then
          call xios_send_field("vrtfast",vrtfast(:,:))
      endif
#  endif
#  endif


#  ifdef DIAGNOSTICS_EK

      if (xios_field_is_active("ekrate")) then
          call xios_send_field("ekrate", ekrate(:,:))
      endif
      if (xios_field_is_active("ekadv")) then
          call xios_send_field("ekadv", ekHadv(:,:) + ekVadv(:,:))
      endif
      if (xios_field_is_active("ekCor")) then
          call xios_send_field("ekCor",ekCor(:,:))
      endif
      if (xios_field_is_active("ekPrsgrd")) then
          call xios_send_field("ekPrsgrd",ekPrsgrd(:,:))
      endif
      if (xios_field_is_active("ekHmix")) then
          call xios_send_field("ekHmix",ekHmix(:,:))
      endif
      if (xios_field_is_active("ekHdiff")) then
          call xios_send_field("ekHdiff",ekHdiff(:,:))
      endif
      if (xios_field_is_active("ekVmix")) then
          call xios_send_field("ekVmix",ekVmix(:,:))
      endif
      if (xios_field_is_active("ekVmix2")) then
          call xios_send_field("ekVmix2",ekVmix2(:,:))
      endif
#  ifdef DIAGNOSTICS_BARO
      if (xios_field_is_active("ekBaro")) then
          call xios_send_field("ekBaro",ekBaro(:,:))
      endif
#  endif
      if (xios_field_is_active("ekvol")) then
          call xios_send_field("ekvol",ekvol(:,:))
      endif
      if (xios_field_is_active("ekDrag")) then
          call xios_send_field("ekDrag",ekDrag(:,:))
      endif
      if (xios_field_is_active("ekWind")) then
          call xios_send_field("ekWind",ekWind(:,:))
      endif
#  ifdef M3FAST
      if (xios_field_is_active("ekfast")) then
          call xios_send_field("ekfast",ekfast(:,:))
      endif
#  endif
#  endif


# ifdef DIAGNOSTICS_PV
      if (xios_field_is_active("u_rhs")) then
          call xios_send_field("u_rhs", Mrhs(:,:,:,1))
      endif
      if (xios_field_is_active("v_rhs")) then
          call xios_send_field("v_rhs", Mrhs(:,:,:,2))
      endif
# ifdef TEMPERATURE          
      if (xios_field_is_active("T_rhs")) then
          call xios_send_field("T_rhs", Trhs(:,:,:,1))
      endif
#  endif
#  ifdef SALINITY
      if (xios_field_is_active("S_rhs")) then
          call xios_send_field("S_rhs", Trhs(:,:,:,2))
      endif
#  endif
# ifdef DIAGNOSTICS_PV_FULL
      if (xios_field_is_active("u_vmix_trans")) then
          call xios_send_field("u_vmix_trans", pv(:,:,:))
      endif
      if (xios_field_is_active("v_vmix_trans")) then
          call xios_send_field("v_vmix_trans", pvd(:,:,:))
      endif
#  endif
#  endif

#  if defined BIOLOGY
#   if defined PISCES
      CALL trc_wri_pisces
#   else
      if (xios_field_is_active("hel")) then
      call xios_send_field("hel",hel)     
      endif
#   if (defined BIO_NChlPZD ||  defined BIO_N2ChlPZD2)
      if (xios_field_is_active("theta")) then
      call xios_send_field("theta",theta)
      endif
#    ifdef OXYGEN
      if (xios_field_is_active("U10")) then
      call xios_send_field("U10",u10)
      endif
      if (xios_field_is_active("KvO2")) then
      call xios_send_field("KvO2",Kv_O2)
      endif
      if (xios_field_is_active("O2sat")) then
      call xios_send_field("O2sat",O2satu)
      endif
#    endif /* OXYGEN */
#   elif defined BIO_BioEBUS 
      if (xios_field_is_active("AOU")) then
      call xios_send_field("AOU",AOU)
      endif
      if (xios_field_is_active("wind10")) then
      call xios_send_field("wind10",wind10)
      endif
#   endif
#   endif
#  endif /* BIOLOGY */


#  ifdef SEDIMENT
      do i=1,NGRAV
         write(nametrc,'(A,I1)') 'gravel_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    t(:,:,:,nstp,1+ntrc_salt+ntrc_pas+ntrc_bio+i))
      endif
      enddo      
      do i=1,NSAND
         write(nametrc,'(A,I1)') 'sand_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    t(:,:,:,nstp,1+ntrc_salt+ntrc_pas+ntrc_bio+ngrav+i)) 
      endif
      enddo      
      do i=1,NMUD
         write(nametrc,'(A,I1)') 'mud_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    t(:,:,:,nstp,1+ntrc_salt+ntrc_pas+ntrc_bio+ngrav+nsand+i)) 
      endif
      enddo      

      if (xios_field_is_active("bed_thick")) then
      call xios_send_field("bed_thick",bed_thick) 
      endif
      if (xios_field_is_active("bed_poros")) then
      call xios_send_field("bed_poros",bed_poros) 
      endif

      do i=1,NGRAV
         write(nametrc,'(A,I1)') 'bed_frac_grav_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bed_frac(:,:,:,i))
      endif
         write(nametrc,'(A,I1)') 'bed_mass_grav_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bed_mass(:,:,:,nstp,i))
      endif
      enddo      
      do i=1,NSAND
         write(nametrc,'(A,I1)') 'bed_frac_sand_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bed_frac(:,:,:,NGRAV+i))
      endif
         write(nametrc,'(A,I1)') 'bed_mass_sand_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bed_mass(:,:,:,nstp,NGRAV+i))
      endif
      enddo      
      do i=1,NMUD
         write(nametrc,'(A,I1)') 'bed_frac_mud_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bed_frac(:,:,:,NGRAV+NSAND+i))
      endif
         write(nametrc,'(A,I1)') 'bed_mass_mud_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bed_mass(:,:,:,nstp,NGRAV+NSAND+i))
      endif
      enddo      

#   ifdef SUSPLOAD
      do i=1,NGRAV
         write(nametrc,'(A,I1)') 'dflux_grav_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    settling_flux(:,:,i))
      endif
         write(nametrc,'(A,I1)') 'eflux_grav_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    ero_flux(:,:,i))
      endif
      enddo      
      do i=1,NSAND
         write(nametrc,'(A,I1)') 'dflux_sand_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    settling_flux(:,:,NGRAV+i))
      endif
         write(nametrc,'(A,I1)') 'eflux_sand_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    ero_flux(:,:,NGRAV+i))
      endif
      enddo      
      do i=1,NMUD
         write(nametrc,'(A,I1)') 'dflux_mud_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
      endif
     &    settling_flux(:,:,NGRAV+NSAND+i))
         write(nametrc,'(A,I1)') 'edflux_mud_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    ero_flux(:,:,NGRAV+NSAND+i))
      endif
      enddo      
#   endif

#   ifdef BEDLOAD
      do i=1,NGRAV
         write(nametrc,'(A,I1)') 'bdlu_gravel_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bedldu(:,:,i))
      endif
         write(nametrc,'(A,I1)') 'bdlv_gravel_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bedldv(:,:,i))
      endif
      enddo      
      do i=1,NSAND
         write(nametrc,'(A,I1)') 'bdlu_sand_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bedldu(:,:,NGRAV+i))
      endif
         write(nametrc,'(A,I1)') 'bdlv_sand_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bedldv(:,:,NGRAV+i))
      endif
      enddo      
      do i=1,NMUD
         write(nametrc,'(A,I1)') 'bdlu_mud_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bedldu(:,:,NGRAV+NSAND+i))
      endif
         write(nametrc,'(A,I1)') 'bdlv_mud_',i
      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),
     &    bedldv(:,:,NGRAV+NSAND+i))
      endif
      enddo      
#   endif
#  endif /* SEDIMENT */



# ifdef MUSTANG
      do itrc=3,NT
         write(nametrc,'(a)') TRIM(name_var(itrc-2))
#if defined key_sand2D
           if (itrc>=itsubs1+nv_grav .AND.
     &         itrc<=itsubs1+nv_grav+nv_sand-1) THEN   ! isand1,isand2
           if (l_subs2D(itrc-(itsubs1+nv_grav)+1)
     &         .AND.l_outsandrouse(itrc-(itsubs1+nv_grav)+1)) then
            workr=t(:,:,:,nstp,itrc)
            do j=Jstr,Jend
              do i=Istr,Iend
                DO k=1,N
           IF( sum_tmp(itrc-itsubs1-nv_grav+1,i,j) > 0.0_rsh) THEN
           workr(i,j,k) =t(i,j,1,nstp,itrc)*(z_w(i,j,1)-z_w(i,j,0))*
     &       (  (z_w(i,j,N)-z_r(i,j,k))
     &          /(z_r(i,j,k)-z_w(i,j,0))
     &       ) **rouse2D(itrc-itsubs1-nv_grav+1,i,j)
           workr(i,j,k)=workr(i,j,k)/sum_tmp(itrc-itsubs1-nv_grav+1,i,j)
           ELSE
           workr(i,j,k)=0.
          ENDIF
          ENDDO
          enddo
          enddo
           else
            workr=t(:,:,:,nstp,itrc)
           endif
           else
            workr=t(:,:,:,nstp,itrc)
           endif
#else
            workr=t(:,:,:,nstp,itrc)
#endif

      if (xios_field_is_active(TRIM(nametrc))) then
         call xios_send_field(TRIM(nametrc),workr)
      endif
      enddo

      CALL sed_MUSTANG_outres(Istr,Iend,Jstr,Jend,ntrc_subs,h)
      if (xios_field_is_active("tauskin")) then
        call xios_send_field("tauskin",tauskin)
      endif
      if (xios_field_is_active("tauskin_c")) then
        call xios_send_field("tauskin_c",tauskin_c)
      endif
      if (xios_field_is_active("tauskin_w")) then
        call xios_send_field("tauskin_w",tauskin_w)
      endif

      if (xios_field_is_active("ksma")) then
        call xios_send_field("ksma",REAL(ksma))
      endif
      IF (l_outsed_hsed) then
            if (xios_field_is_active("hsed")) then
            call xios_send_field("hsed",var2D_hsed)
            endif
      ENDIF

      IF (l_outsed_cv_sed) THEN
        do itrc=1,ntrc_subs
         do k=1,nk_nivsed_out
           workrsed(:,:,k)=var3D_cvsed(k,:,:,itrc)
         enddo
            write (nametrc, '(a,"_sed")') TRIM(name_var(itrc))
            if (xios_field_is_active(TRIM(nametrc))) then
            call xios_send_field(TRIM(nametrc),workrsed)
            endif
        enddo
      ENDIF

      do itrc=1,ntrc_subs
            IF (l_outsed_ws .and. ((itrc .GE. imud1) .AND. (itrc .LE. imud2))) THEN
                  write(nametrc,'(a,"_ws")') TRIM(name_var(itrc))
                  workr=ws_part(:,:,:,itemp+ntrc_salt+itrc)
                  if (xios_field_is_active(TRIM(nametrc))) then
                        call xios_send_field(TRIM(nametrc),workr)
                  endif
            endif
      enddo
      
      IF (l_outsed_dzs) then
            if (xios_field_is_active('dzs')) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3D_dzs(k,:,:)
            enddo
            call xios_send_field('dzs',workrsed)
            endif
      ENDIF

      IF (l_outsed_poro) then
            if (xios_field_is_active('poro')) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3D_poro(k,:,:)
            enddo
            call xios_send_field('poro',workrsed)
            endif
      ENDIF

      IF (l_outsed_temp_sed) then
            if (xios_field_is_active('temp_sed')) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3D_TEMP(k,:,:)
            enddo
            call xios_send_field('temp_sed',workrsed)
            endif
      ENDIF

      IF (l_outsed_salt_sed) then
            if (xios_field_is_active('salt_sed')) then
            do k=1,nk_nivsed_out
            workrsed(:,:,k)=var3D_SAL(k,:,:)
            enddo
            call xios_send_field('salt_sed',workrsed)
            endif
      ENDIF

      DO itrc = 1,nvpc
         IF (l_outsed_toce) then
            write (nametrc, '(a,"_toce")') TRIM(name_var(itrc))
            if (xios_field_is_active(TRIM(nametrc))) then
              call xios_send_field(TRIM(nametrc),var2D_toce(itrc,:,:))
            endif
         ENDIF
         IF (l_outsed_flx_s2w_w2s) then
            write (nametrc, '(a,"_flx_s2w")') TRIM(name_var(itrc))
            if (xios_field_is_active(TRIM(nametrc))) then
              call xios_send_field(TRIM(nametrc),var2D_flx_s2w(itrc,:,:))
            endif
            write (nametrc, '(a,"_flx_w2s")') TRIM(name_var(itrc))
            if (xios_field_is_active(TRIM(nametrc))) then
              call xios_send_field(TRIM(nametrc),var2D_flx_w2s(itrc,:,:))
            endif
         ENDIF
#ifdef  key_MUSTANG_V2
         IF (l_outsed_pephm_fcor) then
            write (nametrc, '(a,"_pephm_fcor")') TRIM(name_var(itrc))
            if (xios_field_is_active(TRIM(nametrc))) then
              call xios_send_field(TRIM(nametrc),
     &                             var2D_pephm_fcor(itrc,:,:))
            endif
         ENDIF
#ifdef key_MUSTANG_bedload
         IF (l_outsed_bedload) then
            write (nametrc, '(a,"_flx_bx")') TRIM(name_var(itrc))
            if (xios_field_is_active(TRIM(nametrc))) then
              call xios_send_field(TRIM(nametrc),
     &                             var2D_flx_bx(itrc,:,:))
            endif
            write (nametrc, '(a,"_flx_by")') TRIM(name_var(itrc))
            if (xios_field_is_active(TRIM(nametrc))) then
              call xios_send_field(TRIM(nametrc),
     &                             var2D_flx_by(itrc,:,:))
            endif
            write (nametrc, '(a,"_bil_bedload")') TRIM(name_var(itrc))
            if (xios_field_is_active(TRIM(nametrc))) then
              call xios_send_field(TRIM(nametrc),
     &                             var2D_bil_bedload(itrc,:,:))
            endif
         ENDIF
         IF (l_outsed_fsusp) then
            write (nametrc, '(a,"_fsusp")') TRIM(name_var(itrc))
            if (xios_field_is_active(TRIM(nametrc))) then
              call xios_send_field(TRIM(nametrc),
     &                             var2D_fsusp(itrc,:,:))
            endif
         ENDIF
#endif /* key_MUSTANG_bedload */
#endif /* key_MUSTANG_V2 */
      ENDDO

      IF (l_outsed_frmudsup) then
            if (xios_field_is_active("frmudsup")) then
            call xios_send_field("frmudsup",var2D_frmudsup)
            endif
      ENDIF
      IF (l_outsed_dzs_ksmax) then
            if (xios_field_is_active("dzs_ksmax")) then
            call xios_send_field("dzs_ksmax",var2D_dzs_ksmax)
            endif
      ENDIF
      if (xios_field_is_active("z0sed")) then
            call xios_send_field("z0sed",z0sed)
      endif
      if (xios_field_is_active("z0hydro")) then
            call xios_send_field("z0hydro",zob)
      endif
      IF (l_outsed_flx_s2w_w2s) then
            if (xios_field_is_active("flx_s2w_noncoh")) then
            call xios_send_field("flx_s2w_noncoh",var2D_flx_s2w_noncoh)
            endif
            if (xios_field_is_active("flx_s2w_coh")) then
            call xios_send_field("flx_s2w_coh",var2D_flx_s2w_coh)
            endif
            if (xios_field_is_active("flx_w2s_noncoh")) then
            call xios_send_field("flx_w2s_noncoh",var2D_flx_w2s_noncoh)
            endif
            if (xios_field_is_active("flx_w2s_coh")) then
            call xios_send_field("flx_w2s_coh",var2D_flx_w2s_coh)
            endif
      ENDIF
#ifdef  key_MUSTANG_V2
      IF (l_outsed_theoric_active_layer) then
            if (xios_field_is_active("theoric_active_layer")) then
            call xios_send_field("theoric_active_layer",
     &                        var2D_theoric_active_layer)
            endif
      ENDIF
      IF (l_outsed_ero_details) then
            if (xios_field_is_active("tero_noncoh")) then
            call xios_send_field("tero_noncoh",
     &                        var2D_tero_noncoh)
            endif
            if (xios_field_is_active("tero_coh")) then
            call xios_send_field("tero_coh",
     &                        var2D_tero_coh)
            endif
            if (xios_field_is_active("pct_iter_noncoh")) then
            call xios_send_field("pct_iter_noncoh",
     &                        var2D_pct_iter_noncoh)
            endif
            if (xios_field_is_active("pct_iter_coh")) then
                  call xios_send_field("pct_iter_coh",
     &                        var2D_pct_iter_coh)
            endif
            if (xios_field_is_active("niter_ero")) then
                  call xios_send_field("niter_ero",
     &                        var2D_niter_ero)
            endif
      ENDIF
#ifdef key_MUSTANG_bedload
      IF (l_outsed_bedload) then
            if (xios_field_is_active("flx_bx_int")) then
            call xios_send_field("flx_bx_int",
     &                        var2D_flx_bx_int)
            endif
            if (xios_field_is_active("flx_by_int")) then
                  call xios_send_field("flx_by_int",
     &                        var2D_flx_by_int)
            endif
            if (xios_field_is_active("bil_bedload_int")) then
            call xios_send_field("bil_bedload_int",
     &                        var2D_bil_bedload_int)
            endif
      ENDIF
#endif /* key_MUSTANG_bedload */
#endif /* key_MUSTANG_V2 */

      IF (l_dyn_insed) THEN
            IF (l_outsed_consolidation) then
                  if (xios_field_is_active('dtsdzs')) then
                        do k=1,nk_nivsed_out
                              workrsed(:,:,k)=var3Dksed_dtsdzs(k,:,:)
                        enddo
                        call xios_send_field('dtsdzs',workrsed)
                  endif
                  if (xios_field_is_active('loadograv')) then
                        do k=1,nk_nivsed_out
                              workrsed(:,:,k)=var3Dksed_loadograv(k,:,:)
                        enddo
                        call xios_send_field('loadograv',workrsed)
                  endif
                  if (xios_field_is_active('permeab')) then
                        do k=1,nk_nivsed_out
                              workrsed(:,:,k)=var3Dksed_permeab(k,:,:)
                        enddo
                        call xios_send_field('permeab',workrsed)
                  endif
                  if (xios_field_is_active('sed_rate')) then
                        do k=1,nk_nivsed_out
                              workrsed(:,:,k)=var3Dksed_sed_rate(k,:,:)
                        enddo
                        call xios_send_field('sed_rate',workrsed)
                  endif
                  if (xios_field_is_active('sigmadjge')) then
                        do k=1,nk_nivsed_out
                              workrsed(:,:,k)=var3Dksed_sigmadjge(k,:,:)
                        enddo
                        call xios_send_field('sigmadjge',workrsed)
                  endif
                  if (xios_field_is_active('sigmapsg')) then
                        do k=1,nk_nivsed_out
                              workrsed(:,:,k)=var3Dksed_sigmapsg(k,:,:)
                        enddo
                        call xios_send_field('sigmapsg',workrsed)
                  endif
                  if (xios_field_is_active('hinder')) then
                        do k=1,nk_nivsed_out
                              workrsed(:,:,k)=var3Dksed_hinder(k,:,:)
                        enddo
                        call xios_send_field('hinder',workrsed)
                  endif
                  if (xios_field_is_active('stateconsol')) then
                        do k=1,nk_nivsed_out
                        workrsed(:,:,k)=var3Dksed_stateconsol(k,:,:)
                        enddo
                        call xios_send_field('stateconsol',workrsed)
                  endif
            ENDIF
      ENDIF
# endif /* MUSTANG */     

#  ifdef BBL
      if (xios_field_is_active("Abed")) then
         call xios_send_field("Abed",Abed) 
      endif
      if (xios_field_is_active("Hripple")) then
         call xios_send_field("Hripple",Hripple)
      endif
      if (xios_field_is_active("Lripple")) then
         call xios_send_field("Lripple",Lripple) 
      endif
      if (xios_field_is_active("Zbnot")) then
         call xios_send_field("Zbnot",Zbnot)
      endif
      if (xios_field_is_active("Zbapp")) then
         call xios_send_field("Zbapp",Zbapp)
      endif

      if (xios_field_is_active("bostrw")) then
         work2d=sqrt(bustrw**2
     &                             +bvstrw**2)
     &                                   *rho0
         call xios_send_field("bostrw",work2d)
      endif
      if (xios_field_is_active("bustrc")) then
!         call xios_send_field("bustrc",bustrc*rho0)
      endif
      if (xios_field_is_active("bvstrc")) then
!         call xios_send_field("bvstrc",bvstrc*rho0)
      endif
      if (xios_field_is_active("bustrw")) then
         call xios_send_field("bustrw",bustrw*rho0)
      endif
      if (xios_field_is_active("bvstrw")) then
         call xios_send_field("bvstrw",bvstrw*rho0)
      endif
      if (xios_field_is_active("bustrcwmax")) then
!         call xios_send_field("bustrcwmax",bustrcwmax*rho0)
      endif
      if (xios_field_is_active("bvstrcwmax")) then
!         call xios_send_field("bvstrcwmax",bvstrcwmax*rho0)
      endif
!     
#  endif /* BBL */


#  ifdef OBSTRUCTION
      do iv=1,obst_nbvar
         call xios_send_field(
     &      TRIM(obst_nout_pos)//'_'//TRIM(obst_varname(iv)),
     &      obst_position(iv,:,:))

         call xios_send_field(
     &      TRIM(obst_nout_height_f)//'_'//TRIM(obst_varname(iv)),
     &      obst_height_inst(iv,:,:))

         call xios_send_field(
     &      TRIM(obst_nout_height_e)//'_'//TRIM(obst_varname(iv)),
     &      obst_height(iv,:,:))

         call xios_send_field(
     &      TRIM(obst_nout_dens_f)//'_'//TRIM(obst_varname(iv)),
     &      obst_dens_inst(iv,:,:))
         
         ! obst_dens3d is not allocated if not wanted in output
         IF (l_obstout_dens_e) then
            do k=1,N
                  workr(:,:,k)=obst_dens3d(iv,k,:,:)
            enddo
            call xios_send_field(
     &         TRIM(obst_nout_dens_e)//'_'//TRIM(obst_varname(iv)),
     &         workr)
         END IF

         call xios_send_field(
     &      TRIM(obst_nout_width_f)//'_'//TRIM(obst_varname(iv)),
     &      obst_width_inst(iv,:,:))

         ! obst_width3d is not allocated if not wanted in output
         IF (l_obstout_width_e) then
            do k=1,N
                  workr(:,:,k)=obst_width3d(iv,k,:,:)
            enddo
            call xios_send_field(
     &         TRIM(obst_nout_width_e)//'_'//TRIM(obst_varname(iv)),
     &         workr)
         END IF

         call xios_send_field(
     &      TRIM(obst_nout_thick_f)//'_'//TRIM(obst_varname(iv)),
     &      obst_thick_inst(iv,:,:))

         ! obst_thick3d is not allocated if not wanted in output
         IF (l_obstout_thick_e) then
            do k=1,N
                  workr(:,:,k)=obst_thick3d(iv,k,:,:)
            enddo
            call xios_send_field(
     &         TRIM(obst_nout_thick_e)//'_'//TRIM(obst_varname(iv)),
     &         workr)
         END IF

         ! obst_theta3d is not allocated if not wanted in output
         IF (l_obstout_theta) then
            do k=1,N
                  workr(:,:,k)=obst_theta3d(iv,k,:,:)*180.0/3.1415
            enddo
            call xios_send_field(
     &         TRIM(obst_nout_theta)//'_'//TRIM(obst_varname(iv)),
     &         workr)
         END IF

         ! obst_fracxy is not allocated if not wanted in output
         IF (l_obstout_frac_xy) then
            call xios_send_field(
     &         TRIM(obst_nout_frac_xy)//'_'//TRIM(obst_varname(iv)),
     &         obst_fracxy(iv,:,:))
         END IF

         ! obst_fracz3d is not allocated if not wanted in output
         IF (l_obstout_frac_z) then
            do k=1,N
                  workr(:,:,k)=obst_fracz3d(iv,k,:,:)
            enddo
            call xios_send_field(
     &         TRIM(obst_nout_frac_z)//'_'//TRIM(obst_varname(iv)),
     &         workr)
         END IF

         ! obst_drag3d is not allocated if not wanted in output
         IF (l_obstout_drag) then
            do k=1,N
                  workr(:,:,k)=obst_drag3d(iv,k,:,:)
            enddo
            call xios_send_field(
     &         TRIM(obst_nout_drag)//'_'//TRIM(obst_varname(iv)),
     &         workr)
         END IF

         ! obst_a2d is not allocated if not wanted in output
         IF (l_obstout_a2d) then
            call xios_send_field(
     &      TRIM(obst_nout_a2d)//'_'//TRIM(obst_varname(iv)),
     &      obst_a2d(iv,:,:))
         END IF

         do k=1,N
            workr(:,:,k)=obst_a3d(iv,k,:,:)
         enddo
         call xios_send_field(
     &      TRIM(obst_nout_a3d)//'_'//TRIM(obst_varname(iv)),
     &      workr)

         ! obst_s2d is not allocated if not wanted in output
         IF (l_obstout_s2d) then
            call xios_send_field(
     &      TRIM(obst_nout_s2d)//'_'//TRIM(obst_varname(iv)),
     &      obst_s2d(iv,:,:))
         END IF

         ! obst_s3d is not allocated if not wanted in output
         IF (l_obstout_s3d) then
            do k=1,N
                  workr(:,:,k)=obst_s3d(iv,k,:,:)
            enddo
            call xios_send_field(
     &      TRIM(obst_nout_s3d)//'_'//TRIM(obst_varname(iv)),
     &      workr)
         END IF

      enddo

      ! obst_a2d is not allocated if not wanted in output
      IF (l_obstout_a2d) then
         call xios_send_field(
     &         TRIM(obst_nout_a2d)//'_NoTurb',
     &         obst_a2d(obst_nbvar + 1,:,:))

         call xios_send_field(
     &         TRIM(obst_nout_a2d)//'_Turb',
     &         obst_a2d(obst_nbvar + 2,:,:))
         call xios_send_field(
     &         TRIM(obst_nout_a2d)//'_All',
     &         obst_a2d(obst_nbvar + 3,:,:))
      END IF   

         do k=1,N
               workr(:,:,k)=obst_a3d(obst_nbvar + 1,k,:,:)
         enddo
         call xios_send_field(
     &         TRIM(obst_nout_a3d)//'_NoTurb',
     &         workr)
         do k=1,N
               workr(:,:,k)=obst_a3d(obst_nbvar + 2,k,:,:)
         enddo
         call xios_send_field(
     &         TRIM(obst_nout_a3d)//'_Turb',
     &         workr)
         do k=1,N
               workr(:,:,k)=obst_a3d(obst_nbvar + 3,k,:,:)
         enddo
         call xios_send_field(
     &         TRIM(obst_nout_a3d)//'_All',
     &         workr)

      ! obst_s2d is not allocated if not wanted in output
      IF (l_obstout_s2d) then
         call xios_send_field(
     &         TRIM(obst_nout_s2d)//'_NoTurb',
     &         obst_s2d(obst_nbvar + 1,:,:))
         call xios_send_field(
     &         TRIM(obst_nout_s2d)//'_Turb',
     &         obst_s2d(obst_nbvar + 2,:,:))
         call xios_send_field(
     &         TRIM(obst_nout_s2d)//'_All',
     &         obst_s2d(obst_nbvar + 3,:,:))
      END IF

      ! obst_s3d is not allocated if not wanted in output
      IF (l_obstout_s3d) then
         do k=1,N
               workr(:,:,k)=obst_s3d(obst_nbvar + 1,k,:,:)
         enddo
         call xios_send_field(
     &         TRIM(obst_nout_s3d)//'_NoTurb',
     &         workr)
         do k=1,N
               workr(:,:,k)=obst_s3d(obst_nbvar + 2,k,:,:)
         enddo
         call xios_send_field(
     &         TRIM(obst_nout_s3d)//'_Turb',
     &         workr)
         do k=1,N
               workr(:,:,k)=obst_s3d(obst_nbvar + 3,k,:,:)
         enddo
         call xios_send_field(
     &         TRIM(obst_nout_s3d)//'_All',
     &         workr)
      END IF
         call xios_send_field(TRIM(obst_nout_tau), obst_tau)
         call xios_send_field(TRIM(obst_nout_fuzvz)//'_uz', obst_fuz)
         call xios_send_field(TRIM(obst_nout_fuzvz)//'_vz', obst_fvz)

#  endif /* OBSTRUCTION */


# endif /* SOLVE3D */
  
      return
      end

#else   /* XIOS */
 
      subroutine send_xios_empty
      return
      end

#endif /* XIOS */
