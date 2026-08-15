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
!----------------------------------------------------------------
! Add attribute if any, begining at index indx = 5, to 11 for now
! but it can be updated easily
!====
!  5 :  the standard_name
!  6 :  coordinates
!  7 :  axis
!  8 :  _FillValue
!  9 :  positive for vertical axis
! 10 :  time origin
! 11 :  calendar for time
!==
! 12 :
! 13 :
!---------------------------------------------------------------
#include "cppdefs.h"
      subroutine nf_add_attribute(ncid, vid, vnameindx, i, nf_type, ierr)
      implicit none
      integer ncid,vid,lvar,ierr,lenstr,i, nf_type
      character*75  vnameindx(20)
#include "param.h"
#include "scalars.h"
#include "ncscrum.h"
#include "netcdf.inc"

      character*70 text

      text=vnameindx(i)
      lvar=lenstr(text)
      if (LEN_TRIM(text).gt.0) then
         ierr=nf_put_att_text (ncid, vid, 'standard_name',
     &        lvar, text(1:lvar))
      endif

      text=vnameindx(i+1)
      lvar=lenstr(text)
      if (LEN_TRIM(text).gt.0) then
         ierr=nf_put_att_text (ncid, vid, 'coordinates',
     &        lvar, text(1:lvar))
      endif

      text=vnameindx(i+2)
      lvar=lenstr(text)
      if (LEN_TRIM(text).gt.0) then
         ierr=nf_put_att_text (ncid, vid, 'axis',
     &        lvar, text(1:lvar))
      endif

      text=vnameindx(i+4)
      lvar=lenstr(text)
      if (LEN_TRIM(text).gt.0) then
         ierr=nf_put_att_text (ncid, vid, 'positive',
     &        lvar, text(1:lvar))
      endif
      text=vnameindx(i+5)
      lvar=lenstr(text)
      if (LEN_TRIM(text).gt.0) then
         ierr=nf_put_att_text (ncid, vid, 'time_origin',
     &        lvar, text(1:lvar))
      endif

      text=vnameindx(i+6)
      lvar=lenstr(text)
      if (LEN_TRIM(text).gt.0) then
         ierr=nf_put_att_text (ncid, vid, 'calendar',
     &        lvar, text(1:lvar))
      endif
#ifdef FILLVAL
      text=vnameindx(i+3)
      if (TRIM(text) /= '.NO.') then
       if (nf_type == NF_REAL) then
        ierr= nf_put_att_real (ncid, vid,'_FillValue', nf_type,
     &     1, spval)
       elseif(nf_type == NF_DOUBLE) then
        ierr= nf_put_att_double (ncid, vid,'_FillValue', nf_type,
     &     1, spval)
       elseif(nf_type == NF_INT) then
        ierr= nf_put_att_int (ncid, vid,'_FillValue', nf_type,
     &     1, int(spval))
       endif
      endif
#endif
      return
      end

