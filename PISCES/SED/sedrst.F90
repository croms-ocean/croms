! $Id: wrt_rst.F 1571 2014-07-01 12:38:05Z gcambon $
!
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

MODULE sedrst

#if defined key_sediment

   !! * Modules used
   USE oce_sed
   USE sed
#ifdef AGRIF
      USE param, ONLY : Lmmpi,Mmmpi
#endif
   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC sed_rst_wri         ! routine called by opa.F90
   PUBLIC sed_rst_read

   !!* Substitution
      !! * Substitutions
#  include "ocean2pisces.h90"
#  include "do_loop_substitute.h90"

CONTAINS

#if defined key_sediment
      SUBROUTINE def_rst_sed( ncid, total_rec, ierr)  ! restart netCDF

# include "netcdf.inc"

      logical :: create_new_file
      integer :: ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim    &
      &      , r2dgrd(3),  auxil(2),  checkdims                           &
#ifdef NC4PAR
      &      , csize,cmode           &
#endif
      &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4),  w3dgrd(4), jn
#ifdef USE_CALENDAR
      CHARACTER (len=19)    :: cdate,tool_sectodat
#endif
      CHARACTER(len=20) :: cltra, cvar

!
! Put time record index into file name. In  the case when model
! output is to be arranged into sequence of named files, the naming
! convention is as follows: 'rst_root.INDEX.[MPI_node.]nc', where
! INDEX is an integer number such that (i) it is divisible by the
! specified number of records per file; and (ii)
!
!      INDEX + record_within_the_file = total_record
!
! where, 1 =< record_within_the_file =< records_per_file, so that
! total_record changes continuously throughout the sequence of files.
!
      ierr=0
      lstr=lenstr(cn_sedrst_out)
      if (nrpfrst.gt.0) then
        lvar=total_rec - (1+mod(total_rec-1, nrpfrst))
        call insert_time_index (cn_sedrst_out, lstr, lvar, ierr)
#ifdef USE_CALENDAR
        if (nrpfrst.eq.1) then
          cdate = tool_sectodat(time)
          cn_sedrst_out=TRIM(cn_sedrst_out(1:lstr-9))//'.'//cdate(7:10)//cdate(4:5)
          cn_sedrst_out=TRIM(cn_sedrst_out)//cdate(1:2)//cdate(12:13)//cdate(15:16)
          cn_sedrst_out=TRIM(cn_sedrst_out)//cdate(18:19)//'.nc'
          lstr=lenstr(cn_sedrst_out)
        end if
#endif
        if (ierr .ne. 0) goto 99
      endif

!
! Decide whether to create a new file, or open existing one.
! Overall the whole code below is organized into 3-way switch,
!
! 10  if (create_new_file) then
!        .... create new file, save netCDF ids for all variables;
!     elseif (ncid.eq.-1) then
!        .... try to open existing file and check its dimensions
!       if (cannot be opened or rejected) then
!         create_new_file=.true.
!         goto 10
!       endif   and prepare
!        .... prepare the file for adding new data,
!        .... find and save netCDF ids for all variables
!     else
!        .... just open, no checking, all ids are assumed to be
!        .... already known (MPI single file output only).
!     endif
!
! which is designed to implement flexible opening policy:
! if ldefhis=.true., it forces creation of a new file [if the
! file already exists, it will be overwritten]; on the other hand,
! ldefhis=.false., it is assumed that the file already exists and
! an attempt to open it is made; if the attempt is successful, the
! file is prepared for appending hew data; if it fails, a new file
! is created.
!
      create_new_file = ldefhis
      IF (ncid .NE. -1) create_new_file = .false.
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      IF (mynode > 0) create_new_file = .false.
#endif
!
! Create new restart file:    Put global attributes
!======= === ======= =====    and define all variables.
!
  10  if (create_new_file) then

#ifndef NC4PAR
        ierr  = nf_create(cn_sedrst_out(1:lstr),NF_CLOBBER, ncid)
#else
        cmode = ior(nf_netcdf4,nf_classic_model)
        cmode = ior(cmode, nf_mpiio)
        csize = xi_rho*eta_rho/NNODES
        WRITE(stdout,*)'CREATE RST NC4 PARALLEL FILE'
        ierr  = nf_create_par(cn_sedrst_out(1:lstr),cmode, &
        &       MPI_COMM_WORLD,MPI_INFO_NULL,ncid)
#endif

        IF (ierr .NE. nf_noerr) THEN
           WRITE(stdout,'(/3(1x,A)/)') 'ERROR in DEF_RST_SED: Cannot',    &
           &             'create restart NetCDF file:', TRIM(cn_sedrst_out)
           GOTO 99                                         !--> ERROR
        ENDIF
        IF (nrpfrst == 0) total_rec = 0
!
! Put global attributes.
! --- ------ -----------
!
        CALL put_global_atts (ncid, ierr)
!
! Define dimensions of staggered fields.
! ------ ---------- -- --------- -------
!
        ierr = nf_def_dim (ncid, 'xi_rho',   xi_rho,  r2dgrd(1))
        ierr = nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr = nf_def_dim (ncid, 'profsed',    jpksed,   r3dgrd(3))
        ierr = nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr = nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        auxil(2)  = timedim

        r2dgrd(3) = timedim           ! Free surface
        r3dgrd(1) = r2dgrd(1)         !
        r3dgrd(2) = r2dgrd(2)         ! 3D RHO-type
        r3dgrd(4) = timedim           !
!
! Define evolving model variables:
! ------ -------- ----- ----------
!
!
! Time step number and time record numbers:
!
        ierr = nf_def_var (ncid, 'time_step', nf_int, 2, auxil,     &
        &       rstsedstep)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstsedstep,nf_collective)
#endif
        ierr = nf_put_att_text (ncid, rstsedstep, 'long_name', 48,    &
        &       'time step and record numbers from initialization')
!
! Time.
!
        lvar = lenstr(vname(1,indxTime))
        ierr = nf_def_var (ncid, vname(1,indxTime)(1:lvar),           &
        &                              NF_DOUBLE, 1, timedim, rstsedtime)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstsedtime,nf_collective)
#endif
        lvar = lenstr(vname(2,indxTime))
        ierr = nf_put_att_text (ncid, rstsedtime, 'long_name', lvar,     &
        &                                  vname(2,indxTime)(1:lvar))
        lvar = lenstr(vname(3,indxTime))
        ierr = nf_put_att_text (ncid, rstsedtime, 'units',     lvar,     &
        &                                  vname(3,indxTime)(1:lvar))
        lvar = lenstr (vname(4,indxTime))
        ierr = nf_put_att_text(ncid, rstsedtime, 'field',     lvar,      &
        &                                  vname(4,indxTime)(1:lvar))

!
! Time2.
!
        lvar = lenstr(vname(1,indxTime2))
        ierr = nf_def_var (ncid, vname(1,indxTime2)(1:lvar),            &
        &                              NF_DOUBLE, 1, timedim, rstsedtime2)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstsedtime2,nf_collective)
#endif
        lvar = lenstr(vname(2,indxTime2))
        ierr = nf_put_att_text (ncid, rstsedtime2, 'long_name', lvar,     &
        &                                  vname(2,indxTime2)(1:lvar))
        lvar = lenstr(vname(3,indxTime2))
        ierr = nf_put_att_text (ncid, rstsedtime2, 'units',     lvar,     &
        &                                  vname(3,indxTime2)(1:lvar))
        lvar = lenstr (vname(4,indxTime2))
        ierr = nf_put_att_text(ncid, rstsedtime2, 'field',     lvar,      &
        &                                  vname(4,indxTime2)(1:lvar))

!
! Tracer variables.
!
        DO jn = 1, jptrased
           cltra = TRIM(sedtrcd(jn))
           ierr  = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstsed(jn))
#ifdef NC4PAR
           ierr = nf_var_par_access(ncid,rstsed(jn),nf_collective)
#endif
           lvar = lenstr(TRIM(sedtrcl(jn)))
           ierr = nf_put_att_text (ncid, rstsed(jn), 'long_name',    &
           &                     lvar, TRIM(sedtrcl(jn)))
           lvar = lenstr(TRIM(sedtrcu(jn)))
           ierr = nf_put_att_text (ncid, rstsed(jn), 'units', lvar, TRIM(sedtrcu(jn)))
        END DO

        cltra = "SedPH"
        ierr  = nf_def_var (ncid, cltra, NF_DOUBLE, 4, r3dgrd, rstph)
#ifdef NC4PAR
        ierr = nf_var_par_access(ncid,rstph,nf_collective)
#endif
        cvar = "PH in sediments"
        lvar = lenstr(TRIM(cvar))
        ierr = nf_put_att_text (ncid, rstph, 'long_name',    &
        &                     lvar, TRIM(cvar) )
        cvar = "-"
        lvar = lenstr(TRIM(cvar))
        ierr = nf_put_att_text (ncid, rstph, 'units', lvar, TRIM(cvar))

        DO jn = 1, jpsol
           cltra = "burial"//TRIM(sedtrcd(jn))
           ierr  = nf_def_var (ncid, cltra, NF_DOUBLE, 3, r2dgrd, rstsol(jn))
#ifdef NC4PAR
           ierr = nf_var_par_access(ncid,rstsol(jn),nf_collective)
#endif
           lvar = lenstr(TRIM(sedtrcl(jn)))
           ierr = nf_put_att_text (ncid, rstsol(jn), 'long_name',    &
           &                     lvar, TRIM(sedtrcl(jn)))
           lvar = lenstr(TRIM(sedtrcu(jn)))
           ierr = nf_put_att_text (ncid, rstsol(jn), 'units', lvar, TRIM(sedtrcu(jn)))
        END DO        
!
! Leave definition mode.                  Also initialize record
! ----- ---------- -----                  dimension size to zero.
!
        ierr = nf_enddef(ncid)
        WRITE(*,'(6x,4A,1x,A,i4)') 'DEF_RST_SED - Created new ',        &
        &              'netCDF file ''', TRIM(cn_sedrst_out), '''.'
!
! Open an existing file and prepare for appending data.
! ==== == ======== ==== === ======= === ========= =====
! Check consistency of the dimensions of fields from the
! file with model dimensions. Determine the current size
! of unlimited dimension and set initial record [in the
! case of MPI serialized output, at this moment the last
! time record is assumed to be **partially** written by
! MPI processes with lower rank. Thus the next write is
! expected to be into the same record rather than next
! one (except MPI-master, who initializes the record).
!
! In the case when file is rejected (whether it cannot
! be opened, or something is wrong with its dimensions, 
! create new file. 
!
      ELSEIF (ncid == -1) THEN
#ifndef NC4PAR
        ierr = nf_open (TRIM(cn_sedrst_out), nf_write, ncid)
#else
        ierr = nf_open_par (TRIM(cn_sedrst_out), IOR(nf_write, nf_mpiio),   &
        &     MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#endif
        IF (ierr == nf_noerr) THEN
           ierr = checkdims (ncid, cn_sedrst_out, lstr, rec)
           IF (ierr == nf_noerr) THEN
              IF (nrpfrst == 0) THEN
                 ierr = rec+1 - nrecsedrst
              ELSE
                 ierr = rec+1 - (1+mod(nrecsedrst-1, abs(nrpfrst)))
              ENDIF
              IF (ierr > 0) THEN
                 MPI_master_only write( stdout,                              &
        &                 '(/1x,A,I5,1x,A/8x,3A,I5/8x,A,I5,1x,A/)'         &
        &           ) 'DEF_RST_SED WARNING: Actual number of records', rec,    &
        &             'in netCDF file',  '''',  TRIM(cn_sedrst_out),       &
        &             ''' exceeds the record number from restart data',    &
        &             rec+1-ierr,'/', total_rec,', restart is assumed.'
                 rec = rec-ierr
              ELSEIF (nrpfrst == 0) THEN
                 total_rec = rec+1           ! <-- set to the next record
#if defined MPI & !defined PARALLEL_FILES
                 IF (mynode > 0) total_rec = total_rec-1
#endif
              ENDIF
              ierr = nf_noerr
           ENDIF
        ENDIF

        IF (ierr .NE. nf_noerr) THEN
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
           IF (mynode == 0) THEN
              create_new_file = .true.
              GOTO 10
           ELSE
              WRITE(stdout,'(/1x,4A, 1x,A,I4/)')     'DEF_RST_SED ERROR: ',    &
              &     'Cannot open restart netCDF file ''',TRIM(cn_sedrst_out),'''.'
              GOTO 99                                     !--> ERROR 
           ENDIF
#else
           create_new_file=.true.
           GOTO 10
#endif
        ENDIF
!
! Find netCDF IDs of evolving model variables:
! ---- ------ --- -- -------- ----- ----------
!
! Time step indices:
!
        ierr = nf_inq_varid (ncid, 'time_step', rstsedstep)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) 'time_step', TRIM(cn_sedrst_out)
          GOTO 99                                         !--> ERROR
        ENDIF
!
! Time.
!
        lvar = lenstr(vname(1,indxTime))
        ierr = nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), rstsedtime)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) vname(1,indxTime)(1:lvar), TRIM(cn_sedrst_out)
          GOTO 99                                         !--> ERROR
        ENDIF
!
! Time2.
!
        lvar = lenstr(vname(1,indxTime2))
        ierr = nf_inq_varid (ncid, vname(1,indxTime2)(1:lvar), rstsedtime2)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) vname(1,indxTime2)(1:lvar), TRIM(cn_sedrst_out)
          GOTO 99                                         !--> ERROR
        ENDIF
!
! Tracer variables.
!
       DO jn = 1, jptrased
          cltra = TRIM(sedtrcd(jn))
          ierr = nf_inq_varid (ncid, cltra, rstsed(jn))
          IF (ierr .NE. nf_noerr) THEN
             WRITE(stdout,1) cltra, cn_sedrst_out(1:lstr)
             GOTO 99                                       !--> ERROR
          ENDIF
       END DO
!
       cltra = "SedPH"
       ierr = nf_inq_varid (ncid, cltra, rstph)
       IF (ierr .NE. nf_noerr) THEN
           WRITE(stdout,1) cltra, cn_sedrst_out(1:lstr)
           GOTO 99                                       !--> ERROR
       ENDIF
!
       DO jn = 1, jpsol
          cltra = "burial"//TRIM(sedtrcd(jn))
          ierr = nf_inq_varid (ncid, cltra, rstsol(jn))
          IF (ierr .NE. nf_noerr) THEN
             WRITE(stdout,1) cltra, cn_sedrst_out(1:lstr)
             GOTO 99                                       !--> ERROR
          ENDIF
       END DO

!
        MPI_master_only WRITE(*,'(6x,2A,i4,1x,A,i4)')              &
        &             'DEF_RST_SED -- Opened ',                        &
        &             'existing restart file,  record =', rec

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      ELSE
         ierr = nf_open (cn_sedrst_out(1:lstr), nf_write, ncid)
         IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,'(/1x,4A, 1x,A,I4/)')        &
            &          'DEF_RST_SED ERROR: Cannot',                       &
            &          'open restart netCDF file ''', cn_sedrst_out(1:lstr), '''.'
            GOTO 99                                         !--> ERROR
         ENDIF
#endif
      ENDIF              !<-- create_new_file
   1  FORMAT(/1x,'DEF_RST_SED ERROR: Cannot find variable ''',        &
      &               A, ''' in netCDF file ''', A, '''.'/)
  99  RETURN                                              !--> ERROR
      END SUBROUTINE def_rst_sed

      SUBROUTINE sed_rst_wri      ! variables into restart
                                  ! netCDF file.
# include "netcdf.inc"

      INTEGER :: ierr, record, lstr, lvar, lenstr   &
      &  , start(2), count(2), ibuff(2), nf_fwrite
      INTEGER :: ji, jj, jk, jn
      REAL(wp), DIMENSION(jpoce,jpksed)            :: zdta
      REAL(wp), DIMENSION(GLOBAL_2D_ARRAY,jpksed)  :: ztratmp
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY,jpksed) :: ztra
      REAL(wp), DIMENSION(GLOBAL_2D_ARRAY)         :: ztrabtmp
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY,jpsol)  :: zburial
      CHARACTER(len=20) :: cltra


#ifdef MPI
#define LOCALLM Lmmpi
#define LOCALMM Mmmpi
#else
#define LOCALLM Lm
#define LOCALMM Mm
#endif


#if defined MPI & !defined PARALLEL_FILES
      INCLUDE 'mpif.h'
      INTEGER status(MPI_STATUS_SIZE), blank
#endif

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      IF (mynode > 0) THEN
         call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,       &
         &                 1, MPI_COMM_WORLD, status, ierr)
      ENDIF
#endif
!
! Create/open restart file; write grid arrays, if so needed.
!
      CALL def_rst_sed (ncidsedrst, nrecsedrst, ierr)
      IF (ierr .NE. nf_noerr) GOTO 99
      lstr = lenstr(cn_sedrst_out)
!                                            !!! WARNING: Here it is
! Set record within the file.                !!! assumed that global
!                                            !!! restart record index 
      nrecsedrst = max(nrecsedrst,1)                 !!! nrecrst is already
      IF (nrpfrst == 0) THEN                 !!! advanced by main.
         record = nrecsedrst
      ELSE
         record = 1+mod(nrecsedrst-1, abs(nrpfrst))
      ENDIF

!
! Write out evolving model variables:
! ----- --- -------- ----- ----------
!
! Time step number and record indices. 
!
      ibuff(1) = iic
      ibuff(2) = nrecsedrst
      start(1) = 1
      start(2) = record
      count(1) = 2
      count(2) = 1
      ierr = nf_put_vara_int (ncidsedrst, rstsedstep, start, count, ibuff)
      IF (ierr .NE. nf_noerr) THEN
         WRITE(stdout,1) 'time_step', record, ierr      
         GOTO 99                                           !--> ERROR
      ENDIF
!
! Time.
!
      ierr = nf_put_var1_FTYPE (ncidsedrst, rstsedtime, record, time)
      IF (ierr .NE. nf_noerr) THEN
         lvar = lenstr(vname(1,indxTime))
         WRITE(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
         GOTO 99                                           !--> ERROR
      ENDIF
!
! Tracer variables.
!

      ! Back to 2D geometry
      DO jn = 1, jpsol
         DO jk = 1, jpksed
            ztra(:,:,jk) = UNPACK( solcp(:,jk,jn), sedmask == 1.0, 0.0 )
         END DO
         ztratmp(A2D(0),:) = ztra(A2D(0),:)
         !
         cltra = TRIM(sedtrcd(jn))
         ierr = nf_fwrite(ztratmp(START_2D_ARRAY,1), ncidsedrst,   &
         &                             rstsed(jn), record, r3dsed)
         IF(ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) cltra, record, ierr
            GOTO 99                                         !--> ERROR
         ENDIF
      END DO
      !
      DO jn = 1, jpwat
         DO jk = 1, jpksed
            ztra(:,:,jk) = UNPACK( pwcp(:,jk,jn), sedmask == 1.0, 0.0 )
         END DO
         ztratmp(A2D(0),:) = ztra(A2D(0),:)
         !
         cltra = TRIM(sedtrcd(jpsol+jn))
         ierr = nf_fwrite(ztratmp(START_2D_ARRAY,1), ncidsedrst,   &
         &                             rstsed(jpsol+jn), record, r3dsed)
         IF(ierr .NE. nf_noerr) THEN
            WRITE(stdout,1) cltra, record, ierr
            GOTO 99                                         !--> ERROR
         ENDIF
      END DO
      ! pH
      DO jk = 1, jpksed
         DO ji = 1, jpoce
            zdta(ji,jk) = -LOG10( hipor(ji,jk) / ( densSW(ji) + rtrn ) + rtrn )
         ENDDO
         ztra(:,:,jk) = UNPACK( zdta(:,jk), sedmask == 1.0, 0.0 )
      ENDDO
      ztratmp(A2D(0),:) = ztra(A2D(0),:)
      cltra = "SedPH"
      ierr = nf_fwrite(ztratmp(START_2D_ARRAY,1), ncidsedrst,   &
      &                             rstph, record, r3dsed)
      IF (ierr .NE. nf_noerr) THEN
         WRITE(stdout,1) cltra, record, ierr
         GOTO 99                                         !--> ERROR
      ENDIF

      ! prognostic variables
      ! --------------------
      DO jn = 1, jpsol
         zburial(:,:,jn)   = UNPACK( burial(:,jn), sedmask == 1.0, 0.0 )
         ztrabtmp(A2D(0)) = zburial(A2D(0),jn)
         !
         cltra = "burial"//TRIM(sedtrcd(jn))
         ierr = nf_fwrite(ztrabtmp(START_2D_ARRAY), ncidsedrst,   &
         &                             rstsol(jn), record, r2dvar)
        IF (ierr .NE. nf_noerr) THEN
          WRITE(stdout,1) cltra, record, ierr
          GOTO 99                                         !--> ERROR
        ENDIF
      END DO

  1   FORMAT(/1x, 'SED_RST_WRI ERROR while writing variable ''', A,   &
       &           ''' into restart file.', /11x, 'Time record:', &
       &               i6, 3x, 'netCDF error code', i4, 3x, A,i4) 
      GOTO 100 
  99  may_day_flag=3
 100  CONTINUE

!
! Synchronize restart netCDF file to disk to allow other
! processes to access data immediately after it is written.
!
#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      ierr = nf_close (ncidsedrst)
      IF (nrpfrst > 0 .AND. record >= nrpfrst) ncidsedrst = -1
#else
      IF (nrpfrst > 0 .AND. record >= nrpfrst) THEN
        ierr = nf_close (ncidsedrst)
        ncidsedrst = -1
      ELSE
        ierr = nf_sync(ncidsedrst)
      ENDIF
#endif
      IF (ierr == nf_noerr) THEN
         MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')    & 
         &            'WRT_RST_SED -- wrote ',                          &
         &            'restart fields into time record =', record, '/',  &
         &             nrecrst  
      ELSE
         MPI_master_only  write(stdout,'(/1x,2A/)')     & 
         &             'WRT_RST_SED ERROR: Cannot ',        &
         &             'synchronize/close restart netCDF file.'
         may_day_flag = 3
      ENDIF

#if defined MPI & !defined PARALLEL_FILES  & !defined NC4PAR
      IF (mynode < NNODES-1) THEN
         CALL MPI_Send (blank, 1, MPI_INTEGER, mynode+1, 1, MPI_COMM_WORLD,  ierr)
      ENDIF
#endif
      RETURN
      END SUBROUTINE sed_rst_wri 


                              ! Read initial conditions for the
      SUBROUTINE sed_rst_read ! primitive variables from NetCDF
                              ! initialization file.

!======================================================
!
!======================================================

# include "netcdf.inc"

      real(wp) :: time_scale
      integer  :: ji, jj, jk, jn
      integer  :: ncid, indx, varid,  ierr, lstr, lvar, latt, lenstr, max_rec,  &
      &        start(2), count(2), ibuff(2), nf_fread, checkdims
      character :: units*180
      REAL(wp), DIMENSION(jpoce,jpksed)            :: zhipor
      REAL(wp), DIMENSION(GLOBAL_2D_ARRAY,jpksed)  :: ztratmp
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY,jpksed) :: ztra
      REAL(wp), DIMENSION(GLOBAL_2D_ARRAY)         :: ztrabtmp
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY)        :: ztrab
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY,jpsol)  :: zburial
      CHARACTER(len=20) :: cltra


#ifdef MPI
#define LOCALLM Lmmpi
#define LOCALMM Mmmpi
#else
#define LOCALLM Lm
#define LOCALMM Mm
#endif

!
! Open initial conditions netCDF file for reading. Check that all
! spatial dimensions in that file are consistent with the model
! arrays, determine how many time records are available in the file
! and set record from which the dada will be read.
!
! The record is set as follows: (1) if only one time record is
! available in the file, then that record is used REGARDLESS of
! value of nrrec supplied from the parameter file; (2) if the
! file has multiple records and nrrec is positive, then nrrec is
! available record is used.
!
      if (may_day_flag.ne.0) return      !-->  EXIT
      lstr = lenstr(cn_sedrst_in)
      ierr = nf_open(TRIM(cn_sedrst_in), nf_nowrite, ncid)
      if (ierr .eq. nf_noerr) then
        ierr=checkdims (ncid, cn_sedrst_in, lstr, max_rec)
        if (ierr .eq. nf_noerr) then
          if (max_rec.gt.0) then
            if (nrrec.gt.0) then
              if (nrrec.le.max_rec) then
                indx=nrrec
              else
       MPI_master_only write(stdout,'(/1x,2A,I4,1x,A/12x,A,I4,1x,3A/)')  &
     &           'SED_RST_READ ERROR :: requested restart time ', &
     &           'record',  nrrec, 'exceeds number',  'of records', &
     &                       max_rec,  'available in netCDF file ''',  &
     &                                        TRIM(cn_sedrst_in), '''.'

                goto 99         !--> ERROR
              endif
            else
              indx=max_rec
            endif
          else
            indx=1
          endif
        else
          goto 99
        endif
      else
        MPI_master_only write(stdout,'(/1x,4A/12x,A/)')  &
     &               'SED_RST_READ WARNING :: ',    &
     &               'Cannot open netCDF file ''',   TRIM(cn_sedrst_in),  &
     &                                  '''.',     nf_strerror(ierr)
        goto 99                                           !--> ERROR
      endif
!
! Read in evolving model variables:
! ---- -- -------- ----- ----------
!
! Time: find netCDF id, read value, read attribute 'units'
! and set starting time index and time clock in days.
!
      lvar = lenstr(vname(1,indxTime))
      ierr = nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), varid)
      IF (ierr == nf_noerr) THEN
        ierr = nf_get_var1_FTYPE (ncid, varid, indx, time)
        IF (ierr == nf_noerr) THEN
          ierr = nf_get_att_text(ncid, varid, 'units', units)
          IF (ierr == nf_noerr) THEN
            latt = lenstr(units)
            IF (units(1:6) == 'second') THEN
               time_scale = 1.
            ELSEIF (units(1:3) == 'day') THEN
              time_scale = day2sec
            ELSE
              WRITE (stdout,'(/1x,4A/8x,3A/)') 'SED_RST_READ ',      &
       &              'ERROR: unknown units of for variable ''',      &
       &               vname(1,indxTime)(1:lvar), '''',               &
       &              'in netCDF file ''', TRIM(cn_sedrst_in),'''.'
              GOTO 99                                    !--> ERROR
            ENDIF
          ELSE
            WRITE (stdout,'(/1x,2A/8x,5A/)') 'SED_RST_READ ERROR: ',   &
       &             'cannot read attribute ''units'' for variable',  &
       &             '''', vname(1,indxTime)(1:lvar),                 &
       &             ''' in netCDF file ''',  TRIM(cn_sedrst_in), '''.'
            GOTO 99                                       !--> ERROR
          ENDIF
        ELSE
          MPI_master_only write(stdout,2) vname(1,indxTime)(1:lvar)  &
          &                                , indx, TRIM(cn_sedrst_in)
          GOTO 99                                         !--> ERROR
        ENDIF
      ELSE
        MPI_master_only write(stdout,1) vname(1,indxTime)(1:lvar), TRIM(cn_sedrst_in)
        GOTO 99                                           !--> ERROR
      ENDIF

!      time = time*time_scale
!      tdays = time*sec2day

      ierr = nf_inq_varid (ncid, 'time_step', varid)
      IF (ierr == nf_noerr) THEN
         start(1) = 1
         start(2) = indx
         count(1) = 2
         count(2) = 1
         ierr = nf_get_vara_int (ncid, varid, start, count, ibuff)
         IF (ierr == nf_noerr) THEN
!            ntstart = ibuff(1)
            nrecsedrst = ibuff(2)

            MPI_master_only WRITE(stdout,                            &
            &     '(6x,A,G12.4,A,I2,A,I6,A,I3,A)')              &
            &     'SED_RST_READ: Restarted from day =', tdays, ' rec =',   &
            &      indx, '(', ntstart, ',', nrecsedrst, ').'

         ELSE
            MPI_master_only write(stdout,'(/1x,2A/)')                     &
            &                            'SED_RST_READ ERROR: Cannot ',    &
            &                            'read time and record indices.'
            GOTO 99                                         !--> ERROR
         ENDIF
      ELSE
 !        ntstart = 1
         nrecsedrst = 0
         MPI_master_only WRITE(stdout,'(6x,2A,G12.4,1x,A,I4)')      &
         &          'SED_RST_READ -- ',                              &
         &          'Processing data for time =', tdays, 'record =', indx
      ENDIF
!      IF (ntstart < 1) ntstart = 1
!      ntimes = ntstart+ntimes-1
!
! Tracer variables.
!
      DO jn = 1, jptrased
        IF( jn <= jpsol ) THEN
           cltra = TRIM(sedtrcd(jn))
        ELSE
           cltra = TRIM(sedtrcd(jpsol+jn))
        ENDIF
        ierr = nf_inq_varid (ncid, cltra, varid)
        IF (ierr == nf_noerr) THEN
          ierr = nf_fread (ztratmp(START_2D_ARRAY,1), ncid,  varid, indx, r3dsed)
          IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_sedrst_in)
            GOTO 99                                       !--> ERROR
          ENDIF
          MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_sedrst_in)
        ENDIF
        !
        ztra(A2D(0),:) = ztratmp(A2D(0),:)
        !
        IF( jn <= jpsol ) THEN
          DO jk = 1, jpksed
             solcp(:,jk,jn) = PACK( ztra(:,:,jk), sedmask == 1.0 )
          END DO
        ELSE
          DO jk = 1, jpksed
             pwcp(:,jk,jn-jpsol) = PACK( ztra(:,:,jk), sedmask == 1.0 )
          END DO
        ENDIF   
      END DO
      ! Initialization of sediment composant only ie jk=2 to jk=jpksed 
      ! ( nothing in jk=1)
      solcp(1:jpoce,1,:) = 0.
      pwcp (1:jpoce,1,:) = 0.

      cltra = "SedPH"
      ierr = nf_inq_varid (ncid, cltra, varid)
      IF(ierr == nf_noerr) THEN
         ierr = nf_fread (ztratmp(START_2D_ARRAY,1), ncid,  varid, indx, r3dsed)
         IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_sedrst_in)
            GOTO 99                                       !--> ERROR
            MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_sedrst_in)
         ENDIF
        !
        ztra(A2D(0),:) = ztratmp(A2D(0),:)
        !
      ENDIF

      zhipor(:,:) = 0.
      DO jk = 1, jpksed
         zhipor(:,jk) = PACK(ztra(:,:,jk), sedmask == 1.0 )
      END DO

      ! Initialization of [h+] in mol/kg
      DO jk = 1, jpksed
         DO ji = 1, jpoce
            hipor(ji,jk) = 10.**( -1. * zhipor(ji,jk) )
         ENDDO
      ENDDO
!
! Additional variables.

      DO jn = 1, jpsol
        cltra = "burial"//TRIM(sedtrcd(jn))
        ierr = nf_inq_varid (ncid, cltra, varid)
        IF (ierr == nf_noerr) THEN
          ierr = nf_fread (ztrabtmp(START_2D_ARRAY), ncid,  varid, indx, r2dvar)
          IF (ierr .NE. nf_noerr) THEN
            MPI_master_only WRITE(stdout,2) cltra, indx, TRIM(cn_sedrst_in)
            GOTO 99                                       !--> ERROR
          ENDIF
          MPI_master_only WRITE(stdout,3) cltra, TRIM(cn_sedrst_in)
        ENDIF
        !
        ztrab(A2D(0)) = ztrabtmp(A2D(0))
        !
        burial(:,jn) = PACK( ztrab(:,:), sedmask == 1.0 )
        !
      END DO

!======================================================
! END MODIF_JG_2
!======================================================

!
!  Close input NetCDF file.
!
      ierr = nf_close(ncid)

  1   FORMAT(/1x,'SED_RST_READ - unable to find variable:',    1x,A,    &
      &                            /15x,'in input NetCDF file:',1x,A/)
  2   FORMAT(/1x,'SED_RST_READ - error while reading variable:',1x, A,  &
      &    2x,'at time record =',i4/15x,'in input NetCDF file:',1x,A/)
  3   format(/1x,'SED_RST_READ : Read variable ',    1x,A,    &
     &                            /15x,'in input NetCDF file:',1x,A/)
      RETURN
  99  may_day_flag = 2
      RETURN
      END SUBROUTINE sed_rst_read
#else
      SUBROUTINE sed_rst_wri 
      END SUBROUTINE sed_rst_wri
      SUBROUTINE sed_rst_read
      END SUBROUTINE sed_rst_read
#endif      

#endif

END MODULE sedrst
