#include "cppdefs.h"

MODULE lib_mpp

#if defined key_pisces

   use scalars

   IMPLICIT NONE
   PUBLIC

# if defined MPI
  include 'mpif.h'
# include "mpi_cpl.h"
# endif

  PUBLIC  mpp_sum, mpp_max, mpp_min
  PUBLIC  ctl_opn, ctl_nam, ctl_stop, ctl_warn
  PUBLIC  timing_start, timing_stop
  PUBLIC  load_nml

  INTERFACE mpp_sum
     MODULE PROCEDURE mppsum_int, mppsum_real
  END INTERFACE
  INTERFACE mpp_max
     MODULE PROCEDURE mppmax_int, mppmax_real
  END INTERFACE
  INTERFACE mpp_min
     MODULE PROCEDURE mppmin_int, mppmin_real
  END INTERFACE

#  include "ocean2pisces.h90"
#  include "do_loop_substitute.h90"

CONTAINS

  FUNCTION getunit()
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  getunit  ***
      !!
      !! ** Purpose :   return the index of an unused logical unit
      !!----------------------------------------------------------------------
      INTEGER :: getunit
      LOGICAL :: llopn
      !!----------------------------------------------------------------------
      !
      getunit = 15   ! choose a unit that is big enough then it is not already used in NEMO
      llopn = .TRUE.
      DO WHILE( (getunit < 998) .AND. llopn )
         getunit = getunit + 1
         INQUIRE( unit = getunit, opened = llopn )
      END DO
      IF( (getunit == 999) .AND. llopn ) THEN
         CALL ctl_stop( 'STOP', 'getunit: All logical units until 999 are used...' )
         getunit = -1
      ENDIF
      !
   END FUNCTION getunit

   SUBROUTINE ctl_opn ( knum, cdfile, cdstat, cdform, cdacce, klengh, kout, ldwp, karea )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctl_opn  ***
      !!
      !! ** Purpose :   Open file and check if required file is available.
      !!
      !! ** Method  :   Fortan open
      !!
      !! History :
      !!        !  1995-12  (G. Madec)  Original code
      !!   8.5  !  2002-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------

      INTEGER          , INTENT(  out) ::   knum      ! logical unit to open
      CHARACTER(len=*) , INTENT(in   ) ::   cdfile    ! file name to open
      CHARACTER(len=*) , INTENT(in   ) ::   cdstat    ! disposition specifier
      CHARACTER(len=*) , INTENT(in   ) ::   cdform    ! formatting specifier
      CHARACTER(len=*) , INTENT(in   ) ::   cdacce    ! access specifier
      INTEGER          , INTENT(in   ) ::   klengh    ! record length
      INTEGER          , INTENT(in   ) ::   kout      ! number of logical units for write
      LOGICAL          , INTENT(in   ) ::   ldwp      ! boolean term for print
      INTEGER, OPTIONAL, INTENT(in   ) ::   karea     ! proc number
      !!
      CHARACTER(len=80) ::   clfile
      INTEGER           ::   iost

      ! adapt filename
      ! ----------------
      clfile = TRIM(cdfile)
      IF( PRESENT( karea ) ) THEN
         IF( karea > 1 )   WRITE(clfile, "(a,'_',i4.4)") TRIM(clfile), karea-1
      ENDIF
      knum=getunit()

      iost=0
      IF( cdacce(1:6) == 'DIRECT' )  THEN
         OPEN( UNIT=knum, FILE=clfile, FORM=cdform, ACCESS=cdacce, STATUS=cdstat, RECL=klengh, ERR=100, IOSTAT=iost )
      ELSE
         OPEN( UNIT=knum, FILE=clfile, FORM=cdform, ACCESS=cdacce, STATUS=cdstat             , ERR=100, IOSTAT=iost )
      ENDIF
      IF( iost == 0 ) THEN
         IF(ldwp) THEN
            WRITE(kout,*)
            WRITE(kout,*) '     file   : ', clfile,' open ok'
            WRITE(kout,*) '     unit   = ', knum
            WRITE(kout,*) '     status = ', cdstat
            WRITE(kout,*) '     form   = ', cdform
            WRITE(kout,*) '     access = ', cdacce
            WRITE(kout,*)
         ENDIF
      ENDIF
100   CONTINUE
      IF( iost /= 0 ) THEN
         IF(ldwp) THEN
            WRITE(kout,*)
            WRITE(kout,*) ' ===>>>> : bad opening file: ', clfile
            WRITE(kout,*) ' =======   ===  '
            WRITE(kout,*) '           unit   = ', knum
            WRITE(kout,*) '           status = ', cdstat
            WRITE(kout,*) '           form   = ', cdform
            WRITE(kout,*) '           access = ', cdacce
            WRITE(kout,*) '           iostat = ', iost
            WRITE(kout,*) '           we stop. verify the file '
            WRITE(kout,*)
         ENDIF
         STOP 'ctl_opn bad opening'
      ENDIF
      
   END SUBROUTINE ctl_opn


   SUBROUTINE ctl_nam ( kios, cdnam, ldwarn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctl_nam  ***
      !!
      !! ** Purpose :   Informations when error while reading a namelist
      !!
      !! ** Method  :   Fortan open
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(inout) ::   kios     ! IO status after reading the namelist
      CHARACTER(len=*) , INTENT(in   ) ::   cdnam    ! group name of namelist for which error occurs
      LOGICAL, OPTIONAL, INTENT(in   ) ::   ldwarn   ! if absent or .TRUE., missing namelist-group records trigger a warning
      !
      CHARACTER(len=5) ::   clios    ! string to convert iostat in character for print
      LOGICAL          ::   llwarn   ! auxiliary variable
      !!----------------------------------------------------------------------
      !
      WRITE (clios, '(I5.0)')   kios
      llwarn = .TRUE.
      IF ( PRESENT( ldwarn ) ) llwarn = ldwarn
      IF( llwarn .AND. kios < 0 ) THEN
         CALL ctl_warn( 'end of record or file while reading namelist ' // TRIM(cdnam) // ' iostat = ' // TRIM(clios) )
      ENDIF
      !
      IF ( kios > 0 ) THEN
         CALL ctl_stop( 'STOP: ','misspelled variable in namelist ' // TRIM(cdnam) // ' iostat = ' // TRIM(clios) )
      ENDIF
      kios = 0
      !
   END SUBROUTINE ctl_nam

   INTEGER FUNCTION get_unit()
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  get_unit  ***
      !!
      !! ** Purpose :   return the index of an unused logical unit
      !!----------------------------------------------------------------------
      LOGICAL :: llopn
      !!----------------------------------------------------------------------
      !
      get_unit = 15   ! choose a unit that is big enough then it is not already used in NEMO
      llopn = .TRUE.
      DO WHILE( (get_unit < 9999) .AND. llopn )
         get_unit = get_unit + 1
         INQUIRE( unit = get_unit, opened = llopn )
      END DO
      IF( (get_unit == 9999) .AND. llopn ) THEN
         CALL ctl_stop( 'STOP',  'get_unit: All logical units until 9999 are used...' )
      ENDIF
      !
   END FUNCTION get_unit

   SUBROUTINE ctl_warn(clname)
      CHARACTER(len=*), INTENT(in) :: clname
      WRITE(numout,"(/,' ===>>> : W A R N I N G', /,'         ===============',/)") 
      IF(mynode .eq. 0 )  WRITE(numout,*) clname
   END SUBROUTINE

   SUBROUTINE ctl_stop(cd1, clname)
      CHARACTER(len=*), INTENT(in   )           ::   cd1
      CHARACTER(len=*), INTENT(in) :: clname
      WRITE(numout,"(/,' ===>>> : E R R O R',     /,'         ===========',/)") 
      IF(mynode .eq. 0 )  WRITE(numout,*) clname
      STOP
   END SUBROUTINE


   SUBROUTINE mppsum_int(ktab)
      INTEGER, INTENT(inout) :: ktab
#ifdef MPI
      INTEGER :: ierror, localcomm     
      INTEGER :: iwork     

      localcomm=MPI_COMM_WORLD
      CALL mpi_allreduce( ktab, iwork, 1, mpi_integer, mpi_sum, localcomm, ierror)
      ktab = iwork
#endif
   END SUBROUTINE mppsum_int

   SUBROUTINE mppsum_real(ptab)
      REAL, INTENT(inout) :: ptab
#ifdef MPI
      INTEGER :: ierror, localcomm     
      REAL    :: zwork     

      localcomm=MPI_COMM_WORLD
      CALL mpi_allreduce( ptab, zwork, 1, mpi_double_precision, mpi_sum, localcomm, ierror )
      ptab = zwork
#endif
   END SUBROUTINE mppsum_real

   SUBROUTINE mppmax_int(ktab)
      INTEGER, INTENT(inout) :: ktab
#ifdef MPI
      INTEGER :: ierror, localcomm     
      INTEGER :: iwork     

      localcomm=MPI_COMM_WORLD
      CALL mpi_allreduce( ktab, iwork, 1, mpi_integer, mpi_max, localcomm, ierror)
      ktab = iwork
#endif
   END SUBROUTINE mppmax_int

   SUBROUTINE mppmax_real(ptab)
      REAL, INTENT(inout) :: ptab
#ifdef MPI
      INTEGER :: ierror, localcomm     
      REAL    :: zwork     

      localcomm=MPI_COMM_WORLD
      CALL mpi_allreduce( ptab, zwork, 1, mpi_double_precision, mpi_max, localcomm, ierror )
      ptab = zwork
#endif
   END SUBROUTINE mppmax_real

    SUBROUTINE mppmin_int(ktab)
      INTEGER, INTENT(inout) :: ktab
#ifdef MPI
      INTEGER :: ierror, localcomm     
      INTEGER :: iwork     

      localcomm=MPI_COMM_WORLD
      CALL mpi_allreduce( ktab, iwork, 1, mpi_integer, mpi_min, localcomm, ierror)
      ktab = iwork
#endif
   END SUBROUTINE mppmin_int

   SUBROUTINE mppmin_real(ptab)
      REAL, INTENT(inout) :: ptab
#ifdef MPI
      INTEGER :: ierror, localcomm     
      REAL    :: zwork     

      localcomm=MPI_COMM_WORLD
      CALL mpi_allreduce( ptab, zwork, 1, mpi_double_precision, mpi_min, localcomm, ierror )
      ptab = zwork
#endif
   END SUBROUTINE mppmin_real


   SUBROUTINE timing_start( cdinfo, ldstatplot )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_start  ***
      !! ** Purpose :   collect execution time
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) :: cdinfo
      LOGICAL, OPTIONAL, INTENT(in) :: ldstatplot   ! .true. if you want to call gnuplot analyses on this timing
      !
   END SUBROUTINE timing_start


   SUBROUTINE timing_stop( cdinfo, kt, ld_finalize )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_stop  ***
      !! ** Purpose :   stop timing window
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdinfo
      INTEGER, OPTIONAL, INTENT(in) ::   kt
      LOGICAL, OPTIONAL, INTENT(in) ::   ld_finalize
      !
   END SUBROUTINE timing_stop

   SUBROUTINE load_nml( cdnambuff , cdnamfile, kout, ldwp)
      CHARACTER(LEN=:)    , ALLOCATABLE, INTENT(INOUT) :: cdnambuff
      CHARACTER(LEN=*), INTENT(IN )                :: cdnamfile
      CHARACTER(LEN=256)                           :: chline
      CHARACTER(LEN=1)                             :: csp
      INTEGER, INTENT(IN)                          :: kout
      LOGICAL, INTENT(IN)                          :: ldwp  !: .true. only for the root broadcaster
      INTEGER                                      :: itot, iun, iltc, inl, ios, itotsav
      !
      !csp = NEW_LINE('A')
      ! a new line character is the best seperator but some systems (e.g.Cray)
      ! seem to terminate namelist reads from internal files early if they
      ! encounter new-lines. Use a single space for safety.
      csp = ' '
      !
      ! Check if the namelist buffer has already been allocated. Return if it has.
      !
      IF ( ALLOCATED( cdnambuff ) ) RETURN
      IF( ldwp ) THEN
         !
         ! Open namelist file
         !
         CALL ctl_opn( iun, cdnamfile, 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, kout, ldwp )
         !
         ! First pass: count characters excluding comments and trimable white space
         !
         itot=0
     10  READ(iun,'(A256)',END=20,ERR=20) chline
         iltc = LEN_TRIM(chline)
         IF ( iltc.GT.0 ) THEN
          inl = INDEX(chline, '!')
          IF( inl.eq.0 ) THEN
           itot = itot + iltc + 1                                ! +1 for the newline character
          ELSEIF( inl.GT.0 .AND. LEN_TRIM( chline(1:inl-1) ).GT.0 ) THEN
           itot = itot + inl                                  !  includes +1 for the newline character
          ENDIF
         ENDIF
         GOTO 10
     20  CONTINUE
         !
         ! Allocate text cdnambuff for condensed namelist
         !
!$AGRIF_DO_NOT_TREAT
         ALLOCATE( CHARACTER(LEN=itot) :: cdnambuff )
!$AGRIF_END_DO_NOT_TREAT
         itotsav = itot
         !
         ! Second pass: read and transfer pruned characters into cdnambuff
         !
         REWIND(iun)
         itot=1
     30  READ(iun,'(A256)',END=40,ERR=40) chline
         iltc = LEN_TRIM(chline)
         IF ( iltc.GT.0 ) THEN
          inl = INDEX(chline, '!')
          IF( inl.eq.0 ) THEN
           inl = iltc
          ELSE
           inl = inl - 1
          ENDIF
          IF( inl.GT.0 .AND. LEN_TRIM( chline(1:inl) ).GT.0 ) THEN
             cdnambuff(itot:itot+inl-1) = chline(1:inl)
             WRITE( cdnambuff(itot+inl:itot+inl), '(a)' ) csp
             itot = itot + inl + 1
          ENDIF
         ENDIF
         GOTO 30
     40  CONTINUE
         itot = itot - 1
         IF( itotsav .NE. itot ) WRITE(*,*) 'WARNING in load_nml. Allocated ',itotsav,' for read buffer; but used ',itot
         !
         ! Close namelist file
         !
         CLOSE(iun)
         !write(*,'(32A)') cdnambuff
      ENDIF
#if ! defined key_mpi_off
      CALL mpp_bcast_nml( cdnambuff, itot )
#endif
   END SUBROUTINE load_nml

   SUBROUTINE mpp_bcast_nml( cdnambuff , kleng )
      CHARACTER(LEN=:)    , ALLOCATABLE, INTENT(INOUT) :: cdnambuff
      INTEGER                          , INTENT(INOUT) :: kleng
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_bcast_nml  ***
      !!
      !! ** Purpose :   broadcast namelist character buffer
      !!
      !!----------------------------------------------------------------------
      !!
#ifdef MPI
      INTEGER ::   iflag, localcomm

      localcomm=MPI_COMM_WORLD
      call MPI_BCAST(kleng, 1, MPI_INT, 0, localcomm, iflag)
      call MPI_BARRIER(localcomm, iflag)
!$AGRIF_DO_NOT_TREAT
      IF ( .NOT. ALLOCATED(cdnambuff) ) ALLOCATE( CHARACTER(LEN=kleng) :: cdnambuff )
!$AGRIF_END_DO_NOT_TREAT
      call MPI_BCAST(cdnambuff, kleng, MPI_CHARACTER, 0, localcomm, iflag)
      call MPI_BARRIER(localcomm, iflag)
#endif
      !
   END SUBROUTINE mpp_bcast_nml

#endif

END MODULE lib_mpp
