#include "cppdefs.h"


MODULE trcnam_pisces
   !!======================================================================
   !!                      ***  MODULE trcnam_pisces  ***
   !! TOP :   initialisation of some run parameters for PISCES bio-model
   !!======================================================================
   !! History :    -   !  1999-10 (M.A. Foujols, M. Levy) original code
   !!              -   !  2000-01 (L. Bopp) hamocc3, p3zd
   !!             1.0  !  2003-08 (C. Ethe)  module F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.pisces.h90
#if defined key_pisces
   !!----------------------------------------------------------------------
   !! trc_nam_pisces   : PISCES model namelist read
   !!----------------------------------------------------------------------
   USE oce_trc
   USE sms_pisces
   USE trc
   USE sed

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_pisces   ! called by trcnam.F90 module

   !! * Substitutions
#  include "ocean2pisces.h90"
#  include "read_nml_substitute.h90"
#  include "do_loop_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcnam_pisces.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_nam_pisces
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_pisces  ***  
      !!
      !! ** Purpose :   read PISCES namelist
      !!
      !! ** input   :   file 'namelist.trc.sms' containing the following
      !!             namelist: natext, natbio, natsms
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      INTEGER  :: jn, ierr, ioptio
      TYPE(PTRACER), DIMENSION(jptra) :: tracer
      CHARACTER(LEN=20)::   clname

      NAMELIST/nampistrc/ tracer
      !!
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_pisces : read PISCES namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      !                               ! Open the namelist file
      !                               ! ----------------------
      clname = 'namelist_pisces'
      CALL ctl_opn( numnatp_ref, TRIM( clname )//'_ref', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwm )
      CALL ctl_opn( numnatp_cfg, TRIM( clname )//'_cfg', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwm )
      
!      CALL load_nml( numnatp_ref, TRIM( clname )//'_ref', numout, lwm )
!      CALL load_nml( numnatp_cfg, TRIM( clname )//'_cfg', numout, lwm )
      IF(lwm) CALL ctl_opn( numonp, 'output.namelist.pis', &
              &            'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

      IF(lwp) WRITE(numout,*) 'number of tracer : ', jptra

      ALLOCATE( ctrcnm(jptra), ctrcnl(jptra), ctrcnu(jptra), STAT = ierr  )  
      IF( ierr /= 0 )   CALL ctl_stop('STOP','trc_nam_pisces: failed to allocate arrays')

      DO jn = 1, jptra
         WRITE( ctrcnm(jn),'("TR_",I2)'           ) jn
         WRITE( ctrcnl(jn),'("TRACER NUMBER ",I2)') jn
         ctrcnu(jn) = 'mmole/m3'
      END DO

      READ_NML_REF(numnatp,nampistrc)
      READ_NML_CFG(numnatp,nampistrc)
      IF(lwm) WRITE( numonp, nampistrc )


#if defined key_pisces_quota
      ln_p5z = .true.
      ln_p2z = .false.
      ln_p4z = .false.
#elif defined key_pisces_npzd
      ln_p4z = .false.
      ln_p2z = .true.
      ln_p5z = .false.
#else
      ln_p2z = .false.
      ln_p4z = .true.
      ln_p5z = .false.
#endif
#if defined key_ligand
      ln_ligand = .true.
#else
      ln_ligand = .false.
#endif
#if defined key_sediment
      ln_sediment = .true.
#else
      ln_sediment = .false.
!      ln_sed_2way = .false.
#endif
#if defined key_trc_diaadd
      l_diaadd = .true.
#else
      l_diaadd = .false.
#endif
      ln_top_euler = .TRUE.
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*) '   Namelist : nampismod '
         WRITE(numout,*) '      Flag to use PISCES reduced model     ln_p2z      = ', ln_p2z
         WRITE(numout,*) '      Flag to use PISCES standard model    ln_p4z      = ', ln_p4z
         WRITE(numout,*) '      Flag to use PISCES quota    model    ln_p5z      = ', ln_p5z
         WRITE(numout,*) '      Flag to ligand                       ln_ligand   = ', ln_ligand
         WRITE(numout,*) '      Flag to use sediment                 ln_sediment = ', ln_sediment
      ENDIF
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         IF( ln_p5z      )  WRITE(numout,*) '   ==>>>   PISCES QUOTA model is used'
         IF( ln_p4z      )  WRITE(numout,*) '   ==>>>   PISCES STANDARD model is used'
         IF( ln_p2z      )  WRITE(numout,*) '   ==>>>   PISCES REDUCED model is used'
         IF( ln_ligand   )  WRITE(numout,*) '   ==>>>   Compute remineralization/dissolution of organic ligands'
         IF( ln_sediment )  WRITE(numout,*) '   ==>>>   Sediment module is used'
      ENDIF
    
      ioptio = 0
      IF( ln_p2z )    ioptio = ioptio + 1
      IF( ln_p4z )    ioptio = ioptio + 1
      IF( ln_p5z )    ioptio = ioptio + 1
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'STOP','Choose ONE PISCES model namelist nampismod' )

      DO jn = 1, jptra
         ctrcnm(jn) = tracer(jn)%clsname
         ctrcnl(jn) = tracer(jn)%cllname
         ctrcnu(jn) = tracer(jn)%clunit
      END DO


      IF(lwp) THEN                   ! control print
         DO jn = 1, jptra
            WRITE(numout,*) '   tracer nb             : ', jn 
            WRITE(numout,*) '   short name            : ', TRIM(ctrcnm(jn))
            WRITE(numout,*) '   long name             : ', TRIM(ctrcnl(jn))
            WRITE(numout,*) '   unit                  : ', TRIM(ctrcnu(jn))
            WRITE(numout,*) ' '
         END DO
      ENDIF
      !
   END SUBROUTINE trc_nam_pisces

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                            No PISCES biochemical model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_pisces             ! Empty routine
   END SUBROUTINE trc_nam_pisces
#endif
   !!======================================================================
END MODULE trcnam_pisces
