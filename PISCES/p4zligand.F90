#include "cppdefs.h"

MODULE p4zligand
   !!======================================================================
   !!                         ***  MODULE p4zligand  ***
   !! TOP :   PISCES Compute remineralization/dissolution of organic ligands
   !!=========================================================================
   !! History :   3.6  !  2016-03  (O. Aumont, A. Tagliabue) Quota model and reorganization
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_ligand     :  Compute remineralization/dissolution of organic ligands
   !!   p4z_ligand_init:  Initialisation of parameters for remineralisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE prtctl          ! print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_ligand         ! called in p4zbio.F90
   PUBLIC   p4z_ligand_init    ! called in trcsms_pisces.F90

   REAL(wp), PUBLIC ::  rlgw     !: lifetime (years) of weak ligands
   REAL(wp), PUBLIC ::  rlgs     !: lifetime (years) of strong ligands
   REAL(wp), PUBLIC ::  rlig     !: Remin ligand production
   REAL(wp), PUBLIC ::  prlgw    !: Photochemical of weak ligand
   REAL(wp), PUBLIC ::  xklig    !: 1/2 saturation constant of photolysis

   REAL(wp) ::  xklig2    !: xklig * xklig
   LOGICAL  ::  l_dia_ligand

   !! * Substitutions
#  include "ocean2pisces.h90"   
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zligand.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_ligand( kt, knt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_ligand  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic ligands
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      INTEGER, INTENT(in)  ::  Kbb, Krhs ! time level indices
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zlig, zlig2, zlig3
      REAL(wp) :: zlgwp, zlgwpr, zlgwr, zlablgw 
      REAL(wp) :: zfecoag, zaggliga, zligco
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      CHARACTER (len=25) ::   charout
      !!---------------------------------------------------------------------
      !
      IF( kt == nittrc000 )  &
          &  l_dia_ligand = iom_use( "LIGREM" ) .OR. iom_use( "LIGPR" ) &
          &            .OR. iom_use( "LPRODR" ) .OR. iom_use( "LGWCOLL" )

      IF( ln_timing )   CALL timing_start('p4z_ligand')
      !
      ! production from remineralisation of organic matter
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         zlig  = tr(ji,jj,jk,jplgw,Kbb)
         zlig2 = zlig  * zlig
         zlig3 = zlig2 * zlig
         zfecoag = xfecolagg(ji,jj,jk) * 1.0E-9
         zligco  = MAX(0., zlig -  zfecoag)
         !
         ! production from remineralisation of organic matter
         zlgwp = orem(ji,jj,jk) * rlig
         ! Decay of weak ligand
         ! This is based on the idea that as LGW is lower
         ! there is a larger fraction of refractory OM
         !
         zlgwr = ( 1.0 / rlgs * zligco + 1.0 / rlgw * zfecoag ) / ( rtrn + zlig )
         zlgwr = zlgwr * tgfunc(ji,jj,jk) * ( xstep / nyear_len(1) ) * blim(ji,jj,jk) * zlig 
         !
         ! photochem loss of weak ligand
         zlgwpr = prlgw * xstep * etot(ji,jj,jk) * zlig3 * (1. - fr_i(ji,jj) ) / ( zlig2 + xklig2 )
         !
         ! Coagulation of ligands due to various processes (Brownian, shear, diff. sedimentation
         ! xcoagfe is computed in p4zfechem
         ! 50% of the ligands are supposed to be in the colloidal size fraction as for FeL
         zaggliga = xcoagfe(ji,jj,jk) * xstep * 0.5 * zligco
         !
         tr(ji,jj,jk,jplgw,Krhs) = tr(ji,jj,jk,jplgw,Krhs) &
                 &                 + zlgwp  - zlgwr - zlgwpr - zaggliga
         !
      END_3D
      
      IF( l_dia_ligand .AND. ( lk_iomput .AND. knt == nrdttrc ) ) THEN
         ALLOCATE( zw3d(GLOBAL_2D_ARRAY,jpk) )  ;  zw3d(:,:,:) = 0._wp     
         !
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jkR) =  orem(ji,jj,jk) * rlig * 1e9 * 1.e+3 * rfact2r * tmask(ji,jj,jk)
         END_3D   
         CALL iom_put( "LPRODR", zw3d )
         !
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zlig = tr(ji,jj,jk,jplgw,Kbb)
            zw3d(ji,jj,jkR)  = ( 1.0 / rlgs * MAX(0., zlig - xfecolagg(ji,jj,jk) * 1.0E-9 )    &
              &                  + 1.0 / rlgw * xfecolagg(ji,jj,jk) * 1.0E-9 ) / ( rtrn + zlig ) &
              &                 * tgfunc(ji,jj,jk) * ( xstep / nyear_len(1) ) * blim(ji,jj,jk) * zlig 
         END_3D   
         CALL iom_put( "LIGREM", zw3d )
         !
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zlig = tr(ji,jj,jk,jplgw,Kbb)
            zw3d(ji,jj,jkR)  = prlgw * xstep * etot(ji,jj,jk) * zlig*zlig*zlig &
              &            * fr_i(ji,jj) / ( zlig*zlig + xklig2 )
         END_3D   
         CALL iom_put( "LIGPR", zw3d )
         !
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zlig = tr(ji,jj,jk,jplgw,Kbb)
            zw3d(ji,jj,jkR)  =  xcoagfe(ji,jj,jk) * xstep &
              &        * 0.5 * MAX(0., zlig - xfecolagg(ji,jj,jk) * 1.0E-9 )
         END_3D   
         CALL iom_put( "LGWCOLL", zw3d )
         !
         DEALLOCATE( zw3d ) 
      ENDIF

      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('ligand1')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
 !        CALL prt_ctl(tab4d_1=tr(:,:,:,:,Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_ligand')
      !
   END SUBROUTINE p4z_ligand


   SUBROUTINE p4z_ligand_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_ligand_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampislig namelist and check the parameters
      !!
      !! ** input   :   Namelist nampislig
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer 
      !
      NAMELIST/nampislig/ rlgw, prlgw, rlgs, rlig, xklig
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_ligand_init : remineralization/scavenging of organic ligands'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF
      READ_NML_REF(numnatp,nampislig)
      READ_NML_CFG(numnatp,nampislig)
      IF(lwm) WRITE ( numonp, nampislig )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : nampislig'
         WRITE(numout,*) '      Lifetime (years) of weak ligands             rlgw  =', rlgw
         WRITE(numout,*) '      Remin ligand production per unit C           rlig  =', rlig
         WRITE(numout,*) '      Photolysis of weak ligand                    prlgw =', prlgw
         WRITE(numout,*) '      Lifetime (years) of strong ligands           rlgs  =', rlgs
         WRITE(numout,*) '      1/2 saturation for photolysis                xklig =', xklig
      ENDIF
      !
      xklig2 = xklig * xklig
      !
   END SUBROUTINE p4z_ligand_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_ligand                    ! Empty routine
   END SUBROUTINE p4z_ligand
#endif

   !!======================================================================
END MODULE p4zligand
