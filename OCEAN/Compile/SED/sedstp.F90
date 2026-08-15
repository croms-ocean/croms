#include "cppdefs.h"


MODULE sedstp
   !!======================================================================
   !!                       ***  MODULE sedstp   ***
   !!   Sediment model : Sediment model time-stepping
   !!======================================================================
#if defined key_sediment
   USE sed      ! sediment global variables
   USE seddta   ! data read
   USE sedchem  ! chemical constant
   USE sedsol   ! Organic reactions and diffusion
   USE sedadv   ! vertical advection
   USE sedsfc   ! sediment surface data
   USE sedrst   ! restart
   USE sedwri   ! outputs
   USE sedini
   USE setavg_sed
!   USE trcdmp_sed
   USE lib_mpp         ! distribued memory computing library
   USE iom

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC sed_stp  ! called by step.F90

   !! * Substitutions
#  include "ocean2pisces.h90"
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"

   !! $Id: sedstp.F90 15450 2021-10-27 14:32:08Z cetlod $
CONTAINS

   SUBROUTINE sed_stp ( kt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sed_stp  ***
      !!
      !! ** Purpose :   Sediment time stepping
      !!                Simulation of pore water chemistry
      !!
      !! ** Action  :
      !!
      !!
      !!   History :
      !!        !  98-08 (E. Maier-Reimer, Christoph Heinze )  Original code
      !!        !  04-10 (N. Emprin, M. Gehlen ) coupled with PISCES
      !!        !  06-04 (C. Ethe)  Re-organization
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! number of iteration
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      INTEGER  :: ilc
      !!----------------------------------------------------------------------
      IF( ln_timing )           CALL timing_start('sed_stp')
        !
#if defined NEMO        
                                CALL sed_rst_opn  ( kt )       ! Open tracer restart file 
      IF( lrst_sed )            CALL sed_rst_cal  ( kt, 'WRITE' )   ! calenda
#endif

      IF (kt /= nitsed000)      CALL sed_dta( kt, Kbb, Kmm )    ! Load  Data for bot. wat. Chem and fluxes

      IF( kt /= nitsed000 )  &
        &  CALL sed_chem( kt )      ! update of chemical constant to account for salinity, temperature changes
           CALL sed_sol( kt )       ! Solute diffusion and reactions 
           CALL sed_adv( kt )       ! advection

      IF (ln_sed_2way) CALL sed_sfc( kt, Kbb )   ! Give back new bottom wat chem to tracer model

#if ! defined NEMO        
#if ! defined XIOS  && defined AVERAGES
      CALL set_avg_sed
      ilc = 1+iic-nit000 ! number of time step since restart
      IF ( iic > nit000 .AND. mod(ilc-1,nwrtsedpis_avg) == 0 .AND. wrtavg(indxTime) ) THEN
         nrecsedpis_avg=nrecsedpis_avg+1
         CALL sed_wri
      ENDIF
#else
      CALL sed_wri
#endif

      ilc = 1+iic-nit000 ! number of time step since restart
      IF( iic > nit000 ) THEN
         IF( MOD( ilc-1, nitrst ) == 0  &
#ifdef EXACT_RESTART
     &                      .OR. MOD(ilc,nitrst) == 0  &
#endif
     &                      )  THEN
#if defined key_sediment      
! need the CPP key to avoid compilation error
            nrecsedrst = nrecsedrst + 1 
#endif            
            CALL sed_rst_wri
        ENDIF
      ENDIF
#endif
      IF( kt == nitsedend )  CLOSE( numsed )

      IF( ln_timing )           CALL timing_stop('sed_stp')

   END SUBROUTINE sed_stp

#endif

END MODULE sedstp
