#include "cppdefs.h"

MODULE seddta
   !!======================================================================
   !!                     ***  MODULE  seddta  ***
   !! Sediment data  :  read sediment input data from a file
   !!=====================================================================
#if defined key_sediment
   !! * Modules used
   USE trc
   USE sed
   USE sedini
   USE iom
   USE lib_mpp         ! distribued memory computing library

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC sed_dta   ! 

   !! *  Module variables
   REAL(wp) ::  conv2    ! [kg/m2/month]-->[g/cm2/s] ( 1 month has 30 days )

   !! * Substitutions
      !!* Substitution
#  include "ocean2pisces.h90"
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"

CONTAINS

   !!---------------------------------------------------------------------------
   !!   sed_dta  : read the NetCDF data file in online version using module iom
   !!---------------------------------------------------------------------------

   SUBROUTINE sed_dta( kt, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE sed_dta  ***
      !!                    
      !! ** Purpose :   Reads data from a netcdf file and 
      !!                initialization of rain and pore water (k=1) components
      !! 
      !!
      !!   History :
      !!        !  04-10  (N. Emprin, M. Gehlen )  Original code
      !!        !  06-04  (C. Ethe)  Re-organization ; Use of iom
      !!----------------------------------------------------------------------

      !! Arguments
      INTEGER, INTENT( in ) ::   kt    ! time-step
      INTEGER, INTENT( in ) ::   Kbb, Kmm ! time level indices

      !! * Local declarations
      INTEGER  ::  ji, jj, js, jw, ikt

      REAL(wp), DIMENSION(jpoce) :: zdtap, zdtag
      REAL(wp), DIMENSION(A2D(0)) :: zwsbio4, zwsbio3, zddust
      REAL(wp) :: zfact, zdep
      REAL(wp) :: zzf0, zzf1, zzf2, zzf3, zzf4, zzf5, zzf6

      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zw2d

      !----------------------------------------------------------------------

      ! Initialization of sediment variable 
      ! Spatial dimension is merged, and unity converted if needed
      !-------------------------------------------------------------

      IF( ln_timing )  CALL timing_start('sed_dta')

      IF (lwp) THEN
         WRITE(numsed,*)
         WRITE(numsed,*) ' sed_dta : Bottom layer fields'
         WRITE(numsed,*) ' ~~~~~~'
         WRITE(numsed,*) ' Data from SMS model'
         WRITE(numsed,*)
      ENDIF


      ! open file
      IF( kt == nitsed000 ) THEN
         IF (lwp) WRITE(numsed,*) ' sed_dta : Sediment fields'
         conv2 = 1.0e+3 /  1.0e+4 
      ENDIF

      ! reading variables
      IF (lwp) WRITE(numsed,*)
      IF (lwp) WRITE(numsed,*) ' sed_dta : Bottom layer fields at time  kt = ', kt
      ! reading variables
      !
      !    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      !    -----------------------------------------------------------
      DO_2D( 0, 0, 0, 0 )
         ikt = mbkt(ji,jj)
         zdep = e3t(ji,jj,ikt,Kmm) / rDt_trc
         zwsbio4(ji,jj) = MIN( 0.99 * zdep, wsbio4(ji,jj,ikt) / rday  * 1E3 * 1E-4 )
         zwsbio3(ji,jj) = MIN( 0.99 * zdep, wsbio3(ji,jj,ikt) / rday  * 1E3 * 1E-4 )
      END_2D
      !
      trc_data(:,:,:) = 0.
      DO_2D( 0, 0, 0, 0 )
         ikt = mbkt(ji,jj)
         IF ( sedmask(ji,jj) == 1.0 ) THEN
            trc_data(ji,jj,jwsil) = tr(ji,jj,ikt,jpsil,Kbb)
            trc_data(ji,jj,jwoxy) = tr(ji,jj,ikt,jpoxy,Kbb)
            trc_data(ji,jj,jwdic) = tr(ji,jj,ikt,jpdic,Kbb)
            trc_data(ji,jj,jwno3) = tr(ji,jj,ikt,jpno3,Kbb) * redNo3 / redC
            trc_data(ji,jj,jwpo4) = tr(ji,jj,ikt,jppo4,Kbb) * redPo4 / redC
            trc_data(ji,jj,jwalk) = tr(ji,jj,ikt,jptal,Kbb) 
            trc_data(ji,jj,jwnh4) = tr(ji,jj,ikt,jpnh4,Kbb) * redNo3 / redC 
            trc_data(ji,jj,jwh2s) = 0.0
            trc_data(ji,jj,jwso4) = 0.14 * ts(ji,jj,ikt,jp_sal,Kmm) / 1.80655 / 96.062
            trc_data(ji,jj,jwfe2) = tr(ji,jj,ikt,jpfer,Kbb)
            trc_data(ji,jj,jwlgw) = 1.E-9
            trc_data(ji,jj,12)    = MIN(tr(ji,jj,ikt,jpgsi,Kbb), 1E-4) * zwsbio4(ji,jj)
            trc_data(ji,jj,13)    = MIN(tr(ji,jj,ikt,jppoc,Kbb), 1E-4) * zwsbio3(ji,jj) &
            &                       + MIN(tr(ji,jj,ikt,jpgoc,Kbb), 1E-4) * zwsbio4(ji,jj)
            trc_data(ji,jj,14)    = MIN(tr(ji,jj,ikt,jpcal,Kbb), 1E-4) * zwsbio4(ji,jj)
            trc_data(ji,jj,15)    = ts(ji,jj,ikt,jp_tem,Kmm)
            trc_data(ji,jj,16)    = ts(ji,jj,ikt,jp_sal,Kmm)
            trc_data(ji,jj,17)    = ( tr(ji,jj,ikt,jpsfe,Kbb) * zwsbio3(ji,jj) &
            &                        + tr(ji,jj,ikt,jpbfe,Kbb) * zwsbio4(ji,jj) ) &
            &                        / ( trc_data(ji,jj,13) + rtrn )
            trc_data(ji,jj,17)    = MIN(1E-3, trc_data(ji,jj,17) )
         ENDIF
      END_2D

      ! Pore water initial concentration [mol/l] in  k=1
      !-------------------------------------------------
      DO jw = 1, jpwat
          pwcp_dta(:,jw) = PACK( trc_data(:,:,jw), sedmask == 1.0 )
      END DO

      !  Solid components : 
      !-----------------------
      !  Sinking fluxes for OPAL in mol.m-2.s-1 ; conversion in mol.cm-2.s-1
      rainrg(:,jsopal) = PACK( trc_data(:,:,12), sedmask == 1.0 ) 

      !  Sinking fluxes for POC in mol.m-2.s-1 ; conversion in mol.cm-2.s-1
      zdtap = PACK( trc_data(:,:,13), sedmask == 1.0 )
      zzf1 = 0.5187
      zzf2 = 0.3525
      zzf3 = 0.0963
      zzf4 = 0.0243
      zzf5 = 0.0061
      zzf6 = 0.0022
      DO ji = 1, jpoce
         rainrg(ji,jspoc1) = zzf1 * zdtap(ji)
         rainrg(ji,jspoc2) = zzf2 * zdtap(ji)
         rainrg(ji,jspoc3) = zzf3 * zdtap(ji)
         rainrg(ji,jspoc4) = zzf4 * zdtap(ji)
         rainrg(ji,jspoc5) = zzf5 * zdtap(ji)
         rainrg(ji,jspoc6) = zzf6 * zdtap(ji)
      END DO

      !  Sinking fluxes for Calcite in mol.m-2.s-1 ; conversion in mol.cm-2.s-1
      rainrg(:,jscal) = PACK( trc_data(:,:,14), sedmask == 1.0 )

      ! vector temperature [°C] and salinity 
      temp(:) = PACK( trc_data(:,:,15), sedmask == 1.0 )
      salt(:) = PACK( trc_data(:,:,16), sedmask == 1.0 )
      
      ! Clay rain rate in [mol/(cm**2.s)] 
      ! inputs data in [kg.m-2.sec-1] ---> 1e+3/(1e+4) [g.cm-2.s-1]   
      ! divided after by molecular weight g.mol-1      
      rainrg(:,jsclay) = PACK( dust(:,:), sedmask == 1.0 )
      
      zfact = por1(jpksed) * dens_sol(jsclay) / ryear
      DO ji = 1, jpoce
         rainrg(ji,jsclay) = rainrg(ji,jsclay) * conv2 + wacc(ji) * zfact 
         rainrg(ji,jsfeo)  = rainrg(ji,jsclay) / mol_wgt(jsfeo) * 0.035 * 0.5
         rainrg(ji,jsclay) = rainrg(ji,jsclay) / mol_wgt(jsclay) * ( 1.0 - 0.035 * 0.5 ) 
      END DO

      ! Iron monosulphide rain rates. Set to 0
      rainrg(1:jpoce,jsfes)  = 0. 

      ! Fe/C ratio in sinking particles that fall to the sediments
      fecratio(:) = PACK( trc_data(:,:,17), sedmask == 1.0 )

      ! sediment pore water at 1st layer (k=1)
      DO js = 1, jpwat
         pwcp(1:jpoce,1,js) = pwcp_dta(1:jpoce,js)
      END DO

      ! Calculation of raintg of each sol. comp.: rainrg in [g/(cm**2.s)]
      ! computation of dzdep = total thickness of solid material rained [cm] in each cell
      dzdep(:) = 0.
      DO js = 1, jpsol
         zfact = dtsed / ( dens_sol(js) * por1(2) )
         DO ji = 1, jpoce
            rainrg(ji,js) = rainrg(ji,js) * mol_wgt(js)
            dzdep(ji) = dzdep(ji) + rainrg(ji,js) * zfact
         END DO
      END DO

      IF( lk_iomput ) THEN
          ALLOCATE( zw2d(GLOBAL_2D_ARRAY) )  ;  zw2d(:,:) = 0._wp
          IF( iom_use("sflxclay" ) ) THEN 
              zddust(:,:) = UNPACK( wacc(:), sedmask == 1.0, 0.0 )
              zddust(:,:) = dust(:,:) + zddust(:,:) / ( rday * 365.0 ) * por1(jpksed) * dens_sol(jsclay) / conv2
              zw2d(A2D(0)) = zddust(A2D(0)) * 1E3 / 1.E4
              CALL iom_put( "sflxclay", zw2d )
          ENDIF
          IF( iom_use("sflxcal" ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),14) 
              CALL iom_put( "sflxcal", zw2d )
          ENDIF
          IF( iom_use("sflxbsi" ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),12) 
              CALL iom_put( "sflxbsi", zw2d )
          ENDIF
          IF( iom_use("sflxpoc" ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),13) 
              CALL iom_put( "sflxpoc", zw2d )
          ENDIF
          IF( iom_use("OceNO3"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),jwno3) 
              CALL iom_put( "OceNO3", zw2d )
          ENDIF
          IF( iom_use("OceDIC"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),jwdic) 
              CALL iom_put( "OceDIC", zw2d )
          ENDIF
          IF( iom_use("OceALK"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),jwalk) 
              CALL iom_put( "OceALK", zw2d )
          ENDIF
          IF( iom_use("OceOXY"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),jwoxy) 
              CALL iom_put( "OceOXY", zw2d )
          ENDIF
          IF( iom_use("OceSIL"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),jwsil) 
              CALL iom_put( "OceSIL", zw2d )
          ENDIF
          IF( iom_use("OcePO4"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),jwpo4) 
              CALL iom_put( "OcePO4", zw2d )
          ENDIF
          IF( iom_use("OceNH4"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),jwnh4) 
              CALL iom_put( "OceNH4", zw2d )
          ENDIF
          IF( iom_use("OceSO4"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),jwso4) 
              CALL iom_put( "OceSO4", zw2d )
          ENDIF
          IF( iom_use("OceFE2"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),jwfe2) 
              CALL iom_put( "OceFE2", zw2d )
          ENDIF
          IF( iom_use("OceTEM"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),15) 
              CALL iom_put( "OceTEM", zw2d )
          ENDIF
          IF( iom_use("OceSAL"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),16) 
              CALL iom_put( "OceSAL", zw2d )
          ENDIF
          IF( iom_use("OceSFE"  ) )  THEN
              zw2d(A2D(0)) = trc_data(A2D(0),17) 
              CALL iom_put( "OceSFE", zw2d )
          ENDIF
          DEALLOCATE( zw2d )
      ENDIF

      IF( ln_timing )  CALL timing_stop('sed_dta')
      
   END SUBROUTINE sed_dta

#endif

END MODULE seddta
