#include "cppdefs.h"

MODULE par_pisces
   !!======================================================================
   !!                        ***  par_pisces  ***
   !! TOP :   set the PISCES parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   IMPLICIT NONE
   PUBLIC

   !                                                                !!** Floating point **
   INTEGER, PUBLIC, PARAMETER ::   sp = SELECTED_REAL_KIND( 6, 37)   !: single precision (real 4)
   INTEGER, PUBLIC, PARAMETER ::   dp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
# if defined key_single
   INTEGER, PUBLIC, PARAMETER ::   wp = sp                              !: working precision
# else
   INTEGER, PUBLIC, PARAMETER ::   wp = dp                              !: working precision
# endif
   !                                                                !!** Integer **
   INTEGER, PUBLIC, PARAMETER ::   i4 = SELECTED_INT_KIND( 9)        !: single precision (integer 4)
   INTEGER, PUBLIC, PARAMETER ::   i8 = SELECTED_INT_KIND(14)        !: double precision (integer 8)

   !                                                                !!** Integer **
   INTEGER, PUBLIC, PARAMETER ::   lc  = 256                          !: Lenght of Character strings
   INTEGER, PUBLIC, PARAMETER ::   lca = 400                          !: Lenght of Character arrays

#if defined key_pisces
   !!---------------------------------------------------------------------
   !!   'key_pisces'   :                         standard PISCES bio-model
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .TRUE.  !: PISCES flag 

#if defined key_pisces_npzd
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 9      !: number of PISCES passive tracers
#elif defined key_pisces_quota
#   if defined key_ligand
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 40      !: number of PISCES passive tracers
#   else
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 39      !: number of PISCES passive tracers
#   endif
#else
#   if defined key_ligand
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 25      !: number of PISCES passive tracers
#   else
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 24      !: number of PISCES passive tracers
#   endif
#endif
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  = 11      !: additional 2d output ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  = 16      !: additional 3d output ('key_trc_diaadd')

   ! assign an index in trc arrays for each LOBSTER prognostic variables
   !    WARNING: be carefull about the order when reading the restart
        !   !!gm  this warning should be obsolet with IOM
#if defined key_pisces_npzd
   INTEGER, PUBLIC, PARAMETER ::   jpdic =  1    !: dissolved inoganic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jptal =  2    !: total alkalinity 
   INTEGER, PUBLIC, PARAMETER ::   jpoxy =  3    !: oxygen carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppoc =  4    !: small particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphy =  5    !: phytoplancton concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpzoo =  6    !: zooplancton concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdoc =  7   !: dissolved organic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpno3 =  8   !: Nitrates Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpfer =  9   !: Iron Concentration
   INTEGER, PUBLIC ::   jpcal     !: calcite  concentration 
   INTEGER, PUBLIC ::   jppo4     !: phosphate concentration 
   INTEGER, PUBLIC ::   jpsil     !: silicate concentration
   INTEGER, PUBLIC ::   jpdia     !: Diatoms Concentration
   INTEGER, PUBLIC ::   jpmes     !: Mesozooplankton Concentration
   INTEGER, PUBLIC ::   jpgsi     !: (big) Silicate Concentration
   INTEGER, PUBLIC ::   jpbfe     !: Big iron particles Concentration
   INTEGER, PUBLIC ::   jpgoc     !: big particulate organic phosphate concentration
   INTEGER, PUBLIC ::   jpsfe     !: Small iron particles Concentration
   INTEGER, PUBLIC ::   jpdfe     !: Diatoms iron Concentration
   INTEGER, PUBLIC ::   jpdsi     !: Diatoms Silicate Concentration
   INTEGER, PUBLIC ::   jpnfe     !: Nano iron Concentration
   INTEGER, PUBLIC ::   jpnch     !: Nano Chlorophyll Concentration
   INTEGER, PUBLIC ::   jpdch     !: Diatoms Chlorophyll Concentration
   INTEGER, PUBLIC ::   jpnh4     !: Ammonium Concentration
   INTEGER, PUBLIC ::   jplgw    !: Ammonium Concentration
   INTEGER, PUBLIC ::   jpdon    !: DON concentration 
   INTEGER, PUBLIC ::   jpdop    !: DOP concentration 
   INTEGER, PUBLIC ::   jppon    !: PON concentration
   INTEGER, PUBLIC ::   jppop    !: POP concentration
   INTEGER, PUBLIC ::   jpnph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpndi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppdi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppic     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpnpi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpppi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppfe     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppch     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpgon    !: GON concentration
   INTEGER, PUBLIC ::   jpgop    !: GOP concentration
#   else
   INTEGER, PUBLIC, PARAMETER ::   jpdic =  1    !: dissolved inoganic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jptal =  2    !: total alkalinity 
   INTEGER, PUBLIC, PARAMETER ::   jpoxy =  3    !: oxygen carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpcal =  4    !: calcite  concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppo4 =  5    !: phosphate concentration 
   INTEGER, PUBLIC, PARAMETER ::   jppoc =  6    !: small particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsil =  7    !: silicate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpphy =  8    !: phytoplancton concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpzoo =  9    !: zooplancton concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdoc = 10    !: dissolved organic carbon concentration 
   INTEGER, PUBLIC, PARAMETER ::   jpdia = 11    !: Diatoms Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpmes = 12    !: Mesozooplankton Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdsi = 13    !: Diatoms Silicate Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpfer = 14    !: Iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpbfe = 15    !: Big iron particles Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgoc = 16    !: big particulate organic phosphate concentration
   INTEGER, PUBLIC, PARAMETER ::   jpsfe = 17    !: Small iron particles Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdfe = 18    !: Diatoms iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgsi = 19    !: Sinking Biogenic Silicate Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnfe = 20    !: Nano iron Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnch = 21    !: Nano Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdch = 22    !: Diatoms Chlorophyll Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpno3 = 23    !: Nitrates Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnh4 = 24    !: Ammonium Concentration
#if defined key_ligand
   INTEGER, PUBLIC, PARAMETER ::   jplgw = 25    !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpp4z = 25    !: Number of P4Z tracers
#else
   INTEGER, PUBLIC ::   jplgw    !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpp4z = 24    !: Number of p4z tracers
#endif
#if defined key_pisces_quota
   INTEGER, PUBLIC, PARAMETER ::   jpdon = jpp4z + 1   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpdop = jpp4z + 2   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppon = jpp4z + 3    !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppop = jpp4z + 4   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnph = jpp4z + 5   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppph = jpp4z + 6   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpndi = jpp4z + 7   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppdi = jpp4z + 8   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppic = jpp4z + 9   !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpnpi = jpp4z + 10  !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpppi = jpp4z + 11  !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppfe = jpp4z + 12  !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jppch = jpp4z + 13  !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgon = jpp4z + 14  !: Ammonium Concentration
   INTEGER, PUBLIC, PARAMETER ::   jpgop = jpp4z + 15  !: Ammonium Concentration
#else
   INTEGER, PUBLIC ::   jpdon    !: DON concentration 
   INTEGER, PUBLIC ::   jpdop    !: DOP concentration 
   INTEGER, PUBLIC ::   jppon    !: PON concentration
   INTEGER, PUBLIC ::   jppop    !: POP concentration
   INTEGER, PUBLIC ::   jpnph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpndi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppdi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppic     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpnpi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpppi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppfe     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppch     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpgon    !: GON concentration
   INTEGER, PUBLIC ::   jpgop    !: GOP concentration
#endif
#endif

#else

   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .FALSE.  !: PISCES flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 0      !: number of PISCES passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  = 0      !: additional 2d output ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  = 0      !: additional 3d output ('key_trc_diaadd')
   !
   INTEGER, PUBLIC ::   jpdic     !: dissolved inoganic carbon concentration 
   INTEGER, PUBLIC ::   jptal     !: total alkalinity 
   INTEGER, PUBLIC ::   jpoxy     !: oxygen carbon concentration 
   INTEGER, PUBLIC ::   jppoc     !: small particulate organic phosphate concentration
   INTEGER, PUBLIC ::   jpphy     !: phytoplancton concentration 
   INTEGER, PUBLIC ::   jpzoo     !: zooplancton concentration
   INTEGER, PUBLIC ::   jpdoc    !: dissolved organic carbon concentration 
   INTEGER, PUBLIC ::   jpno3    !: Nitrates Concentration
   INTEGER, PUBLIC ::   jpfer    !: Iron Concentration
   INTEGER, PUBLIC ::   jpcal     !: calcite  concentration 
   INTEGER, PUBLIC ::   jppo4     !: phosphate concentration 
   INTEGER, PUBLIC ::   jpsil     !: silicate concentration
   INTEGER, PUBLIC ::   jpdia     !: Diatoms Concentration
   INTEGER, PUBLIC ::   jpmes     !: Mesozooplankton Concentration
   INTEGER, PUBLIC ::   jpgsi     !: (big) Silicate Concentration
   INTEGER, PUBLIC ::   jpbfe     !: Big iron particles Concentration
   INTEGER, PUBLIC ::   jpgoc     !: big particulate organic phosphate concentration
   INTEGER, PUBLIC ::   jpsfe     !: Small iron particles Concentration
   INTEGER, PUBLIC ::   jpdfe     !: Diatoms iron Concentration
   INTEGER, PUBLIC ::   jpdsi     !: Diatoms Silicate Concentration
   INTEGER, PUBLIC ::   jpnfe     !: Nano iron Concentration
   INTEGER, PUBLIC ::   jpnch     !: Nano Chlorophyll Concentration
   INTEGER, PUBLIC ::   jpdch     !: Diatoms Chlorophyll Concentration
   INTEGER, PUBLIC ::   jpnh4     !: Ammonium Concentration
   INTEGER, PUBLIC ::   jplgw    !: Ammonium Concentration
   INTEGER, PUBLIC ::   jpdon    !: DON concentration 
   INTEGER, PUBLIC ::   jpdop    !: DOP concentration 
   INTEGER, PUBLIC ::   jppon    !: PON concentration
   INTEGER, PUBLIC ::   jppop    !: POP concentration
   INTEGER, PUBLIC ::   jpnph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpndi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppdi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppic     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpnpi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpppi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppfe     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppch     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpgon    !: GON concentration
   INTEGER, PUBLIC ::   jpgop    !: GOP concentration
#endif

   INTEGER, PUBLIC ::   jp_flxco2  
   INTEGER, PUBLIC ::   jp_flxo2   
   INTEGER, PUBLIC ::   jp_kgco2   
   INTEGER, PUBLIC ::   jp_dpco2   
   INTEGER, PUBLIC ::   jp_sinkco2 
   INTEGER, PUBLIC ::   jp_sinkfer 
   INTEGER, PUBLIC ::   jp_sinksil 
   INTEGER, PUBLIC ::   jp_sinkcal 
   INTEGER, PUBLIC ::   jp_heup    
   INTEGER, PUBLIC ::   jp_sildep   
   INTEGER, PUBLIC ::   jp_po4dep   
   INTEGER, PUBLIC ::   jp_no3dep   
   INTEGER, PUBLIC ::   jp_nh4dep   
   INTEGER, PUBLIC ::   jp_nfix 

   INTEGER, PUBLIC ::   jp_hi      
   INTEGER, PUBLIC ::   jp_co3     
   INTEGER, PUBLIC ::   jp_co3sat  
   INTEGER, PUBLIC ::   jp_etot    
   INTEGER, PUBLIC ::   jp_pphy    
   INTEGER, PUBLIC ::   jp_pphy2   
   INTEGER, PUBLIC ::   jp_pnew    
   INTEGER, PUBLIC ::   jp_pnew2   
   INTEGER, PUBLIC ::   jp_pbsi    
   INTEGER, PUBLIC ::   jp_pfed    
   INTEGER, PUBLIC ::   jp_pfen    
   INTEGER, PUBLIC ::   jp_pnewo2  
   INTEGER, PUBLIC ::   jp_prego2  
   INTEGER, PUBLIC ::   jp_grapoc   
   INTEGER, PUBLIC ::   jp_grapoc2   
   INTEGER, PUBLIC ::   jp_mico2  
   INTEGER, PUBLIC ::   jp_meso2  
   INTEGER, PUBLIC ::   jp_nitrifo2 
   INTEGER, PUBLIC ::   jp_remino2 
   INTEGER, PUBLIC ::   jp_nfixo2  
   INTEGER, PUBLIC ::   jp_irondep  
   INTEGER, PUBLIC ::   jp_ironsed  

   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jptra       = jp_pisces                  !: First index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0     = 1                  !: First index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1     = jp_pisces          !: Last  index of PISCES tracers

   REAL(wp), PUBLIC ::  mMass_C  = 12.00      ! Molar mass of carbon
   REAL(wp), PUBLIC ::  mMass_N  = 14.00      ! Molar mass of nitrogen
   REAL(wp), PUBLIC ::  mMass_P  = 31.00      ! Molar mass of phosphorus
   REAL(wp), PUBLIC ::  mMass_Fe = 55.85      ! Molar mass of iron
   REAL(wp), PUBLIC ::  mMass_Si = 28.00      ! Molar mass of silver

   !!======================================================================
END MODULE par_pisces
