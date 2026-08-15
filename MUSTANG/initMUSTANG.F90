!------------------------------------------------------------------------------
MODULE initMUSTANG
!------------------------------------------------------------------------------

#include "cppdefs.h"

#ifdef MUSTANG

!&E============================================================================
!&E                   ***  MODULE  initMUSTANG  ***
!&E
!&E ** Purpose : concerns all subroutines related to module MUSTANG
!&E              initialization
!&E 
!&E ** Description :
!&E     subroutine MUSTANG_init           ! initialize the 
!&E                                         parameters of sediment module
!&E                                       ! and the sediment layers
!&E     subroutine MUSTANG_readnml        ! read namelists
!&E     subroutine MUSTANG_param_log      ! write parameters for information
!&E     subroutine MUSTANG_compatibility  ! check compatibility 
!&E                                         between the various parameters
!&E     subroutine MUSTANG_sedinit        ! initialize the sediment if not 
!&E                                         from file
!&E     subroutine MUSTANG_init_hsed      ! initialize the sediment thickness 
!&E                                         function sediment parameters 
!&E                                         (MUSTANG_mixsed)
!&E     subroutine MUSTANG_init_param     ! initialize the settling velocity 
!&E                                         of sand variables and erosion 
!&E                                         parameters
!&E     subroutine MUSTANG_init_output    ! initialize dimensions of output 
!&E                                         tables
!&E     subroutine MUSTANG_morphoinit      ! initialize 
!&E
!&E============================================================================
    !! * Modules used
    USE comMUSTANG
    USE sed_MUSTANG,  ONLY : MUSTANG_E0sand
    USE comsubstance
    USE module_substance
#ifdef key_MUSTANG_flocmod
    USE flocmod, ONLY : flocmod_alloc, flocmod_init
    USE flocmod, ONLY : f_ws, f_diam, f_vol, f_rho, f_mass
#endif

    IMPLICIT NONE

    !! * Accessibility
    PUBLIC MUSTANG_init

    PRIVATE

    !! * Local declarations
    ! declaration of namelists, variables described in comMUSTANG
    namelist /namsedim_init/ filrepsed, l_unised,                   &
                             date_start_morpho, date_start_dyninsed,          &
                             hseduni, cseduni, ksmiuni, ksmauni,              &
                             sini_sed, tini_sed,                              &
                             l_unised_adjust_hsed, csed_mud_ini,              &
                             l_initsed_vardiss, poro_mud_ini

    namelist /namsedim_layer/ l_dzsminuni, dzsminuni,                         &
                              l_dzsmaxuni, dzsmaxuni,                         &
                              dzsmax_bottom, dzsmin,                          &
                              nlayer_surf_sed,                                &
                              k1HW97, k2HW97,                                 &
                              fusion_para_activlayer

    namelist /namsedim_erosion/ activlayer, frmudcr2, coef_frmudcr1,          &
                                x1toce_mud, x2toce_mud,                       &
                                E0_sand_option, E0_sand_para, n_eros_sand,    &
                                E0_mud, n_eros_mud,                           &
                                ero_option, xexp_ero,                         &
                                E0_sand_Cst,                    &
                                tau_cri_option,                               &
                                tau_cri_mud_option_eroindep,                  &
                                l_peph_suspension, l_xexp_ero_cst,            &
                                l_eroindep_mud, l_eroindep_noncoh,            &
                                E0_mud_para_indep

    namelist /namsedim_bottomstress/ l_z0seduni,                              &
                                     z0seduni, z0sedmud, z0sedbedrock,        &
                                     l_fricwave, fricwav,                     &
                                     l_z0hydro_coupl_init,                    & 
                                     l_z0hydro_coupl,                         &
                                     coef_z0_coupl,                           &
                                     z0_hydro_mud,                            &
                                     z0_hydro_bed

    namelist /namsedim_deposition/ cfreshmud, csedmin, cmudcr, aref_sand,     &
                                   cvolmaxsort, cvolmaxmel, slopefac

    namelist /namsedim_lateral_erosion/ l_erolat, coef_erolat,&
                                        coef_tauskin_lat, l_erolat_wet_cell, &
                                        htncrit_eros 

    namelist /namsedim_consolidation/ l_consolid, xperm1, xperm2, xsigma1,    &
                                      xsigma2, csegreg, csandseg,             &
                                      dt_consolid, subdt_consol

    namelist /namsedim_diffusion/ l_diffused, choice_flxdiss_diffsed,         &
                                  xdifs1, xdifsi1,           &
                                  epdifi, fexcs, dt_diffused

    namelist /namsedim_bioturb/ l_bioturb, l_biodiffs,                        &
                                xbioturbmax_part, xbioturbk_part,             &
                                dbiotu0_part, dbiotum_part,                   &
                                xbioturbmax_diss, xbioturbk_diss,             &
                                dbiotu0_diss, dbiotum_diss,                   &
                                frmud_db_min, frmud_db_max,                   &
                                dt_bioturb, subdt_bioturb

    namelist /namsedim_morpho/ l_morphocoupl, MF, l_MF_dhsed, dt_morpho                                  

    namelist /namsedoutput/ choice_nivsed_out,                                &
                            nk_nivsed_out, ep_nivsed_out, epmax_nivsed_out,   &
                            l_outsed_nb_lay_sed, &
                            l_outsed_hsed, &
                            l_outsed_tauskin, &
                            l_outsed_tauskin_cw, &
                            l_outsed_poro, &
                            l_outsed_dzs, &
                            l_outsed_temp_sed, &
                            l_outsed_salt_sed, &
                            l_outsed_cv_sed, &
                            l_outsed_ws, &
                            l_outsed_toce, &
                            l_outsed_flx_s2w_w2s, &
                            l_outsed_pephm_fcor, &
                            l_outsed_bedload, &
                            l_outsed_fsusp, &
                            l_outsed_frmudsup, &
                            l_outsed_dzs_ksmax, &
                            l_outsed_theoric_active_layer, &
                            l_outsed_ero_details, &
                            l_outsed_z0sed, &
                            l_outsed_z0hydro, &
                            l_outsed_consolidation

    namelist /namdredging/ dredging_location_file, dredging_settings_file,    &
            dredging_out_file, dredging_dumping_layer, &
            dredging_dt, dredging_dt_out

    namelist /namsedim_poro/ poro_option, poro_min,                           &
                             Awooster, Bwooster, Bmax_wu 
#ifdef key_MUSTANG_V2
#ifdef key_MUSTANG_bedload
    namelist /namsedim_bedload/ l_peph_bedload, l_slope_effect_bedload,       &
                                alphabs, alphabn, hmin_bedload, l_fsusp
#endif
#endif

#ifdef key_MUSTANG_flocmod
    namelist /namflocmod/ l_ADS, l_ASH, l_COLLFRAG,                           &
                          f_dp0, f_nf, f_nb_frag, f_alpha, f_beta, f_ater,    &
                          f_ero_frac, f_ero_nbfrag, f_ero_iv, f_mneg_param,   &
                          f_collfragparam, f_dmin_frag, f_cfcst, f_fp, f_fy,  &
                          f_clim
#endif
#if !defined key_noTSdiss_insed
    namelist /namtempsed/ mu_tempsed1, mu_tempsed2, mu_tempsed3,              &
                          epsedmin_tempsed,                                   &
                          epsedmax_tempsed
#endif

    ! logging IO
    INTEGER :: ierrorlog, iwarnlog, iscreenlog


CONTAINS
  
!!=============================================================================
    SUBROUTINE MUSTANG_init(ifirst, ilast, jfirst, jlast,                     &
#if defined MORPHODYN  
            dhsed,                                                            &
#endif
            h0fondin, z0hydro, WATER_CONCENTRATION )
    !&E------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_init  ***
    !&E
    !&E ** Purpose : initialize the sediment parameters and variables 
    !&E
    !&E ** Description : called at the beginning of the simulation
    !&E
    !&E ** Called by : mustang_init_main
    !&E
    !&E------------------------------------------------------------------------
   !! * Modules used
    USE dredging, ONLY : dredging_init_param
    USE coupler_MUSTANG,  ONLY : coupl_conv2MUSTANG
    USE sed_MUSTANG,  ONLY : sed_MUSTANG_comp_z0hydro
#ifdef key_MUSTANG_splitlayersurf
    USE sed_MUSTANG,  ONLY : sed_MUSTANG_split_surflayer
#endif
#ifdef key_MUSTANG_V2
    USE sed_MUSTANG,  ONLY : MUSTANGV2_comp_poro_mixsed
#ifdef key_MUSTANG_bedload
    USE sed_MUSTANG_CROCO,  ONLY : sed_bottom_slope
#endif
#endif

    !! * Arguments
    INTEGER, INTENT(IN)                    :: ifirst, ilast, jfirst, jlast
    REAL(KIND=rsh),INTENT(IN)              :: h0fondin
    REAL(KIND=rsh),DIMENSION(GLOBAL_2D_ARRAY),INTENT(INOUT) :: z0hydro                      
    REAL(KIND=rsh),DIMENSION(GLOBAL_2D_ARRAY,N,3,NT), INTENT(IN) :: WATER_CONCENTRATION  
#if defined MORPHODYN 
    REAL(KIND=rsh),DIMENSION(GLOBAL_2D_ARRAY),INTENT(INOUT)          :: dhsed                       
#endif

   !! * Local declarations
    INTEGER   :: i, j, k, iv, isplit
    CHARACTER(len=lchain) :: filepc
    CHARACTER(len=lchain) :: filepc_user
    INTEGER               :: lstr
    INTEGER               :: lenstr
#ifdef key_MUSTANG_V2
    REAL(KIND=rsh)                  :: mass_tot, poro_kij, crel_mud_ini
    REAL(KIND=rsh),DIMENSION(1:nvpc):: frac_sed
    REAL(KIND=rsh),DIMENSION(1:nvp) :: mass_sed
#else
    REAL(KIND=rsh)                  :: somalp
#endif
#ifdef key_MUSTANG_flocmod
    LOGICAL :: l_0Dcase
#endif

    ierrorlog = stdout
    iwarnlog = stdout
    iscreenlog = stdout
    !! * Executable part
    h0fond = h0fondin

    ! Read default and then user namelist
    lstr = lenstr(sedname_must)
    filepc_user = sedname_must(1:lstr)
    filepc = "MUSTANG_NAMELIST/paraMUSTANG_default.txt"
    CALL MUSTANG_readnml(filepc) ! read default namelist
    CALL MUSTANG_readnml(filepc_user) ! read user namelist (replace default param)

    ! Allocate tables 
    CALL MUSTANG_alloc()

    ! Initialize MUSTANG parameters from namelist
    CALL MUSTANG_init_param()

    ! Floculation module
#ifdef key_MUSTANG_flocmod
    CALL flocmod_alloc(nv_mud)
#ifdef SED_TOY_FLOC_0D
    ! for 0D test case, we need to suppress all settling process
    l_0Dcase = .true.
#else
    l_0Dcase = .false.
#endif
    CALL flocmod_init(l_ADS, l_ASH, l_COLLFRAG,             &
        f_dp0, f_nf, f_nb_frag, f_alpha, f_beta, f_ater,    &
        f_ero_frac, f_ero_nbfrag, f_ero_iv, f_mneg_param,   &
        f_collfragparam, f_dmin_frag, f_cfcst, f_fp, f_fy,  &
        f_clim, diam_sed(imud1:nvpc), ros(imud1:nvpc),      &
        rho0, l_0Dcase, ierrorlog)
#endif

    CALL dredging_init_param(ifirst,ilast,jfirst,jlast)


    ! Information on screen
    CALL MUSTANG_param_log()  
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Recovery of concentrations at the bottom layer
    ! which could be used to initiate concentrations in the sediment
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! here iappel = 0 
    ! only WATER_CONCENTRATION, SALINITY_MOD, TEMPERATURE_MOD are used     
    CALL coupl_conv2MUSTANG(ifirst,ilast,jfirst,jlast,0,WATER_CONCENTRATION )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Evaluation of slope for bedload
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
#if defined MORPHODYN
    it_morphoYes=0
#endif
    CALL sed_bottom_slope(ifirst, ilast, jfirst, jlast, h)
#endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Definition of initial conditions in sediment
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (l_unised) THEN
        CALL MUSTANG_sedinit(ifirst, ilast, jfirst, jlast, h)
    ELSE
        CALL MUSTANG_init_fromfile
    END IF
      
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Porosity estimation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO j = jfirst, jlast
        DO i = ifirst, ilast
          IF (ksma(i,j) .NE. 0) THEN

#ifdef key_MUSTANG_splitlayersurf
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! Splitting surface layers if too thick
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            isplit = 0
            DO k = ksma(i,j), ksma(i,j)-nlayer_surf_sed+1, -1
                IF (k > ksmi(i,j)) THEN
                    IF (dzs(k,i,j) > dzsmax(i,j) + 5.0_rsh* dzsmin) isplit=1
                ENDIF
            ENDDO
            IF(isplit == 1) THEN            
                CALL sed_MUSTANG_split_surflayer(i,j,ksma(i,j))
            ENDIF
#endif

            DO k=ksmi(i,j),ksma(i,j)
#ifdef key_MUSTANG_V2
             mass_tot = 0.0_rsh
             DO iv=igrav1,imud2
               mass_sed(iv)=cv_sed(iv,k,i,j)*dzs(k,i,j) 
               mass_tot=mass_tot+mass_sed(iv)
             END DO
             DO iv=igrav1,imud2
               frac_sed(iv)=mass_sed(iv)/mass_tot
             END DO
             crel_mud_ini=(1.0_rsh-poro_mud_ini)*ros(1)
             CALL MUSTANGV2_comp_poro_mixsed(frac_sed,poro_mud_ini,  &
                                crel_mud_ini,poro_kij)
             poro(k,i,j)=poro_kij
             crel_mud(k,i,j)=crel_mud_ini
             poro_mud(k,i,j)=poro_mud_ini
             c_sedtot(k,i,j)=(1.0_rsh-poro(k,i,j))*ros(1)
             DO iv=igrav1,imud2
               cv_sed(iv,k,i,j)=c_sedtot(k,i,j)*frac_sed(iv)
             END DO
#else
             somalp=0.0_rsh
             DO iv=1,nvpc
                somalp=somalp+cv_sed(iv,k,i,j)/ros(iv)
             ENDDO
             poro(k,i,j)=1.0_rsh-somalp
#endif

            ENDDO
          ELSE
            poro(:,i,j)=0.0_rsh
#ifdef key_MUSTANG_V2
            crel_mud(:,i,j)=0.0_rsh
            poro_mud(:,i,j)=0.0_rsh
#endif 
          ENDIF
        ENDDO
      ENDDO
      
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Estimation of sediment heights and bathy bedrock for morphodynamic
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    ! call even if no morpho, initialization of hsed in all cases, 
    ! even if initfromfile
    CALL MUSTANG_morphoinit(ifirst, ilast, jfirst, jlast   &
#if defined MORPHODYN  
                                                ,dhsed               &
#endif
            )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Evaluation of Z0_hydro if l_z0_coupl_init
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    IF(l_z0hydro_coupl_init) THEN
        CALL sed_MUSTANG_comp_z0hydro(ifirst, ilast, jfirst, jlast, z0hydro)
    ENDIF


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialization of output tables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    CALL MUSTANG_init_output

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check CPP keys and parameters for compatibility
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    CALL MUSTANG_compatibility(ifirst, ilast, jfirst, jlast)

    END SUBROUTINE MUSTANG_init
!!===========================================================================
 
    SUBROUTINE MUSTANG_readnml(filein)
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_readnml  ***
    !&E
    !&E ** Purpose : reads namelist file from filein
    !&E
    !&E ** Description : read namelist file paraMUSTANG
    !&E
    !&E ** Called by :  MUSTANG_init
    !&E
    !&E ** External calls : tool_datosec()
    !&E
    !&E--------------------------------------------------------------------------
    !! * Arguments
    CHARACTER(len=lchain), INTENT(IN) :: filein

    !! * Local declarations
    REAL(KIND=rlg)        :: tool_datosec

    !! * Executable part
    OPEN(unit = 50, file = filein, status = 'old', action = 'read')
 
    MPI_master_only write(*,*) '*****************************************************'
    MPI_master_only write(*,*) 'READING MUSTANG input file'
    MPI_master_only write(*,*) TRIM(filein)
    MPI_master_only write(*,*) '*****************************************************'

    READ(50, namsedim_init); rewind(50)
    READ(50, namsedim_layer); rewind(50)
    READ(50, namsedim_bottomstress); rewind(50)
    READ(50, namsedim_deposition); rewind(50)
    READ(50, namsedim_erosion); rewind(50)
    READ(50, namsedim_poro); rewind(50)
#ifdef key_MUSTANG_V2
#ifdef key_MUSTANG_bedload
    READ(50, namsedim_bedload); rewind(50)
#endif
#endif
    READ(50, namsedim_lateral_erosion); rewind(50)
    READ(50, namsedim_consolidation); rewind(50)
    READ(50, namsedim_diffusion); rewind(50)
    READ(50, namsedim_bioturb); rewind(50)
    READ(50, namsedim_morpho); rewind(50)
#if !defined key_noTSdiss_insed
    READ(50, namtempsed); rewind(50)
#endif
    READ(50, namsedoutput); rewind(50)
#ifdef key_MUSTANG_flocmod
    ! module FLOCULATION
    READ(50, namflocmod); rewind(50)
#endif
    READ(50, namdredging); rewind(50)

    CLOSE(50) 
   
    END SUBROUTINE MUSTANG_readnml
!!===========================================================================

    SUBROUTINE MUSTANG_param_log()
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_param_log  ***
    !&E
    !&E ** Purpose : write parameters in log file
    !&E
    !&E ** Description : namelists writing for information
    !&E
    !&E ** Called by :  MUSTANG_init
    !&E
    !&E--------------------------------------------------------------------------
    !! * Local declarations
    INTEGER   :: iv

    MPI_master_only WRITE(iscreenlog, namsedim_init)
    MPI_master_only WRITE(iscreenlog, namsedim_layer)
    MPI_master_only WRITE(iscreenlog, namsedim_bottomstress)
    MPI_master_only WRITE(iscreenlog, namsedim_deposition)
    MPI_master_only WRITE(iscreenlog, namsedim_poro)
#ifdef key_MUSTANG_V2
#ifdef key_MUSTANG_bedload
    MPI_master_only WRITE(iscreenlog, namsedim_bedload)
#endif  
#else
    MPI_master_only WRITE(iscreenlog, namsedim_erosion)
#endif
    MPI_master_only WRITE(iscreenlog, namsedim_lateral_erosion)
    MPI_master_only WRITE(iscreenlog, namsedim_consolidation)
    MPI_master_only write(iscreenlog, *) 'MUSTANG_param_log dt_consolid=', dt_consolid
    MPI_master_only WRITE(iscreenlog, namsedim_diffusion)
    MPI_master_only WRITE(iscreenlog, namsedim_bioturb)
    MPI_master_only WRITE(iscreenlog, namsedim_morpho)
#if !defined key_noTSdiss_insed
    MPI_master_only WRITE(iscreenlog, namtempsed)
#endif

#ifdef key_MUSTANG_flocmod
    !! module floculation
    MPI_master_only WRITE(iscreenlog, *) ' '
    MPI_master_only WRITE(iscreenlog, *) '    FLOCMOD'
    MPI_master_only WRITE(iscreenlog, *) '***********************'
    MPI_master_only WRITE(iscreenlog, *) 'class  diameter  volume  density  mass Ws'
    DO iv = 1, nv_mud
        MPI_master_only WRITE(iscreenlog, *) iv, f_diam(iv), f_vol(iv), f_rho(iv), f_mass(iv), f_ws(iv)
    ENDDO
    MPI_master_only WRITE(iscreenlog, *) ' '
    MPI_master_only WRITE(iscreenlog, *) ' *** PARAMETERS ***'
    MPI_master_only WRITE(iscreenlog, *) &
        'Primary particle size (f_dp0)                                : ', f_dp0
    MPI_master_only WRITE(iscreenlog, *) &
        'Fractal dimension (f_nf)                                     : ', f_nf
    MPI_master_only WRITE(iscreenlog, *) &
        'Flocculation efficiency (f_alpha)                            : ', f_alpha
    MPI_master_only WRITE(iscreenlog, *) &
        'Floc break up parameter (f_beta)                             : ', f_beta
    MPI_master_only WRITE(iscreenlog, *) &
        'Nb of fragments (f_nb_frag)                                  : ', f_nb_frag
    MPI_master_only WRITE(iscreenlog, *) &
        'Ternary fragmentation (f_ater)                               : ', f_ater
    MPI_master_only WRITE(iscreenlog, *) &
        'Floc erosion (% of mass) (f_ero_frac)                        : ', f_ero_frac
    MPI_master_only WRITE(iscreenlog, *) &
        'Nb of fragments by erosion (f_ero_nbfrag)                    : ', f_ero_nbfrag
    MPI_master_only WRITE(iscreenlog, *) &
        'fragment class (f_ero_iv)                                    : ', f_ero_iv
    MPI_master_only WRITE(iscreenlog, *) &
        'negative mass tolerated before redistribution (f_mneg_param) : ', f_mneg_param
    MPI_master_only WRITE(iscreenlog, *) &
        'Boolean for differential settling aggregation (L_ADS)        : ', l_ADS
    MPI_master_only WRITE(iscreenlog, *) &
        'Boolean for shear aggregation (L_ASH)                        : ', l_ASH
    MPI_master_only WRITE(iscreenlog, *) &
        'Boolean for collision fragmenation (L_COLLFRAG)              : ', l_COLLFRAG
    MPI_master_only WRITE(iscreenlog, *) &
        'Collision fragmentation parameter (f_collfragparam)          : ', f_collfragparam
    MPI_master_only WRITE(iscreenlog, *) &
        'Min concentration below which flocculation is not calculated : ', f_clim
    MPI_master_only WRITE(iscreenlog, *) ' '
    MPI_master_only WRITE(iscreenlog, *) '*** END FLOCMOD INIT *** '    
#endif
   
    END SUBROUTINE MUSTANG_param_log
!!===========================================================================
 
    SUBROUTINE MUSTANG_compatibility(ifirst, ilast, jfirst, jlast)
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_compatibility  ***
    !&E
    !&E ** Purpose : verifies compatibility between the various parameters (sedim)
    !&E
    !&E ** Description :
    !&E
    !&E ** Called by : MUSTANG_init
    !&E
    !&E ** External calls : tool_datosec()
    !&E
    !&E--------------------------------------------------------------------------
    !! * Arguments
    INTEGER, INTENT(IN) :: ifirst, ilast, jfirst, jlast
   
    !! * Local declarations
    INTEGER   :: iv, iv2, isubs, i, j 
    REAL(rsh) :: rohomog
 
    ! test dissolved variables in sediment must be used without  key_nofluxwat_IWS
#if defined key_nofluxwat_IWS && ! defined key_noTSdiss_insed
    MPI_master_only WRITE(iwarnlog, *)
    MPI_master_only WRITE(iwarnlog, *) '****************************'
    MPI_master_only WRITE(iwarnlog, *) '        WARNING           ' 
    MPI_master_only WRITE(iwarnlog, *) '****************************'
    MPI_master_only WRITE(iwarnlog, *) &
        ' You are using CPP key : key_nofluxwat_IWS ' 
    MPI_master_only WRITE(iwarnlog, *) &
        ' So you are not taking into account water fluxes threw water/sediment interface ' 
    MPI_master_only WRITE(iwarnlog, *) &
        ' which are generated by erosion/ settling/ consolidation..'
    MPI_master_only WRITE(iwarnlog, *) &
        ' AND you are NOT using CPP key : key_noTSdiss_insed ' 
    MPI_master_only WRITE(iwarnlog, *) &
        ' so you want to simulate dissolved variables in sediment'
    MPI_master_only WRITE(iwarnlog, *) &
        ' You will have probably conservativity problems '
    MPI_master_only WRITE(iwarnlog, *) &
        ' If you simulate dissolved variables in water AND in sediment, don t use key_nofluxwat_IWS '
#endif

#ifdef key_MUSTANG_flocmod
    IF (.not.l_ADS .and. .not.l_ASH) THEN
        MPI_master_only write(ierrorlog, *) 'CAUTION : incompatible flocculation kernel options : '
        MPI_master_only write(ierrorlog, *) '*****************************************************'
        MPI_master_only write(ierrorlog, *) 'l_ADS=', l_ADS
        MPI_master_only write(ierrorlog, *) 'l_ASH=', l_ASH
        MPI_master_only write(ierrorlog, *) 'simulation stopped'
        STOP
    ENDIF
#endif

    ! test compatibility tocd /consolidation
    IF (l_consolid) THEN
        MPI_master_only WRITE(iwarnlog, *)
        MPI_master_only WRITE(iwarnlog, *) '****************************'
        MPI_master_only WRITE(iwarnlog, *) ' You are taking into account CONSOLIDATION PROCESS'
        MPI_master_only WRITE(iwarnlog, *) ' and fresh deposit concentration (cfreshmud) = ', cfreshmud
        MPI_master_only WRITE(iwarnlog, *) &
            ' and critical stress of deposition (tocd) of each sand and mud variable are equal to :'
        DO isubs = imud1, imud2
            MPI_master_only WRITE(iwarnlog, *) &
                ' variable', name_var(isubs),'  tocd=', tocd(isubs)
        ENDDO
        MPI_master_only WRITE(iwarnlog, *) &
            ' If cfreshmud is big (> 100), you have to choose a tocd between 1 and 10'
        MPI_master_only WRITE(iwarnlog, *) &
            ' otherwise you have to choose a tocd > 10 or 20'
       
    ELSE
       DO isubs = isand2+1, nvp
         IF (tocd(isubs) > 1.5 ) THEN
           MPI_master_only WRITE(iwarnlog, *) '****************************'
           MPI_master_only WRITE(iwarnlog, *) '  WARNING    '
           MPI_master_only WRITE(iwarnlog, *) ' You are not taking into account CONSOLIDATION PROCESS'
           MPI_master_only WRITE(iwarnlog, *) ' fresh deposit concentration (cfreshmud) = ', cfreshmud
           MPI_master_only WRITE(iwarnlog, *) ' and critical stress of deposition (tocd) of isusb=', isubs, &
                                              ' exceed 1.5 Pa :', tocd(isubs)
           MPI_master_only WRITE(iwarnlog, *) ' it should be smaller (0.5 to 1 ?) '
         ENDIF
       ENDDO
    ENDIF
    
    ! test ros homogen for sand in all cases
    IF (l_consolid) THEN
       rohomog = ros(irk_fil(1))
       DO iv = 2, nvpc
        IF(ros(irk_fil(iv)) .NE. rohomog) THEN
            MPI_master_only WRITE(ierrorlog, *)
            MPI_master_only WRITE(ierrorlog, *) '****************************'
            MPI_master_only WRITE(ierrorlog, *) ' You are taking into account CONSOLIDATION PROCESS'
            MPI_master_only WRITE(ierrorlog, *) ' and you must have same density for all constitutive particles'
            MPI_master_only WRITE(ierrorlog, *) ' It is not the case : (see in substance file) '
            MPI_master_only WRITE(ierrorlog, *) ' ros(1) = ', ros(irk_fil(1))
            MPI_master_only WRITE(ierrorlog, *) ' ros(', irk_fil(iv), ') = ', ros(irk_fil(iv))
            STOP
         ENDIF
       ENDDO
     ELSE
       IF(isand2 > isand1) THEN
        rohomog = ros(irk_fil(isand1))
        DO iv = isand1+1, isand2
         IF(ros(irk_fil(iv)) .NE. rohomog) THEN
            MPI_master_only WRITE(ierrorlog, *)
            MPI_master_only WRITE(ierrorlog, *) '****************************'
            MPI_master_only WRITE(ierrorlog, *) ' you must have same density for all sand particles'
            MPI_master_only WRITE(ierrorlog, *) ' It is not the case : (see in substance file) '
            MPI_master_only WRITE(ierrorlog, *) ' ros(1) = ', ros(irk_fil(isand1))
            MPI_master_only WRITE(ierrorlog, *) ' ros(', irk_fil(iv), ') = ', ros(irk_fil(iv))
            STOP
          ENDIF
          DO iv2 = isand1, iv
           IF(diam_sed(irk_fil(iv)) > diam_sed(irk_fil(iv2))) THEN
             MPI_master_only WRITE(ierrorlog, *)
             MPI_master_only WRITE(ierrorlog, *) '****************************'
             MPI_master_only WRITE(ierrorlog, *) ' The "SAND" variables should be stored in substance file '
             MPI_master_only WRITE(ierrorlog, *) ' in order of decreasing diameters (the coarsest sand to the finest sand)'
             MPI_master_only WRITE(ierrorlog, *) ' It is not the case : (see in substance file) '
             STOP
           ENDIF
          ENDDO
        ENDDO
       ENDIF
     ENDIF

! test compatibility keys sedimento+biolo
#if defined key_biolo && defined key_benthos && ! defined key_Pconstitonly_insed
    MPI_master_only WRITE(ierrorlog, *)
    MPI_master_only WRITE(ierrorlog, *) '****************************'
    MPI_master_only WRITE(ierrorlog, *) ' You are using key_MUSTANG + key_biolo + key_benthos'
    MPI_master_only WRITE(ierrorlog, *) ' It is not compatible'
    MPI_master_only WRITE(ierrorlog, *) ' You must add the key key_Pconstitonly_insed and biological variables '
    MPI_master_only WRITE(ierrorlog, *) ' will not be included in sedim module but in benthos module -in biolo)'
    MPI_master_only WRITE(ierrorlog, *)
    MPI_master_only WRITE(ierrorlog, *) ' OR you want simulate biological variables in sediment module'
    MPI_master_only WRITE(ierrorlog, *) ' and you have to use : key_MUSTANG + key_biolo without key_benthos'
    STOP
#endif

! test compatibility key_BLOOM_insed
#if defined key_BLOOM_insed && ! defined key_oxygen
    MPI_master_only WRITE(ierrorlog, *)
    MPI_master_only WRITE(ierrorlog, *) '****************************'
    MPI_master_only WRITE(ierrorlog, *) ' It is not possible to use key_BLOOM_insed  without the key_oxygen'
    MPI_master_only WRITE(ierrorlog, *) ' WARNING : use CPP key : key_oxygen in Makefile'
    STOP
#endif

   ! test compatibility bedload with MUSTANG Version 2 
#if defined key_MUSTANG_bedload && ! defined key_MUSTANG_V2
    MPI_master_only WRITE(ierrorlog, *)
    MPI_master_only WRITE(ierrorlog, *) '****************************'
    MPI_master_only WRITE(ierrorlog, *) ' You are using key_MUSTANG_bedload but '
    MPI_master_only WRITE(ierrorlog, *) ' You must add the key key_MUSTANG_V2 '
    MPI_master_only WRITE(ierrorlog, *) ' in order to use the new version V2 of MUSTANG which include bedload behavior'
    MPI_master_only WRITE(ierrorlog, *)
    STOP
#endif

#ifdef key_MUSTANG_V2
! test version V2 consolidation/poro_option
    IF (l_consolid .AND. poro_option .NE. 2) THEN
        MPI_master_only WRITE(iwarnlog, *) '****************************'
        MPI_master_only WRITE(iwarnlog, *) ' WARNING    '
        MPI_master_only WRITE(iwarnlog, *) ' You are  taking into account CONSOLIDATION PROCESS'
        MPI_master_only WRITE(iwarnlog, *) ' and have choose poro_optio n= 1'
        MPI_master_only WRITE(iwarnlog, *) ' poro_option must be = 2. it has been changed automatically'
        MPI_master_only WRITE(iwarnlog, *)
        poro_option = 2
    ENDIF
#endif
   
    END SUBROUTINE MUSTANG_compatibility
!!===========================================================================
 
    SUBROUTINE MUSTANG_init_param()
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_init_param  ***
    !&E
    !&E ** Purpose : initialize MUSTANG parameters  
    !&E
    !&E ** Description : 
    !&E    - compute for each sand and gravel (for 1 to isand2)  
    !&E       * diamstar
    !&E       * ws_sand
    !&E       * tetacri0
    !&E       * stresscri0
    !&E       * psi_sed (in V2)
    !&E    - compute for each sand (for isand1 to isand2) 
    !&E       * E0_sand (in V2)
    !&E    - compute for each gravel, dans and mud (for 1 to imud2) 
    !&E       * rosmrowsros
    !&E    - ros_sand_homogen = ros(isand1) if there is sand (isand1 > 0)
    !&E
    !&E ** Note : g and rho0 must be known 
    !&E
    !&E ** Called by :  MUSTANG_init
    !&E
    !&E--------------------------------------------------------------------------

    !! * Local declarations
    INTEGER        :: iv
    REAL(KIND=rlg) :: dtsedc, dtsedd, dtsedb  ! to compute consolidation subdt_consol
    REAL(KIND=rlg) :: tool_datosec
    REAL(KIND=rsh),PARAMETER :: shield_cri_wu = 0.03_rsh          

    !!--------------------------------------------------------------------------
    !! * Executable part

    IF (csed_mud_ini == 0.0_rsh ) csed_mud_ini = cfreshmud
    tstart_dyninsed = max(tool_datosec(date_start_dyninsed), time)
    t_dyninsed = tstart_dyninsed
    tstart_morpho = tool_datosec(date_start_morpho)
    dtsedc = 365._rlg * 86400._rlg
    MPI_master_only write(*,*) 'init dtsedc = ', dtsedc
    dtsedd = dtsedc
    dtsedb = dtsedc
    IF (l_consolid) THEN
        dtsedc = dt_consolid
        MPI_master_only write(*,*) 'init dt_consolid = ',dt_consolid
    ELSE
        dt_consolid = dtsedc
        MPI_master_only write(*,*) 'init dt_consolid = ',dt_consolid
    ENDIF
    IF (l_diffused) THEN
        dtsedd = dt_diffused
    ELSE
        dt_diffused = dtsedd
    ENDIF
    IF (l_bioturb .OR. l_biodiffs) THEN
        dtsedb = dt_bioturb
    ELSE
        dt_bioturb = dtsedb
    ENDIF
    dt_dyninsed = MIN(dtsedc,dtsedd,dtsedb)
    dtsedc = 365._rlg*86400._rlg
    dtsedb = dtsedc
    IF (l_consolid) THEN
        dtsedc = subdt_consol
    ELSE
        subdt_consol = dtsedc
    ENDIF
    IF (l_bioturb) THEN
        dtsedb = subdt_bioturb
    ELSE
        subdt_bioturb = dtsedb
    ENDIF    
    subdt_consol = MIN(dtsedb, dtsedc, dt_dyninsed)
#ifdef key_MUSTANG_V2    
    IF (l_dzsminuni) THEN 
        coeff_dzsmin = 0._rsh        
    ELSE           
        coeff_dzsmin = 1._rsh  
    ENDIF   
#endif  
    IF (l_morphocoupl) THEN
    t_morpho = tstart_morpho
        IF (l_MF_dhsed) THEN
            MF_dhsed = MF
            MF = 1._rsh
        ELSE
            MF_dhsed = 0._rsh  ! unused
            ! MF used, read from namelist
        ENDIF
    ELSE
        MF = 1.0_rsh
    ENDIF
    cexcs = 1.0_rsh - fexcs
    fws2 = fricwav * 0.5_rsh
    l_dyn_insed = .FALSE.
    IF (l_consolid .OR. l_bioturb .OR. l_diffused .OR. l_biodiffs) l_dyn_insed = .TRUE.



    ! conversion of thickness reading in mm to m for computation
    epmax_nivsed_out = epmax_nivsed_out / 1000.0_rsh
    ep_nivsed_out(:) = ep_nivsed_out(:) / 1000.0_rsh
    IF (choice_nivsed_out == 1) nk_nivsed_out = ksdmax
                               
    MPI_master_only write(*,*) 'MUSTANG_init_param dt_consolid=', dt_consolid

    DO iv = 1, isand2

      ! here, for simplification, water density = rho0 
      ! here, for simplification, viscosity = 1.e-6 constant
      ! 10000 = (1/viscosite**2)**(1/3)
      rosmrowsros(iv) = (ros(iv) - rho0) / ros(iv)
      diamstar(iv) = diam_sed(iv) * 10000.0_rsh * (g * (ros(iv) / rho0 - 1.0_rsh))**0.33_rsh
      ! according to Soulsby, 1997, and if viscosity = 10-6 m/s :
      ws_sand(iv)=.000001_rsh*((107.33_rsh+1.049_rsh*diamstar(iv)**3)**0.5_rsh-10.36_rsh)/diam_sed(iv)


       ! Critical shear stress in erosion law (N/m2)
      IF (tau_cri_option == 0) THEN
        tetacri0(iv)=0.3_rsh/(1.0_rsh+1.2_rsh*diamstar(iv))+0.055_rsh*(1.0_rsh-EXP(-0.02_rsh*diamstar(iv)))
        stresscri0(iv)=tetacri0(iv)*g*(ros(iv)-rho0)*diam_sed(iv)
      ELSE IF (tau_cri_option == 1) THEN ! Wu and Lin (2014)
        stresscri0(iv)=g*(ros(iv)-rho0)*diam_sed(iv)*shield_cri_wu
      END IF

#ifdef key_MUSTANG_V2
      psi_sed(iv)=LOG(diam_sed(iv)*1000.0_rsh)/LOG(2.0_rsh)
      IF (iv .ge. isand1) THEN
            E0_sand(iv) = MUSTANG_E0sand(diam_sed(iv), stresscri0(iv), ros(iv), ws_sand(iv))
      ENDIF
#endif
    ENDDO ! iv = 1, isand2

    DO iv=imud1,imud2
      rosmrowsros(iv)=(ros(iv)-rho0)/ros(iv)
    ENDDO

    ! homogeneous sand density
    IF (isand1 > 0) THEN
       iv=isand1
       ros_sand_homogen=ros(iv)
    ENDIF

    END SUBROUTINE MUSTANG_init_param
!!===========================================================================

    SUBROUTINE MUSTANG_init_hsed(cv_sedini)
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_init_hsed  ***
    !&E
    !&E ** Purpose : initialize sediment thickness, in order to be coherent with sediment parameters
    !&E     function of cseduni, cvolmax, csed_ini of each sediment,
    !&E      valid only if cv_sed  and dzs not variable
    !&E
    !&E ** Description : 
    !&E
    !&E ** Called by :  MUSTANG_init
    !&E
    !&E--------------------------------------------------------------------------
    !! * Arguments
    REAL(KIND=rsh), INTENT(IN),DIMENSION(nvpc) :: cv_sedini

    !! * Local declarations
    INTEGER        :: iv                
    REAL(KIND=rsh) :: dzsgrv,dzssan,dzsmud,voldepgrv,voldepsan,masdepmud,     &
                        dvolsan,dmasmud,dmasmudgrav,dmasmudsand

    !! * Executable part

    ! **TODO** translate all comments in english
      !!! Volume a deposer
    voldepgrv=0.0_rsh
    voldepsan=0.0_rsh        
    DO iv=igrav1,igrav2
        voldepgrv=voldepgrv+cv_sedini(iv)/ros(iv)*hseduni  ! kg/m3 * m /(kg/m3)  (vol_g)
    ENDDO
    DO iv=isand1,isand2
        voldepsan=voldepsan+cv_sedini(iv)/ros(iv)*hseduni  ! kg/m3 * m /(kg/m3)  (vol_s)
    ENDDO

    !!! Masse a deposer
    masdepmud=0.0_rsh
    DO iv=imud1,imud2
        masdepmud=masdepmud+cv_sedini(iv)*hseduni  ! kg/m2  (m_v)
    ENDDO

    !!! nouvelle Couche de gravier
    IF(voldepgrv.GE.0.0_rsh)THEN
        dzsgrv=voldepgrv/cvolmaxsort
    ENDIF

    !!! Couche de sable
    IF(voldepsan.GE.0.0_rsh)THEN
        dvolsan=MIN(voldepsan,dzsgrv*(cvolmaxmel-cvolmaxsort)) ! volume de sable qui s insere dans le gravier (cvolmaxmel-cvolmaxsort)
        dzssan=MAX(0.,(voldepsan-dvolsan)/cvolmaxsort)
    ENDIF
    
    !!! Couche de vase
    IF(masdepmud.GE.0.0_rsh) THEN
        IF(dzsgrv >0.0_rsh .and. dzssan > 0.0_rsh) THEN
            ! melange sable+gravier  
            dmasmudgrav=dzsgrv*csed_mud_ini*(1.0_rsh-cvolmaxmel)  ! masse de vase qui s insere dans le gravier (1-cvolmaxmel)
            IF (isand2>isand1 ) THEN
                dmasmudsand=dzssan*csed_mud_ini*(1.0_rsh-cvolmaxmel) ! masse de vase qui s insere dans le sable (1-cvolmaxsort)
            ELSE
                dmasmudsand=dzssan*csed_mud_ini*(1.0_rsh-cvolmaxsort) ! masse de vase qui s insere dans le sable (1-cvolmaxsort)
            ENDIF
        ELSE IF (dzssan > 0.0_rsh .AND. isand2>isand1 ) THEN
            ! melange sables  
            dmasmudgrav=0.
            dmasmudsand=dzssan*csed_mud_ini*(1.0_rsh-cvolmaxmel) ! masse de vase qui s insere dans le sable (1-cvolmaxsort)
        ELSE IF(dzsgrv >0.0_rsh .AND. igrav2>igrav1 ) THEN
            ! melange graviers  
            dmasmudgrav=dzsgrv*csed_mud_ini*(1.0_rsh-cvolmaxmel)  ! masse de vase qui s insere dans le gravier (1-cvolmaxmel)
            dmasmudsand=0.
        ELSE
            ! un seul type de sediment 
            dmasmudgrav=dzsgrv*csed_mud_ini*(1.0_rsh-cvolmaxsort)  ! masse de vase qui s insere dans le gravier ou sable (1-cvolmaxsort)
            dmasmudsand=dzssan*csed_mud_ini*(1.0_rsh-cvolmaxsort) ! masse de vase qui s insere dans le sable (1-cvolmaxsort)
        ENDIF
        ! nouvelle couche de vase
        dzsmud=(masdepmud-dmasmudgrav-dmasmudsand)/csed_mud_ini
        IF(dzsmud <0.0_rsh) THEN
            dzsmud=0.0_rsh
            IF (masdepmud-dmasmudgrav < 0.0_rsh) THEN
                dmasmudgrav = masdepmud
                dmasmudsand = 0.0_rsh
            ELSE
                dmasmudsand = masdepmud-dmasmudgrav
            ENDIF
        ENDIF     
    dmasmud = dmasmudgrav +dmasmudsand
    ENDIF

    !!! Nouvelle Hauteur totale
    hsed_new = dzsgrv +dzssan +dzsmud

    END SUBROUTINE MUSTANG_init_hsed

!!=============================================================================
  SUBROUTINE MUSTANG_init_fromfile
 
    !&E-------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_init_fromfile  ***
    !&E
    !&E ** Purpose : manages the fields to be re-initialized in case of the 
    !&E              continuation of a previous run
    !&E
    !&E ** Description : open and read a netcdf file, written during a previous run 
    !&E
    !&E ** Called by :  MUSTANG_init 
    !&E
    !&E-------------------------------------------------------------------------
    !! * Modules used
     USE dredging, ONLY : l_dredging, dredg_hsed_init
     
     implicit none
                       
 
# include "netcdf.inc"
     real time_scale
     integer iv, k, itrc, indWrk
     integer ncid, indx, varid,  ierr, lstr, lvar, latt, lenstr,       &
     start(2), count(2), ibuff(6), nf_fread, checkdims
     character units*180, nomcv*30
     character inised_name*180
     real tmp(GLOBAL_2D_ARRAY)
     real tmp3d(GLOBAL_2D_ARRAY, ksdmin:ksdmax)
 
     tmp(GLOBAL_2D_ARRAY) = 0
     tmp3d(GLOBAL_2D_ARRAY,ksdmin:ksdmax) = 0
     z0sed(GLOBAL_2D_ARRAY) = z0seduni
     dzsmax(GLOBAL_2D_ARRAY) = dzsmaxuni
     ksmi(GLOBAL_2D_ARRAY) = ksmiuni
     ksma(GLOBAL_2D_ARRAY) = 0
     hsed(GLOBAL_2D_ARRAY) = 0.0_rsh
     dzs(ksdmin:ksdmax,GLOBAL_2D_ARRAY) = 0.0_rsh
     cv_sed(-1:nv_tot,ksdmin:ksdmax,GLOBAL_2D_ARRAY) = 0.0_rsh
     c_sedtot(ksdmin:ksdmax,GLOBAL_2D_ARRAY) = 0.0_rsh
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
 ! used, provided that nrrec is within the available records; 
 ! (3) if the file has multiple records and nrrec<0, then THE LAST 
 ! available record is used.
 
 !      if (may_day_flag.ne.0) return      !-->  EXIT
       inised_name = filrepsed
       lstr = lenstr(inised_name)
       ierr = nf_open(inised_name(1:lstr), nf_nowrite, ncid)
       if (ierr.eq.nf_noerr) then
         ierr = checkdims(ncid, inised_name, lstr, indx)
 
         if (ierr.ne.nf_noerr) then
          goto 99
         elseif (indx.eq.0) then
           indx = 1
         elseif (indx.gt.0 .and. nrrec.gt.0 .and. nrrec.le.indx) then
           indx = nrrec
         elseif (indx.gt.0 .and. nrrec.gt.indx) then
           write(stdout,'(/1x,A,I4,A/16x,A,I4,A/16x,3A/)')   &
                'MUSTANG_init_fromfile ERROR: requested restart time record', &
                nrrec, ' exceeds',  'number of available records', &
                indx,'in netCDF file', '''',inised_name(1:lstr),'''.'
           goto 99                                        !--> ERROR
         endif
       else
         write(stdout,'(/1x,2A/15x,3A)') &
                'MUSTANG_init_fromfile ERROR: Cannot ', &
                'open netCDF file', '''', inised_name(1:lstr) ,'''.'
         goto 99                                           !--> ERROR
       endif
 
 
 ! ksmi, ksma
       ierr=nf_inq_varid (ncid,'ksmi', varid)
       if (ierr .eq. nf_noerr) then
         ierr=nf_fread (tmp, ncid, varid, indx, 0)
         ksmi(GLOBAL_2D_ARRAY)=INT(tmp(GLOBAL_2D_ARRAY))
         if (ierr .ne. nf_noerr) then
           MPI_master_only write(stdout,2) 'ksmi', indx, inised_name(1:lstr)
           goto 99                                         !--> ERROR
         endif
       else
         MPI_master_only  write(stdout,1) 'ksmi', inised_name(1:lstr)
         goto 99                                           !--> ERROR
       endif
       
        ierr=nf_inq_varid (ncid,'ksma', varid)
       if (ierr .eq. nf_noerr) then
         ierr=nf_fread (tmp, ncid, varid, indx, 0)
         ksma(GLOBAL_2D_ARRAY)=INT(tmp(GLOBAL_2D_ARRAY))
 
         if (ierr .ne. nf_noerr) then
           MPI_master_only write(stdout,2) 'ksma', indx, inised_name(1:lstr)
           goto 99                                         !--> ERROR
         endif
       else
         MPI_master_only  write(stdout,1) 'ksma', inised_name(1:lstr)
         goto 99                                           !--> ERROR
       endif
 
      WHERE (ksmi(GLOBAL_2D_ARRAY) < ksdmin ) 
         ksmi(GLOBAL_2D_ARRAY) = 1
         ksma(GLOBAL_2D_ARRAY) = 0
      END WHERE
 !  DZS
       ierr=nf_inq_varid (ncid,'DZS', varid)
       if (ierr .eq. nf_noerr) then
         ierr=nf_fread (tmp3d, ncid, varid, indx, 12)
          do k=ksdmin,ksdmax
             dzs(k,:,:)=tmp3d(:,:,k)
          enddo
         
         if (ierr .ne. nf_noerr) then
           MPI_master_only write(stdout,2) 'DZS', indx, inised_name(1:lstr)
           goto 99                                         !--> ERROR
         endif
       else
         MPI_master_only  write(stdout,1) 'DZS', inised_name(1:lstr)
         goto 99                                           !--> ERROR
       endif
 
 
 ! CVSED
  
       c_sedtot(:,:,:)=0.0_rsh
       do iv=-1,nv_tot 
 
        if (iv == -1) then
           nomcv='temp_sed'
        elseif (iv==0) then
           nomcv='salt_sed'
        else
           nomcv=TRIM(name_var(iv))//'_sed'
        endif
 
        ierr=nf_inq_varid (ncid,nomcv, varid)
        if (ierr .eq. nf_noerr) then
         ierr=nf_fread (tmp3d, ncid, varid, indx, 12)
 
         do k=ksdmin,ksdmax
             cv_sed(iv,k,:,:)=tmp3d(:,:,k)
             c_sedtot(k,:,:)=c_sedtot(k,:,:)+cv_sed(iv,k,:,:)*typart(iv)  
         enddo
 
         if (ierr .ne. nf_noerr) then
           MPI_master_only write(stdout,2) nomcv, indx, inised_name(1:lstr)
           goto 99                                         !--> ERROR
         endif
        
        else
         if (iv > nvpc .and. iv < nvp+1 ) then
          ! not constitutive particulate  variables (Initial in M/Msed converted to M/m3 sed)
          IF (irkm_var_assoc(iv) >0) THEN
            do k=ksdmin,ksdmax
              cv_sed(iv,k,:,:)=cini_sed(iv)*cv_sed(irkm_var_assoc(iv),k,:,:)
            enddo
          ELSE
            do k=ksdmin,ksdmax
             cv_sed(iv,k,:,:)=cini_sed(iv)*c_sedtot(k,:,:)
            enddo
          END IF
          MPI_master_only  write(stdout,3) nomcv, inised_name(1:lstr)
         else if (iv > nvp) then
           ! dissolved variables (M/m3 EI)
           do k=ksdmin,ksdmax
             cv_sed(iv,k,:,:)=cini_sed(iv)             
           enddo
          MPI_master_only  write(stdout,3) nomcv,name_var(iv), &
                                           inised_name(1:lstr)
         else
          MPI_master_only  write(stdout,1) nomcv, inised_name(1:lstr)
          goto 99                          !--> ERROR
         endif
        endif
       enddo
      
      do k=ksdmin,ksdmax
      WHERE (c_sedtot(k,GLOBAL_2D_ARRAY) == -valmanq) c_sedtot(k,GLOBAL_2D_ARRAY)=0.0_rsh
      enddo
 
 
 ! dredging
      if (l_dredging) then
         ierr=nf_inq_varid (ncid,'dredg_hsed_init', varid)
         if (ierr .eq. nf_noerr) then
           ierr=nf_fread (tmp, ncid, varid, indx, 0)
           dredg_hsed_init(GLOBAL_2D_ARRAY)=tmp(GLOBAL_2D_ARRAY)
           if (ierr .ne. nf_noerr) then
             MPI_master_only write(stdout,2) 'dredg_hsed_init', indx, inised_name(1:lstr)
             goto 99                                         !--> ERROR
           endif
         else
           MPI_master_only  write(stdout,1) 'dredg_hsed_init', inised_name(1:lstr)
           goto 99                                           !--> ERROR
         endif
      endif
 
 ! Close input NetCDF file.
 !
       ierr=nf_close(ncid)
 
   1   format(/1x,'MUSTANG_init_fromfile - unable to find variable:', 1x,A, &
                                  /15x,'in input NetCDF file:',1x,A/)
   2   format(/1x,'MUSTANG_init_fromfile - error while reading variable:',1x, A, &
          2x,'at time record =',i4/15x,'in input NetCDF file:',1x,A/)
   3   format(/1x,'MUSTANG_init_fromfile - unable to find variable:', 1x,A,/10x,A, &
      &                            /15x,'in input NetCDF file:',1x,A,  &
          1x,'-> analytical value (cini_sed)'/)
       return
   99  may_day_flag=2
       return
       
   END SUBROUTINE MUSTANG_init_fromfile

!!===========================================================================================
    
    SUBROUTINE MUSTANG_sedinit(ifirst, ilast, jfirst, jlast, h)
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_sedinit  ***
    !&E
    !&E ** Purpose : initialize the sediment geometry & content
    !&E
    !&E ** Description : 
    !&E
    !&E ** Called by :  MUSTANG_init_sediment only if l_unised
    !&E
    !&E--------------------------------------------------------------------------
    USE dredging, ONLY : l_dredging, dredg_hsed_init
         
    implicit none

    !! * Arguments
    INTEGER,INTENT(IN)                                        :: ifirst, ilast, jfirst, jlast
    REAL(KIND=rsh),DIMENSION(GLOBAL_2D_ARRAY),INTENT(IN)       :: h      

    !! * Local declarations
    CHARACTER*50   :: filesand, filemud
    INTEGER        :: i, j, k, iv, ivg, imud_nobio, IERR_MPI
    REAL(KIND=rsh) :: somalp,dzsf,dzsgrv,dzssan,dzsmud,voldepgrv,voldepsan,masdepmud,h0max_def,     &
                        dvolsan,dmasmud,dmasmudgrav,dmasmudsand,cini_modif,tot_csed_init,dzsmaxmin

    !! * Executable part

    dzsmax(GLOBAL_2D_ARRAY)=dzsmaxuni
    z0sed(GLOBAL_2D_ARRAY)=z0seduni
    ksmi(GLOBAL_2D_ARRAY)=ksmiuni
    ksma(GLOBAL_2D_ARRAY)=0
    hsed(GLOBAL_2D_ARRAY)=-valmanq
    dzs(ksdmin:ksdmax,GLOBAL_2D_ARRAY)=-valmanq
    cv_sed(-1:nv_tot,ksdmin:ksdmax,GLOBAL_2D_ARRAY)=-valmanq
    c_sedtot(ksdmin:ksdmax,GLOBAL_2D_ARRAY)=-valmanq


        ! uniform bed :
        WHERE (h(GLOBAL_2D_ARRAY) /= -valmanq) ksmi(GLOBAL_2D_ARRAY)=ksmiuni
        WHERE (h(GLOBAL_2D_ARRAY) /= -valmanq) ksma(GLOBAL_2D_ARRAY)=ksmauni
        DO k=ksdmin,ksdmax
            WHERE (h(GLOBAL_2D_ARRAY) /= -valmanq) dzs(k,GLOBAL_2D_ARRAY)=0.0_rsh
            WHERE (h(GLOBAL_2D_ARRAY) /= -valmanq) c_sedtot(k,GLOBAL_2D_ARRAY)=0.0_rsh
            DO iv=-1,nv_tot
                WHERE (h(GLOBAL_2D_ARRAY) /= -valmanq) cv_sed(iv,k,GLOBAL_2D_ARRAY)=0.0_rsh
            ENDDO
        ENDDO
 
 ! if biological variable (detrital organic matter) are constitutive particulaite variable
 !     readjustment of cini_sed (fraction of csedtot : somme=1)
        cini_modif=0.0_rsh
        imud_nobio=1
        DO iv=1,isand2
            cini_modif=cini_modif+cini_sed(iv)
        ENDDO
        DO iv=imud1,imud2
            IF(unit_modif_mudbio_N2dw(irk_fil(iv)) .NE. 1.0_rsh) THEN
                ! variable bio (Ndet) mis en constitutive mais exprime dans la donnee initiale en mmole/kg de sediment total
                !cini_sed(iv)=cini_sed(iv)/cseduni
                cini_modif=cini_modif+cini_sed(iv)
            ELSE 
                ! vase non detritique
                IF(cini_sed(iv) > (1.0_rsh-cini_modif)) imud_nobio=iv
                cini_modif=cini_modif+cini_sed(iv)
            ENDIF    
        ENDDO
        IF(cini_modif > 1.0_rsh) cini_sed(imud_nobio)=cini_sed(imud_nobio)-cini_modif+1.0_rsh
   
        tot_csed_init=0.0_rsh
        DO iv=1,nvpc
            ! constitutive variables : cini_sed in fraction of cseduni
            tot_csed_init=tot_csed_init+cini_sed(iv)              
        ENDDO

        IF( abs(tot_csed_init - 1.0_rsh) > epsilon_MUSTANG ) THEN
            MPI_master_only WRITE(ierrorlog,*)' '
            MPI_master_only WRITE(ierrorlog,*)' inital concentrations of constitutive particulate matter'
            MPI_master_only WRITE(ierrorlog,*)' have not been well defined'
            MPI_master_only WRITE(ierrorlog,*)' See variable.dat to give fraction of total concentration '
            MPI_master_only WRITE(ierrorlog,*)' for each constitutive variable (GRAV, SAND, MUD) (but not MUDB) '
            MPI_master_only WRITE(ierrorlog,*)' : the sum of all must be =1 and it is :'
            MPI_master_only WRITE(ierrorlog,*)  tot_csed_init
            MPI_master_only WRITE(ierrorlog,*)' difference with 1 is :'
            MPI_master_only WRITE(ierrorlog,*)  abs(tot_csed_init - 1.0_rsh) 
            STOP   
        ELSE
            DO iv=1,nvpc
                ! force fraction sum to equal 1.
                cini_sed(iv) = cini_sed(iv)/tot_csed_init            
            ENDDO    
        ENDIF
      
        ! ajustement du hsed et cv_sed
        ! ==============================
        IF(l_unised_adjust_hsed .AND. hseduni > 0.0_rsh) THEN
            ALLOCATE(cv_sedini(1:nvpc))
            !!! Concentration massique
            DO iv=1,nvpc
                ! constitutive variables : cini_sed in % of cseduni
                cv_sedini(iv)=cseduni*cini_sed(iv)
            ENDDO

            CALL MUSTANG_init_hsed(cv_sedini)
            MPI_master_only WRITE(iscreenlog,*)'NEW TOTAL SEDIMENT THICKNESS AFTER ADJUSTEMENT = ',hsed_new
            cseduni=cseduni*hseduni/hsed_new         
            hseduni=hsed_new
            WHERE (h(GLOBAL_2D_ARRAY) /= -valmanq) hsed(GLOBAL_2D_ARRAY)=hseduni
            DEALLOCATE(cv_sedini)
        ELSE
            WHERE (h(GLOBAL_2D_ARRAY) /= -valmanq) hsed(GLOBAL_2D_ARRAY)=hseduni
        ENDIF

        DO j=jfirst,jlast
        DO i=ifirst,ilast
            IF(h(i,j) /= -valmanq) THEN
                DO k=ksmi(i,j),ksma(i,j)
                    c_sedtot(k,i,j)=0.0_rsh
                    DO iv=1,nvpc
                        ! constitutive variables : cini_sed in % of cseduni
                        cv_sed(iv,k,i,j)=cseduni*cini_sed(iv)        
                        c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
                    ENDDO
                    ! porosity estimation moved after dzs evaluation below
#ifdef key_Pconstitonly_insed
                    cv_sed(nvpc+1:nv_state,k,i,j)=0.0_rsh 
#else
                    DO iv=nvpc+1,nvp
                        ! non constitutive variables : cini_sed in conc/kg of sediment
                        IF (irkm_var_assoc(iv) >0) THEN
                            cv_sed(iv,k,i,j)=cini_sed(iv)*cv_sed(irkm_var_assoc(iv),k,i,j)
                        ELSE
                            cv_sed(iv,k,i,j)=cini_sed(iv)*c_sedtot(k,i,j)
                        END IF
                    ENDDO
                    IF(l_initsed_vardiss) THEN
                    ! initalisation avec cv_wat (CW) du fond
                        DO iv=nvp+1,nv_adv
                            ! dissolved variables : cini_sed in conc/m3 of pore water
                            cv_sed(iv,k,i,j)=cw_bottom_MUSTANG(iv,i,j)
                        ENDDO
                    ELSE
                        ! initalisation avec cini (variable.dat)
                        DO iv=nvp+1,nv_adv
                            ! dissolved variables : cini_sed in conc/m3 of pore water
                            cv_sed(iv,k,i,j)=cini_sed(iv)
                        ENDDO
                    ENDIF
                    DO iv=nv_adv+1,nv_state
                        ! + fixed variables : cini_sed in conc/m3 of sediment
                        cv_sed(iv,k,i,j)=cini_sed(iv)
                    ENDDO
#endif
                    IF(l_initsed_vardiss) THEN
                    ! initalisation avec valeurs du fond
                        cv_sed(-1,k,i,j)=temp_bottom_MUSTANG(i,j)
                        cv_sed( 0,k,i,j)=sal_bottom_MUSTANG(i,j)
                    ELSE
                        cv_sed(-1,k,i,j)=tini_sed
                        cv_sed( 0,k,i,j)=sini_sed
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        ENDDO


!   regular vertical mesh
!   ---------------------

    DO j=jfirst,jlast
        DO i=ifirst,ilast
            IF(ksma(i,j).NE.0)THEN
                dzsf=hsed(i,j)/float(ksma(i,j)-ksmi(i,j)+1)
                DO k=ksmi(i,j),ksma(i,j)
                    dzs(k,i,j)=dzsf
                ENDDO
            ELSE
                dzs(:,i,j) = 0.0_rsh
                hsed(i,j) = 0.0_rsh
            ENDIF
        ENDDO
    ENDDO

    IF (l_dredging) THEN
        dredg_hsed_init(GLOBAL_2D_ARRAY) = hsed(GLOBAL_2D_ARRAY) 
    ENDIF

! DZSmax non uniform
    IF(.NOT. l_dzsmaxuni) THEN
        IF(dzsmaxuni == 0.0_rsh) THEN
            MPI_master_only write(ierrorlog,*)' ERROR in paraMUSTANG : dzsmaxuni =0. and l_unised=.TRUE.'
            MPI_master_only write(ierrorlog,*)' if dzsmaxuni =0. dzsmax must be read in an initial file '
            MPI_master_only write(ierrorlog,*)'  dzsmax must be evaluated from dzsmaxuni not null'
            STOP
        ENDIF

        ! definition of dzsmax
        ! estimated as a function of depth and dzsmaxuni considered as the max value of dzsmax
        ! here minimum value of dzsmax in depth = dzsmaxuni / 100
        ! max depth giving the lowest value of dzsma to inform, setting by default = 50m
       dzsmaxmin=dzsmaxuni/100._rsh
       h0max_def=50.0_rsh
       MPI_master_only write(iscreenlog,*)
       MPI_master_only write(iscreenlog,*)'DZSMAX NON UNIFORM'
       MPI_master_only write(iscreenlog,*) &
       'by default linearly computed in MUSTANG_sedinit '
       MPI_master_only write(iscreenlog,*) &
       'from dzsmaxuni (',dzsmaxuni,'m) to dzsmaxuni/100 (',dzsmaxmin,'m) depending on water depth'
       MPI_master_only write(iscreenlog,*) &
       'dzsmax minimum for depth > ',h0max_def
       MPI_master_only write(iscreenlog,*) &
       'PROGRAM your own initialization of dzsmax (or READ a file) in MUSTANG_sedinit'
       dzsmax(GLOBAL_2D_ARRAY)=(h0max_def-h(GLOBAL_2D_ARRAY))/h0max_def*dzsmaxuni+dzsmaxmin
       WHERE (dzsmax(GLOBAL_2D_ARRAY) > dzsmaxuni ) dzsmax(GLOBAL_2D_ARRAY)=dzsmaxuni
       WHERE (dzsmax(GLOBAL_2D_ARRAY) < dzsmaxmin ) dzsmax(GLOBAL_2D_ARRAY)=dzsmaxmin

       ! or to be read in a data file (to be programmed) 
        
    ENDIF 



    END SUBROUTINE MUSTANG_sedinit

!!===========================================================================

    SUBROUTINE MUSTANG_check_output_param
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_check_output_param  ***
    !&E
    !&E ** Purpose : check input options and input parameters compatibilities
    !&E--------------------------------------------------------------------------
        !! * Local declarations
        INTEGER        :: k, err, nk_nivsed_outlu
        REAL(KIND=rsh) :: dzs_estim

        nk_nivsed_outlu = nk_nivsed_out
        IF(choice_nivsed_out == 1 ) THEN
            MPI_master_only WRITE(iscreenlog,*)'results in sediment will be save on all the sediment layers (ksdmax)'
            nk_nivsed_out = ksdmax
        ELSE IF(choice_nivsed_out == 2 ) THEN
            IF(nk_nivsed_outlu > ksdmax) THEN
                nk_nivsed_out = ksdmax
            ELSE IF (nk_nivsed_outlu < 1 ) THEN
                MPI_master_only WRITE(iscreenlog,*) &
        'Wrong nk_nivsed_out value, should be between 1 and',ksdmax
                STOP
            ENDIF 
            MPI_master_only WRITE(iscreenlog,*) &
        'results in sediment will be save only for the first ', nk_nivsed_out, &
                                    ' layers from sediment surface, one average for the last'
        ELSE IF(choice_nivsed_out == 3 ) THEN
            dzs_estim=MIN(dzsmaxuni,hseduni/ksmauni)
            nk_nivsed_out = MIN(ksdmax, INT(epmax_nivsed_out / dzs_estim) +3) 
            MPI_master_only WRITE(iscreenlog,*) &
        'results in sediment will be save from sediment surface till ',  &
                                epmax_nivsed_out*1000,'mmeters (integration below)'
            MPI_master_only WRITE(iscreenlog,*) &
        'the last layer saved in the sediment output file will be that  ',  &
                                'the bottom of which exceeds the desired maximum thickness '          
            IF (.NOT. l_dzsmaxuni) THEN
                MPI_master_only WRITE(iscreenlog,*)
                MPI_master_only WRITE(iscreenlog,*) &
        '            if dzsmax is not uniform, it could be smaller than dzsmaxuni and  &
                            then there may be points where '
                MPI_master_only WRITE(iscreenlog,*) &
        '             this number is not sufficient to describe the maximum thickness'
                MPI_master_only WRITE(iscreenlog,*)
                MPI_master_only WRITE(iscreenlog,*) &
        'This option is not recommended with l_dzsmaxuni=.False.'
            ENDIF
    
        ELSE IF(choice_nivsed_out == 4 ) THEN
            IF(nk_nivsed_outlu > min(5,ksdmax-1) ) THEN
                MPI_master_only WRITE(iscreenlog,*)'ERROR in namelist namsedoutput (paraMUSTANG.txt) '
                MPI_master_only WRITE(iscreenlog,*)'nk_nivsed_out must be <=min(5,',ksdmax-1,')'
                MPI_master_only WRITE(iscreenlog,*)'nk_nivsed_out is automatically set to 5 '
                MPI_master_only WRITE(iscreenlog,*) &
        'and an latter layer (6th) is added to integrate till the bottom sediment  '
                nk_nivsed_outlu = min(5,ksdmax-1)
            ELSE IF (nk_nivsed_outlu < 1 ) THEN
                MPI_master_only WRITE(iscreenlog,*) &
        'Wrong nk_nivsed_out value, should be between 1 and min(5,',ksdmax-1,')'
                STOP 1
            ENDIF   
            err = 0
            DO k=1,nk_nivsed_outlu
                IF (ep_nivsed_out(k) == 0.) THEN
                    MPI_master_only WRITE(iscreenlog,*) &
        'Wrong ep_nivsed_out(', k ,') value, should not be zero'
                    err = err +1 
                ENDIF
            ENDDO
            IF (err>0) THEN
                STOP 1
            ENDIF     
            
            nk_nivsed_out =  nk_nivsed_outlu+1   
            ALLOCATE(ep_nivsed_outp1(nk_nivsed_out))
            ep_nivsed_outp1(1:nk_nivsed_outlu)=ep_nivsed_out(1:nk_nivsed_outlu)
            ep_nivsed_outp1(nk_nivsed_out)=10.0_rsh  
    
            MPI_master_only WRITE(iscreenlog,*) &
        'results in sediment will be save on ',nk_nivsed_out-1, &
                'integrated layers whom thickness are constant and given by user - first is  sediment surface'
            MPI_master_only WRITE(iscreenlog,*) &
        'the ',nk_nivsed_out,'layer will be an integrated layer till the bottom' 
    
        ELSE
            MPI_master_only WRITE(iscreenlog,*) &
        'Wrong choice_nivsed_out value, should be 1,2,3 or 4'
            STOP
        ENDIF !choice_nivsed_out

    END SUBROUTINE MUSTANG_check_output_param

!!===========================================================================

    SUBROUTINE MUSTANG_allocate_output_arrays
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_allocate_output_arrays  ***
    !&E
    !&E ** Purpose : prepare needed output arrays
    !&E--------------------------------------------------------------------------

#if defined key_BLOOM_insed
        USE bioloinit,  ONLY : ndiag_tot, ndiag_3d_sed, ndiag_2d_sed, ndiag_1d, ndiag_2d

        ALLOCATE(var2D_diagsed(GLOBAL_2D_ARRAY,ndiag_1d+ndiag_2d-ndiag_2d_sed+1:ndiag_1d+ndiag_2d))
        ALLOCATE(var3D_diagsed(nk_nivsed_out,GLOBAL_2D_ARRAY,ndiag_tot-ndiag_3d_sed+1:ndiag_tot))
#endif

        IF (l_outsed_hsed) THEN
            ALLOCATE(var2D_hsed(GLOBAL_2D_ARRAY))
            var2D_hsed(GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
        IF (l_outsed_dzs) THEN
            ALLOCATE(var3D_dzs(nk_nivsed_out,GLOBAL_2D_ARRAY))
            var3D_dzs(:,GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
        IF (l_outsed_poro) THEN
            ALLOCATE(var3D_poro(nk_nivsed_out,GLOBAL_2D_ARRAY))
            var3D_poro(:,GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
        IF (l_outsed_temp_sed) THEN
            ALLOCATE(var3D_TEMP(nk_nivsed_out,GLOBAL_2D_ARRAY))
            var3D_TEMP(:,GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
        IF (l_outsed_salt_sed) THEN
            ALLOCATE(var3D_SAL(nk_nivsed_out,GLOBAL_2D_ARRAY))
            var3D_SAL(:,GLOBAL_2D_ARRAY)= 0.0_rsh
        ENDIF
        IF (l_outsed_cv_sed) THEN
            ALLOCATE(var3D_cvsed(nk_nivsed_out,GLOBAL_2D_ARRAY,ntrc_subs))
            var3D_cvsed(:,GLOBAL_2D_ARRAY,ntrc_subs) = 0.0_rsh
        ENDIF
        IF (l_outsed_toce) THEN
            ALLOCATE(var2D_toce(nvpc,GLOBAL_2D_ARRAY))
            var2D_toce(:,GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
        IF (l_outsed_flx_s2w_w2s) THEN
            ALLOCATE(var2D_flx_s2w(nvpc,GLOBAL_2D_ARRAY))
            var2D_flx_s2w(:,GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
        IF (l_outsed_flx_s2w_w2s) THEN
            ALLOCATE(var2D_flx_w2s(nvpc,GLOBAL_2D_ARRAY))
            var2D_flx_w2s(:,GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF

        IF (l_outsed_frmudsup) THEN
            ALLOCATE(var2D_frmudsup(GLOBAL_2D_ARRAY))
            var2D_frmudsup(GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
        IF (l_outsed_dzs_ksmax) THEN
            ALLOCATE(var2D_dzs_ksmax(GLOBAL_2D_ARRAY))
            var2D_dzs_ksmax(GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
        IF (l_outsed_flx_s2w_w2s) THEN
            ALLOCATE(var2D_flx_s2w_coh(GLOBAL_2D_ARRAY))
            var2D_flx_s2w_coh(GLOBAL_2D_ARRAY) = 0.0_rsh
            ALLOCATE(var2D_flx_w2s_coh(GLOBAL_2D_ARRAY))
            var2D_flx_w2s_coh(GLOBAL_2D_ARRAY) = 0.0_rsh
            ALLOCATE(var2D_flx_s2w_noncoh(GLOBAL_2D_ARRAY))
            var2D_flx_s2w_noncoh(GLOBAL_2D_ARRAY) = 0.0_rsh
            ALLOCATE(var2D_flx_w2s_noncoh(GLOBAL_2D_ARRAY))
            var2D_flx_w2s_noncoh(GLOBAL_2D_ARRAY)= 0.0_rsh
        ENDIF

#ifdef key_MUSTANG_V2
        IF (l_outsed_pephm_fcor) THEN
            ALLOCATE(var2D_pephm_fcor(nvpc,GLOBAL_2D_ARRAY))
            var2D_pephm_fcor(:,GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
        IF (l_outsed_theoric_active_layer) THEN
            ALLOCATE(var2D_theoric_active_layer(GLOBAL_2D_ARRAY))
            var2D_theoric_active_layer(GLOBAL_2D_ARRAY)= 0.0_rsh
        ENDIF
        IF (l_outsed_ero_details) THEN
            ALLOCATE(var2D_tero_noncoh(GLOBAL_2D_ARRAY))
            var2D_tero_noncoh(GLOBAL_2D_ARRAY)= 0.0_rsh
            ALLOCATE(var2D_tero_coh(GLOBAL_2D_ARRAY))
            var2D_tero_coh(GLOBAL_2D_ARRAY) = 0.0_rsh
            ALLOCATE(var2D_pct_iter_noncoh(GLOBAL_2D_ARRAY))
            var2D_pct_iter_noncoh(GLOBAL_2D_ARRAY) = 0.0_rsh
            ALLOCATE(var2D_pct_iter_coh(GLOBAL_2D_ARRAY))
            var2D_pct_iter_coh(GLOBAL_2D_ARRAY) = 0.0_rsh
            ALLOCATE(var2D_niter_ero(GLOBAL_2D_ARRAY))
            var2D_niter_ero(GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
#ifdef key_MUSTANG_bedload
        IF (l_outsed_bedload) THEN
            ALLOCATE(var2D_flx_bx(nvpc,GLOBAL_2D_ARRAY))
            var2D_flx_bx(:,GLOBAL_2D_ARRAY) = 0.0_rsh
            ALLOCATE(var2D_flx_by(nvpc,GLOBAL_2D_ARRAY))
            var2D_flx_by(:,GLOBAL_2D_ARRAY) = 0.0_rsh
            ALLOCATE(var2D_bil_bedload(nvpc,GLOBAL_2D_ARRAY))
            var2D_bil_bedload(:,GLOBAL_2D_ARRAY) = 0.0_rsh
            ALLOCATE(var2D_flx_bx_int(GLOBAL_2D_ARRAY))
            var2D_flx_bx_int(GLOBAL_2D_ARRAY) = 0.0_rsh
            ALLOCATE(var2D_flx_by_int(GLOBAL_2D_ARRAY))
            var2D_flx_by_int(GLOBAL_2D_ARRAY) = 0.0_rsh
            ALLOCATE(var2D_bil_bedload_int(GLOBAL_2D_ARRAY)) 
            var2D_bil_bedload_int(GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
        IF (l_outsed_fsusp) THEN
            ALLOCATE(var2D_fsusp(nvpc,GLOBAL_2D_ARRAY))
            var2D_fsusp(:,GLOBAL_2D_ARRAY) = 0.0_rsh
        ENDIF
#endif
#endif
        IF (l_dyn_insed) THEN
            IF (l_outsed_consolidation) THEN
                IF (choice_nivsed_out == 1) then
                    ALLOCATE(var3Dksed_loadograv(nk_nivsed_out,GLOBAL_2D_ARRAY)) 
                    var3Dksed_loadograv(:,GLOBAL_2D_ARRAY) = 0.0_rsh
                    ALLOCATE(var3Dksed_permeab(nk_nivsed_out,GLOBAL_2D_ARRAY)) 
                    var3Dksed_permeab(:,GLOBAL_2D_ARRAY) = 0.0_rsh
                    ALLOCATE(var3Dksed_sigmapsg(nk_nivsed_out,GLOBAL_2D_ARRAY)) 
                    var3Dksed_sigmapsg(:,GLOBAL_2D_ARRAY) = 0.0_rsh
                    ALLOCATE(var3Dksed_dtsdzs(nk_nivsed_out,GLOBAL_2D_ARRAY)) 
                    var3Dksed_dtsdzs(:,GLOBAL_2D_ARRAY) = 0.0_rsh
                    ALLOCATE(var3Dksed_hinder(nk_nivsed_out,GLOBAL_2D_ARRAY)) 
                    var3Dksed_hinder(:,GLOBAL_2D_ARRAY) = 0.0_rsh
                    ALLOCATE(var3Dksed_sed_rate(nk_nivsed_out,GLOBAL_2D_ARRAY)) 
                    var3Dksed_sed_rate(:,GLOBAL_2D_ARRAY) = 0.0_rsh
                    ALLOCATE(var3Dksed_sigmadjge(nk_nivsed_out,GLOBAL_2D_ARRAY)) 
                    var3Dksed_sigmadjge(:,GLOBAL_2D_ARRAY) = 0.0_rsh
                    ALLOCATE(var3Dksed_stateconsol(nk_nivsed_out,GLOBAL_2D_ARRAY)) 
                    var3Dksed_stateconsol(:,GLOBAL_2D_ARRAY) = 0.0_rsh
                ELSE
                    l_outsed_consolidation = .FALSE. ! no output if not all sediment in output
                ENDIF
            ENDIF
                    
        ENDIF

    END SUBROUTINE MUSTANG_allocate_output_arrays

!!===========================================================================

    SUBROUTINE MUSTANG_init_output
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_init_output  ***
    !&E
    !&E ** Purpose : check parameters and prepare needed arrays
    !&E--------------------------------------------------------------------------

    MPI_master_only WRITE(iscreenlog, *)
    MPI_master_only WRITE(iscreenlog, *) '***************************************************************'
    MPI_master_only WRITE(iscreenlog, *) '*****************     SEDIMENT INIT OUTPUT  *******************'
    MPI_master_only WRITE(iscreenlog, *) '***************************************************************'
    MPI_master_only WRITE(iscreenlog, *)

    CALL MUSTANG_check_output_param

    CALL MUSTANG_allocate_output_arrays

    CALL MUSTANG_init_vname

    MPI_master_only WRITE(iscreenlog, *)
    MPI_master_only WRITE(iscreenlog, *) '***************************************************************'
    MPI_master_only WRITE(iscreenlog, *) '***************   SEDIMENT END INIT OUTPUT  *******************'
    MPI_master_only WRITE(iscreenlog, *) '***************************************************************'
    MPI_master_only WRITE(iscreenlog, *)


    END SUBROUTINE MUSTANG_init_output
!!===========================================================================


    SUBROUTINE MUSTANG_init_vname
      !!---------------------------------------------------------------------
      !!                 *** MUSTANG_init_vname  ***
      !!
      !! ** Purpose : initialize MUSTANG output vname and variables
      !!
      !!---------------------------------------------------------------------
        USE dredging, ONLY : l_dredging

        IMPLICIT NONE

        INTEGER indx, isubs

        if (l_dredging) then
            rstMust_nbvar = NT+4
        else
            rstMust_nbvar = NT+3
        endif

        ALLOCATE (rstMust(1:rstMust_nbvar))
        ALLOCATE (rstoutintegerMust(1:rstMust_nbvar))
        ALLOCATE (rstout2DMust(1:rstMust_nbvar))
        ALLOCATE (rstout3DsedMust(1:rstMust_nbvar))
        ALLOCATE (vname_rstMust(20, 1:rstMust_nbvar))

        vname_rstMust(:,:) = ""
        rstoutintegerMust(:) = .FALSE.
        rstout2DMust(:) = .FALSE.
        rstout3DsedMust(:) = .FALSE.

        outMust_nbvar = 2*ntrc_subs + 3*nvpc + 17
#ifdef  key_MUSTANG_V2
        outMust_nbvar = outMust_nbvar + nvpc + 6
#ifdef  key_MUSTANG_bedload
        outMust_nbvar = outMust_nbvar + 4*nvpc + 3
# endif
# endif
        if (l_dyn_insed) outMust_nbvar = outMust_nbvar + 8

        ALLOCATE (hisMust(1:outMust_nbvar))
        ALLOCATE (avgMust(1:outMust_nbvar))
        ALLOCATE (outMust(1:outMust_nbvar))
        ALLOCATE (out2DMust(1:outMust_nbvar))
        ALLOCATE (out3DsedMust(1:outMust_nbvar))
        ALLOCATE (vname_Must(20, 1:outMust_nbvar))

        vname_Must(:,:) = ""
        outMust(:) = .FALSE.
        out2DMust(:) = .FALSE.
        out3DsedMust(:) = .FALSE.


        ! RESTART file
        indx = 1
        vname_rstMust(1,indx) = 'ksmi'
        vname_rstMust(2,indx) = 'lower sediment layer index'
        vname_rstMust(3,indx) = 'no units'
        vname_rstMust(4,indx) = ' '
        vname_rstMust(5,indx) = ' '
        vname_rstMust(6,indx) = ' '
        vname_rstMust(7,indx) = ' '
        rstoutintegerMust(indx) = .TRUE.
        rstout2DMust(indx) = .FALSE.
        rstout3DsedMust(indx) = .FALSE.
        
        indx = 2
        vname_rstMust(1,indx) = 'ksma'
        vname_rstMust(2,indx) = 'upper sediment layer index'
        vname_rstMust(3,indx) = 'no units'
        vname_rstMust(4,indx) = ' '
        vname_rstMust(5,indx) = ' '
        vname_rstMust(6,indx) = ' '
        vname_rstMust(7,indx) = ' '
        rstoutintegerMust(indx) = .TRUE.
        rstout2DMust(indx) = .FALSE.
        rstout3DsedMust(indx) = .FALSE.

        indx = 3
        vname_rstMust(1,indx) = 'DZS'
        vname_rstMust(2,indx) = 'thickness of sediment layer'
        vname_rstMust(3,indx) = 'meter'
        vname_rstMust(4,indx) = ' '
        vname_rstMust(5,indx) = ' '
        vname_rstMust(6,indx) = ' '
        vname_rstMust(7,indx) = ' '
        rstoutintegerMust(indx) = .FALSE.
        rstout2DMust(indx) = .FALSE.
        rstout3DsedMust(indx) = .TRUE.

        indx = 4
        vname_rstMust(1,indx) = 'temp_sed'
        vname_rstMust(2,indx) = 'sediment temperature'
        vname_rstMust(3,indx) = 'Celsius'
        vname_rstMust(4,indx) = ' '
        vname_rstMust(5,indx) = ' '
        vname_rstMust(6,indx) = ' '
        vname_rstMust(7,indx) = ' '
        rstoutintegerMust(indx) = .FALSE.
        rstout2DMust(indx) = .FALSE.
        rstout3DsedMust(indx) = .TRUE.

        indx = 5
        vname_rstMust(1,indx) = 'salt_sed'
        vname_rstMust(2,indx) = 'sediment salinity'
        vname_rstMust(3,indx) = 'PSU'
        vname_rstMust(4,indx) = ' '
        vname_rstMust(5,indx) = ' '
        vname_rstMust(6,indx) = ' '
        vname_rstMust(7,indx) = ' '
        rstoutintegerMust(indx) = .FALSE.
        rstout2DMust(indx) = .FALSE.
        rstout3DsedMust(indx) = .TRUE.

        DO isubs = 1, ntrc_subs
            indx = indx + 1
            vname_rstMust(1,indx) = TRIM(name_var(isubs))//'_sed'
            vname_rstMust(2,indx) = TRIM(long_name_var(isubs))//'_sed'
            vname_rstMust(3,indx) = unit_var(isubs)
            vname_rstMust(4,indx) = TRIM(ADJUSTL(ADJUSTR(standard_name_var(isubs))))//', scalar, series'
            vname_rstMust(5,indx) = ' '
            vname_rstMust(6,indx) = ' '
            vname_rstMust(7,indx) = ' '
            rstoutintegerMust(indx) = .FALSE.
            rstout2DMust(indx) = .FALSE.
            rstout3DsedMust(indx) = .TRUE.
        ENDDO

        if (l_dredging) then
            indx = indx + 1
            vname_rstMust(1,indx) = 'dredg_hsed_init'
            vname_rstMust(2,indx) = 'initial sediment height'
            vname_rstMust(3,indx) = 'm'
            vname_rstMust(4,indx) = ' '
            vname_rstMust(5,indx) = ' '
            vname_rstMust(6,indx) = ' '
            vname_rstMust(7,indx) = ' '
            rstoutintegerMust(indx) = .FALSE.
            rstout2DMust(indx) = .TRUE.
            rstout3DsedMust(indx) = .FALSE.
        endif

        ! end for RESTART file

        ! HIS file

        indx = 1
        vname_Must(1,indx) = 'NB_LAY_SED'
        vname_Must(2,indx) = 'number of sediment layers'
        vname_Must(3,indx) = 'no units'
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_nb_lay_sed) outMust(indx) = .TRUE.
        out2DMust(indx) = .TRUE.
        out3DsedMust(indx) = .FALSE.

        indx = indx+1
        vname_Must(1,indx) = 'HSED'
        vname_Must(2,indx) = 'total thickness of sediment'
        vname_Must(3,indx) = 'meter'
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_hsed) outMust(indx) = .TRUE.
        out2DMust(indx) = .TRUE.
        out3DsedMust(indx) = .FALSE.

        indx = indx+1
        vname_Must(1,indx) = 'TAUSKIN'
        vname_Must(2,indx) = 'total bottom shear stress for erosion'
        vname_Must(3,indx) = 'N/m2'
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_tauskin) outMust(indx) = .TRUE.
        out2DMust(indx) = .TRUE.
        out3DsedMust(indx) = .FALSE.

        indx = indx+1
        vname_Must(1,indx) = 'TAUSKIN_C'
        vname_Must(2,indx) = 'current bottom shear stress for erosion'
        vname_Must(3,indx) = 'N/m2'
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_tauskin_cw) outMust(indx) = .TRUE.
        out2DMust(indx) = .TRUE.
        out3DsedMust(indx) = .FALSE.

        indx = indx+1
        vname_Must(1,indx) = 'TAUSKIN_W'
        vname_Must(2,indx) = 'wave bottom shear stress for erosion'
        vname_Must(3,indx) = 'N/m2'
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_tauskin_cw) outMust(indx) = .TRUE.
        out2DMust(indx) = .TRUE.
        out3DsedMust(indx) = .FALSE.

        indx = indx+1
        vname_Must(1,indx) = 'DZS'
        vname_Must(2,indx) = 'thickness of sediment layer'
        vname_Must(3,indx) = 'meter'
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_dzs) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.

        indx = indx+1
        vname_Must(1,indx) = 'poro'
        vname_Must(2,indx) = 'porosity'
        vname_Must(3,indx) = 'no unit'
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_poro) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.

        indx = indx+1
        vname_Must(1,indx) = 'temp_sed'
        vname_Must(2,indx) = 'sediment temperature'
        vname_Must(3,indx) = 'Celsius'
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_temp_sed) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.


        indx = indx+1
        vname_Must(1,indx) = 'salt_sed'
        vname_Must(2,indx) = 'sediment salinity'
        vname_Must(3,indx) = 'PSU'
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_salt_sed) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.

        DO isubs = 1,ntrc_subs
            indx = indx + 1
            vname_Must(1,indx) = TRIM(name_var(isubs))//'_sed'
            vname_Must(2,indx) = TRIM(long_name_var(isubs))//'_sed'
            vname_Must(3,indx) = unit_var(isubs)
            vname_Must(4,indx) = TRIM(ADJUSTL(ADJUSTR(standard_name_var(isubs))))//', scalar, series'
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_cv_sed .and. l_out_subs(isubs)) outMust(indx) = .TRUE.
            out2DMust(indx) = .FALSE.
            out3DsedMust(indx) = .TRUE.
        ENDDO

        DO isubs = 1,ntrc_subs
            indx = indx + 1
            vname_Must(1,indx) = TRIM(name_var(isubs))//'_ws'
            vname_Must(2,indx) = TRIM(long_name_var(isubs))//' settling velocity'
            vname_Must(3,indx) = "m/s"
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_ws .and. l_out_subs(isubs) .and. &
                ((isubs .GE. imud1) .AND. (isubs .LE. imud2))) THEN
                outMust(indx) = .TRUE.
            ENDIF
            out2DMust(indx) = .FALSE.
            out3DsedMust(indx) = .FALSE.
        ENDDO
    
        DO isubs = 1,nvpc
            indx = indx + 1
            vname_Must(1,indx) = TRIM(name_var(isubs))//'_toce'
            vname_Must(2,indx) = 'critical shear stress'
            vname_Must(3,indx) = 'N/m2'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_toce .and. l_out_subs(isubs)) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

            indx = indx + 1
            vname_Must(1,indx) = TRIM(name_var(isubs))//'_flx_s2w'
            vname_Must(2,indx) = 'erosion flux'
            vname_Must(3,indx) = 'kg.m-2'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_flx_s2w_w2s .and. l_out_subs(isubs)) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

            indx = indx + 1
            vname_Must(1,indx) = TRIM(name_var(isubs))//'_flx_w2s'
            vname_Must(2,indx) = 'deposition flux'
            vname_Must(3,indx) = 'kg.m-2'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_flx_s2w_w2s .and. l_out_subs(isubs)) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

#ifdef  key_MUSTANG_V2
            indx = indx + 1
            vname_Must(1,indx) = TRIM(name_var(isubs))//'_pephm_fcor'
            vname_Must(2,indx) = 'Hindering exposure factor on toce'
            vname_Must(3,indx) = 'no units'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_pephm_fcor .and. l_out_subs(isubs)) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

    
#ifdef key_MUSTANG_bedload
            indx = indx + 1
            vname_Must(1,indx) = TRIM(name_var(isubs))//'_flx_bx'
            vname_Must(2,indx) = 'bedload flux along x-axis'
            vname_Must(3,indx) = 'kg/m/s'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_bedload .and. l_out_subs(isubs) .and. &
                ((isubs .GE. ibedload1) .AND. (isubs .LE. ibedload2))) THEN
                outMust(indx) = .TRUE.
            ENDIF
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

            indx = indx + 1
            vname_Must(1,indx) = TRIM(name_var(isubs))//'_flx_by'
            vname_Must(2,indx) = 'bedload flux along y-axis'
            vname_Must(3,indx) = 'kg/m/s'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_bedload .and. l_out_subs(isubs) .and. &
                ((isubs .GE. ibedload1) .AND. (isubs .LE. ibedload2))) THEN
                outMust(indx) = .TRUE.
            ENDIF
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

            indx = indx + 1
            vname_Must(1,indx) = TRIM(name_var(isubs))//'_bil_bedload'
            vname_Must(2,indx) = 'divergence of bedload flux'
            vname_Must(3,indx) = 'kg.m-2'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_bedload .and. l_out_subs(isubs) .and. &
                ((isubs .GE. ibedload1) .AND. (isubs .LE. ibedload2))) THEN
                outMust(indx) = .TRUE.
            ENDIF
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

            indx = indx + 1
            vname_Must(1,indx) = TRIM(name_var(isubs))//'_fsusp'
            vname_Must(2,indx) = 'fraction of transport in suspension'
            vname_Must(3,indx) = 'no units'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_fsusp .and. l_out_subs(isubs) .and. &
                ((isubs .GE. ibedload1) .AND. (isubs .LE. ibedload2))) THEN
                outMust(indx) = .TRUE.
            ENDIF
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
#endif /* key_MUSTANG_bedload */
#endif /* key_MUSTANG_V2 */
    ENDDO
            indx = indx + 1
            vname_Must(1,indx) = 'frmudsup'
            vname_Must(2,indx) = 'mud fraction in the ksmax layer'
            vname_Must(3,indx) = 'no units'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_frmudsup) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

            indx = indx + 1
            vname_Must(1,indx) = 'dzs_ksmax'
            vname_Must(2,indx) = 'layer thickness at sediment surface'
            vname_Must(3,indx) = 'meter'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_dzs_ksmax) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

            indx = indx + 1
            vname_Must(1,indx) = 'z0sed'
            vname_Must(2,indx) = 'Skin roughness length'
            vname_Must(3,indx) = 'meter'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_z0sed) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

            indx = indx + 1
            vname_Must(1,indx) = 'z0hydro'
            vname_Must(2,indx) = 'hydrodynamic roughness length'
            vname_Must(3,indx) = 'm'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_z0hydro) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

                           
            ! 11 : flx_s2w_noncoh
            indx = indx + 1
            vname_Must(1,indx) = 'flx_s2w_noncoh'
            vname_Must(2,indx) = 'erosion flux of non-cohesive sediments (sum: isand1 to isand2)' 
            vname_Must(3,indx) = 'kg.m-2'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_flx_s2w_w2s) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
            
            ! 12 : flx_w2s_noncoh
            indx = indx + 1
            vname_Must(1,indx) = 'flx_w2s_noncoh'
            vname_Must(2,indx) = 'deposition flux of non-cohesive sediments (sum: isand1 to isand2)'
            vname_Must(3,indx) = 'kg.m-2'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_flx_s2w_w2s) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
            
            ! 13 : flx_s2w_coh
            indx = indx + 1
            vname_Must(1,indx) = 'flx_s2w_coh'
            vname_Must(2,indx) = 'erosion flux of cohesive sediments (sum: imud1 to imud2)'
            vname_Must(3,indx) = 'kg.m-2'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_flx_s2w_w2s) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
            
            ! 14 : flx_w2s_coh
            indx = indx + 1
            vname_Must(1,indx) = 'flx_w2s_coh'
            vname_Must(2,indx) = 'deposition flux of cohesive sediments (sum: imud1 to imud2)'
            vname_Must(3,indx) = 'kg.m-2'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_flx_s2w_w2s) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

#ifdef key_MUSTANG_V2
            ! 3 : dzs_aclay_comp_save
            indx = indx + 1
            vname_Must(1,indx) = 'dzs_aclay_comp_save'
            vname_Must(2,indx) = 'Theoretical active layer thickness Harris and Wiberg 1997'
            vname_Must(3,indx) = 'meter'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_theoric_active_layer) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.

            ! 5 : tero_noncoh (cumulated time (in hours) elapsed in non cohesive regime)
            indx = indx + 1
            vname_Must(1,indx) = 'tero_noncoh'
            vname_Must(2,indx) = 'time elapsed in the non-cohesive erosion regime'
            vname_Must(3,indx) = 'hours'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_ero_details) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
            
            ! 6 : tero_coh (cumulated time (in hours) elapsed in cohesive regime)
            indx = indx + 1
            vname_Must(1,indx) = 'tero_coh'
            vname_Must(2,indx) = 'time elapsed in the cohesive erosion regime'
            vname_Must(3,indx) = 'hours'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_ero_details) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
            
            ! 7 : pct_iter_noncoh
            indx = indx + 1
            vname_Must(1,indx) = 'pct_iter_noncoh'
            vname_Must(2,indx) = 'part of erosion iterations in the non-cohesive regime'
            vname_Must(3,indx) = 'percent'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_ero_details) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
            
            ! 8 : pct_iter_coh
            indx = indx + 1
            vname_Must(1,indx) = 'pct_iter_coh'
            vname_Must(2,indx) = 'part of erosion iterations in the cohesive regime'
            vname_Must(3,indx) = 'percent'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_ero_details) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
            
            ! 9 : niter_ero
            indx = indx + 1
            vname_Must(1,indx) = 'niter_ero'
            vname_Must(2,indx) = 'Number of iterations in sed_erosion during time step'
            vname_Must(3,indx) = 'no units'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_ero_details) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
      
            
            
#ifdef key_MUSTANG_bedload
            ! 16 : flx_bx_int
            indx = indx + 1
            vname_Must(1,indx) = 'flx_bx_int'
            vname_Must(2,indx) = 'total bedload flux along  x-axis (sum: igrav1 to isand2)'
            vname_Must(3,indx) = 'kg/m/s'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_bedload) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
            
            ! 17 : flx_by_int
            indx = indx + 1
            vname_Must(1,indx) = 'flx_by_int'
            vname_Must(2,indx) = 'total bedload flux along y-axis (sum: igrav1 to isand2)'
            vname_Must(3,indx) = 'kg/m/s'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_bedload) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
            
            ! 18 : bil_bedload_int
            indx = indx + 1
            vname_Must(1,indx) = 'bil_bedload_int'
            vname_Must(2,indx) = 'divergence of total bedload flux (sum: igrav1 to isand2)'
            vname_Must(3,indx) = 'kg/m2'
            vname_Must(4,indx) = ' '
            vname_Must(5,indx) = ' '
            vname_Must(6,indx) = ' '
            vname_Must(7,indx) = ' '
            IF (l_outsed_bedload) outMust(indx) = .TRUE.
            out2DMust(indx) = .TRUE.
            out3DsedMust(indx) = .FALSE.
            
#endif /* key_MUSTANG_bedload */
#endif /* key_MUSTANG_V2 */

    IF (l_dyn_insed) THEN
        indx = indx + 1
        vname_Must(1,indx) = 'loadograv'
        vname_Must(2,indx) = 'excess of interstitial water pressure in the middle of the layer'
        vname_Must(3,indx) = ' '
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_consolidation) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.

        indx = indx + 1
        vname_Must(1,indx) = 'sigmadjge'
        vname_Must(2,indx) = 'sigma unseparated (without the share of water)'
        vname_Must(3,indx) = ' '
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_consolidation) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.

        indx = indx + 1
        vname_Must(1,indx) = 'sigmapsg'
        vname_Must(2,indx) = 'effective stress (transmitted from grain to grain)'
        vname_Must(3,indx) = ' '
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_consolidation) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.

        indx = indx + 1
        vname_Must(1,indx) = 'stateconsol'
        vname_Must(2,indx) = 'state of consolidation indicator'
        vname_Must(3,indx) = ' '
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_consolidation) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.

        indx = indx + 1
        vname_Must(1,indx) = 'permeab'
        vname_Must(2,indx) = 'permeability'
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_consolidation) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.

        indx = indx + 1
        vname_Must(1,indx) = 'hinder'
        vname_Must(2,indx) = 'shackling sand / gravel between 0 and 1'
        vname_Must(3,indx) = 'no unit'
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_consolidation) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.

        indx = indx + 1
        vname_Must(1,indx) = 'sed_rate'
        vname_Must(2,indx) = 'advection speed of mud particles'
        vname_Must(3,indx) = ' '
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_consolidation) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.

        indx = indx + 1
        vname_Must(1,indx) = 'dtsdzs'
        vname_Must(2,indx) = 'dtsdzs'
        vname_Must(3,indx) = ' '
        vname_Must(4,indx) = ' '
        vname_Must(5,indx) = ' '
        vname_Must(6,indx) = ' '
        vname_Must(7,indx) = ' '
        IF (l_outsed_consolidation) outMust(indx) = .TRUE.
        out2DMust(indx) = .FALSE.
        out3DsedMust(indx) = .TRUE.
    ENDIF

    END SUBROUTINE MUSTANG_init_vname

!!==========================================================================================================


    SUBROUTINE MUSTANG_morphoinit(ifirst, ilast, jfirst, jlast   &
#if defined MORPHODYN 
                  , dhsed                                                &
#endif
                  )
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE morphoinit  ***
    !&E
    !&E ** Purpose : initialization of hsed, hsed_previous
    !&E ** Called by :  MUSTANG_init_sediment 
    !&E--------------------------------------------------------------------------

    !! * Arguments 
    INTEGER, INTENT(IN)                    :: ifirst, ilast, jfirst, jlast
#if defined MORPHODYN
    REAL(KIND=rsh),DIMENSION(GLOBAL_2D_ARRAY),INTENT(INOUT)           :: dhsed                       
#endif
    !! * Local declarations
    INTEGER :: i, j, k
 
    !!  initialisation of hsed  inside the domain  !!!
    !!  eliminate the meshes at open boundaries ??? not reuse from MARS? TODO or not ?
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO j=jfirst,jlast
    DO i=ifirst,ilast
        !  WARNING : not in first and last mesh on i and j axis (boundaries) ! ??? not reuse from MARS? TODO or not ?
        hsed(i,j)=0.0_rsh
        DO k=ksmi(i,j),ksma(i,j)
            hsed(i,j)=hsed(i,j)+dzs(k,i,j)
        ENDDO   

        IF(l_morphocoupl)THEN
            !  WARNING: not in first and last mesh on i and j axis (boundaries)! ??? not reuse from MARS? TODO or not ?
            hsed_previous(i,j)=hsed(i,j)      
#if defined MORPHODYN
            dhsed(i,j)=hsed_previous(i,j)-hsed(i,j)
#endif
        ENDIF        ! endif l_morphocoupl
    ENDDO
    ENDDO

    END SUBROUTINE MUSTANG_morphoinit
!!===========================================================================

    SUBROUTINE MUSTANG_alloc()
    !&E--------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_alloc  ***
    !&E
    !&E ** Purpose : allocation of arrays relative to sediment
    !&E
    !&E ** Called by : MUSTANG_init
    !&E
    !&E--------------------------------------------------------------------------
    
    !! * Local declarations
    INTEGER :: iv

    !! * Executable part
    ALLOCATE(ws_sand(nvp))        
    ALLOCATE(diamstar(nvp))
    ALLOCATE(rosmrowsros(nvp))
    ALLOCATE(stresscri0(nvp))
    ALLOCATE(tetacri0(nvp))
    ALLOCATE (typart(-1:nv_adv))
    typart(-1:0)=0.0_rsh
    typart(1:nvpc)=1.0_rsh
    typart(nvpc+1:nv_adv)=0.0_rsh
#ifdef key_MUSTANG_V2
    ALLOCATE(E0_sand(nvp))
    E0_sand(1:nvp)=0.0_rsh
#endif   
#if  ! defined key_noTSdiss_insed
    ALLOCATE(ivdiss(-1:nv_adv-nvp))
#endif
 
 
    !  allocation of spatial variables  
    ALLOCATE(ksmi(GLOBAL_2D_ARRAY))
    ALLOCATE(ksma(GLOBAL_2D_ARRAY))
    ALLOCATE(hsed(GLOBAL_2D_ARRAY))  
    ALLOCATE(z0sed(GLOBAL_2D_ARRAY))
    ALLOCATE(tauskin(GLOBAL_2D_ARRAY))
    ALLOCATE(tauskin_c(GLOBAL_2D_ARRAY))
    ALLOCATE(tauskin_w(GLOBAL_2D_ARRAY))
    ALLOCATE(ustarbot(GLOBAL_2D_ARRAY))
    ALLOCATE(dzsmax(GLOBAL_2D_ARRAY))
    ALLOCATE(htot(GLOBAL_2D_ARRAY))
    ALLOCATE(alt_cw1(GLOBAL_2D_ARRAY))
    ALLOCATE(epn_bottom_MUSTANG(GLOBAL_2D_ARRAY))  
    ALLOCATE(sal_bottom_MUSTANG(GLOBAL_2D_ARRAY))  
    ALLOCATE(temp_bottom_MUSTANG(GLOBAL_2D_ARRAY))
    ALLOCATE(cw_bottom_MUSTANG(nv_tot,GLOBAL_2D_ARRAY))
    ALLOCATE(ws3_bottom_MUSTANG(nvp,GLOBAL_2D_ARRAY))  
    ALLOCATE(roswat_bot(GLOBAL_2D_ARRAY))  
    ksmi(GLOBAL_2D_ARRAY) = 0
    ksma(GLOBAL_2D_ARRAY) = 0
    hsed(GLOBAL_2D_ARRAY) = 0.0_rsh
    z0sed(GLOBAL_2D_ARRAY) = 0.0_rsh
    tauskin(GLOBAL_2D_ARRAY) = 0.0_rsh
    tauskin_c(GLOBAL_2D_ARRAY) = 0.0_rsh
    tauskin_w(GLOBAL_2D_ARRAY) = 0.0_rsh
    ustarbot(GLOBAL_2D_ARRAY) = 0.0_rsh
    dzsmax(GLOBAL_2D_ARRAY) = 0.0_rsh
    htot(GLOBAL_2D_ARRAY) = 0.0_rsh
    alt_cw1(GLOBAL_2D_ARRAY) = 0.0_rsh
    epn_bottom_MUSTANG(GLOBAL_2D_ARRAY) = 0.0_rsh
    sal_bottom_MUSTANG(GLOBAL_2D_ARRAY) = 0.0_rsh
    temp_bottom_MUSTANG(GLOBAL_2D_ARRAY) = 0.0_rsh
    cw_bottom_MUSTANG(nv_tot,GLOBAL_2D_ARRAY) = 0.0_rsh
    ws3_bottom_MUSTANG(nvp,GLOBAL_2D_ARRAY) = 0.0_rsh
    roswat_bot(GLOBAL_2D_ARRAY) = 0.0_rsh
#if ! defined key_noTSdiss_insed
    ALLOCATE(phitemp_s(GLOBAL_2D_ARRAY))
    phitemp_s(GLOBAL_2D_ARRAY) = 0.0_rsh
#endif
#ifdef key_MUSTANG_V2
    ALLOCATE(sigmapsg(ksdmin:ksdmax))
    ALLOCATE(stateconsol(ksdmin:ksdmax))
    ALLOCATE(permeab(ksdmin:ksdmax))
    ALLOCATE(psi_sed(nvp))
    psi_sed(1:nvp) = 0.0_rsh
#endif
#ifdef key_sand2D
    ALLOCATE(rouse2D(nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(sum_tmp(nv_adv,GLOBAL_2D_ARRAY))
    rouse2D(1:nv_adv,GLOBAL_2D_ARRAY) = 0.0_rsh
    sum_tmp(1:nv_adv,GLOBAL_2D_ARRAY) = 0.0_rsh
#endif
    
    ALLOCATE(cv_sed(-1:nv_tot,ksdmin:ksdmax,GLOBAL_2D_ARRAY))
    ALLOCATE(c_sedtot(ksdmin:ksdmax,GLOBAL_2D_ARRAY))
    ALLOCATE(poro(ksdmin:ksdmax,GLOBAL_2D_ARRAY))
    ALLOCATE(dzs(ksdmin:ksdmax,GLOBAL_2D_ARRAY))
    ALLOCATE(flx_s2w(-1:nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(flx_w2s(-1:nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(flx_w2s_sum(-1:nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(corflux(nv_adv, GLOBAL_2D_ARRAY))
    ALLOCATE(corfluy(nv_adv, GLOBAL_2D_ARRAY))
    ALLOCATE(fludif(-1:nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(fluconsol(-1:nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(fluconsol_drycell(-1:nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(flu_dyninsed(-1:nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(gradvit(N,GLOBAL_2D_ARRAY))

    cv_sed(-1:nv_tot,ksdmin:ksdmax,GLOBAL_2D_ARRAY) = 0.0_rsh
    c_sedtot(ksdmin:ksdmax,GLOBAL_2D_ARRAY) = 0.0_rsh
    poro(ksdmin:ksdmax,GLOBAL_2D_ARRAY) = 0.0_rsh
    dzs(ksdmin:ksdmax,GLOBAL_2D_ARRAY) = 0.0_rsh
    flx_s2w(-1:nv_adv,GLOBAL_2D_ARRAY) = 0.0_rsh
    flx_w2s(-1:nv_adv,GLOBAL_2D_ARRAY) = 0.0_rsh
    flx_w2s_sum(-1:nv_adv,GLOBAL_2D_ARRAY) = 0.0_rsh
    corflux(1:nv_adv, GLOBAL_2D_ARRAY) = 1.0_rsh
    corfluy(1:nv_adv, GLOBAL_2D_ARRAY) = 1.0_rsh
    fludif(-1:nv_adv,GLOBAL_2D_ARRAY) = 0.0_rsh
    fluconsol(-1:nv_adv,GLOBAL_2D_ARRAY) = 0.0_rsh
    fluconsol_drycell(-1:nv_adv,GLOBAL_2D_ARRAY) = 0.0_rsh
    flu_dyninsed(-1:nv_adv,GLOBAL_2D_ARRAY) = 0.0_rsh
    gradvit(N,GLOBAL_2D_ARRAY) = 0.0_rsh


#ifdef key_MUSTANG_V2
    ALLOCATE(poro_mud(ksdmin:ksdmax,GLOBAL_2D_ARRAY))
    ALLOCATE(crel_mud(ksdmin:ksdmax,GLOBAL_2D_ARRAY))
    ALLOCATE(l_isitcohesive(GLOBAL_2D_ARRAY))
    poro_mud(ksdmin:ksdmax,GLOBAL_2D_ARRAY) = 0.0_rsh
    crel_mud(ksdmin:ksdmax,GLOBAL_2D_ARRAY) = 0.0_rsh
    l_isitcohesive(GLOBAL_2D_ARRAY) = .FALSE.

#ifdef key_MUSTANG_bedload
    ALLOCATE( flx_bx(1:nvp,GLOBAL_2D_ARRAY)) ! Warning /Baptiste : m1p1 sur les 2 indices au lieu d1
    ALLOCATE( flx_by(1:nvp,GLOBAL_2D_ARRAY)) ! Warning /Baptiste : m1p1 sur les 2 indices au lieu d1
    ALLOCATE( slope_dhdx(GLOBAL_2D_ARRAY))
    ALLOCATE( slope_dhdy(GLOBAL_2D_ARRAY))
    ALLOCATE( sedimask_h0plusxe(GLOBAL_2D_ARRAY)) ! Warning /Baptiste : m1p1 sur les 2 indices au lieu d1
    flx_bx(1:nvp,GLOBAL_2D_ARRAY)=0.0_rsh 
    flx_by(1:nvp,GLOBAL_2D_ARRAY)=0.0_rsh
    slope_dhdx(GLOBAL_2D_ARRAY)=0.0_rsh
    slope_dhdy(GLOBAL_2D_ARRAY)=0.0_rsh
    sedimask_h0plusxe(GLOBAL_2D_ARRAY)=0.0_rsh
#endif
#endif

    ALLOCATE(phieau_s2w(GLOBAL_2D_ARRAY))
    ALLOCATE(phieau_s2w_drycell(GLOBAL_2D_ARRAY))
    ALLOCATE(phieau_s2w_consol(GLOBAL_2D_ARRAY))
    phieau_s2w(GLOBAL_2D_ARRAY)=0.0_rlg
    phieau_s2w_drycell(GLOBAL_2D_ARRAY)=0.0_rlg
    phieau_s2w_consol(GLOBAL_2D_ARRAY)=0.0_rlg

    ALLOCATE( raphbx(GLOBAL_2D_ARRAY), raphby(GLOBAL_2D_ARRAY) )
    ALLOCATE( tauskin_x(GLOBAL_2D_ARRAY), tauskin_y(GLOBAL_2D_ARRAY) )
    raphbx(GLOBAL_2D_ARRAY)=0.0_rsh
    raphby(GLOBAL_2D_ARRAY)=0.0_rsh

    ALLOCATE(flx_s2w_corim1(-1:nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(flx_s2w_corip1(-1:nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(flx_s2w_corjm1(-1:nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(flx_s2w_corjp1(-1:nv_adv,GLOBAL_2D_ARRAY))
    ALLOCATE(flx_w2s_corin(nvp,GLOBAL_2D_ARRAY))
    ALLOCATE(flx_w2s_corim1(nvp,GLOBAL_2D_ARRAY))
    ALLOCATE(flx_w2s_corip1(nvp,GLOBAL_2D_ARRAY))
    ALLOCATE(flx_w2s_corjm1(nvp,GLOBAL_2D_ARRAY))
    ALLOCATE(flx_w2s_corjp1(nvp,GLOBAL_2D_ARRAY))
#if ! defined key_nofluxwat_IWS
    ALLOCATE(phieau_s2w_corim1(GLOBAL_2D_ARRAY))
    ALLOCATE(phieau_s2w_corip1(GLOBAL_2D_ARRAY))
    ALLOCATE(phieau_s2w_corjm1(GLOBAL_2D_ARRAY))
    ALLOCATE(phieau_s2w_corjp1(GLOBAL_2D_ARRAY))
#endif
  
    ! Initialization
    ALLOCATE(cini_sed(nv_state))
    DO iv = 1, nv_adv
        cini_sed(iv) = cini_sed_r(irk_fil(iv)) * unit_modif_mudbio_N2dw(irk_fil(iv))
    ENDDO   

#ifdef key_Pconstitonly_insed 
    nv_use = nvpc
#else
    nv_use = nvp
#endif  

#if ! defined key_noTSdiss_insed
    ! counting of dissolved variables for diffusion in the sediment
    ivdiss(:) = 0
    ivdiss(-1) = -1
    ivdiss(0) = 0
#if ! defined key_Pconstitonly_insed
    DO iv = 1, nv_adv-nvp
        ivdiss(iv) = iv + nvp
    ENDDO  
#endif
#endif

    !  option morpho
    !!!!!!!!!!!!!!!!!!!!
    IF(l_morphocoupl) THEN
    !! Warning : no hx, hy in other model than MARS 
        ALLOCATE(hsed_previous(GLOBAL_2D_ARRAY))
        hsed_previous(GLOBAL_2D_ARRAY) = 0.0_rsh
    ENDIF

    !! declaration of the MUSTANG variables needed in the hydro model 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(flx_s2w_CROCO(GLOBAL_2D_ARRAY,1:NT))
    ALLOCATE(flx_w2s_CROCO(GLOBAL_2D_ARRAY,1:NT))
    ALLOCATE(flx_w2s_sum_CROCO(GLOBAL_2D_ARRAY,1:NT))
    flx_s2w_CROCO(GLOBAL_2D_ARRAY,1:NT) = 0.0_rsh
    flx_w2s_CROCO(GLOBAL_2D_ARRAY,1:NT) = 0.0_rsh
    flx_w2s_sum_CROCO(GLOBAL_2D_ARRAY,1:NT) = 0.0_rsh
#if ! defined key_nofluxwat_IWS && ! defined key_noTSdiss_insed
    ALLOCATE(WATER_FLUX_INPUTS(GLOBAL_2D_ARRAY,1:N)) ! not operationnal, stil to code **TODO**
    WATER_FLUX_INPUTS(GLOBAL_2D_ARRAY,1:N) = 0.0_rsh
#endif    
    
    ALLOCATE(fwet(GLOBAL_2D_ARRAY))
    fwet(:,:) = 1.0_rsh
  

    END SUBROUTINE MUSTANG_alloc



#endif /* MUSTANG */

END MODULE initMUSTANG
