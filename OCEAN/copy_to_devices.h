!$acc enter data if(compute_on_device) copyin(

!ocean2d.h
!$acc& zeta
!$acc&, ubar
!$acc&, vbar
#if !defined SOLVE3D && defined M2_HADV_UP3
!$acc&, urhs
!$acc&, vrhs
!$acc&, Duon
!$acc&, DVom
#endif

!scalars.h
#ifdef SOLVE3D
!$acc&, sc_w, Cs_w, sc_r, Cs_r
# ifdef TRACERS
!$acc&, tnu2, tnu4
# endif
!$acc&, weight
#endif
#if defined SPONGE || \
     defined TNUDGING   || defined M2NUDGING  || \
     defined M3NUDGING  || defined ZNUDGING
#endif
#if defined T_FRC_BRY     || defined M2_FRC_BRY    || \
     defined M3_FRC_BRY    || defined Z_FRC_BRY     || \
     defined W_FRC_BRY     || defined NBQ_FRC_BRY   || \
     defined TCLIMATOLOGY  || defined M2CLIMATOLOGY || \
     defined M3CLIMATOLOGY || defined ZCLIMATOLOGY  || \
     defined WCLIMATOLOGY  || defined NBQCLIMATOLOGY
#endif
#ifdef ABL1D
#endif
#if defined DIAGNOSTICS_TS
#endif
#if defined DIAGNOSTICS_UV
#endif
#ifdef DIAGNOSTICS_VRT
#endif
#ifdef DIAGNOSTICS_EK
#endif
#ifdef DIAGNOSTICS_PV
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
#endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
#endif
#ifdef DIAGNOSTICS_BIO
#endif
#if defined SOLVE3D && defined TRACERS
!$acc&, got_tini
#endif
#ifdef SEDIMENT
!$acc&, got_inised
#endif
#ifdef BBL
!$acc&, got_inibed
#endif
#if defined DIAGNOSTICS_TS
#endif
#if defined DIAGNOSTICS_UV
#endif
#if defined DIAGNOSTICS_VRT
#endif
#if defined DIAGNOSTICS_EK
#endif
#if defined DIAGNOSTICS_PV
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
#endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
#endif
#ifdef DIAGNOSTICS_BIO
#endif
#ifdef SOLVE3D
#endif
#if defined SPONGE || \
     defined TNUDGING   || defined M2NUDGING  || \
     defined M3NUDGING  || defined ZNUDGING
#endif
#if defined T_FRC_BRY     || defined M2_FRC_BRY    || \
     defined M3_FRC_BRY    || defined Z_FRC_BRY     || \
     defined W_FRC_BRY     ||                          \
     defined TCLIMATOLOGY  || defined M2CLIMATOLOGY || \
     defined M3CLIMATOLOGY || defined ZCLIMATOLOGY  || \
     defined WCLIMATOLOGY
#endif
#if defined DIAGNOSTICS_TS
#endif
#if defined DIAGNOSTICS_UV
#endif
#ifdef DIAGNOSTICS_VRT
#endif
#ifdef DIAGNOSTICS_EK
#endif
#ifdef DIAGNOSTICS_PV
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
#endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
#endif
#ifdef DIAGNOSTICS_BIO
#endif
#ifdef ABL1D
#endif
#if defined SOLVE3D  && !defined LMD_MIXING
# ifdef TRACERS
!$acc&, Akt_bak
# endif
#endif
#if defined BIOLOGY && defined TRACERS
!$acc&, global_sum
#endif
#ifdef SOLITON
#else
#endif

!ocean3d.h
#ifdef SOLVE3D
!$acc&, u
!$acc&, v
!$acc&, t
!$acc&, Hz
!$acc&, Hz_bak
!$acc&, z_r
!$acc&, z_w
!$acc&, Huon
!$acc&, Hvom
!$acc&, We
# ifdef VADV_ADAPT_IMP
!$acc&, Wi
# endif
# ifdef NBQ
!$acc&, wz
#  ifdef NBQ_MASS
!$acc&, Hzr
#  else
#   define Hzr Hz
#  endif
# endif
# if defined UV_VIS4 && defined UV_MIX_GEO
!$acc&, z_u
!$acc&, z_v
!$acc&, dz_u
!$acc&, dz_v
# endif
!$acc&, rho1
!$acc&, rho
# if defined NONLIN_EOS && defined SPLIT_EOS
!$acc&, qp1
# endif
# ifdef BIOLOGY
#  ifdef BIO_NChlPZD
!$acc&, theta
#  elif defined BIO_BioEBUS
!$acc&, AOU
!$acc&, wind10
#  endif
# endif  /* BIOLOGY */
# ifdef OXYGEN
!$acc&, u10
!$acc&, Kv_O2
!$acc&, O2satu
# endif /* OXYGEN */
#endif  /* SOLVE3D */

!grid.h
!$acc&, h
!$acc&, hinv
!$acc&, f
!$acc&, fomn
#ifdef MORPHODYN
!$acc&, dh
#endif
#ifdef CURVGRID
!$acc&, angler
#endif
#ifdef SPHERICAL
!$acc&, latr
!$acc&, lonr
!$acc&, latu
!$acc&, lonu
!$acc&, latv
!$acc&, lonv
#else
!$acc&, xp
!$acc&, xr
!$acc&, yp
!$acc&, yr
#endif
!$acc&, pm
!$acc&, pn
!$acc&, om_r
!$acc&, on_r
!$acc&, om_u
!$acc&, on_u
!$acc&, om_v
!$acc&, on_v
!$acc&, om_p
!$acc&, on_p
!$acc&, pn_u
!$acc&, pm_v
!$acc&, pm_u
!$acc&, pn_v
#if (defined CURVGRID && defined UV_ADV)
!$acc&, dmde
!$acc&, dndx
#endif
!$acc&, pmon_p
!$acc&, pmon_r
!$acc&, pmon_u
!$acc&, pnom_p
!$acc&, pnom_r
!$acc&, pnom_v
!$acc&, grdscl
#ifdef MASKING
!$acc&, rmask
!$acc&, pmask
!$acc&, umask
!$acc&, vmask
!$acc&, pmask2
#endif
#ifdef WET_DRY
!$acc&, rmask_wet
!$acc&, pmask_wet
!$acc&, umask_wet
!$acc&, vmask_wet
!$acc&, rmask_wet_avg
!$acc&, Dcrit
!$acc&, wetdry
#endif
#ifdef WAVE_DRY
!$acc&, rmask_wavewet
!$acc&, pmask_wavewet
!$acc&, umask_wavewet
!$acc&, vmask_wavewet
!$acc&, Dcrit_wave
!$acc&, wave_wetdry
#endif
#ifdef REDUC_SECTION
!$acc&, ureduc
!$acc&, vreduc
#endif
!$acc&, zob
#if defined UV_COR_NT || defined CROCO_QH
!$acc&, e
!$acc&, eomn
!$acc&, cosa
!$acc&, sina
#endif

!coupling.h
#ifdef SOLVE3D
# ifdef VAR_RHO_2D
!$acc&, rhoA
!$acc&, rhoS
# endif
!$acc&, rufrc
!$acc&, rvfrc
!$acc&, rufrc_bak
!$acc&, rvfrc_bak
!$acc&, Zt_avg1
!$acc&, DU_avg1
!$acc&, DV_avg1
!$acc&, DU_avg2
!$acc&, DV_avg2
#endif

!private_scratch.h
!$acc&, A2d, A3d
#if defined SEDIMENT || defined LMD_MIXING
!$acc&, B2d
#endif
#if defined ABL1D
!$acc&, T1d
!$acc&, T2d, T3d
#endif

!mixing.h
#if defined UV_VIS2 || defined SPONGE_VIS2
!$acc&, visc2_r
!$acc&, visc2_p
!$acc&, visc2_sponge_r
!$acc&, visc2_sponge_p
#endif
#if defined UV_VIS4
# if !defined SPONGE_VIS2
!$acc&, visc2_sponge_r
!$acc&, visc2_sponge_p
# endif
!$acc&, visc4_sponge_r
!$acc&, visc4_sponge_p
!$acc&, visc4_r
!$acc&, visc4_p
#endif
#if defined TS_DIF2 || defined SPONGE_DIF2
!$acc&, diff2_sponge
!$acc&, diff2
#endif
#if defined TS_DIF4
# if !defined SPONGE_DIF2
!$acc&, diff2_sponge
# endif
!$acc&, diff4_sponge
!$acc&, diff4
#endif
#ifdef VIS_COEF_3D
!$acc&, visc3d_r
!$acc&, visc3d_p
#endif
#ifdef DIF_COEF_3D
!$acc&, diff3d_u
!$acc&, diff3d_v
# if defined TS_DIF_SMAGO || defined GLS_MIXING_3D
!$acc&, diff3d_r
# endif
#endif
#if defined TS_MIX_ISO || defined TS_MIX_GEO
!$acc&, dRdx
!$acc&, dRde
!$acc&, idRz
#endif /*  TS_MIX_ISO || TS_MIX_GEO */
#ifdef SPONGE_SED
!$acc&, sed_sponge
#endif
#ifdef SOLVE3D
!$acc&, Akv
!$acc&, Akt
# ifdef GLS_MIXING
!$acc&, Akv_old
!$acc&, Akt_old
# endif
# if defined ANA_VMIX   || defined LMD_MIXING \
  || defined LMD_SKPP   || defined LMD_BKPP \
  || defined GLS_MIXING || defined UV_VIS_SMAGO_3D
!$acc&, bvf
# endif
# ifdef BIOLOGY
!$acc&, hel
# endif
# ifdef LMD_MIXING
!$acc&, ustar
#  ifdef LMD_LANGMUIR
!$acc&, Langmuir
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP
!$acc&, kbl
!$acc&, kbbl
!$acc&, hbbl
#   ifdef LMD_SKPP2005
!$acc&, hbls
#   else
!$acc&, hbl
#   endif
#   ifdef LMD_NONLOCAL
!$acc&, ghats
#   endif
#  endif /* LMD_SKPP || LMD_BKPP */
# elif defined GLS_MIXING
!$acc&, trb
!$acc&, Lscale
!$acc&, Eps_gls
!$acc&, kbl
!$acc&, hbl
# endif /* GLS_MIXING */
#else
# define u(i,j,k,nrhs) ubar(i,j,kstp)
# define v(i,j,k,nrhs) vbar(i,j,kstp)
# define visc3d_r(i,j,k)  visc2_r(i,j)
# define visc3d_p(i,j,k)  visc2_p(i,j)
# define exchange_p3d_tile(a,b,c,d,visc3d_p) exchange_p2d_tile(a,b,c,d,visc2_p)
#endif /* SOLVE3D */

!forces.h
!$acc&, sustr
!$acc&, svstr
#if defined OA_COUPLING || defined OW_COUPLING
!$acc&, smstr
# ifdef READ_PATM
!$acc&, patm2d
# endif
#endif
#ifdef OW_COUPLING
!$acc&, twox
!$acc&, twoy
!$acc&, foc
!$acc&, tawx
!$acc&, tawy
#endif
#ifndef ANA_SMFLUX
!$acc&, sustrg
!$acc&, svstrg
!$acc&, sustrp, svstrp, sms_time
# if defined SFLUX_CFB && !defined BULK_FLUX
!$acc&, wspd
# endif
# ifdef CFB_WIND_TRA
!$acc&, wspd_cfb
# endif
#endif /* !ANA_SMFLUX */
!$acc&, bustr
!$acc&, bvstr
#ifndef ANA_BMFLUX
!$acc&, bustrg
!$acc&, bvstrg
!$acc&, bms_tintrp, bustrp, bvstrp, tbms
# undef BMFLUX_DATA
#endif /* !ANA_BMFLUX */
#ifdef SOLVE3D
!$acc&, stflx
# if defined BULK_FLUX
!$acc&, shflx_rsw
!$acc&, shflx_rlw
!$acc&, shflx_lat
!$acc&, shflx_sen
# endif /* BULK_FLUX */
# if defined SST_SKIN && defined TEMPERATURE
!$acc&, sst_skin
!$acc&, dT_skin
# endif/* SST_SKIN && TEMPERATURE */
# if defined TRACERS
#  if !defined ANA_STFLUX || !defined ANA_SSFLUX
!$acc&, stflxg
!$acc&, stflxp, stf_time
!$acc&, stf_cycle, stf_scale
!$acc&, itstf, stf_ncycle, stf_rec
!$acc&, lstfgrd, stf_tid, stf_id
#   undef STFLUX_DATA
#  endif /* !ANA_STFLUX || !ANA_SSFLUX */
!$acc&, btflx
#  if defined BHFLUX || defined BWFLUX
!$acc&, btflxg
!$acc&, btflxp, btf_time
!$acc&, btf_cycle, btf_scale
!$acc&, itbtf, btf_ncycle, btf_rec
!$acc&, lbtfgrd, btf_tid, btf_id
#   undef BTFLUX_DATA
#  endif /*  BHFLUX */
# endif /* TRACERS */
# if defined QCORRECTION && (defined TEMPERATURE || defined SALINITY)
!$acc&, dqdt
!$acc&, sst
#  ifndef ANA_SST
!$acc&, dqdtg
!$acc&, sstg
!$acc&, sstp, dqdtp, sst_time
#   undef SST_DATA
#  endif /* !ANA_SST */
# endif /* QCORRECTION && TEMPERATURE */
# if defined SALINITY && defined SFLX_CORR
!$acc&, sss
#  if !defined QCORRECTION
!$acc&, dqdt
#  endif
#  ifndef ANA_SSS
!$acc&, sssg
!$acc&, sssp, sss_time
#   if !defined QCORRECTION
!$acc&, dqdtg
!$acc&, dqdtp
#   endif
#   undef SSS_DATA
#  endif /* !ANA_SSS */
# endif /* SALINITY && SFLX_CORR */
# if defined BULK_FLUX
!$acc&, tair
!$acc&, rhum
!$acc&, prate
!$acc&, radlw
!$acc&, radsw
!$acc&, wspd
#  ifdef READ_PATM
!$acc&, patm2d
#  endif
!$acc&, uwnd
!$acc&, vwnd
#  ifdef DIURNAL_INPUT_SRFLX
!$acc&, radswbio
#  endif
!$acc&, tairg
!$acc&, rhumg
!$acc&, prateg
!$acc&, radlwg
!$acc&, radswg
#  ifdef READ_PATM
!$acc&, patmg
#  endif
#  ifdef ONLINE
!$acc&, uwndg_norot
!$acc&, radswg_down
#  endif
!$acc&, uwndg
!$acc&, vwndg
!$acc&, wspdg
#  ifdef DIURNAL_INPUT_SRFLX
!$acc&, radswbiog
#  endif
!$acc&, tairp, rhump, pratep, radlwp, radswp
#  ifdef READ_PATM
!$acc&, patmp
#  endif
!$acc&, uwndp, vwndp
#  ifdef DIURNAL_INPUT_SRFLX
!$acc&, radswbiop
#  endif
!$acc&, bulk_time
# endif /* BULK_FLUX */
!$acc&, srflx
# ifdef ANA_DIURNAL_SW
!$acc&, sin_phi
!$acc&, cos_phi
!$acc&, tan_phi
# endif
# ifdef DIURNAL_INPUT_SRFLX
!$acc&, srflxbio
# endif
# ifndef ANA_SRFLUX
!$acc&, srflxg
!$acc&, srflxp, srf_time
#  ifdef DIURNAL_INPUT_SRFLX
!$acc&, srflxbiog
!$acc&, srflxbiop
#  endif /*  DIURNAL_INPUT_SRFLX   */
#  undef SRFLUX_DATA
# endif /* !ANA_SRFLUX */
#endif /* SOLVE3D */
#if defined BBL || defined MRL_WCI || defined OW_COUPLING
!$acc&, wfrq
!$acc&, wwkx
!$acc&, wwke
!$acc&, ubr
#endif
#ifdef BBL
!$acc&, uorb
!$acc&, vorb
#endif   /* BBL */
#if defined MRL_WCI || defined OW_COUPLING
!$acc&, whrm
!$acc&, wepb
!$acc&, wepd
!$acc&, wdrx
!$acc&, wdre
!$acc&, wlm
# ifdef WAVE_ROLLER
!$acc&, wepr
# endif
!$acc&, brk2dx
!$acc&, brk2de
!$acc&, ust2d
!$acc&, vst2d
!$acc&, frc2dx
!$acc&, frc2de
!$acc&, sup
!$acc&, ust_ext
# ifdef SOLVE3D
!$acc&, calP
!$acc&, Kapsrf
!$acc&, bhd
#  ifndef WAVE_SFC_BREAK
!$acc&, brk3dx
!$acc&, brk3de
#  endif
#  ifdef WAVE_STREAMING
!$acc&, frc3dx
!$acc&, frc3de
#  endif
!$acc&, ust
!$acc&, vst
!$acc&, wst
!$acc&, kvf
!$acc&, Akb
!$acc&, Akw
!$acc&, E_pre
# endif  /* SOLVE3D */
#endif   /* MRL_WCI */
#if defined BBL || defined MRL_WCI \
     ||  (defined MUSTANG && defined WAVE_OFFLINE)
!$acc&, Awave
!$acc&, Dwave
!$acc&, Pwave
# ifdef WAVE_OFFLINE
!$acc&, wwag
!$acc&, wwdg
!$acc&, wwpg
#  ifdef MUSTANG
!$acc&, Uwave
!$acc&, wwug
#  endif
!$acc&, wwfrq
#  ifdef BBL_OFFLINE
!$acc&, wwub
!$acc&, wwuob
!$acc&, wwvob
#  endif /* BBL_OFFLINE */
#  ifdef MRL_WCI
#   ifdef WAVE_OFFLINE_BREAKING
!$acc&, wweb
#   endif
#   ifdef WAVE_OFFLINE_FRICTION
!$acc&, wwed
#   endif
#   ifdef WAVE_OFFLINE_ROLLER
!$acc&, wwer
#   endif
#  endif /* MRL_WCI */
!$acc&, wwv_time
!$acc&, wwap, wwdp, wwpp
!$acc&, wwebp, wwedp, wwerp
#  ifdef MUSTANG
!$acc&, wwup
#  endif
# endif /* WAVE_OFFLINE */
#endif /* BBL || MRL_WCI */
#ifdef WAVE_MAKER
!$acc&, wf_bry, wk_bry, wa_bry
!$acc&, wd_bry, wa_bry_d
# ifdef WAVE_MAKER_DSPREAD
!$acc&, wpha_bry
# else
!$acc&, wpha_bry
# endif
#endif

!work.h
#ifdef SOLVE3D
!$acc&, work
!$acc&, workr
#endif
#ifdef ABL1D
!$acc&, work3dabl
#endif
!$acc&, work2d
!$acc&, work2d2

!ncscrum.h
#ifdef SOLVE3D
# ifdef TRACERS
#  ifdef PASSIVE_TRACER
!$acc&, indxTPAS
#  endif
# endif
# ifdef BIOLOGY
#  ifdef PISCES
#   ifdef key_pisces_quota
#    ifdef key_ligand
#    else
#    endif
#   endif
#  elif defined BIO_NChlPZD
#  elif defined BIO_N2ChlPZD2
#  elif defined BIO_BioEBUS
#  endif
# endif /* BIOLOGY */
# ifdef SEDIMENT
!$acc&, indxGRAV
!$acc&, indxSAND
!$acc&, indxMUD
# endif
# ifdef DIAGNOSTICS_TS
#  ifdef DIAGNOSTICS_TSVAR
#  else
#  endif
# endif
# ifdef DIAGNOSTICS_UV
# endif
# ifdef DIAGNOSTICS_VRT
# endif
# ifdef DIAGNOSTICS_EK
#  ifdef DIAGNOSTICS_EK_MLD
#  endif
# endif
# ifdef DIAGNOSTICS_PV
# endif
# if defined BIOLOGY && defined DIAGNOSTICS_BIO
# endif /* BIOLOGY && DIAGNOSTICS_BIO */
# ifdef MUSTANG
# endif
# ifdef BIOLOGY
#  ifdef BIO_BioEBUS
#  endif
# endif
#endif
#if defined BIOLOGY && !defined PISCES
# if (defined BIO_NChlPZD || defined BIO_N2ChlPZD2)
# endif
#endif /* BIOLOGY*/
#ifdef SOLVE3D
# if defined BIOLOGY && !defined PISCES
#  if (defined BIO_NChlPZD || defined BIO_N2ChlPZD2)
#   ifdef OXYGEN
#   else
#   endif
#  else
#  endif
# else
# endif
#else
# if defined BIOLOGY && !defined PISCES
#  ifdef BIO_NChlPZD
#   ifdef OXYGEN
#   else
#   endif
#  endif
# else
# endif
#endif /* SOLVE3D */
#if defined ANA_VMIX  || defined LMD_MIXING \
  || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
#endif
#ifdef EXACT_RESTART
# ifdef M3FAST
#  ifdef TS_MIX_ISO_FILT
#  else
#  endif
# else
#  ifdef TS_MIX_ISO_FILT
#  else
#  endif
# endif
#else
#endif
#ifdef SOLVE3D
# ifdef SALINITY
# else
# endif
# ifdef BULK_FLUX
# endif
#endif /* SOLVE3D */
#ifdef SOLVE3D
# ifdef SEDIMENT
!$acc&, indxBFRA
#  ifdef SUSPLOAD
!$acc&, indxDFLX
!$acc&, indxEFLX
#   ifdef BEDLOAD
!$acc&, indxBDLU
!$acc&, indxBDLV
#   endif
#  else
#   ifdef BEDLOAD
!$acc&, indxBDLU
!$acc&, indxBDLV
#   endif
#  endif
#  if defined MIXED_BED || defined COHESIVE_BED
#   if defined BEDLOAD && defined SUSPLOAD
#   else
#   endif
#  endif
# endif
#endif /* SOLVE3D */
#ifdef BBL
# ifdef SEDIMENT
# else
# endif
# ifndef ANA_WWAVE
# endif /* !ANA_WWAVE */
#else /* BBL */
#endif  /* BBL */
#if defined MRL_WCI || defined OW_COUPLING
# ifdef SEDIMENT
# else
# endif
#endif  /* MRL_WCI */
#ifdef PSOURCE_NCFILE
#endif /* PSOURCE_NCFILE */
#ifndef AGRIF
# if defined MPI && defined PARALLEL_FILES
# else
# endif
#else
# if defined MPI && defined PARALLEL_FILES
# else
# endif
#endif /* AGRIF */
#if defined SOLVE3D && defined TRACERS
!$acc&, nttclm, ntstf, nttsrc
!$acc&, ntbtf
#endif
#ifdef SOLVE3D
# if defined TRACERS
!$acc&, rstT
# endif
# if defined GLS_MIXING
# endif
# ifdef M3FAST
# endif
# ifdef SEDIMENT
!$acc&, rstSed
# endif
# ifdef MUSTANG
!$acc&, rstMUS
# endif
#endif
#ifdef EXACT_RESTART
#endif
#ifdef BBL
!$acc&, rstBBL
#endif
#ifdef WAVE_IO
!$acc&, rstWAVE, hisWAVE
#endif
#ifdef MRL_WCI
#endif
#ifdef BBL
!$acc&, hisBBL
#endif
#ifdef SOLVE3D
# if defined ANA_VMIX || defined LMD_MIXING \
  || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
# endif
# ifdef BIOLOGY
#  ifdef BIO_NChlPZD
#  elif defined BIO_BioEBUS
#  endif
# endif  /* BIOLOGY */
!$acc&, hisT
# ifdef SEDIMENT
# endif /* SEDIMENT */
# ifdef MUSTANG
# endif /* MUSTANG */
# if defined DIAGNOSTICS_TS
!$acc&, diaTXadv, diaTYadv, diaTVadv
!$acc&, diaTHmix, diaTVmix
#  ifdef DIAGNOSTICS_TSVAR
!$acc&, diaTVmixt
#  endif
!$acc&, diaTForc, diaTrate
#  if defined DIAGNOSTICS_TS_MLD
!$acc&, diaTXadv_mld, diaTYadv_mld, diaTVadv_mld
!$acc&, diaTHmix_mld, diaTVmix_mld
!$acc&, diaTForc_mld, diaTrate_mld, diaTentr_mld
#  endif
# endif
# ifdef DIAGNOSTICS_UV
!$acc&, diaMXadv, diaMYadv, diaMVadv
!$acc&, diaMCor, diaMPrsgrd, diaMHmix
!$acc&, diaMHdiff
!$acc&, diaMVmix, diaMVmix2, diaMrate
#  ifdef DIAGNOSTICS_BARO
!$acc&, diaMBaro
#  endif
#  ifdef M3FAST
!$acc&, diaMfast
#  endif
#  ifdef MRL_WCI
!$acc&, diaMvf, diaMbrk, diaMStCo
!$acc&, diaMVvf, diaMPrscrt, diaMsbk
!$acc&, diaMbwf, diaMfrc
#  endif
# endif
# ifdef DIAGNOSTICS_VRT
!$acc&, diags_vrtXadv, diags_vrtYadv, diags_vrtHdiff
!$acc&, diags_vrtCor, diags_vrtPrsgrd, diags_vrtHmix
!$acc&, diags_vrtVmix, diags_vrtrate
!$acc&, diags_vrtVmix2, diags_vrtWind, diags_vrtDrag
#  ifdef DIAGNOSTICS_BARO
!$acc&, diags_vrtBaro
#  endif
#  ifdef M3FAST
!$acc&, diags_vrtfast
#  endif
# endif
# ifdef DIAGNOSTICS_EK
!$acc&, diags_ekHadv, diags_ekHdiff, diags_ekVadv
!$acc&, diags_ekCor, diags_ekPrsgrd, diags_ekHmix
!$acc&, diags_ekVmix, diags_ekrate, diags_ekvol
!$acc&, diags_ekVmix2, diags_ekWind, diags_ekDrag
#  ifdef DIAGNOSTICS_BARO
!$acc&, diags_ekBaro
#  endif
#  ifdef M3FAST
!$acc&, diags_ekfast
#  endif
#  ifdef DIAGNOSTICS_EK_MLD
!$acc&, diags_ekHadv_mld, diags_ekHdiff_mld
!$acc&, diags_ekVadv_mld, diags_ekCor_mld
!$acc&, diags_ekPrsgrd_mld, diags_ekHmix_mld
!$acc&, diags_ekVmix_mld, diags_ekrate_mld
!$acc&, diags_ekvol_mld, diags_ekVmix2_mld
#  endif
#  ifdef DIAGNOSTICS_BARO
!$acc&, diags_ekBaro_mld
#  endif
# endif
# ifdef DIAGNOSTICS_PV
#  ifdef DIAGNOSTICS_PV_FULL
!$acc&, diags_pvpv, diags_pvpvd
#  endif
!$acc&, diags_pvMrhs, diags_pvTrhs
# endif
# if defined DIAGNOSTICS_EDDY && ! defined XIOS
!$acc&, diags_eddyzz
!$acc&, diags_eddyuu, diags_eddyvv, diags_eddyuv
!$acc&, diags_eddyub, diags_eddyvb, diags_eddywb
!$acc&, diags_eddyuw, diags_eddyvw
!$acc&, diags_eddyubu, diags_eddyvbv
!$acc&, diags_eddyusu, diags_eddyvsv
!$acc&, diags_eddyugsu, diags_eddyvgsv
# endif
# if defined OUTPUTS_SURFACE && ! defined XIOS
!$acc&, surf_surft, surf_surfs, surf_surfz
!$acc&, surf_surfu, surf_surfv
# endif
# ifdef DIAGNOSTICS_BIO
!$acc&, diabioFlux
!$acc&, diabioVSink
!$acc&, diabioGasExc
# endif
#elif defined DIAGNOSTICS_UV && defined MRL_WCI
!$acc&, diaMvf, diaMbrk, diaMStCo
!$acc&, diaMVvf, diaMPrscrt, diaMsbk
!$acc&, diaMbwf, diaMfrc
#endif /* SOLVE3D */
#ifdef AVERAGES
# ifdef SOLVE3D
#  if defined ANA_VMIX || defined LMD_MIXING \
   || defined LMD_SKPP || defined LMD_BKPP \
   || defined GLS_MIXING
#  endif
#  ifdef BIOLOGY
#   ifdef BIO_NChlPZD
#   elif defined BIO_BioEBUS
#   endif
#  endif  /* BIOLOGY */
#  if defined TRACERS
!$acc&, avgT
#  endif
#  ifdef SEDIMENT
#  endif
#  ifdef MUSTANG
!$acc&, avgMust
#  endif
# endif /* SOLVE3D */
# ifdef BBL
!$acc&, avgBBL
# endif
# ifdef WAVE_IO
!$acc&, avgWAVE
# endif
# ifdef MRL_WCI
# endif
# if defined SOLVE3D && defined TRACERS
#  if defined DIAGNOSTICS_TS && defined TRACERS
!$acc&, diaTXadv_avg, diaTYadv_avg, diaTVadv_avg
!$acc&, diaTHmix_avg, diaTVmix_avg
#   ifdef DIAGNOSTICS_TSVAR
!$acc&, diaTVmixt_avg
#   endif
!$acc&, diaTForc_avg, diaTrate_avg
#   ifdef DIAGNOSTICS_TS_MLD
!$acc&, diaTXadv_mld_avg, diaTYadv_mld_avg
!$acc&, diaTVadv_mld_avg
!$acc&, diaTHmix_mld_avg, diaTVmix_mld_avg
!$acc&, diaTForc_mld_avg, diaTrate_mld_avg
!$acc&, diaTentr_mld_avg
#   endif
#  endif
#  ifdef DIAGNOSTICS_UV
!$acc&, diaMXadv_avg, diaMYadv_avg, diaMVadv_avg
!$acc&, diaMCor_avg, diaMPrsgrd_avg, diaMHmix_avg
!$acc&, diaMHdiff_avg
!$acc&, diaMVmix_avg, diaMVmix2_avg, diaMrate_avg
#   ifdef DIAGNOSTICS_BARO
!$acc&, diaMBaro_avg
#   endif
#   ifdef M3FAST
!$acc&, diaMfast_avg
#   endif
#  endif
#  ifdef DIAGNOSTICS_VRT
!$acc&, diags_vrtXadv_avg, diags_vrtYadv_avg, diags_vrtHdiff_avg
!$acc&, diags_vrtCor_avg, diags_vrtPrsgrd_avg, diags_vrtHmix_avg
!$acc&, diags_vrtVmix_avg, diags_vrtrate_avg
!$acc&, diags_vrtVmix2_avg, diags_vrtWind_avg, diags_vrtDrag_avg
#   ifdef DIAGNOSTICS_BARO
!$acc&, diags_vrtBaro_avg
#   endif
#   ifdef M3FAST
!$acc&, diags_vrtfast_avg
#   endif
#  endif
#  ifdef DIAGNOSTICS_EK
!$acc&, diags_ekHadv_avg, diags_ekHdiff_avg, diags_ekVadv_avg
!$acc&, diags_ekCor_avg, diags_ekPrsgrd_avg, diags_ekHmix_avg
!$acc&, diags_ekVmix_avg, diags_ekrate_avg, diags_ekvol_avg
!$acc&, diags_ekVmix2_avg, diags_ekWind_avg, diags_ekDrag_avg
#   ifdef DIAGNOSTICS_BARO
!$acc&, diags_ekBaro_avg
#   endif
#   ifdef M3FAST
!$acc&, diags_ekfast_avg
#   endif
#   ifdef DIAGNOSTICS_EK_MLD
!$acc&, diags_ekHadv_mld_avg, diags_ekHdiff_mld_avg
!$acc&, diags_ekVadv_mld_avg, diags_ekCor_mld_avg
!$acc&, diags_ekPrsgrd_mld_avg, diags_ekHmix_mld_avg
!$acc&, diags_ekVmix_mld_avg, diags_ekrate_mld_avg
!$acc&, diags_ekvol_mld_avg, diags_ekVmix2_mld_avg
#   endif
#   ifdef DIAGNOSTICS_BARO
!$acc&, diags_ekBaro_mld_avg
#   endif
#  endif
#  ifdef DIAGNOSTICS_PV
#   ifdef DIAGNOSTICS_PV_FULL
!$acc&, diags_pvpv_avg, diags_pvpvd_avg
#   endif
!$acc&, diags_pvMrhs_avg, diags_pvTrhs_avg
#  endif
#  if defined DIAGNOSTICS_EDDY && ! defined XIOS
!$acc&, diags_eddyzz_avg
!$acc&, diags_eddyuu_avg, diags_eddyvv_avg, diags_eddyuv_avg
!$acc&, diags_eddyub_avg, diags_eddyvb_avg, diags_eddywb_avg
!$acc&, diags_eddyuw_avg, diags_eddyvw_avg
!$acc&, diags_eddyubu_avg, diags_eddyvbv_avg
!$acc&, diags_eddyusu_avg, diags_eddyvsv_avg
!$acc&, diags_eddyugsu_avg, diags_eddyvgsv_avg
#  endif
#  if defined OUTPUTS_SURFACE && ! defined XIOS
!$acc&, surf_surft_avg, surf_surfs_avg, surf_surfz_avg
!$acc&, surf_surfu_avg, surf_surfv_avg
#  endif
#  ifdef DIAGNOSTICS_BIO
!$acc&, diabioFlux_avg
!$acc&, diabioVSink_avg
!$acc&, diabioGasExc_avg
#  endif
# endif /* SOLVE3D */
# ifdef MRL_WCI
!$acc&, diaMvf_avg, diaMbrk_avg, diaMStCo_avg
!$acc&, diaMVvf_avg, diaMPrscrt_avg, diaMsbk_avg
!$acc&, diaMbwf_avg, diaMfrc_avg
# endif
#endif /* AVERAGES */
#ifdef SOLVE3D
# define NWRTHIS 1000+NT
#else
# define NWRTHIS 500
#endif
!$acc&, wrthis
#ifdef AVERAGES
!$acc&, wrtavg
#endif
#ifdef ABL1D
!$acc&, wrtabl
# ifdef AVERAGES
!$acc&, wrtabl_avg
# endif
#endif
#ifdef DIAGNOSTICS_TS
!$acc&, wrtdia3D
!$acc&, wrtdia2D
# ifdef AVERAGES
!$acc&, wrtdia3D_avg
!$acc&, wrtdia2D_avg
# endif
#endif
#ifdef DIAGNOSTICS_UV
!$acc&, wrtdiaM
# ifdef AVERAGES
!$acc&, wrtdiaM_avg
# endif
#endif
#ifdef DIAGNOSTICS_VRT
!$acc&, wrtdiags_vrt
# ifdef AVERAGES
!$acc&, wrtdiags_vrt_avg
# endif
#endif
#ifdef DIAGNOSTICS_EK
!$acc&, wrtdiags_ek
# ifdef AVERAGES
!$acc&, wrtdiags_ek_avg
# endif
#endif
#ifdef DIAGNOSTICS_PV
!$acc&, wrtdiags_pv
# ifdef AVERAGES
!$acc&, wrtdiags_pv_avg
# endif
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
!$acc&, wrtdiags_eddy
# ifdef AVERAGES
!$acc&, wrtdiags_eddy_avg
# endif
#endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
!$acc&, wrtsurf
# ifdef AVERAGES
!$acc&, wrtsurf_avg
# endif
#endif
#ifdef DIAGNOSTICS_BIO
!$acc&, wrtdiabioFlux
!$acc&, wrtdiabioVSink
!$acc&, wrtdiabioGasExc
# ifdef AVERAGES
!$acc&, wrtdiabioFlux_avg
!$acc&, wrtdiabioVSink_avg
!$acc&, wrtdiabioGasExc_avg
# endif
#endif
#ifdef SOLVE3D
# if defined GLS_MIXING
# endif
# ifdef M3FAST
# endif
# ifdef EXACT_RESTART
# endif
#endif
#ifdef SOLVE3D
# if defined ANA_VMIX || defined LMD_MIXING \
  || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
# endif
# ifdef BIOLOGY
#  ifdef BIO_NChlPZD
#  elif defined BIO_BioEBUS
#  endif
# endif  /* BIOLOGY */
#endif
#ifdef DIAGNOSTICS_TS
# ifdef AVERAGES
# endif
#endif
#ifdef DIAGNOSTICS_UV
# ifdef AVERAGES
# endif
#endif
#ifdef DIAGNOSTICS_VRT
# ifdef AVERAGES
# endif
#endif
#ifdef DIAGNOSTICS_EK
# ifdef AVERAGES
# endif
# ifdef DIAGNOSTICS_EK_MLD
#  ifdef AVERAGES
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_PV
# ifdef AVERAGES
# endif
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
#endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
#endif
#ifdef DIAGNOSTICS_BIO
#endif
#ifdef AVERAGES
# ifdef SOLVE3D
#  if defined ANA_VMIX || defined LMD_MIXING \
  || defined LMD_SKPP  || defined LMD_BKPP \
  || defined GLS_MIXING
#  endif
#  ifdef BIOLOGY
#   ifdef BIO_NChlPZD
#   elif defined BIO_BioEBUS
#   endif
#  endif  /* BIOLOGY */
# endif /* SOLVE3D */
# ifdef MRL_WCI
# endif
#endif /* AVERAGES */
#ifdef DIAGNOSTICS_TS
#endif
#ifdef DIAGNOSTICS_UV
#endif
#ifdef DIAGNOSTICS_VRT
#endif
#ifdef DIAGNOSTICS_EK
#endif
#ifdef DIAGNOSTICS_PV
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
#endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
#endif
#ifdef DIAGNOSTICS_BIO
#endif
#ifdef ABL1D
#endif
#ifdef DIAGNOSTICS_TS
#endif
#ifdef DIAGNOSTICS_UV
#endif
#ifdef DIAGNOSTICS_VRT
#endif
#ifdef DIAGNOSTICS_EK
#endif
#ifdef DIAGNOSTICS_PV
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
#endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
#endif
#ifdef DIAGNOSTICS_BIO
#endif
#if (defined TCLIMATOLOGY  && !defined ANA_TCLIMA)\
 || (defined ZCLIMATOLOGY  && !defined ANA_SSH)\
 || (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA)\
 || (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
#endif
#ifdef SEDIMENT
#elif defined MUSTANG
#endif
#ifdef SOLVE3D
!$acc&, vname
#else
!$acc&, vname
#endif
#ifdef ABL1D
!$acc&, vname_abl1d
#endif
#ifdef ABL1D
#endif
#ifdef DIAGNOSTICS_TS
#endif
#ifdef DIAGNOSTICS_UV
#endif
#ifdef DIAGNOSTICS_VRT
#endif
#ifdef DIAGNOSTICS_EK
#endif
#ifdef DIAGNOSTICS_PV
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
#endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
#endif
#ifdef DIAGNOSTICS_BIO
#endif
#if (defined TCLIMATOLOGY  && !defined ANA_TCLIMA)\
 || (defined ZCLIMATOLOGY  && !defined ANA_SSH)\
 || (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA)\
 || (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
#endif
#ifdef SEDIMENT
#elif defined MUSTANG
#endif

!averages.h
#ifdef AVERAGES
!$acc&, zeta_avg
!$acc&, ubar_avg
!$acc&, vbar_avg
!$acc&, bostr_avg
!$acc&, bustr_avg
!$acc&, bvstr_avg
!$acc&, wstr_avg
!$acc&, sustr_avg
!$acc&, svstr_avg
!$acc&, srflx_avg
# ifdef MORPHODYN
!$acc&, h_avg
# endif
# ifdef SOLVE3D
!$acc&, u_avg
!$acc&, v_avg
!$acc&, t_avg
!$acc&, rho_avg
#  if defined ANA_VMIX || defined LMD_MIXING \
   || defined LMD_SKPP || defined LMD_BKPP \
   || defined GLS_MIXING
!$acc&, bvf_avg
#  endif
!$acc&, omega_avg
#  ifdef NBQ
!$acc&, w_avg
#  else
!$acc&, w_avg
#  endif
#  if defined ANA_VMIX || defined LMD_MIXING \
   || defined LMD_SKPP || defined LMD_BKPP \
   || defined GLS_MIXING
#  endif
!$acc&, stflx_avg
!$acc&, btflx_avg
#  if defined LMD_SKPP || defined GLS_MIXING
!$acc&, hbl_avg
#  endif
#  ifdef LMD_BKPP
!$acc&, hbbl_avg
#  endif
#  ifdef GLS_MIXING
!$acc&, tke_avg
!$acc&, gls_avg
!$acc&, Lscale_avg
#  endif
#  ifdef BULK_FLUX
!$acc&, shflx_rsw_avg
!$acc&, shflx_rlw_avg
!$acc&, shflx_lat_avg
!$acc&, shflx_sen_avg
#  endif
#  ifdef SST_SKIN
!$acc&, sst_skin_avg
#  endif
#  ifdef BIOLOGY
!$acc&, hel_avg
#   ifdef BIO_NChlPZD
!$acc&, theta_avg
#    ifdef OXYGEN
!$acc&, u10_avg
!$acc&, Kv_O2_avg
!$acc&, O2satu_avg
#    endif /* OXYGEN */
#   elif defined BIO_BioEBUS
!$acc&, AOU_avg
!$acc&, wind10_avg
#   endif
#  endif /* BIOLOGY */
#  ifdef VIS_COEF_3D
!$acc&, visc3d_avg
#  endif
#  ifdef DIF_COEF_3D
!$acc&, diff3d_avg
#  endif
#  ifdef AVERAGES_K
!$acc&, Akv_avg
!$acc&, Akt_avg
#  endif
# endif
# ifdef WAVE_IO
!$acc&, whrm_avg
!$acc&, wfrq_avg
!$acc&, wac_avg
!$acc&, wkx_avg
!$acc&, wke_avg
!$acc&, wepb_avg
!$acc&, wepd_avg
#  ifdef WAVE_ROLLER
!$acc&, war_avg
!$acc&, wepr_avg
#  endif
# endif
# ifdef MRL_WCI
!$acc&, sup_avg
!$acc&, ust2d_avg
!$acc&, vst2d_avg
#  ifdef SOLVE3D
!$acc&, ust_avg
!$acc&, vst_avg
!$acc&, wst_avg
!$acc&, akb_avg
!$acc&, akw_avg
!$acc&, kvf_avg
!$acc&, calp_avg
!$acc&, kaps_avg
#  endif
# endif
#endif /* AVERAGES */

!lmd_kpp.h
#if defined LMD_SKPP || defined LMD_BKPP || defined GLS_MIXING
!$acc&, Jwtype
#endif

!climat.h
#if defined ZCLIMATOLOGY || defined AGRIF
!$acc&, ssh
#endif
#ifdef ZCLIMATOLOGY
# ifdef ZNUDGING
!$acc&, Znudgcof
# endif
# ifndef ANA_SSH
!$acc&, sshg
!$acc&, ssh_time
#  undef SSH_DATA
# endif /* !ANA_SSH */
#endif
#ifdef SOLVE3D
# if defined TRACERS && (defined TCLIMATOLOGY || (defined AGRIF && !defined T_FRC_BRY))
!$acc&, tclm
# endif
# if defined TRACERS && defined TCLIMATOLOGY
#  ifdef TNUDGING
!$acc&, Tnudgcof
#  endif
#  ifndef ANA_TCLIMA
!$acc&, tclima
!$acc&, tclm_time
!$acc&, tclm_cycle
!$acc&, ittclm, tclm_ncycle, tclm_rec
!$acc&, tclm_tid, tclm_id
!$acc&, got_tclm
#   undef TCLIMA_DATA
#  endif /* !ANA_TCLIMA */
# endif /* TCLIMATOLOGY */
#endif /* SOLVE3D */
#if defined M2CLIMATOLOGY || (defined AGRIF && !defined M2_FRC_BRY)
!$acc&, ubclm
!$acc&, vbclm
#endif
#if defined SOLVE3D && (defined M3CLIMATOLOGY || \
                        (defined AGRIF && !defined M3_FRC_BRY))
!$acc&, uclm
!$acc&, vclm
#endif
#ifdef M2CLIMATOLOGY
# ifdef M2NUDGING
!$acc&, M2nudgcof
# endif
# ifndef ANA_M2CLIMA
!$acc&, ubclima
!$acc&, vbclima
# endif
#endif
#if defined SOLVE3D && defined M3CLIMATOLOGY
# ifdef M3NUDGING
!$acc&, M3nudgcof
# endif
# ifndef ANA_M3CLIMA
!$acc&, uclima
!$acc&, vclima
# endif
#endif
#if defined M2CLIMATOLOGY || defined M3CLIMATOLOGY
!$acc&, uclm_time
#endif
#ifdef ZONAL_NUDGING
# define GLOBAL_1D_ETA 0:Mm+1
!$acc&, zetazon
!$acc&, ubzon
!$acc&, vbzon
!$acc&, uzon
!$acc&, vzon
!$acc&, sshzon
!$acc&, ubclmzon
!$acc&, vbclmzon
!$acc&, uclmzon
!$acc&, vclmzon
# ifdef TRACERS
!$acc&, tzon
!$acc&, tclmzon
# endif
# undef GLOBAL_1D_ETA
#endif
#if defined M3FAST && (defined NBQCLIMATOLOGY || \
                   (defined AGRIF && !defined NBQ_FRC_BRY))
!$acc&, unbqclm
!$acc&, vnbqclm
# ifdef NBQ
!$acc&, wnbqclm
!$acc&, rnbqclm
# endif
#endif

!nbq.h
#ifdef M3FAST
!$acc&, soundspeed_nbq
!$acc&, soundspeed2_nbq
!$acc&, qdmu_nbq
!$acc&, qdmv_nbq
# ifdef NBQ
!$acc&, qdmw_nbq
# endif
# ifdef NBQ
!$acc&, thetadiv_nbq
!$acc&, thetadiv2_nbq
!$acc&, thetadiv3_nbq
# endif
!$acc&, ru_int_nbq
!$acc&, rv_int_nbq
!$acc&, ru_nbq
!$acc&, rv_nbq
!$acc&, ru_nbq_avg2
!$acc&, rv_nbq_avg2
# ifdef NBQ
!$acc&, rw_int_nbq
!$acc&, rw_nbq
!$acc&, rw_nbq_avg2
!$acc&, rho_nbq
# endif
!$acc&, DU_nbq
!$acc&, DV_nbq
!$acc&, ru_int_nbq_2d
!$acc&, rv_int_nbq_2d
# ifdef NBQ
!$acc&, rho_grd
!$acc&, rho_bak
#  ifdef NBQ_MASS
!$acc&, rho_nbq_avg1
!$acc&, rhobar_nbq
!$acc&, rhobar_nbq_avg1
#  endif
# else
!$acc&, rubar_nbq
!$acc&, rvbar_nbq
!$acc&, rubar_sum
!$acc&, rvbar_sum
# endif
!$acc&, Hzw_half_nbq
# ifdef NBQ
!$acc&, zw_nbq
# endif
# ifdef NBQ_HZCORRECT
!$acc&, Hz_correct
# endif
# ifdef NBQ_HZ_PROGNOSTIC
!$acc&, Hz_bak2
# endif
# ifdef NBQ
#  ifdef NBQ_GRID_SLOW
!$acc&, dthetadiv_nbqdz
!$acc&, dZdxq_w
!$acc&, dZdyq_w
#  else
!$acc&, dthetadiv_nbqdz
!$acc&, dZdxq_w
!$acc&, dZdyq_w
#  endif /* NBQ_GRID_SLOW */
# endif
# ifdef NBQ
!$acc&, wsurf_nbq
!$acc&, usurf_nbq
!$acc&, vsurf_nbq
# endif
# if defined OBC_NBQ && defined OBC_NBQORLANSKI
#  ifdef OBC_COM_WEST
!$acc&, qdmu_nbq_west
!$acc&, qdmv_nbq_west
#   ifdef NBQ
!$acc&, qdmw_nbq_west
!$acc&, rho_nbq_west
#   endif
#  endif
#  ifdef OBC_COM_EAST
!$acc&, qdmu_nbq_east
!$acc&, qdmv_nbq_east
#   ifdef NBQ
!$acc&, qdmw_nbq_east
!$acc&, rho_nbq_east
#   endif
#  endif
#  ifdef OBC_COM_SOUTH
!$acc&, qdmu_nbq_south
!$acc&, qdmv_nbq_south
#   ifdef NBQ
!$acc&, qdmw_nbq_south
!$acc&, rho_nbq_south
#   endif
#  endif
#  ifdef OBC_COM_NORTH
!$acc&, qdmu_nbq_north
!$acc&, qdmv_nbq_north
#   ifdef NBQ
!$acc&, qdmw_nbq_north
!$acc&, rho_nbq_north
#   endif
#  endif
# endif
# ifdef NBQ_NUDGING
!$acc&, NBQnudgcof
# endif
#endif /* M3FAST */

!sources.h
#if defined PSOURCE || defined PSOURCE_MASS || defined PSOURCE_NCFILE
!$acc&, Qbar0
!$acc&, Qbar
!$acc&, Qsrc
!$acc&, Qshape
# if defined TRACERS
!$acc&, Tsrc
!$acc&, Tsrc0
# endif
!$acc&, lasrc
!$acc&, losrc
!$acc&, Dsrc
!$acc&, Isrc
!$acc&, Jsrc
# if defined TRACERS
!$acc&, Lsrc
# endif
# ifdef PSOURCE_NCFILE
!$acc&, qbarg
!$acc&, qbar_time
!$acc&, qbardir
#  if defined PSOURCE_NCFILE_TS && defined TRACERS
!$acc&, tsrcg
!$acc&, tsrc_time
!$acc&, tsrc_cycle
!$acc&, ittsrc, tsrc_ncycle, tsrc_rec, tsrc_tid
!$acc&, tsrc_id
!$acc&, got_tsrc
#  endif
# endif /* PSOURCE_NCFILE */
# ifdef MPI
!$acc&, Isrc_mpi
!$acc&, Jsrc_mpi
# endif
#endif

!wkb_wwave.h
#ifdef WKB_WWAVE
!$acc&, wkx
!$acc&, wke
!$acc&, wac
!$acc&, hrm
!$acc&, frq
!$acc&, wcg
!$acc&, wsb
!$acc&, wvn
!$acc&, wfc
# ifdef WAVE_ROLLER
!$acc&, war
!$acc&, wcr
!$acc&, wsr
# endif
# if defined MRL_CEW || !defined WKB_STEADY
#  ifdef WKB_TIME_FILTER
!$acc&, uwave
!$acc&, vwave
!$acc&, zwave
#  else
!$acc&, uwave
!$acc&, vwave
!$acc&, zwave
#  endif
# endif
#endif /* WKB_WWAVE */

!boundary.h
#ifdef T_FRC_BRY
!$acc&, got_tbry
#endif
#ifndef ANA_BRY
!$acc&, bry_time
# if defined BIOLOGY || defined PISCES
!$acc&, bry_time1
!$acc&, bry_cycle1
!$acc&, bry_tid, bry_ncycle1
!$acc&, bry_rec1, itbry1, ntbry1
# endif
# if defined OBC_WEST || defined AGRIF_OBC_WEST
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY
!$acc&, tbry_west_id
#   endif
#  endif
# endif
# if defined OBC_EAST || defined AGRIF_OBC_EAST
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY
!$acc&, tbry_east_id
#   endif
#  endif
# endif
# if defined OBC_SOUTH || defined AGRIF_OBC_SOUTH
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY
!$acc&, tbry_south_id
#   endif
#  endif
# endif
# if defined OBC_NORTH || defined AGRIF_OBC_NORTH
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY
!$acc&, tbry_north_id
#   endif
#  endif
# endif
#endif  /* ANA_BRY */
#if defined OBC_WEST || defined AGRIF_OBC_WEST
# ifdef Z_FRC_BRY
!$acc&, zetabry_west
!$acc&, zetabry_west_dt
# endif
# ifdef M2_FRC_BRY
!$acc&, ubarbry_west
!$acc&, ubarbry_west_dt
!$acc&, vbarbry_west
!$acc&, vbarbry_west_dt
# endif
# ifdef SOLVE3D
#  ifdef M3_FRC_BRY
!$acc&, ubry_west
!$acc&, ubry_west_dt
!$acc&, vbry_west
!$acc&, vbry_west_dt
#  endif
#  ifdef T_FRC_BRY
!$acc&, tbry_west
!$acc&, tbry_west_dt
#   ifdef WAVE_MAKER_INTERNAL
!$acc&, t0bry_west
!$acc&, bvf0bry_west
#   endif
#  endif
# endif
#endif
#if defined OBC_EAST || defined AGRIF_OBC_EAST
# ifdef Z_FRC_BRY
!$acc&, zetabry_east
!$acc&, zetabry_east_dt
# endif
# ifdef M2_FRC_BRY
!$acc&, ubarbry_east
!$acc&, ubarbry_east_dt
!$acc&, vbarbry_east
!$acc&, vbarbry_east_dt
# endif
# ifdef SOLVE3D
#  ifdef M3_FRC_BRY
!$acc&, ubry_east
!$acc&, ubry_east_dt
!$acc&, vbry_east
!$acc&, vbry_east_dt
#  endif
#  ifdef T_FRC_BRY
!$acc&, tbry_east
!$acc&, tbry_east_dt
#  endif
# endif
#endif
#if defined OBC_SOUTH || defined AGRIF_OBC_SOUTH
# ifdef Z_FRC_BRY
!$acc&, zetabry_south
!$acc&, zetabry_south_dt
# endif
# ifdef M2_FRC_BRY
!$acc&, ubarbry_south
!$acc&, ubarbry_south_dt
!$acc&, vbarbry_south
!$acc&, vbarbry_south_dt
# endif
# ifdef SOLVE3D
#  ifdef M3_FRC_BRY
!$acc&, ubry_south
!$acc&, ubry_south_dt
!$acc&, vbry_south
!$acc&, vbry_south_dt
#  endif
#  ifdef T_FRC_BRY
!$acc&, tbry_south
!$acc&, tbry_south_dt
#  endif
# endif
#endif
#if defined OBC_NORTH || defined AGRIF_OBC_NORTH
# ifdef Z_FRC_BRY
!$acc&, zetabry_north
!$acc&, zetabry_north_dt
# endif
# ifdef M2_FRC_BRY
!$acc&, ubarbry_north
!$acc&, ubarbry_north_dt
!$acc&, vbarbry_north
!$acc&, vbarbry_north_dt
# endif
# ifdef SOLVE3D
#  ifdef M3_FRC_BRY
!$acc&, ubry_north
!$acc&, ubry_north_dt
!$acc&, vbry_north
!$acc&, vbry_north_dt
#  endif
#  ifdef T_FRC_BRY
!$acc&, tbry_north
!$acc&, tbry_north_dt
#  endif
# endif
#endif
#ifdef WKB_WWAVE
# ifndef ANA_BRY_WKB
!$acc&, brywkb_time
# endif  /* ANA_BRY_WKB */
# if defined WKB_OBC_WEST || defined AGRIF_OBC_WEST
!$acc&, wacbry_west
!$acc&, wacbry_west_dt
!$acc&, wkxbry_west
!$acc&, wkxbry_west_dt
!$acc&, wkebry_west
!$acc&, wkebry_west_dt
!$acc&, warbry_west
!$acc&, warbry_west_dt
# endif
# if defined WKB_OBC_EAST || defined AGRIF_OBC_EAST
!$acc&, wacbry_east
!$acc&, wacbry_east_dt
!$acc&, wkxbry_east
!$acc&, wkxbry_east_dt
!$acc&, wkebry_east
!$acc&, wkebry_east_dt
!$acc&, warbry_east
!$acc&, warbry_east_dt
# endif
# if defined WKB_OBC_SOUTH || defined AGRIF_OBC_SOUTH
!$acc&, wacbry_south
!$acc&, wacbry_south_dt
!$acc&, wkxbry_south
!$acc&, wkxbry_south_dt
!$acc&, wkebry_south
!$acc&, wkebry_south_dt
!$acc&, warbry_south
!$acc&, warbry_south_dt
# endif
# if defined WKB_OBC_NORTH || defined AGRIF_OBC_NORTH
!$acc&, wacbry_north
!$acc&, wacbry_north_dt
!$acc&, wkxbry_north
!$acc&, wkxbry_north_dt
!$acc&, wkebry_north
!$acc&, wkebry_north_dt
!$acc&, warbry_north
!$acc&, warbry_north_dt
# endif
#endif  /* WKB_WWAVE */
#if defined M3FAST && defined NBQ_FRC_BRY
# if defined OBC_WEST || defined AGRIF_OBC_WEST
!$acc&, unbqbry_west
!$acc&, vnbqbry_west
#  ifdef NBQ
!$acc&, wnbqbry_west
!$acc&, rnbqbry_west
#  endif
# endif
# if defined OBC_EAST || defined AGRIF_OBC_EAST
!$acc&, unbqbry_east
!$acc&, vnbqbry_east
#  ifdef NBQ
!$acc&, wnbqbry_east
!$acc&, rnbqbry_east
#  endif
# endif
# if defined OBC_SOUTH || defined AGRIF_OBC_SOUTH
!$acc&, unbqbry_south
!$acc&, vnbqbry_south
#  ifdef NBQ
!$acc&, wnbqbry_south
!$acc&, rnbqbry_south
#  endif
# endif
# if defined OBC_NORTH || defined AGRIF_OBC_NORTH
!$acc&, unbqbry_north
!$acc&, vnbqbry_north
#  ifdef NBQ
!$acc&, wnbqbry_north
!$acc&, rnbqbry_north
#  endif
# endif
#endif /* M3FAST */
#if defined NBQ && (defined W_FRC_BRY || defined AGRIF )
# if defined OBC_WEST || defined AGRIF_OBC_WEST
!$acc&, wbry_west
# endif
# if defined OBC_EAST || defined AGRIF_OBC_EAST
!$acc&, wbry_east
# endif
# if defined OBC_SOUTH || defined AGRIF_OBC_SOUTH
!$acc&, wbry_south
# endif
# if defined OBC_NORTH || defined AGRIF_OBC_NORTH
!$acc&, wbry_north
# endif
#endif /* NBQ */

!tides.h
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
!$acc&, Tperiod
#endif /* SSH_TIDES || UV_TIDES */
#if defined SSH_TIDES && (defined ZCLIMATOLOGY || defined Z_FRC_BRY)
!$acc&, SSH_Tamp
!$acc&, SSH_Tphase
#endif
#if defined UV_TIDES && (defined M2CLIMATOLOGY || defined M2_FRC_BRY)
!$acc&, UV_Tangle
!$acc&, UV_Tmajor
!$acc&, UV_Tminor
!$acc&, UV_Tphase
#endif
#ifdef POT_TIDES
!$acc&, POT_Tamp
!$acc&, POT_Tphase
!$acc&, PTide
#endif

!bbl.h
#if defined BBL || defined SEDIMENT
!$acc&, Abed
!$acc&, Hripple
!$acc&, Lripple
!$acc&, w_set
!$acc&, Sdens
!$acc&, Ssize
!$acc&, taucb
!$acc&, Zbnot
!$acc&, Zbapp
!$acc&, bustrw
!$acc&, bvstrw
#endif

!$acc& )
