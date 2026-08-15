#ifdef XIOS
#ifdef key_pisces
#endif
#if defined ONLINE_ANALYSIS
#endif    
# ifdef MPI
# else
# endif
# ifdef THREE_GHOST_POINTS
#  ifdef MPI
#  endif
# else
#  ifdef MPI
#  endif
# endif
#ifdef MASKING
# ifdef MPI
# else
# endif
#endif
#if defined ONLINE_ANALYSIS
# ifndef XIOS2
# else
# endif
#endif  /* ONLINE_ANALYSIS */
#if defined ONLINE_ANALYSIS
# ifndef XIOS2
# else
# endif
#endif  /* ONLINE_ANALYSIS */
#ifndef ANA_GRID
#endif
#ifndef ANA_INITIAL
#endif
#if !defined ANA_SMFLUX || !defined ANA_STFLUX 
#endif
#ifdef PSOURCE_NCFILE
#endif
#ifdef ASSIMILATION
#endif
#ifdef SOLVE3D
# ifdef  NEW_S_COORD
# endif
# ifdef  LMD_SKPP2005
# endif
#endif
#ifdef UV_VIS2
#endif
#ifdef UV_VIS4
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
# endif
# ifdef TS_DIF4
# endif
# if !defined LMD_MIXING
# endif
#endif
#ifdef SOLVE3D
# ifndef NONLIN_EOS
# endif
# ifdef BODYFORCE
# endif
#endif /* SOLVE3D */
# ifdef SPONGE
# endif
# ifdef SEDIMENT
# endif
# ifdef MPI
# endif
# ifdef SPHERICAL
#  ifdef MASKING
#  endif
#  ifdef MASKING
#  endif
#  ifdef MASKING
#  endif
# else
#  ifdef MASKING
#  endif
# ifdef MASKING
# endif
#  ifdef MASKING
#  endif
# endif
# ifdef SPHERICAL
# endif
# ifdef MPI
# else
# endif
# ifdef USE_CALENDAR
# endif
#if defined MUSTANG
# endif
#if defined ONLINE_ANALYSIS
#ifdef MPI
#endif
#endif
#ifdef USE_CALENDAR
#endif
# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
# ifdef AGRIF
# endif
# if defined MPI
# else
# endif
# ifdef USE_CALENDAR
# endif
# ifdef SPHERICAL
# else
# endif
# ifdef MASKING
# endif
# ifdef SEDIMENT
# elif defined MUSTANG
# else
# endif
#ifdef ONLINE_ANALYSIS
#ifdef MPI
#else
#endif /* MPI */
#endif  /* ONLINE_ANALYSIS */
#ifdef key_pisces
#endif
# ifdef AGRIF
# else
# endif
# ifdef SOLVE3D
# endif
#if defined BIOLOGY && defined PISCES
#endif
#if defined SUBSTANCE
#endif
#if defined MUSTANG
# ifdef WAVE_OFFLINE
# endif
#ifdef key_sand2D
#endif
#endif
#if defined MUSTANG
#endif
# ifdef SEDIMENT
# endif
# ifdef BBL
# endif
# ifdef WKB_WWAVE
# endif
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_TS
# endif
# if defined DIAGNOSTICS_VRT
# endif
# ifdef DIAGNOSTICS_EK
# endif
# ifdef DIAGNOSTICS_PV
# endif
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
# ifdef K3FAST
# endif
#ifdef DIAGNOSTICS_EDDY
#endif
#if defined MUSTANG
#endif
#  define I_EXT_RANGE Istr-1,Iend+1
#  define J_EXT_RANGE Jstr-1,Jend+1
#  ifdef EW_PERIODIC 
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
#  else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
#  endif
#  ifdef NS_PERIODIC 
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
#  else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
#  endif
#ifdef OPENACC
#endif
c       if (xios_field_is_active("time")) then
c !$acc update host( time ) 
c       endif
# ifdef SOLVE3D
      if (xios_field_is_active("Cs_r")) then
!$acc update host( Cs_r ) 
      endif
      if (xios_field_is_active("Cs_w")) then
!$acc update host( Cs_w ) 
      endif
      if (xios_field_is_active("sc_r")) then
!$acc update host( sc_r ) 
      endif
      if (xios_field_is_active("sc_w")) then
!$acc update host( sc_w ) 
      endif
c       if (xios_field_is_active("hc")) then
c !$acc update host( hc ) 
c       endif
c       if (xios_field_is_active("theta_s")) then
c !$acc update host( theta_s ) 
c       endif
c       if (xios_field_is_active("theta_b")) then
c !$acc update host( theta_b ) 
c       endif
c       if (xios_field_is_active("Tcline")) then
c !$acc update host( Tcline ) 
c       endif
c       if (xios_field_is_active("Vtransform")) then
c !$acc update host( Vtransform ) 
c       endif
      if (xios_field_is_active("levels_rho")) then
!$acc update host( z_r ) 
      endif
      if (xios_field_is_active("levels_w")) then
!$acc update host( z_w ) 
      endif
c       if (xios_field_is_active("levels_u")) then
c !$acc update host( levels_u ) 
c       endif
c       if (xios_field_is_active("levels_v")) then
c !$acc update host( levels_v ) 
c       endif
# endif
      if (xios_field_is_active("h")) then
!$acc update host( h ) 
      endif
      if (xios_field_is_active("f")) then
!$acc update host( f ) 
      endif
      if (xios_field_is_active("pm")) then
!$acc update host( pm ) 
      endif
      if (xios_field_is_active("pn")) then
!$acc update host( pn ) 
      endif
# ifdef SPHERICAL
      if (xios_field_is_active("lon_rho")) then
!$acc update host( lonr ) 
      endif
      if (xios_field_is_active("lat_rho")) then
!$acc update host( latr ) 
      endif
      if (xios_field_is_active("lon_u")) then
!$acc update host( lonu ) 
      endif
      if (xios_field_is_active("lat_u")) then
!$acc update host( latu ) 
      endif
      if (xios_field_is_active("lon_v")) then
!$acc update host( lonv ) 
      endif
      if (xios_field_is_active("lat_v")) then
!$acc update host( latv ) 
      endif
# else  
      if (xios_field_is_active("x_rho")) then
!$acc update host( xr ) 
      endif
      if (xios_field_is_active("y_rho")) then
!$acc update host( yr ) 
      endif
# endif
# ifdef CURVGRID
      if (xios_field_is_active("angle")) then
!$acc update host( angler ) 
      endif
# endif
# ifdef MASKING
      if (xios_field_is_active("mask_rho")) then
!$acc update host( rmask ) 
      endif
# endif
# ifdef WET_DRY
      if (xios_field_is_active("Dcrit")) then
!$acc update host( Dcrit ) 
      endif
# endif
      if (xios_field_is_active("zeta")) then
!$acc update host( zeta ) 
      endif
      if (xios_field_is_active("ubar")) then
!$acc update host( ubar ) 
      endif
      if (xios_field_is_active("vbar")) then
!$acc update host( vbar ) 
      endif
# ifdef MORPHODYN
      if (xios_field_is_active("hmorph")) then
!$acc update host( hmorph ) 
      endif
# endif
# ifdef WET_DRY
      if (xios_field_is_active("rmask_wet")) then
!$acc update host( rmask_wet ) 
      endif
      if (xios_field_is_active("umask_wet")) then
!$acc update host( umask_wet ) 
      endif
      if (xios_field_is_active("vmask_wet")) then
!$acc update host( vmask_wet ) 
      endif
#endif
      if (xios_field_is_active("bustr")) then
!$acc update host( bustr ) 
      endif
      if (xios_field_is_active("bvstr")) then
!$acc update host( bvstr ) 
      endif
      if (xios_field_is_active("bostr")) then
!$acc update host( bustr, bvstr ) 
      endif
      if (xios_field_is_active("wstr")) then
!$acc update host( sustr, svstr ) 
      endif
      if (xios_field_is_active("sustr")) then
!$acc update host( sustr ) 
      endif
      if (xios_field_is_active("svstr")) then
!$acc update host( svstr ) 
      endif
# ifdef BULK_SM_UPDATE
      if (xios_field_is_active("uwnd")) then
!$acc update host( uwnd ) 
      endif
      if (xios_field_is_active("vwnd")) then
!$acc update host( vwnd ) 
      endif
# endif
# ifdef WKB_WWAVE
      if (xios_field_is_active("hrm")) then
!$acc update host( hrm ) 
      endif
      if (xios_field_is_active("frq")) then
!$acc update host( frq ) 
      endif
      if (xios_field_is_active("wac")) then
!$acc update host( wac ) 
      endif
      if (xios_field_is_active("wkx")) then
!$acc update host( wkx ) 
      endif
      if (xios_field_is_active("wke")) then
!$acc update host( wke ) 
      endif
      if (xios_field_is_active("epb")) then
!$acc update host( epb ) 
      endif
      if (xios_field_is_active("epd")) then
!$acc update host( epd ) 
      endif
#  ifdef WAVE_ROLLER
      if (xios_field_is_active("war")) then
!$acc update host( war ) 
      endif
      if (xios_field_is_active("epr")) then
!$acc update host( epr ) 
      endif
#  endif
# endif
# ifdef MRL_WCI
      if (xios_field_is_active("sup")) then
!$acc update host( sup ) 
      endif
      if (xios_field_is_active("ust2d")) then
!$acc update host( ust2d ) 
      endif
      if (xios_field_is_active("vst2d")) then
!$acc update host( vst2d ) 
      endif
#  ifdef SOLVE3D
      if (xios_field_is_active("ust")) then
!$acc update host( ust ) 
      endif
      if (xios_field_is_active("vst")) then
!$acc update host( vst ) 
      endif
      if (xios_field_is_active("wst")) then
!$acc update host( wst ) 
      endif
      if (xios_field_is_active("Akb")) then
!$acc update host( Akb ) 
      endif
      if (xios_field_is_active("Akw")) then
!$acc update host( Akw ) 
      endif
      if (xios_field_is_active("kvf")) then
!$acc update host( kvf ) 
      endif
      if (xios_field_is_active("calP")) then
!$acc update host( calP ) 
      endif
      if (xios_field_is_active("Kapsrf")) then
!$acc update host( Kapsrf ) 
      endif
#  endif
# endif
# ifdef OW_COUPLING
      if (xios_field_is_active("hrm")) then
!$acc update host( hrm ) 
      endif
      if (xios_field_is_active("frq")) then
!$acc update host( frq ) 
      endif
      if (xios_field_is_active("wkx")) then
!$acc update host( wkx ) 
      endif
      if (xios_field_is_active("wke")) then
!$acc update host( wke ) 
      endif
      if (xios_field_is_active("epb")) then
!$acc update host( epb ) 
      endif
      if (xios_field_is_active("epd")) then
!$acc update host( epd ) 
      endif
#  endif
#  ifdef SOLVE3D
      if (xios_field_is_active("u")) then
!$acc update host( u ) 
      endif
      if (xios_field_is_active("v")) then
!$acc update host( v ) 
      endif
      if (xios_field_is_active("u_surf")) then
!$acc update host( u(:,:,N,nstp) ) 
      endif
      if (xios_field_is_active("v_surf")) then
!$acc update host( v(:,:,N,nstp) ) 
      endif
# ifdef TEMPERATURE
      if (xios_field_is_active("temp")) then
!$acc update host( t(:,:,:,nstp,itemp) ) 
      endif
      if (xios_field_is_active("temp_surf")) then
!$acc update host( t(:,:,N,nstp,itemp) ) 
      endif
# endif
# ifdef SALINITY
      if (xios_field_is_active("salt")) then
!$acc update host( t(:,:,:,nstp,isalt) ) 
      endif
      if (xios_field_is_active("salt_surf")) then
!$acc update host( t(:,:,N,nstp,isalt) ) 
      endif
# endif /* SALINITY */
      if (xios_field_is_active("rho")) then
!$acc update host( rho ) 
      endif
# if defined ANA_VMIX
      if (xios_field_is_active("bvf")) then
!$acc update host( bvf ) 
      endif
# endif
      if (xios_field_is_active("omega")) then
!$acc update host( We
#  ifdef VADV_ADAPT_IMP
!$acc&  , Wi
#  endif
#  ifdef KNBQ_MASS
!$acc&  , /rho_nbq_avg1
#  endif
!$acc& )
!
      endif
#  ifdef VADV_ADAPT_IMP
#  endif                  
#  ifdef KNBQ_MASS
#  endif
#  ifdef K3FAST
      if (xios_field_is_active("u_nbq")) then
!$acc update host( qdmu_nbq ) 
      endif
      if (xios_field_is_active("v_nbq")) then
!$acc update host( qdmv_nbq ) 
      endif
      if (xios_field_is_active("w_nbq")) then
!$acc update host( qdmw_nbq ) 
      endif
#   ifndef M3FAST_SEDLAYERS
#   else
#   endif
      if (xios_field_is_active("rho_nbq")) then
!$acc update host( rho_nbq ) 
      endif
#  endif
#  ifdef KNBQ      
      if (xios_field_is_active("w")) then
!$acc update host( wz ) 
      endif
#  endif      
#  ifdef DIAGNOSTICS_EDDY
      if (xios_field_is_active("uu")) then
!$acc update host( uu ) 
      endif
      if (xios_field_is_active("uv")) then
!$acc update host( uv ) 
      endif
      if (xios_field_is_active("vv")) then
!$acc update host( vv ) 
      endif
      if (xios_field_is_active("ub")) then
!$acc update host( ub ) 
      endif
      if (xios_field_is_active("vb")) then
!$acc update host( vb ) 
      endif
      if (xios_field_is_active("wb")) then
!$acc update host( wb ) 
      endif
      if (xios_field_is_active("uw")) then
!$acc update host( uw ) 
      endif
      if (xios_field_is_active("vw")) then
!$acc update host( vw ) 
      endif
      if (xios_field_is_active("usustr")) then
!$acc update host( usustr ) 
      endif
      if (xios_field_is_active("vsvstr")) then
!$acc update host( vsvstr ) 
      endif
      if (xios_field_is_active("ugsustr")) then
!$acc update host( ugsustr ) 
      endif
      if (xios_field_is_active("vgsvstr")) then
!$acc update host( vgsvstr ) 
      endif
      if (xios_field_is_active("ubustr")) then
!$acc update host( ubustr ) 
      endif
      if (xios_field_is_active("vbvstr")) then
!$acc update host( vbvstr ) 
      endif
      if (xios_field_is_active("uT")) then
!$acc update host( uT ) 
      endif
      if (xios_field_is_active("vT")) then
!$acc update host( vT ) 
      endif
      if (xios_field_is_active("wT")) then
!$acc update host( wT ) 
      endif
      if (xios_field_is_active("uS")) then
!$acc update host( uS ) 
      endif
      if (xios_field_is_active("vS")) then
!$acc update host( vS ) 
      endif
      if (xios_field_is_active("wS")) then
!$acc update host( wS ) 
      endif
#  endif
#  ifdef VIS_COEF_3D
      if (xios_field_is_active("visc3d")) then
!$acc update host( visc3d_r )
!$acc update host( visc3d_p )       	
      endif
#  endif
#  ifdef DIF_COEF_3D
      if (xios_field_is_active("diff3d")) then
!$acc update host( diff3d ) 
      endif
#   ifdef TS_DIF2
#    ifdef TS_DIF_SMAGO
#    endif
#   elif defined TS_DIF4
#    ifdef TS_DIF_SMAGO
#    endif
#   endif
#  endif
      if (xios_field_is_active("AKv")) then
!$acc update host( AKv ) 
      endif
#  ifdef TEMPERATURE
      if (xios_field_is_active("AKt")) then
!$acc update host( Akt(:,:,:,itemp) ) 
      endif
#  endif
#  ifdef SALINITY      
      if (xios_field_is_active("AKs")) then
!$acc update host( Akt(:,:,:,isalt) ) 
      endif
#  endif
#  ifdef GLS_MIXING
c       if (xios_field_is_active("AKk")) then
c !$acc update host( AKk ) 
c       endif
c       if (xios_field_is_active("AKp")) then
c !$acc update host( AKp ) 
c       endif
#  endif
#  ifdef LMD_SKPP
      if (xios_field_is_active("hbl")) then
!$acc update host( hbl ) 
      endif
#   ifdef LMD_SKPP2005
#   else
#   endif
#  elif defined GLS_MIXING
      if (xios_field_is_active("hbl")) then
!$acc update host( hbl ) 
      endif
#  endif
#  ifdef LMD_BKPP
      if (xios_field_is_active("hbbl")) then
!$acc update host( hbbl ) 
      endif
#  endif
#  ifdef GLS_MIXING
      if (xios_field_is_active("tke")) then
!$acc update host(  trb(:,:,:,nstp,itke)) 
      endif
      if (xios_field_is_active("gls")) then
!$acc update host( trb(:,:,:,nstp,igls) ) 
      endif
      if (xios_field_is_active("Lscale")) then
!$acc update host( Lscale ) 
      endif
#  endif
#  ifdef TEMPERATURE
c       if (xios_field_is_active("shflx")) then
c !$acc update host( shflx ) 
c       endif
#  endif      
#  ifdef SALINITY    
c       if (xios_field_is_active("swflx")) then
c !$acc update host( swflx ) 
c       endif
#  endif
#  ifdef BULK_FLUX
      if (xios_field_is_active("radsw")) then
!$acc update host( radsw ) 
      endif
#  else 
      if (xios_field_is_active("swrad")) then
!$acc update host( srflx ) 
      endif
#  endif
#  ifdef BULK_FLUX
      if (xios_field_is_active("shflx_rlw")) then
!$acc update host( shflx_rlw ) 
      endif
      if (xios_field_is_active("shflx_lat")) then
!$acc update host( shflx_lat ) 
      endif
      if (xios_field_is_active("shflx_sen")) then
!$acc update host( shflx_sen ) 
      endif
#  endif
#  ifdef SST_SKIN
      if (xios_field_is_active("sst_skin")) then
!$acc update host( sst_skin ) 
      endif
#  endif
#  ifdef DIAGNOSTICS_UV
      if (xios_field_is_active("u_rate")) then
!$acc update host( u_rate ) 
      endif
      if (xios_field_is_active("u_adv")) then
!$acc update host( u_adv ) 
      endif
      if (xios_field_is_active("u_Cor")) then
!$acc update host( u_Cor ) 
      endif
      if (xios_field_is_active("u_Prsgrd")) then
!$acc update host( u_Prsgrd ) 
      endif
      if (xios_field_is_active("u_Hmix")) then
!$acc update host( u_Hmix ) 
      endif
      if (xios_field_is_active("u_Hdiff")) then
!$acc update host( u_Hdiff ) 
      endif
      if (xios_field_is_active("u_Vmix")) then
!$acc update host( u_Vmix ) 
      endif
      if (xios_field_is_active("u_Vmix2")) then
!$acc update host( u_Vmix2 ) 
      endif
#  ifdef DIAGNOSTICS_BARO
      if (xios_field_is_active("u_Baro")) then
!$acc update host( u_Baro ) 
      endif
#  endif
#  ifdef K3FAST
      if (xios_field_is_active("u_fast")) then
!$acc update host( u_fast ) 
      endif
#  endif
#  ifdef MRL_WCI
      if (xios_field_is_active("u_vf")) then
!$acc update host( u_vf ) 
      endif
      if (xios_field_is_active("u_brk")) then
!$acc update host( u_brk ) 
      endif
      if (xios_field_is_active("u_StCo")) then
!$acc update host( u_StCo ) 
      endif
      if (xios_field_is_active("u_Vvf")) then
!$acc update host( u_Vvf ) 
      endif
      if (xios_field_is_active("u_Prscrt")) then
!$acc update host( u_Prscrt ) 
      endif
      if (xios_field_is_active("u_sbk")) then
!$acc update host( u_sbk ) 
      endif
      if (xios_field_is_active("u_bwf")) then
!$acc update host( u_bwf ) 
      endif
      if (xios_field_is_active("u_frc")) then
!$acc update host( u_frc ) 
      endif
#  endif
      if (xios_field_is_active("v_rate")) then
!$acc update host( v_rate ) 
      endif
      if (xios_field_is_active("v_adv")) then
!$acc update host( v_adv ) 
      endif
      if (xios_field_is_active("v_Cor")) then
!$acc update host( v_Cor ) 
      endif
      if (xios_field_is_active("v_Prsgrd")) then
!$acc update host( v_Prsgrd ) 
      endif
      if (xios_field_is_active("v_Hmix")) then
!$acc update host( v_Hmix ) 
      endif
      if (xios_field_is_active("v_Hdiff")) then
!$acc update host( v_Hdiff ) 
      endif
      if (xios_field_is_active("v_Vmix")) then
!$acc update host( v_Vmix ) 
      endif
      if (xios_field_is_active("v_Vmix2")) then
!$acc update host( v_Vmix2 ) 
      endif
#  ifdef DIAGNOSTICS_BARO
      if (xios_field_is_active("v_Baro")) then
!$acc update host( v_Baro ) 
      endif
#  endif
#  ifdef K3FAST
      if (xios_field_is_active("v_fast")) then
!$acc update host( v_fast ) 
      endif
#  endif
#  ifdef MRL_WCI
      if (xios_field_is_active("v_vf")) then
!$acc update host( v_vf ) 
      endif
      if (xios_field_is_active("v_brk")) then
!$acc update host( v_brk ) 
      endif
      if (xios_field_is_active("v_StCo")) then
!$acc update host( v_StCo ) 
      endif
      if (xios_field_is_active("v_Vvf")) then
!$acc update host( v_Vvf ) 
      endif
      if (xios_field_is_active("v_Prscrt")) then
!$acc update host( v_Prscrt ) 
      endif
      if (xios_field_is_active("v_sbk")) then
!$acc update host( v_sbk ) 
      endif
      if (xios_field_is_active("v_bwf")) then
!$acc update host( v_bwf ) 
      endif
      if (xios_field_is_active("v_frc")) then
!$acc update host( v_frc ) 
      endif
#  endif
#  endif
#  ifdef DIAGNOSTICS_TS
#  ifdef TEMPERATURE
      if (xios_field_is_active("T_rate")) then
!$acc update host( T_rate ) 
      endif
      if (xios_field_is_active("T_adv")) then
!$acc update host( T_adv ) 
      endif
      if (xios_field_is_active("T_Hmix")) then
!$acc update host( T_Hmix ) 
      endif
      if (xios_field_is_active("T_Vmix")) then
!$acc update host( T_Vmix ) 
      endif
      if (xios_field_is_active("T_Forc")) then
!$acc update host( T_Forc ) 
      endif
#  endif          
#  ifdef SALINITY
      if (xios_field_is_active("S_rate")) then
!$acc update host( S_rate ) 
      endif
      if (xios_field_is_active("S_adv")) then
!$acc update host( S_adv ) 
      endif
      if (xios_field_is_active("S_Hmix")) then
!$acc update host( S_Hmix ) 
      endif
      if (xios_field_is_active("S_Vmix")) then
!$acc update host( S_Vmix ) 
      endif
      if (xios_field_is_active("S_Forc")) then
!$acc update host( S_Forc ) 
      endif
#  endif
#  endif
#  ifdef DIAGNOSTICS_VRT
      if (xios_field_is_active("vrtrate")) then
!$acc update host( vrtrate ) 
      endif
      if (xios_field_is_active("vrtadv")) then
!$acc update host( vrtadv ) 
      endif
      if (xios_field_is_active("vrtCor")) then
!$acc update host( vrtCor ) 
      endif
      if (xios_field_is_active("vrtPrsgrd")) then
!$acc update host( vrtPrsgrd ) 
      endif
      if (xios_field_is_active("vrtHmix")) then
!$acc update host( vrtHmix ) 
      endif
      if (xios_field_is_active("vrtHdiff")) then
!$acc update host( vrtHdiff ) 
      endif
      if (xios_field_is_active("vrtVmix")) then
!$acc update host( vrtVmix ) 
      endif
      if (xios_field_is_active("vrtVmix2")) then
!$acc update host( vrtVmix2 ) 
      endif
#  ifdef DIAGNOSTICS_BARO
      if (xios_field_is_active("vrtBaro")) then
!$acc update host( vrtBaro ) 
      endif
#  endif
      if (xios_field_is_active("vrtDrag")) then
!$acc update host( vrtDrag ) 
      endif
      if (xios_field_is_active("vrtWind")) then
!$acc update host( vrtWind ) 
      endif
#  ifdef K3FAST
      if (xios_field_is_active("vrtfast")) then
!$acc update host( vrtfast ) 
      endif
#  endif
#  endif
#  ifdef DIAGNOSTICS_EK
      if (xios_field_is_active("ekrate")) then
!$acc update host( ekrate ) 
      endif
      if (xios_field_is_active("ekadv")) then
!$acc update host( ekadv ) 
      endif
      if (xios_field_is_active("ekCor")) then
!$acc update host( ekCor ) 
      endif
      if (xios_field_is_active("ekPrsgrd")) then
!$acc update host( ekPrsgrd ) 
      endif
      if (xios_field_is_active("ekHmix")) then
!$acc update host( ekHmix ) 
      endif
      if (xios_field_is_active("ekHdiff")) then
!$acc update host( ekHdiff ) 
      endif
      if (xios_field_is_active("ekVmix")) then
!$acc update host( ekVmix ) 
      endif
      if (xios_field_is_active("ekVmix2")) then
!$acc update host( ekVmix2 ) 
      endif
#  ifdef DIAGNOSTICS_BARO
      if (xios_field_is_active("ekBaro")) then
!$acc update host( ekBaro ) 
      endif
#  endif
      if (xios_field_is_active("ekvol")) then
!$acc update host( ekvol ) 
      endif
      if (xios_field_is_active("ekDrag")) then
!$acc update host( ekDrag ) 
      endif
      if (xios_field_is_active("ekWind")) then
!$acc update host( ekWind ) 
      endif
#  ifdef K3FAST
      if (xios_field_is_active("ekfast")) then
!$acc update host( ekfast ) 
      endif
#  endif
#  endif
# ifdef DIAGNOSTICS_PV
      if (xios_field_is_active("u_rhs")) then
!$acc update host( u_rhs ) 
      endif
      if (xios_field_is_active("v_rhs")) then
!$acc update host( v_rhs ) 
      endif
# ifdef TEMPERATURE          
      if (xios_field_is_active("T_rhs")) then
!$acc update host( T_rhs ) 
      endif
#  endif
#  ifdef SALINITY
      if (xios_field_is_active("S_rhs")) then
!$acc update host( S_rhs ) 
      endif
#  endif
# ifdef DIAGNOSTICS_PV_FULL
      if (xios_field_is_active("u_vmix_trans")) then
!$acc update host( u_vmix_trans ) 
      endif
      if (xios_field_is_active("v_vmix_trans")) then
!$acc update host( v_vmix_trans ) 
      endif
#  endif
#  endif
#  if defined BIOLOGY
#   if defined PISCES
#   else
      if (xios_field_is_active("hel")) then
!$acc update host( hel ) 
      endif
#   if (defined BIO_NChlPZD ||  defined BIO_N2ChlPZD2)
      if (xios_field_is_active("theta")) then
!$acc update host( theta ) 
      endif
#    ifdef OXYGEN
      if (xios_field_is_active("U10")) then
!$acc update host( U10 ) 
      endif
      if (xios_field_is_active("KvO2")) then
!$acc update host( KvO2 ) 
      endif
      if (xios_field_is_active("O2sat")) then
!$acc update host( O2sat ) 
      endif
#    endif /* OXYGEN */
#   elif defined BIO_BioEBUS 
      if (xios_field_is_active("AOU")) then
!$acc update host( AOU ) 
      endif
      if (xios_field_is_active("wind10")) then
!$acc update host( wind10 ) 
      endif
#   endif
#   endif
#  endif /* BIOLOGY */
#  ifdef SEDIMENT
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active("bed_thick")) then
!$acc update host( bed_thick ) 
      endif
      if (xios_field_is_active("bed_poros")) then
!$acc update host( bed_poros ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
#   ifdef SUSPLOAD
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
#   endif
#   ifdef BEDLOAD
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
      if (xios_field_is_active(TRIM(nametrc))) then
!$acc update host(  ) 
      endif
#   endif
#  endif /* SEDIMENT */
# ifdef MUSTANG
#if defined key_sand2D
#else
#endif
      if (xios_field_is_active("tauskin")) then
!$acc update host( tauskin ) 
      endif
# ifdef WAVE_OFFLINE
      if (xios_field_is_active("tauskin_c")) then
!$acc update host( tauskin_c ) 
      endif
      if (xios_field_is_active("tauskin_w")) then
!$acc update host( tauskin_w ) 
      endif
# endif
      if (xios_field_is_active("ksma")) then
!$acc update host( ksma ) 
      endif
      if (xios_field_is_active('dzs')) then
!$acc update host(  ) 
      endif
# endif /* MUSTANG */     
#  ifdef BBL
      if (xios_field_is_active("Abed")) then
!$acc update host( Abed ) 
      endif
      if (xios_field_is_active("Hripple")) then
!$acc update host( Hripple ) 
      endif
      if (xios_field_is_active("Lripple")) then
!$acc update host( Lripple ) 
      endif
      if (xios_field_is_active("Zbnot")) then
!$acc update host( Zbnot ) 
      endif
      if (xios_field_is_active("Zbapp")) then
!$acc update host( Zbapp ) 
      endif
      if (xios_field_is_active("bostrw")) then
!$acc update host( bostrw ) 
      endif
      if (xios_field_is_active("bustrc")) then
!$acc update host( bustrc ) 
      endif
      if (xios_field_is_active("bvstrc")) then
!$acc update host( bvstrc ) 
      endif
      if (xios_field_is_active("bustrw")) then
!$acc update host( bustrw ) 
      endif
      if (xios_field_is_active("bvstrw")) then
!$acc update host( bvstrw ) 
      endif
      if (xios_field_is_active("bustrcwmax")) then
!$acc update host( bustrcwmax ) 
      endif
      if (xios_field_is_active("bvstrcwmax")) then
!$acc update host( bvstrcwmax ) 
      endif
#  endif /* BBL */
# endif /* SOLVE3D */
#else   /* XIOS */
#endif /* XIOS */
