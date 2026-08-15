!======================================================================
!
/* This is include file "abl1d.h".
  --------------------------------------------
*/
#ifdef ABL1D
      integer  nt_n,nt_a
      common /time_indices_abl/  nt_n,nt_a

      integer jptq,jp_ta,jp_qa,jptime
      parameter(jptq=2,jp_ta=1,jp_qa=2,jptime=2)

      real u_abl (GLOBAL_2D_ARRAY,N_abl,jptime     )
      real v_abl (GLOBAL_2D_ARRAY,N_abl,jptime     )
      real tq_abl(GLOBAL_2D_ARRAY,N_abl,jptime,jptq)
      common /abl_u/u_abl /abl_v/v_abl /abl_tq/tq_abl

      real zr_abl(N_abl)
      real zw_abl(N_abl)
      real Hzr_abl(N_abl)
      real Hzw_abl(N_abl)
      common /ablgrid_ght/zr_abl /ablgrid_ghw/zw_abl
      common /ablgrid_e3t/Hzr_abl /ablgrid_e3w/Hzw_abl

      real avm_abl (GLOBAL_2D_ARRAY,N_abl,jptime)
      real avt_abl (GLOBAL_2D_ARRAY,N_abl,jptime)
      real mxld_abl(GLOBAL_2D_ARRAY,N_abl,jptime)
      real mxlm_abl(GLOBAL_2D_ARRAY,N_abl,jptime)
      real tke_abl (GLOBAL_2D_ARRAY,N_abl,jptime)
      real ablh    (GLOBAL_2D_ARRAY,      jptime)
      common /ablturb_avm/avm_abl /ablturb_avt/avt_abl
      common /ablturb_mxld/mxld_abl /ablturb_mxlm/mxlm_abl
      common /ablturb_tke/tke_abl /ablturb_ablh/ablh

      real pu_dta  (GLOBAL_2D_ARRAY,N_abl )
      real pv_dta  (GLOBAL_2D_ARRAY,N_abl )
      real pt_dta  (GLOBAL_2D_ARRAY,N_abl )
      real pq_dta  (GLOBAL_2D_ARRAY,N_abl )
      real pgu_dta (GLOBAL_2D_ARRAY,N_abl )
      real pgv_dta (GLOBAL_2D_ARRAY,N_abl )

      real pu_dtag  (GLOBAL_2D_ARRAY,N_abl,2 )
      real pv_dtag  (GLOBAL_2D_ARRAY,N_abl,2 )
      real pt_dtag  (GLOBAL_2D_ARRAY,N_abl,2 )
      real pq_dtag  (GLOBAL_2D_ARRAY,N_abl,2 )
      real pgu_dtag (GLOBAL_2D_ARRAY,N_abl,2 )
      real pgv_dtag (GLOBAL_2D_ARRAY,N_abl,2 )

      common /abl_u_dtag/pu_dtag /abl_v_dtag/pv_dtag
      common /abl_t_dtag/pt_dtag /abl_q_dtag/pq_dtag
      common /abl_pgu_dtag/pgu_dtag /abl_pgv_dtag/pgv_dtag

      common /abl_u_dta/pu_dta /abl_v_dta/pv_dta
      common /abl_t_dta/pt_dta /abl_q_dta/pq_dta
      common /abl_pgu_dta/pgu_dta /abl_pgv_dta/pgv_dta

      real Cd_du  (GLOBAL_2D_ARRAY        )
      real Ch_du  (GLOBAL_2D_ARRAY        )
      real Ce_du  (GLOBAL_2D_ARRAY        )
      real sst_abl(GLOBAL_2D_ARRAY        )
      real ssq_abl(GLOBAL_2D_ARRAY        )
      real rho_abl(GLOBAL_2D_ARRAY        )
      real ustar2 (GLOBAL_2D_ARRAY        )
      real z0_abl (GLOBAL_2D_ARRAY        )
      common /abl_sflux_Cd/Cd_du
      common /abl_sflux_Ch/Ch_du
      common /abl_sflux_Ce/Ce_du
      common /abl_sflux_st/sst_abl
      common /abl_sflux_sq/ssq_abl
      common /abl_sflux_ra/rho_abl
      common /abl_sflux_us/ustar2
      common /abl_sflux_zo/z0_abl

# if defined ABL_NUDGING && defined ABL_DYN_RESTORE_EQ
      real rest_eq    (GLOBAL_2D_ARRAY        )
      common /ablrest_eq/rest_eq
# endif
# if defined ABL_NUDGING
      real bmin, bmax
      parameter( bmin = 0.5, bmax = 1.5 )
      real pblh_min, pblh_max
      common /abl_nudging/pblh_min,pblh_max
#  ifdef ABL_NUDGING_DYN
      real alp0_dyn, alp1_dyn
      real alp2_dyn, alp3_dyn
      real ldyn_min,ldyn_max
      common /abl_nudging_dyn/alp0_dyn, alp1_dyn,
     &                        alp2_dyn, alp3_dyn, ldyn_min,ldyn_max
#  endif
#  ifdef ABL_NUDGING_TRA
      real alp0_tra, alp1_tra
      real alp2_tra, alp3_tra
      real ltra_min,ltra_max
      common /abl_nudging_tra/alp0_tra, alp1_tra,
     &                        alp2_tra, alp3_tra, ltra_min,ltra_max
#  endif
# endif

      real mxl_min,Cm,Ct,Ceps,Sch,lsfc,esfc,rctv0,Ric_abl
      real tke_min,avm_bak,avt_bak,itvref,phimax,Cek,epssfc,grav
      real vkarmn,ctke
      parameter(tke_min=1.e-6,avm_bak=1.e-4,avt_bak=1.e-5,itvref=1./288.)
      parameter(phimax =(1.-2.2)/2.2, Cek = 258)
      parameter(epssfc = 1. /(1.+2.8*2.8 ) )
      parameter(grav=9.81,vkarmn=0.41)
      parameter(Ric_abl = 0.143)
      parameter(ctke    = 0.340)
      parameter(cm      = 0.126)
      parameter(ct      = 0.143)
      parameter(ceps    = 0.845)
      parameter(Sch     = ctke/Cm)
      parameter(mxl_min = (avm_bak / cm) / sqrt( tke_min ))
      parameter(esfc    = 1./SQRT(cm*ceps) )
      parameter(lsfc    = vkarmn * SQRT(SQRT(cm*ceps))/cm)
      parameter(rctv0   = 0.608)

# ifdef AVERAGES
      real timeablavg
      common /timeabl_avg/timeablavg

      real u_abl_avg (GLOBAL_2D_ARRAY,N_abl)
      real v_abl_avg (GLOBAL_2D_ARRAY,N_abl)
      real t_abl_avg (GLOBAL_2D_ARRAY,N_abl)
      real q_abl_avg (GLOBAL_2D_ARRAY,N_abl)
      common /abl_u_avg/u_abl_avg /abl_v_avg/v_abl_avg
      common /abl_t_avg/t_abl_avg /abl_q_avg/q_abl_avg

      real avm_abl_avg (GLOBAL_2D_ARRAY,N_abl)
      real avt_abl_avg (GLOBAL_2D_ARRAY,N_abl)
      real mxld_abl_avg(GLOBAL_2D_ARRAY,N_abl)
      real mxlm_abl_avg(GLOBAL_2D_ARRAY,N_abl)
      real tke_abl_avg (GLOBAL_2D_ARRAY,N_abl)
      real ablh_avg    (GLOBAL_2D_ARRAY      )
      common /ablturb_avm_avg/avm_abl_avg /ablturb_avt_avg/avt_abl_avg
      common /ablturb_mxld_avg/mxld_abl_avg /ablturb_mxlm_avg/mxlm_abl_avg
      common /ablturb_tke_avg/tke_abl_avg /ablturb_ablh_avg/ablh_avg

      real pu_dta_avg  (GLOBAL_2D_ARRAY,N_abl)
      real pv_dta_avg  (GLOBAL_2D_ARRAY,N_abl)
      real pt_dta_avg  (GLOBAL_2D_ARRAY,N_abl)
      real pq_dta_avg  (GLOBAL_2D_ARRAY,N_abl)
      real pgu_dta_avg (GLOBAL_2D_ARRAY,N_abl)
      real pgv_dta_avg (GLOBAL_2D_ARRAY,N_abl)
      common /abl_u_dta_avg/pu_dta_avg /abl_v_dta_avg/pv_dta_avg
      common /abl_t_dta_avg/pt_dta_avg /abl_q_dta_avg/pq_dta_avg
      common /abl_pgu_dta_avg/pgu_dta_avg /abl_pgv_dta_avg/pgv_dta_avg

# endif /* AVERAGES */
#endif /* ABL1D */
