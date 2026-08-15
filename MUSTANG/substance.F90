#include "cppdefs.h"

MODULE substance

!!======================================================================
!! ***  MODULE  substance  ***
!! Substances conservative or not, dissolved or particulate   :  
!! module for tracers/substances defined
!!======================================================================

#if defined SUBSTANCE

   USE module_substance
   USE comsubstance
   USE submassbalance, ONLY :  submassbalance_readdomain

   IMPLICIT NONE
   PRIVATE
    
   PUBLIC   substance_read_alloc   ! called by main.F
   PUBLIC   substance_surfcell     ! called by main.F
   
   ! variables for all ntrc_subs
   CHARACTER(LEN=lchain),DIMENSION(ntrc_subs)      :: name_var_r, long_name_var_r, standard_name_var_r, &
                                                      unit_var_r, init_cv_name_r, obc_cv_name_r
   REAL(KIND=rsh), DIMENSION(ntrc_subs)            :: flx_atm_r, cv_rain_r, cini_wat_r, cini_air_r, cobc_wat_r
   LOGICAL, DIMENSION(ntrc_subs)                   :: l_out_subs_r 

   ! local variable to read each namelist with different length (variables by group of substances)
   CHARACTER(LEN=lchain),DIMENSION(:),ALLOCATABLE  :: name_var_n, long_name_var_n, standard_name_var_n, &
                                                      unit_var_n, init_cv_name_n, obc_cv_name_n
   REAL(KIND=rsh), DIMENSION(:),ALLOCATABLE        :: flx_atm_n, cv_rain_n, cini_wat_n, cini_air_n, cobc_wat_n  
   LOGICAL, DIMENSION(:),ALLOCATABLE               :: l_out_subs_n    
#ifdef MUSTANG
   REAL(KIND=rsh), DIMENSION(:),ALLOCATABLE   :: cini_sed_n
#endif

   !!----------------------------------------------------------------------
   
CONTAINS

   !!======================================================================

  SUBROUTINE substance_read_alloc(may_day_flag, indxT, indxTsrc)
      !!-------------------------------------------------------------------
      !!                    *** ROUTINE substance_read_alloc ***
      !!-------------------------------------------------------------------
      !
   
   !! Argument
   INTEGER,INTENT(INOUT)                     ::  may_day_flag
   INTEGER,INTENT(IN)                        ::  indxT, indxTsrc
   
   !! Local declarations
   LOGICAL                                   :: l_varassoc
   INTEGER                                   :: ivpc, ivp, iv, iv0, indx, ivTS
   INTEGER                                   :: isubs, nballoc, ivr, it, ntypvar
   INTEGER                                   :: lstr, lenstr
   INTEGER :: rc
   CHARACTER(LEN=lchain) :: msg

!! tables (_n) sized to read in namelist by number of substances of such and such a type
!! tables (_r) intermediates sized to the number of substances, will then be copied into the final 
!!             table sized to NT after adding additional variables and reordering as needed 
   INTEGER,DIMENSION(ntrc_subs)            :: itypv_r
   REAL(KIND=rsh), DIMENSION(ntrc_subs)    :: ws_free_min_r, ws_free_max_r
   CHARACTER(LEN=lchain),DIMENSION(ntfix)  :: long_name_var_fix, unit_var_fix, standard_name_var_fix, name_var_fix
   LOGICAL,DIMENSION(ntfix)                :: l_out_subs_fix

   REAL(KIND=rsh), DIMENSION(:),ALLOCATABLE       :: ws_free_min_n,ws_free_max_n
   CHARACTER(LEN=lchain),DIMENSION(ntrc_subs)     :: name_varpc_assoc
   CHARACTER(LEN=lchain),DIMENSION(:),ALLOCATABLE :: name_varpc_assoc_n
#if defined MUSTANG
   REAL(KIND=rsh),DIMENSION(ntrc_subs)        :: diam_r, ros_r, tocd_r 
   REAL(KIND=rsh), DIMENSION(4,ntrc_subs)     :: ws_free_para_r
   REAL(KIND=rsh), DIMENSION(2,ntrc_subs)     :: ws_hind_para_r 
   INTEGER, DIMENSION(ntrc_subs)              :: ws_hind_opt_r, ws_free_opt_r 
   LOGICAL, DIMENSION(ntrc_subs)              :: l_bedload_r
   LOGICAL, DIMENSION(:),ALLOCATABLE          :: l_sand2D_n, l_outsandrouse_n, l_bedload_n
   REAL(KIND=rsh), DIMENSION(:),ALLOCATABLE   :: tocd_n, ros_n, diam_n
   REAL(KIND=rsh), DIMENSION(:),ALLOCATABLE   :: ws_free_opt_n, ws_hind_opt_n
   REAL(KIND=rsh), DIMENSION(:,:),ALLOCATABLE :: ws_free_para_n, ws_hind_para_n
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
   LOGICAL                                    ::  l_ibedload1, l_ibedload2
#endif
#if defined key_sand2D
   LOGICAL, DIMENSION(ntrc_subs)              :: l_outsandrouse_r, l_sand2D_r
#endif
#endif

   REAL(KIND=rlg)        :: tool_datosec
                                    
   !! *  define namelists reading in parasubstance.txt
#ifdef MUSTANG
   NAMELIST/nmlnbvar/ nv_dis, nv_ncp, nv_bent, nv_fix, nv_grav, nv_sand, nv_mud, nv_sorb
   NAMELIST/nmlgravels/name_var_n,long_name_var_n,standard_name_var_n,unit_var_n, &
                       flx_atm_n,cv_rain_n,cini_wat_n,cini_sed_n,l_bedload_n, &
                       cini_air_n,l_out_subs_n,init_cv_name_n,obc_cv_name_n, &
                       tocd_n,ros_n,diam_n
   NAMELIST/nmlsands/name_var_n,long_name_var_n,standard_name_var_n,unit_var_n, &
                       flx_atm_n,cv_rain_n,cini_wat_n,cini_sed_n,l_bedload_n, &
                       cini_air_n,l_out_subs_n,init_cv_name_n,obc_cv_name_n, &
                       tocd_n,ros_n,diam_n,l_sand2D_n,l_outsandrouse_n
   NAMELIST/nmlmuds/name_var_n,long_name_var_n,standard_name_var_n,unit_var_n, &
                       flx_atm_n,cv_rain_n,cini_wat_n,cini_sed_n,cobc_wat_n, &
                       cini_air_n,l_out_subs_n,init_cv_name_n,obc_cv_name_n, &
                       tocd_n,ros_n,diam_n, &
                       ws_free_opt_n,ws_free_min_n,ws_free_max_n,ws_free_para_n, &
                       ws_hind_opt_n,ws_hind_para_n
   NAMELIST/nmlpartnc/name_var_n,long_name_var_n,standard_name_var_n,unit_var_n, &
                       flx_atm_n,cv_rain_n,cini_wat_n,cini_sed_n,cobc_wat_n, &
                       cini_air_n,l_out_subs_n,init_cv_name_n,obc_cv_name_n,tocd_n,ros_n, &
                       ws_free_opt_n,ws_free_min_n,ws_free_max_n,ws_free_para_n, &
                       ws_hind_opt_n,ws_hind_para_n
   NAMELIST/nmlpartsorb/name_var_n,long_name_var_n,standard_name_var_n,unit_var_n, &
                       flx_atm_n,cv_rain_n,cini_wat_n,cini_sed_n,cobc_wat_n, &
                       cini_air_n,l_out_subs_n,init_cv_name_n,obc_cv_name_n, &
                       name_varpc_assoc_n
   NAMELIST/nmlvardiss/name_var_n,long_name_var_n,standard_name_var_n,unit_var_n, &
                       flx_atm_n,cv_rain_n,cini_wat_n,cini_sed_n,cobc_wat_n, &
                       cini_air_n,l_out_subs_n,init_cv_name_n,obc_cv_name_n
#else     
   NAMELIST/nmlnbvar/ nv_dis, nv_ncp, nv_bent, nv_fix, nv_sorb                 
   NAMELIST/nmlpartnc/name_var_n,long_name_var_n,standard_name_var_n,unit_var_n, &
                       flx_atm_n,cv_rain_n,cini_wat_n,cobc_wat_n, &
                       cini_air_n,l_out_subs_n,init_cv_name_n,obc_cv_name_n, &
                       ws_free_min_n,ws_free_max_n
   NAMELIST/nmlpartsorb/name_var_n,long_name_var_n,standard_name_var_n,unit_var_n, &
                       flx_atm_n,cv_rain_n,cini_wat_n,cobc_wat_n,            &
                       cini_air_n,l_out_subs_n,init_cv_name_n,obc_cv_name_n, &
                       name_varpc_assoc_n
   NAMELIST/nmlvardiss/name_var_n,long_name_var_n,standard_name_var_n,unit_var_n, &
                       flx_atm_n,cv_rain_n,cini_wat_n,cobc_wat_n, &
                       cini_air_n,l_out_subs_n,init_cv_name_n,obc_cv_name_n
#endif                      
   NAMELIST/nmlvarfix/name_var_fix,long_name_var_fix,standard_name_var_fix,unit_var_fix, &
                       cini_wat_fix,l_out_subs_fix,init_cv_name_fix  
#ifdef key_benthic 
   NAMELIST/nmlvarbent/name_var_bent, long_name_var_bent, standard_name_var_bent, unit_var_bent, &
                       cini_bent, l_out_subs_bent
#endif 
   NAMELIST/nmlsubmassbalance/submassbalance_l, submassbalance_nb_border, &
                    submassbalance_input_file, submassbalance_output_file, &
                    submassbalance_dtout, submassbalance_date_start

   !!----------------------------------------------------------------------
   !! * Executable part
   MPI_master_only   WRITE(stdout,*) ' '
   MPI_master_only   WRITE(stdout,*) ' '
   MPI_master_only   WRITE(stdout,*) ' '
   MPI_master_only   WRITE(stdout,*) '**************************************************'
   MPI_master_only   WRITE(stdout,*) '**************** MODULE SUBSTANCE     ************'
   MPI_master_only   WRITE(stdout,*) '**************** substance_read_alloc ************'
   MPI_master_only   WRITE(stdout,*) '**************************************************'
   MPI_master_only   WRITE(stdout,*) ' '


   lstr=lenstr(subsfilename)
   MPI_master_only  WRITE(stdout,*),'SUBS:',subsfilename(1:lstr)
   OPEN(500,file=subsfilename(1:lstr),status='old',form='formatted',access='sequential')
   READ(500,nmlnbvar)
#if !defined MUSTANG
   nv_grav=0
   nv_sand=0
   nv_mud=0
#endif


! check coherence between parasubstance and param.h dimensions
#ifdef MUSTANG
   IF(nv_dis+nv_ncp+nv_grav+nv_sand+nv_mud+nv_sorb .NE. ntrc_subs) THEN
     MPI_master_only  WRITE(stdout,*)'WARNING - the total number of substances read from the file'
     MPI_master_only  WRITE(stdout,*)' parasubstance.txt is DIFFERENT from the ntrc_subs parameter'
     MPI_master_only  WRITE(stdout,*)'in param.h  '
     MPI_master_only  WRITE(stdout,*)'ntrc_subs in param.h = ',ntrc_subs
     MPI_master_only  WRITE(stdout,*)'  nv_dis read in parasubstance.txt ',nv_dis
     MPI_master_only  WRITE(stdout,*)'+ nv_ncp read in parasubstance.txt ',nv_ncp
     MPI_master_only  WRITE(stdout,*)'+ nv_grav read in parasubstance.txt ',nv_grav
     MPI_master_only  WRITE(stdout,*)'+ nv_sand read in parasubstance.txt ',nv_sand
     MPI_master_only  WRITE(stdout,*)'+ nv_mud read in parasubstance.txt ',nv_mud
     MPI_master_only  WRITE(stdout,*)'+ nv_sorb read in parasubstance.txt ',nv_sorb
     MPI_master_only  WRITE(stdout,*)'The simulation will stop'
     may_day_flag=77
     goto 99
   END IF
#else
   IF(nv_dis+nv_ncp+nv_sorb .NE. ntrc_subs) THEN
     MPI_master_only  WRITE(stdout,*)'WARNING - the total number of substances read from the file'
     MPI_master_only  WRITE(stdout,*)' parasubstance.txt is DIFFERENT from the ntrc_subs parameter'
     MPI_master_only  WRITE(stdout,*)'in param.h  '
     MPI_master_only  WRITE(stdout,*)'ntrc_subs in param.h = ',ntrc_subs
     MPI_master_only  WRITE(stdout,*)'  nv_dis read in parasubstance.txt ',nv_dis
     MPI_master_only  WRITE(stdout,*)'+ nv_ncp read in parasubstance.txt ',nv_ncp
     MPI_master_only  WRITE(stdout,*)'+ nv_sorb read in parasubstance.txt ',nv_sorb
     MPI_master_only  WRITE(stdout,*)'The simulation will stop'
     may_day_flag=77
     goto 99
   END IF 
#endif

 
   IF(nv_fix .NE. ntfix) THEN
     MPI_master_only  WRITE(stdout,*)'WARNING - the number of FIXED substances read from the file'
     MPI_master_only  WRITE(stdout,*)' parasubstance.txt is DIFFERENT from the ntfix parameter'
     MPI_master_only  WRITE(stdout,*)'in param.h  '
     MPI_master_only  WRITE(stdout,*)'ntfix in param.h = ',ntfix
     MPI_master_only  WRITE(stdout,*)'  nv_fix read in parasubstance.txt ',nv_fix
     MPI_master_only  WRITE(stdout,*)'The simulation will stop'
     may_day_flag=77
     goto 99
   END IF
   ivTS=itsubs1-1
   ivp=0
   iv=0
   ALLOCATE(name_var(ntrc_subs))
   ALLOCATE(long_name_var(ntrc_subs))
   ALLOCATE(standard_name_var(ntrc_subs))
   ALLOCATE(unit_var(ntrc_subs))
   ALLOCATE(l_out_subs(ntrc_subs))
 
    !******************************************
    !    reading NAMELISTS of parasubstance.txt
    !******************************************
 
    ! the global arrays must already be allocated
    ! and we automatically put them in order thanks to the successive namelists
#ifdef MUSTANG
   ALLOCATE(cini_sed_r(ntrc_subs))

   ! reading gravels variables 
   !------------------------------
   IF(nv_grav > 0) THEN
    CALL ALLOC_DEFVAR(nv_grav)   
    ALLOCATE(tocd_n(nv_grav))
    ALLOCATE(diam_n(nv_grav))
    ALLOCATE(ros_n(nv_grav))
    ALLOCATE(l_bedload_n(nv_grav))
    READ(500,nmlgravels)
    iv0=iv
    CALL DEFVAR_DEALLOC(nv_grav,iv)
    DO ivr=1,nv_grav
     ivp=ivp+1
     tocd_r(ivp)=tocd_n(ivr)
     diam_r(ivp)=diam_n(ivr)
     ros_r(ivp)=ros_n(ivr)
     l_bedload_r(ivp)=l_bedload_n(ivr)
     itypv_r(iv0+ivr)=1
    ENDDO
    DEALLOCATE(tocd_n,diam_n,ros_n,l_bedload_n)
   ENDIF
   
  ! reading sand variables
  !------------------------------ 
   IF(nv_sand > 0) THEN
    CALL ALLOC_DEFVAR(nv_sand)   
    ALLOCATE(tocd_n(nv_sand))
    ALLOCATE(diam_n(nv_sand))
    ALLOCATE(ros_n(nv_sand))
    ALLOCATE(l_sand2D_n(nv_sand))
    ALLOCATE(l_outsandrouse_n(nv_sand))
    ALLOCATE(l_bedload_n(nv_sand))
    READ(500,nmlsands)
    iv0=iv
    CALL DEFVAR_DEALLOC(nv_sand,iv)   
    DO ivr=1,nv_sand
     ivp=ivp+1
     tocd_r(ivp)=tocd_n(ivr)
     diam_r(ivp)=diam_n(ivr)
     ros_r(ivp)=ros_n(ivr)
     l_bedload_r(ivp)=l_bedload_n(ivr)
#ifdef key_sand2D
     l_sand2D_r(ivp)=l_sand2D_n(ivr)
     l_outsandrouse_r(ivp)=l_outsandrouse_n(ivr)
#endif
     itypv_r(iv0+ivr)=2
    ENDDO
    DEALLOCATE(tocd_n,diam_n,ros_n,l_sand2D_n,l_outsandrouse_n,l_bedload_n)
   ENDIF

  ! reading mud variables
  !------------------------------ 
   IF(nv_mud > 0) THEN
    CALL ALLOC_DEFVAR(nv_mud)   
    ALLOCATE(tocd_n(nv_mud))
    ALLOCATE(diam_n(nv_mud))
    ALLOCATE(ros_n(nv_mud))
    ALLOCATE(ws_free_opt_n(nv_mud))
    ALLOCATE(ws_free_min_n(nv_mud))
    ALLOCATE(ws_free_max_n(nv_mud))
    ALLOCATE(ws_free_para_n(4,nv_mud))
    ALLOCATE(ws_hind_opt_n(nv_mud))
    ALLOCATE(ws_hind_para_n(2,nv_mud))
    READ(500,nmlmuds)
    iv0=iv   
    CALL DEFVAR_DEALLOC(nv_mud,iv)   
    DO ivr=1,nv_mud
     ivp=ivp+1
     ws_free_opt_r(ivp)=ws_free_opt_n(ivr)
     ws_free_min_r(ivp)=ws_free_min_n(ivr)
     ws_free_max_r(ivp)=ws_free_max_n(ivr)
     ws_free_para_r(1:4,ivp)=ws_free_para_n(1:4,ivr)
     ws_hind_opt_r(ivp)=ws_hind_opt_n(ivr)
     ws_hind_para_r(1:2,ivp)=ws_hind_para_n(1:2,ivr)
     tocd_r(ivp)=tocd_n(ivr)
     diam_r(ivp)=diam_n(ivr)
     ros_r(ivp)=ros_n(ivr)
     itypv_r(iv0+ivr)=3
     
    ENDDO
    DEALLOCATE(ws_free_opt_n,ws_free_min_n,ws_free_max_n,ws_free_para_n, &
                   ws_hind_opt_n,ws_hind_para_n,tocd_n,diam_n,ros_n)
   ENDIF

#else  /* MUSTANG*/
    nv_grav=0
    nv_sand=0
    nv_mud=0
#endif

   ! reading non constitutive particulate variables
   !----------------------------------------------- 
   IF(nv_ncp > 0) THEN
    CALL ALLOC_DEFVAR(nv_ncp)  
#if defined MUSTANG 
    ALLOCATE(tocd_n(nv_ncp))
    ALLOCATE(ros_n(nv_ncp))
    ALLOCATE(ws_free_opt_n(nv_ncp))
    ALLOCATE(ws_free_para_n(4,nv_ncp))
    ALLOCATE(ws_hind_opt_n(nv_ncp))
    ALLOCATE(ws_hind_para_n(2,nv_ncp))
#endif
    ALLOCATE(ws_free_min_n(nv_ncp))
    ALLOCATE(ws_free_max_n(nv_ncp))

    READ(500,nmlpartnc)
    iv0=iv   
    CALL DEFVAR_DEALLOC(nv_ncp,iv)   
    DO ivr=1,nv_ncp
     ivp=ivp+1
#if defined MUSTANG 
     tocd_r(ivp)=tocd_n(ivr)
     ros_r(ivp)=ros_n(ivr)
     ws_free_opt_r(ivp)=ws_free_opt_n(ivr)
     ws_free_para_r(1:4,ivp)=ws_free_para_n(1:4,ivr)
     ws_hind_opt_r(ivp)=ws_hind_opt_n(ivr)
     ws_hind_para_r(1:2,ivp)=ws_hind_para_n(1:2,ivr)
#endif
     ws_free_min_r(ivp)=ws_free_min_n(ivr)
     ws_free_max_r(ivp)=ws_free_max_n(ivr)
     itypv_r(iv0+ivr)=4
    ENDDO
#if defined MUSTANG 
    DEALLOCATE(ws_free_opt_n,ws_free_min_n,ws_free_max_n,ws_free_para_n, &
                   ws_hind_opt_n,ws_hind_para_n,tocd_n,ros_n)
#else
    DEALLOCATE(ws_free_min_n,ws_free_max_n)
#endif
   ENDIF

  ! reading non constitutive SORBED particulate variables 
  !--------------------------------------------------------
   IF(nv_sorb > 0) THEN
    CALL ALLOC_DEFVAR(nv_sorb)   
    ALLOCATE(name_varpc_assoc_n(nv_sorb))
    READ(500,nmlpartsorb)
    iv0=iv   
    CALL DEFVAR_DEALLOC(nv_sorb,iv)   
    DO ivr=1,nv_sorb
      ivp=ivp+1
      name_varpc_assoc(ivp)=name_varpc_assoc_n(ivr)
      itypv_r(iv0+ivr)=5
    ENDDO
   ENDIF

   ! reading dissolved variables characteristics
   !-------------------------------------------
   IF(nv_dis > 0) THEN
    CALL ALLOC_DEFVAR(nv_dis)   
    READ(500,nmlvardiss)
    iv0=iv   
    CALL DEFVAR_DEALLOC(nv_dis,iv)   
    itypv_r(iv0+1:iv)=6
   ENDIF

   ! reading Fixed variables characteristics
   !-----------------------------------------
   IF(nv_fix > 0) THEN
    nballoc=ntfix
    ALLOCATE(cini_wat_fix(nballoc))
    ALLOCATE(init_cv_name_fix(nballoc))
    READ(500,nmlvarfix)
   ENDIF
 
#ifdef key_benthic
   ! reading Benthic variables characteristics
   !-----------------------------------------
   IF(nv_bent > 0) THEN
    nballoc=nv_bent
    ALLOCATE(name_var_bent(nballoc))
    ALLOCATE(long_name_var_bent(nballoc))
    ALLOCATE(standard_name_var_bent(nballoc))
    ALLOCATE(unit_var_bent(nballoc))
    ALLOCATE(cini_bent(nballoc))
    ALLOCATE(l_out_subs_bent(nballoc))
    READ(500,nmlvarbent)
   ENDIF
#endif

   ! initialize the number of variables according to their type
   nvpc=nv_mud+nv_sand+nv_grav
   nvp=nvpc+nv_ncp+nv_sorb
   nv_adv=nvp+nv_dis ! ntrc_subs
   nv_state=nv_adv+nv_fix
   nv_tot=nv_state ! ntrc_substot

    READ (500, nmlsubmassbalance, iostat=rc); REWIND (500)
    IF (rc /= 0) THEN
       msg = "WARNING : nmlsubmassbalance, namelist not found, default values are used :"
       MPI_master_only  WRITE(stdout,*) msg
       MPI_master_only  WRITE(stdout,nml=nmlsubmassbalance)
    END IF

    if (submassbalance_l) then
        submassbalance_tdeb = tool_datosec(submassbalance_date_start)
        CALL submassbalance_readdomain()
    endif

    !******************************************
    !    create new variables 
    !******************************************

#if defined BLOOM && (defined key_N_tracer || defined key_P_tracer)
  CALL bloom_create_vartracer(flx_atm_r,cv_rain_r,cini_wat_r,cini_air_r,          &
                                  l_out_subs_r,init_cv_name_r,obc_cv_name_r,            &
#if defined MUSTANG
                                  ws_free_opt_r,ws_free_para_r,ws_hind_opt_r,           &
                                  ws_hind_para_r,tocd_r,diam_r,ros_r,cini_sed_r,        &
#endif
                                  !itypv_r,name_var,long_name_var,standard_name_var,     &
                                  itypv_r,long_name_var,                                &
                                  name_varpc_assoc,unit_var_r,l_out_subs_fix,           &
                                  long_name_var_fix,unit_var_fix,                       &
                                  ws_free_min_r,ws_free_max_r)

#endif

#if defined PEPTIC
   CALL peptic_create_var(iv,flx_atm_r,cv_rain_r,cini_wat_r,cini_air_r, &
                                l_out_subs_r,ws_free_min_r,ws_free_max_r,cobc_wat_r,   &
                                itypv_r,init_cv_name_r,obc_cv_name_r,unit_var_r,       &
                                long_name_var)
#endif


!#ifdef key_contaminant
!   ! for the variable 'CONTA' : introduction of contaminant species  
!   ! and then possible addition of dissolved or sorbed variables 
!   ! ( a modifier dans conta_read_date la lecture du paraconta (sans fileconta par exemple)
!   IF (l_conta) THEN
!        ALLOCATE(ivdl_conta(nb_conta))
!        CALL conta_read_data(obc_cv_name_r,nb_frac_var_r,                         &
!                             icon_var_r, init_cv_name_r,ispc_var_r,               &
!                             flx_atm_r,cv_rain_r,cini_wat_r,cini_sed_r,           &
!                             cini_air_r,l_out_subs_r,tocd_r,ros_r,diam_r,         &
!                             ws_free_opt_r,ws_free_min_r,ws_free_max_r,           &
!                             ws_free_para_r,ws_hind_opt_r,ws_hind_para_r          &
!
!                                            )
!    END IF   
!#endif

   isubs=0

  ! 

    !******************************************
    !    identification of indices 
    !******************************************

#ifdef MUSTANG
   ! identification of igrav1, igrav2, isand1, isand2, imud1, imud2
   ! ---------------------------------------------------------
   igrav1 = 1
   igrav2 = nv_grav
   isand1 = igrav2 + 1
   isand2 = igrav2 + nv_sand
   imud1 = isand2 + 1
   imud2 = isand2 + nv_mud
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
   ! initialize without any sediment activating bedload
   ibedload1 = isand2 + 1
   ibedload2 = isand2
   l_ibedload1 = .FALSE.
   l_ibedload2 = .FALSE.

   do iv = igrav1,isand2
     if (l_bedload_r(iv) .AND. (.NOT. l_ibedload1)) then
       l_ibedload1 = .TRUE.
       ibedload1 = iv
     endif
   enddo

   do iv = ibedload1+1, isand2
      if ((l_ibedload2) .AND. l_bedload_r(iv)) then
         MPI_master_only  write(stdout,*)'ERROR - for iv = ', iv
         MPI_master_only  write(stdout,*)'ibedload1 = ', ibedload1
         MPI_master_only  write(stdout,*)'ibedload2 = ', ibedload2
         MPI_master_only  write(stdout,*)'iv is greater than ibedload2 and l_bedload is True'
         MPI_master_only  write(stdout,*)'l_bedload should be true for all sediment between ibedload1 and iv'
         MPI_master_only  write(stdout,*)'The simulation is stopped, you should review parasubstance.txt'
         may_day_flag=77
         goto 99
      endif

      if ((.NOT. l_bedload_r(iv)) .AND. (.NOT. l_ibedload2)) then
         l_ibedload2 = .TRUE.
         ibedload2 = iv - 1
      endif 

   enddo

#endif
#endif

   ! initialize the number of variables according to their type
   nvpc = nv_mud + nv_sand + nv_grav
   nvp = nvpc + nv_ncp + nv_sorb
   nv_adv = nvp + nv_dis
   nv_state = nv_adv+nv_fix
   nv_tot = nv_state

   ALLOCATE(irk_mod(nv_tot))
   ALLOCATE(irk_fil(nv_tot))
#if defined PEPTIC || defined key_N_tracer || defined key_P_tracer
   ! ordering the added variables to get the particulate variables first, then the dissolved variables ...
   ! classification of substances by type
   ! there are 6 types of variables arranged in the following order : 
   ! GRAV (1), SAND(2), MUDS+PART(3), NoCP(4), SORB(5), DISS(6)
   ! conversion de l ordre lu dans l ordre modele : irk_mod(ilu)
   ! conversion de l ordre modele dans l ordre lu : irk_fil(imod)
   !  Exemple : lecture des variables dans fichier dans l odre : MOP,Vas1,O2D,Sabl
   !            rangement dans l ordre du modele : Sabl,Vas1,MOP,O2D
   !            on obtient : irk_mod(ilu=1)=3 ; irk_mod(2)=2 ; irk_mod(3)=4 ; irk_mod(4)=1
   !            et en sens inverse : irk_fil(ivmod=1)=4;irk_fil(2)=2;irk_fil(3)=1;irk_fil(4)=3

   it=1
   DO ntypvar=1,6
     DO isubs=1,nv_adv
       IF (itypv_r(isubs) == ntypvar) THEN
         irk_mod(isubs)=it
         irk_fil(it)=isubs
         it=it+1
       END IF
     END DO
   END DO
#else
   DO isubs=1,nv_adv
     irk_fil(isubs)=isubs
     irk_mod(isubs)=isubs
   ENDDO
#endif
   DO isubs=nv_adv+1,nv_tot
     irk_fil(isubs)=isubs
     irk_mod(isubs)=isubs
   ENDDO

   IF (nv_adv /= ntrc_subs) THEN
     MPI_master_only  WRITE(stdout,*)'WARNING - the total number of substances read from the file'
     MPI_master_only  WRITE(stdout,*)' parasubstance.txt is DIFFERENT from the ntrc_subs parameter'
     MPI_master_only  WRITE(stdout,*)'in param.h  '
     MPI_master_only  WRITE(stdout,*)'ntrc_subs in param.h = ',ntrc_subs
     MPI_master_only  WRITE(stdout,*)'nv_adv read in parasubstance.txt ',nv_adv
     MPI_master_only  WRITE(stdout,*)'The simulation is stopped'
     may_day_flag=77
     goto 99
   END IF

    !******************************************
    !    tables allocation 
    !******************************************

   ! allocation des tableaux se rapportant aux variables 
   ! et dont les dimensions dependent des nombres qui viennent d etre lus
   ! ---------------------------------------------------------------------------------
   ! declaration du tableau t par CROCO dans ocean3d.h avec la dimension declaree dans param.h
   ! cette dimension GLOBAL_2D_ARRAY a ete definie dans set_global_definitions.h, lui meme introduit dans cppdefs.h
   
   ! allocation des tableaux supplementaires pour variables fixees et benthiques
   IF(nv_fix > 0) THEN
       ALLOCATE(cvfix_wat(GLOBAL_2D_ARRAY,N,nv_fix))
       cvfix_wat(:,:,:,:)=0.0
   ENDIF
#ifdef key_benthic
   IF(nv_bent > 0)ALLOCATE(cv_bent(GLOBAL_2D_ARRAY,N,nv_bent))
#endif

   ! tables for substances - without temperature and salinity                              
    ALLOCATE(obc_cv_name(itsubs1:itsubs2))
    ALLOCATE(init_cv_name(itsubs1:itsubs2))
    ALLOCATE(sub_flx_atm(itsubs1:itsubs2))
    ALLOCATE(cv_rain(itsubs1:itsubs2))
    ALLOCATE(cini_wat(itsubs1:itsubs2))
    ALLOCATE(cobc_wat(itsubs1:itsubs2))
    ALLOCATE(cini_air(itsubs1:itsubs2))
    ALLOCATE(ws_part(GLOBAL_2D_ARRAY,N,itsubs1:itsubs2))
    ALLOCATE(typdiss(itsubs1:itsubs2)) 
    ws_part(:,:,:,:)=0.0
   
   IF (nvp > 0) THEN
#ifdef MUSTANG
      ALLOCATE(ws_free_opt(imud1:nvp)) 
      ALLOCATE(ws_hind_opt(nvp)) 
      ALLOCATE(ws_free_para(4,nvp)) 
      ALLOCATE(ws_free_min(nvp)) 
      ALLOCATE(ws_free_max(nvp)) 
      ALLOCATE(ws_hind_para(2,nvp)) 
      ALLOCATE(tocd(nvp))
      ALLOCATE(ros(nvp))
      ALLOCATE(diam_sed(nvp))      
      ws_free_min(:)=0.0
      ws_free_max(:)=0.0
      ws_free_para(:,:)=0.0
      ws_hind_para(:,:)=0.0
      ws_free_opt(:)=0.
      ws_hind_opt(:)=0.
      tocd(:)=0.0
      ros(:)=0.0
      diam_sed(:)=0.0
#else
      ALLOCATE(ws_free_min(nvp)) 
      ALLOCATE(ws_free_max(nvp)) 
      ws_free_min(:)=0.0
      ws_free_max(:)=0.0
#endif
   ENDIF
   ALLOCATE(l_subs2D(-1:nvp))
   l_subs2D(:)=.false.
#if defined key_sand2D
   ALLOCATE(l_outsandrouse(nvp))
   l_outsandrouse(:)=.false.
#endif

! a priori si pour l instant on ne rajoute pas des variables supplementaires crees a partir 
! des premieres (varaibles de tracage N ou P avec BIOLO)
!  on a plus besoin de reordonner le tableau puisque qu'on rempli le tableau 
!  au fur et a mesure en lisant les namelist

!  mais avec MUSTANG, quand on aura des variables particulaires 3D uniquement pelagiques ou uniquement benthiques
!  on devra reperer les indices car il y aura des variables dans t qui ne seront pas dans cv_sed (pelagiques)
!              et des variables dans cv_sed qui ne seront pas dans t (benthiques 3D)


   ! reperage de variables particulaires  associees aux variables part. SORB
   ! ------------------------------------------------------------------------------------
   ALLOCATE(irkm_var_assoc(nvp))
   irkm_var_assoc(:)=0
#ifdef MUSTANG
   DO iv=1,nv_sorb
     isubs=nvpc+nv_ncp+iv
     irkm_var_assoc(isubs)=0
     l_varassoc=.FALSE.
     DO ivpc=1,nvpc
       IF ( TRIM(ADJUSTL(ADJUSTR(name_var(irk_fil(ivpc))))) == TRIM(ADJUSTL(ADJUSTR(name_varpc_assoc(irk_fil(isubs))))) )THEN
         irkm_var_assoc(isubs)=ivpc
         MPI_master_only  WRITE(stdout,*) ' '
         MPI_master_only  WRITE(stdout,*)'constitutive part. variable where is sorbed'
         MPI_master_only  WRITE(stdout,*)'the variable :',name_var(irk_fil(isubs))
         MPI_master_only  WRITE(stdout,*)'is the variable ', name_var(irk_fil(ivpc))
         l_varassoc=.TRUE.
         ws_free_min_r(irk_fil(isubs))=ws_free_min_r(irk_fil(ivpc))
         ws_free_max_r(irk_fil(isubs))=ws_free_max_r(irk_fil(ivpc))
         ws_free_para_r(:,irk_fil(isubs))=ws_free_para_r(:,irk_fil(ivpc))
         ws_hind_para_r(:,irk_fil(isubs))=ws_hind_para_r(:,irk_fil(ivpc))
         ws_free_opt_r(irk_fil(isubs))=ws_free_opt_r(irk_fil(ivpc))
         ws_hind_opt_r(irk_fil(isubs))=ws_hind_opt_r(irk_fil(ivpc))
         tocd_r(irk_fil(isubs))=tocd_r(irk_fil(ivpc))
         ros_r(irk_fil(isubs))=ros_r(irk_fil(ivpc))
       END IF
     END DO
     IF (.NOT.l_varassoc) THEN
       MPI_master_only   WRITE(stdout,*)' '
       MPI_master_only   WRITE(stdout,*)'the SORB variable :',name_var(irk_fil(isubs))
       MPI_master_only   WRITE(stdout,*)'does not have associated constitutive part. variable'
!       MPI_master_only   WRITE(stdout,*)'See parasubstance_MUSTANG.txt to give exactly the name of the constitutive associated variable'
       MPI_master_only   WRITE(stdout,*)'otherwise, it is not a SORB variable, but a NoCP variable '
       may_day_flag=78
       goto 99
     END IF
   END DO
#else
   DO iv=1,nv_sorb
     isubs=nv_ncp+iv
     irkm_var_assoc(isubs)=0
     l_varassoc=.FALSE.
     DO ivp=1,nv_ncp
       IF ( TRIM(ADJUSTL(ADJUSTR(name_var(irk_fil(ivp))))) == TRIM(ADJUSTL(ADJUSTR(name_varpc_assoc(irk_fil(isubs))))) )THEN
         irkm_var_assoc(isubs)=ivp
         MPI_master_only  WRITE(stdout,*) ' '
         MPI_master_only  WRITE(stdout,*)'particulate variable where is sorbed'
         MPI_master_only  WRITE(stdout,*)'the variable :',name_var(irk_fil(isubs))
         MPI_master_only  WRITE(stdout,*)'is the variable ', name_var(irk_fil(ivp))
         l_varassoc=.TRUE.
         ws_free_min_r(irk_fil(isubs))=ws_free_min_r(irk_fil(ivp))
         ws_free_max_r(irk_fil(isubs))=ws_free_max_r(irk_fil(ivp))
       END IF
     END DO
     IF (.NOT.l_varassoc) THEN
       MPI_master_only   WRITE(stdout,*)' '
       MPI_master_only   WRITE(stdout,*)'the SORB variable :',name_var(irk_fil(isubs))
       MPI_master_only   WRITE(stdout,*)'does not have associated particulate variable'
       MPI_master_only   WRITE(stdout,*)'See parasubstance.txt to give exactly the name of the associated variable'
       MPI_master_only   WRITE(stdout,*)'otherwise, it is not a SORB variable, but a NoCP variable '
       may_day_flag=78
       goto 99
     END IF
   END DO

#endif

   ! save into simu.log
   !-------------------
   MPI_master_only WRITE(stdout,*) ' '
   MPI_master_only WRITE(stdout,*) 'TRACER-SUBSTANCE NUMBER        NAME             UNIT            TYPE      '
    DO isubs=1,ntrc_subs
      MPI_master_only WRITE(stdout,'(5x,i4,5x,a30,2x,a18,5x,i4)')  &
                 isubs,TRIM(name_var_r(isubs)),TRIM(unit_var_r(isubs)),itypv_r(isubs)
    END DO
    DO isubs=1,ntrc_subs
     MPI_master_only WRITE(stdout,*)' '
     MPI_master_only WRITE(stdout,*)'VARIABLE NAME : ',TRIM(name_var_r(isubs))
     IF (itypv_r(isubs)==3 .OR. itypv_r(isubs)==4) THEN
       
#ifdef MUSTANG
       MPI_master_only WRITE(stdout,*)'Free settling velocity method : ',ws_free_opt_r(isubs)
       MPI_master_only WRITE(stdout,*)'Free settling velocity MIN and MAX : ',ws_free_min_r(isubs),ws_free_max_r(isubs)
       MPI_master_only WRITE(stdout,*)'Free settling velocity parameters : ',ws_free_para_r(1:4,isubs)
       MPI_master_only WRITE(stdout,*)'Hindered settling velocity method : ',ws_hind_opt_r(isubs)
       MPI_master_only WRITE(stdout,*)'Hindered settling velocity parameters : ',ws_hind_para_r(1:2,isubs)
#else
       MPI_master_only WRITE(stdout,*)'Settling velocity MIN and MAX : ',ws_free_min_r(isubs),ws_free_max_r(isubs)
#endif
#ifdef MUSTANG
     ELSE IF (itypv_r(isubs) == 5) THEN
       MPI_master_only WRITE(stdout,*)'Particulate Constitutive associated Variable  : ',name_varpc_assoc(isubs)
#endif
     ENDIF
#ifdef MUSTANG
     IF (itypv_r(isubs)< 6) THEN
       MPI_master_only WRITE(stdout,*)'critical shear stress for deposit : ',tocd_r(isubs)
       MPI_master_only WRITE(stdout,*)'grain density                     : ',ros_r(isubs)
       MPI_master_only WRITE(stdout,*)'grain diameter                    : ',diam_r(isubs)
     END IF
#endif
     MPI_master_only WRITE(stdout,*)'depot atmospherique (masse/m2/seconde) : ',flx_atm_r(isubs)
     MPI_master_only WRITE(stdout,*)'concentration in rain water            : ',cv_rain_r(isubs)
     MPI_master_only WRITE(stdout,*)'uniform initial conc. in water column  : ',cini_wat_r(isubs)
     MPI_master_only WRITE(stdout,*)'OBC uniforme and constant conc.        : ',cobc_wat_r(isubs)
#ifdef MUSTANG
     MPI_master_only WRITE(stdout,*)'uniform initial conc. in sediment      : ',cini_sed_r(isubs)
#endif
     MPI_master_only WRITE(stdout,*)'uniform initial conc. in air           : ',cini_air_r(isubs)
    END DO

    IF (nv_fix > 0) THEN
      DO isubs=1,nv_fix
      MPI_master_only WRITE(stdout,*)' '
      MPI_master_only WRITE(stdout,*)'FIXED VARIABLE NAME : ',TRIM(name_var_fix(isubs))
      MPI_master_only WRITE(stdout,*)'uniform initial conc. in water column  : ',cini_wat_fix(isubs)
      END DO
    END IF
#ifdef key_benthic
    DO isubs=1,nv_bent
     MPI_master_only WRITE(stdout,*)' '
     MPI_master_only WRITE(stdout,*)'BENTHIC VARIABLE NAME : ',TRIM(name_var_bent(isubs))
    MPI_master_only  WRITE(stdout,*)'uniform initial conc.   : ',cini_wat_bent(isubs)
    END DO
#endif

    MPI_master_only WRITE(stdout,*)' '
    MPI_master_only WRITE(stdout,*)' TRACERS_SUBSTANCES'
#ifdef MUSTANG
    MPI_master_only WRITE(stdout,*)'number of GRAV                        : ',nv_grav
    MPI_master_only WRITE(stdout,*)'number of SAND                        : ',nv_sand
    MPI_master_only WRITE(stdout,*)'number of MUDS                        : ',nv_mud
    MPI_master_only WRITE(stdout,*)'number of part. var. constitutive     : ',nvpc
    MPI_master_only WRITE(stdout,*)'number of part. var. SORB             : ',nv_sorb
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
    MPI_master_only WRITE(stdout,*)' ibedload1 = ',ibedload1,' ibedload2 = ',ibedload2
#endif
#endif
    MPI_master_only WRITE(stdout,*)'number of part. var. NO constitutive  : ',nv_ncp
    MPI_master_only WRITE(stdout,*)'number of part. var. SORB             : ',nv_sorb
    MPI_master_only WRITE(stdout,*)'number of dissolved var. (DISS)       : ',nv_dis
    MPI_master_only WRITE(stdout,*)'number of FIXE                        : ',nv_fix
    MPI_master_only WRITE(stdout,*)'number of BENTHIC                     : ',nv_bent
    MPI_master_only WRITE(stdout,*)

   ! -------------------------------------------------------------------
   ! final tables
   ! -------------------------------------------------------------------
    DO iv=1,ntrc_subs
      ivr=iv+ivTS
      name_var(iv)=name_var_r(irk_fil(iv))
      standard_name_var(iv)=standard_name_var_r(irk_fil(iv))
      long_name_var(iv)=long_name_var_r(irk_fil(iv))
      unit_var(iv)=unit_var_r(irk_fil(iv))
      l_out_subs(iv)=l_out_subs_r(irk_fil(iv)) 

      cini_wat(ivr)=cini_wat_r(irk_fil(iv))
      cobc_wat(ivr)=cobc_wat_r(irk_fil(iv))
      cini_air(ivr)=cini_air_r(irk_fil(iv))
      init_cv_name(ivr)=init_cv_name_r(irk_fil(iv))
      obc_cv_name(ivr)=obc_cv_name_r(irk_fil(iv))
      cv_rain(ivr)=cv_rain_r(irk_fil(iv))
      sub_flx_atm(ivr)=flx_atm_r(irk_fil(iv))
    END DO
  ! memorisation des noms des variables
  !  remplissage du tableau vname declare dans ncscrum.h
   !write(*,*)'in substance indxT=',indxT
   DO isubs=1,ntrc_subs
     indx=indxT+ntrc_salt+isubs
     vname(1,indx)=name_var(isubs)
     MPI_master_only write(*,*)'vname(1,',indx,')=', vname(1,indx)
     vname(2,indx)=long_name_var(isubs)
     vname(3,indx)=unit_var(isubs)
     vname(4,indx)=TRIM(ADJUSTL(ADJUSTR(standard_name_var(isubs))))//', scalar, series'
     vname(5,indx)=' '
     vname(6,indx)=' '
     vname(7,indx)=' '
     IF (isubs .GT. igrav2) THEN
         wrthis(indx) = l_out_subs(isubs) 
     ELSE
         wrthis(indx) = .FALSE. !! no output in water for gravel
     ENDIF
   ENDDO
   IF (nv_fix>0) THEN
      DO isubs=1,nv_fix
      indx=indxT+ntrc_salt+ntrc_subs+isubs
      !  write(*,*)'fix, indice vname',indx
      vname(1,indx)=name_var_fix(isubs)
      MPI_master_only write(*,*)'vname(1,',indx,')=', vname(1,indx)
      vname(2,indx)=long_name_var_fix(isubs)
      vname(3,indx)=unit_var_fix(isubs)
      vname(4,indx)=TRIM(ADJUSTL(ADJUSTR(standard_name_var_fix(isubs))))//', scalar, series'
      vname(5,indx)=' '
      vname(6,indx)=' '
      vname(7,indx)=' '
      wrthis(indx)=l_out_subs_fix(isubs) 
      ! MPI_master_only write(*,*)'fix, indice wrthis fixed variables',indx,wrthis(indx),name_var_fix(isubs)
      ENDDO
   END IF
     

#ifdef PSOURCE_NCFILE_TS
   DO isubs=1,ntrc_subs
          indx=indxT+ntrc_salt+isubs
          vname(1,indxTsrc+isubs+1)=trim(ADJUSTL(ADJUSTR(vname(1,indx) )))//'_src        '
          vname(2,indxTsrc+isubs+1)='Tracer source concentration     '
          vname(3,indxTsrc+isubs+1)=trim(ADJUSTL(ADJUSTR(vname(3,indx) )))//''
          vname(4,indxTsrc+isubs+1)='                                '
          vname(5,indxTsrc+isubs+1)='                                '
          vname(6,indxTsrc+isubs+1)='                                '
          vname(7,indxTsrc+isubs+1)='                                '
!     write(*,*)'Traceurs SUBSTANCE:',vname(1,indxTsrc+ntrc_salt+isubs)
    enddo
#endif

#ifdef MUSTANG
   ! ---------------------------------------------------------
   ! SEDIMENTOLOGY : module MUSTANG
   ! ---------------------------------------------------------
   DO iv=1,isand2
     ros(iv)=ros_r(irk_fil(iv))
     diam_sed(iv)=diam_r(irk_fil(iv))
   ENDDO
   DO iv=imud1,nvp
     ws_free_min(iv)=ws_free_min_r(irk_fil(iv))
     ws_free_max(iv)=ws_free_max_r(irk_fil(iv))
     ws_part(:,:,:,ivTS+iv)=ws_free_max(iv)
     ws_free_opt(iv)=ws_free_opt_r(irk_fil(iv))
     ws_free_para(:,iv)=ws_free_para_r(:,irk_fil(iv))
     ws_hind_opt(iv)=ws_hind_opt_r(irk_fil(iv))
     ws_hind_para(:,iv)=ws_hind_para_r(:,irk_fil(iv))
     tocd(iv)=tocd_r(irk_fil(iv))
     ros(iv)=ros_r(irk_fil(iv))
     diam_sed(iv)=diam_r(irk_fil(iv))
   END DO
#ifdef key_sand2D
   DO iv=igrav1,igrav2
     l_subs2D(iv)=.TRUE.
   ENDDO
   DO iv=isand1,isand2
     l_subs2D(iv)=l_sand2D_r(irk_fil(iv))
     l_outsandrouse(iv)=l_outsandrouse_r(irk_fil(iv))
   ENDDO
#endif
#else
   DO iv=1,nvp
     ws_free_min(iv)=ws_free_min_r(irk_fil(iv))
     ws_free_max(iv)=ws_free_max_r(irk_fil(iv))
     ws_part(:,:,:,ivTS+iv)=ws_free_max(iv)
   END DO
#endif

 ! dans CROCO : pas de mise en groupe des variables particulaires

   typdiss(ivTS+nvp+1:ivTS+nv_adv)=1.
   typdiss(ivTS+1:ivTS+nvp)=0.

! pour CROCO a revoir
    ALLOCATE(unit_modif_mudbio_N2dw(nv_tot))
!    ALLOCATE(l_subs2D(nv_adv))
!    l_subs2D(:)=.false.
    DO iv=1,nv_tot
        unit_modif_mudbio_N2dw(iv)=1.0
    ENDDO


!  initialization  of l_sflxsubatm
   l_subflxatm=.false.
   IF(l_subflxatm_xyt) THEN
     l_subflxatm=.true.
   ELSE
     DO iv=ivTS+1,ivTS+nv_adv
       IF (cv_rain(iv) /= 0.0 .OR. sub_flx_atm(iv) /= 0.0) THEN
         l_subflxatm=.true.
       ENDIF
     END DO
   ENDIF
 
 99 CONTINUE

  END SUBROUTINE substance_read_alloc
!!======================================================================

 SUBROUTINE ALLOC_DEFVAR(nballoc)
   INTEGER, INTENT(IN)   :: nballoc
  
  
   ALLOCATE(name_var_n(nballoc))
   ALLOCATE(long_name_var_n(nballoc))
   ALLOCATE(standard_name_var_n(nballoc))
   ALLOCATE(unit_var_n(nballoc))
   ALLOCATE(flx_atm_n(nballoc))
   ALLOCATE(cv_rain_n(nballoc))
   ALLOCATE(cini_wat_n(nballoc))
   ALLOCATE(cini_air_n(nballoc))
   ALLOCATE(cobc_wat_n(nballoc))
   ALLOCATE(l_out_subs_n(nballoc))
   ALLOCATE(init_cv_name_n(nballoc))
   ALLOCATE(obc_cv_name_n(nballoc))
   !initialization
   flx_atm_n(:)=0.
   cv_rain_n(:)=0.
   cini_wat_n(:)=0.
   cobc_wat_n(:)=0.
#ifdef MUSTANG
   ALLOCATE(cini_sed_n(nballoc))
   cini_sed_n(:)=0.
#endif
   cini_air_n(:)=0.
   init_cv_name_n(:)=''
   obc_cv_name_n(:)=''
   name_var_n(:)=''
   long_name_var_n(:)=''
   standard_name_var_n(:)=''
   unit_var_n(:)=''
   l_out_subs_n(:)=.TRUE.

  END SUBROUTINE ALLOC_DEFVAR
!!======================================================================

  SUBROUTINE DEFVAR_DEALLOC(nballoc,iv)

   INTEGER, INTENT(IN)   :: nballoc
   INTEGER, INTENT(INOUT)   :: iv

   INTEGER :: ivr
   
   DO ivr=1,nballoc
     IF(init_cv_name_n(ivr)=='') THEN
        init_cv_name_n(ivr)=name_var_n(ivr)
     ENDIF
     IF(obc_cv_name_n(ivr)=='') THEN
        obc_cv_name_n(ivr)=name_var_n(ivr)
     ENDIF
   ENDDO
   DO ivr=1,nballoc
     iv=iv+1
     name_var_r(iv)=name_var_n(ivr)
     long_name_var_r(iv)=long_name_var_n(ivr)
     standard_name_var_r(iv)=standard_name_var_n(ivr)
     unit_var_r(iv)=unit_var_n(ivr)
     flx_atm_r(iv)=flx_atm_n(ivr)
     cv_rain_r(iv)=cv_rain_n(ivr)
     cini_wat_r(iv)=cini_wat_n(ivr)
     cobc_wat_r(iv)=cobc_wat_n(ivr)
#ifdef MUSTANG
     cini_sed_r(iv)=cini_sed_n(ivr)
#endif
     cini_air_r(iv)=cini_air_n(ivr)
     l_out_subs_r(iv)=l_out_subs_n(ivr)
     init_cv_name_r(iv)=init_cv_name_n(ivr)
     obc_cv_name_r(iv)=obc_cv_name_n(ivr)
   ENDDO
   DEALLOCATE(name_var_n,long_name_var_n,standard_name_var_n,unit_var_n)
   DEALLOCATE(flx_atm_n,cv_rain_n)
   DEALLOCATE(cini_wat_n,cini_air_n,cobc_wat_n)
   DEALLOCATE(l_out_subs_n)
   DEALLOCATE(init_cv_name_n,obc_cv_name_n)
#ifdef MUSTANG
   DEALLOCATE(cini_sed_n)
#endif

  END SUBROUTINE DEFVAR_DEALLOC

!!======================================================================

 SUBROUTINE substance_surfcell
     !!-------------------------------------------------------------------
     !!                    *** ROUTINE substance_surfcell ***
     !!-------------------------------------------------------------------
     !
! evaluation of cell surface if not known in hydro model
    ALLOCATE(surf_cell(GLOBAL_2D_ARRAY))
    surf_cell(:,:)=om_r(:,:)*on_r(:,:)

 END SUBROUTINE substance_surfcell

#endif /* ifdef SUBSTANCE */
!!======================================================================

END MODULE substance
