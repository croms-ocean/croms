!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
!======================================================================
!
!> @brief Croco interface for Online Analysis (OA=OnlineAnalysis)
!! - Prepares spectral and wavelet analysis outputs.
!!
!! @details called in output_oa.F, scalogram only
!!
!> @authors  
!! - B. Lemieux-Dudon
!!  - based on a preliminary Croco-OA interface version : spring 2020.
!! - More history (authors, comments) in source module_interface_oa.F90
!> @todo 
!
!  REFERENCE:
!  
!======================================================================

#include "cppdefs.h"
! BLXD includes :
! - cppdefs_dev.h
! - set_global_definitions.h
!   with GLOBAL_2D_ARRAY
!        PRIVATE_SCRATCH_2D_ARRAY... 
#if defined ONLINE_ANALYSIS && defined SOLVE3D

      subroutine scal0d_oa_out_loc (imm 
     &                             , scal0d_re, scal0d_im
     &                             , vnam_oa_r, vnam_oa_i
     &                             , per0d
     &                             , pnam_oa
     &                             , inam_oa, jnam_oa
     &                             , knam_oa
     &                             , if_imag_part )

     !& ,scal0d_oa, scal0d_cr,verbose_oa, io_hist, directory_out_oa
     !& ,scalogram_analysis

      use module_interface_oa, only : !scal0d_oa, per0d_oa
     &  mdl, model_mode_loc
     & ,tvar_oa , tgv_oa, tgvnam_oa, nzvc_oa, nzupd0d_oa
     & ,nper_sclg_max, pervnam
     & ,nsclg_loc, v2locsclg, verbose_oa
      use module_oa_variables,  only : tv_oa     
     & ,tupd_oa, tv_sclg_oa, cnb_oa
      use module_oa_periode,    only : begc_oa, tvc_oa, tc_oa
     & ,nzc_oa, swt_wfpf_oa, nper_sclg
      !use module_oa_space,      only : ptij_oa
      use scalars
      implicit none
! >> Include with variable declaration
! # include "param.h"
      integer, intent(in) :: imm 
      real(8), intent(inout) :: 
     &        scal0d_re(1:nper_sclg_max(imm),nsclg_loc(imm))
     &       ,scal0d_im(1:nper_sclg_max(imm),nsclg_loc(imm))
     &       ,per0d    (1:nper_sclg_max(imm),nsclg_loc(imm))

      logical, intent(inout) :: if_imag_part
      character(len=9), intent(inout)   :: vnam_oa_r, vnam_oa_i 
      character(len=5), intent(inout)   :: pnam_oa 
      character(len=7), intent(inout)   :: inam_oa, jnam_oa
     &                                    ,knam_oa

! Local variables
      integer :: iv, ic, isclg
      integer :: ip, ivc, io_node
      integer, parameter :: verbl=6

      type(model_mode_loc), pointer :: pdl => null()


      if_imag_part = .false.
      isclg=min(1,nsclg_loc(imm))
      scal0d_re(1:nper_sclg_max(imm),isclg:nsclg_loc(imm)) = 0. 
      scal0d_im(1:nper_sclg_max(imm),isclg:nsclg_loc(imm)) = 0.
      per0d    (1:nper_sclg_max(imm),isclg:nsclg_loc(imm)) = 0.

      oa_config_loop : do ic = 1, nzc_oa

        ! BLXD TODO debugging test if cppkeys      
        ! ==>
        ivc = 0
        do iv = begc_oa(ic), begc_oa(ic+1) -1
            ivc = ivc + 1
        end do
        oa_stop1 : if ( ivc /= nzvc_oa( tc_oa(ic) ) ) then
           MPI_master_only write(stdout,fmt='(a)') '...OA-ERR-iv-ic'
           stop
           !MPI_master_only write(stdout,*) 'Expected total # of var nzvc_oa( tc_oa(ic) ) ',nzvc_oa( tc_oa(ic) )
        end if oa_stop1
        ! <==

        !could enable to retrieve a name, requires to use tupd_oa
        ! # of variable per configuration ic
        ! do iv = 1, nzvc_oa( tc_oa(ic) )
        ivc = 0
        the_oa_var_loop : do iv = begc_oa(ic), begc_oa(ic+1) -1

          is_a_scalogram : if ( tv_sclg_oa(iv) 
     &                  .and. ( cnb_oa(iv)==imm ) ) then 

        !ivc = nzvc_oa( tc_oa(ic) )

           ivc = ivc + 1
 
        ! BLXD TODO debugging test if cppkeys      
        ! ==>
        oa_stop2 : if ( 
     &     tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ) /= tupd_oa(iv) ) 
     &     then
           MPI_master_only write(stdout,fmt='(a)') '...OA-ERR tupd_oa'
           stop
        end if oa_stop2
        ! <==

#ifdef OA_TRACES_XIOS             
        if ( v2locsclg(iv)/=-99 ) then
        write(*,fmt='(a,4(i3))') 
     &      '>><<OA-Scalogram ALL iv node/ic/iv/begc_oa '
     &      ,mynode,ic,iv,begc_oa(ic),(begc_oa(ic+1)-1)
        endif
#endif

#ifdef OA_TRACES
#if defined MILES || defined IGW || defined TANK
        if ( v2locsclg(iv)/=-99 ) call traces_oa1(verbose_oa)
#endif
#endif

                  ! iv <-> ONE SCALOGRAM 
                  ! - complex 2D array scal0d_oa( 1: nper_sclg_max, tupd_oa(iv))
                  ! - dumped in 1D array w_re(1: nper_sclg_max)


        pdl   => mdl(imm) 

#ifdef OA_TRACES_XIOS             
            write(*,fmt='(a,3(i3))') 
     &      '>><<OA-Scalogram IS SCLG node/ic/ip '
     &      ,mynode,ic,iv
#endif

        if_oa_re_im : if ( swt_wfpf_oa(iv)==4 ) then

      ! #BLXD 2020 REAL/DIMAG for real/imag part of dble prec scal0d_oa 
      ! Warning mpc.F parser REAL -> REAL real_type
      ! Warning : different gnu, intel standards 

        ! If one single variable is imag all MUST be
        if_imag_part = .true.

        if ( v2locsclg(iv)/=-99 ) then
        loop_per_re_im : do ip=1,nper_sclg(iv)
#ifdef OA_TRACES_XIOS             
            MPI_master_only write(stdout,fmt='(a,2(i3))') 
     &      '...OA-Scalogram Part/ImPart for period ',ip,iv
            write(*,fmt='(a,4(i3))') 
     &      '>><<OA-Scalogram node/ic/iv/ip '
     &      ,mynode,ic,iv,ip
#endif
                       scal0d_re(ip,v2locsclg(iv)) = 
     & REAL(
     & DBLE( 
     & pdl%scal0d_oa( ip, v2locsclg(iv) )))

                       scal0d_im(ip,v2locsclg(iv)) = 
     & REAL(
     & DIMAG( 
     & pdl%scal0d_oa( ip, v2locsclg(iv) )))


        enddo loop_per_re_im
        endif

            ! scal_ + real or imag
            vnam_oa_r = get_vnam_full( vnamfmt='(a5,a4)'
     &              ,vnams='real',verbose=verbose_oa)
            vnam_oa_i = get_vnam_full( vnamfmt='(a5,a4)'
     &              ,vnams='imag',verbose=verbose_oa)

        else if_oa_re_im

        if ( v2locsclg(iv) /= -99 ) then
        loop_per_re : do ip=1,nper_sclg(iv)
                       scal0d_re(ip,v2locsclg(iv)) = 
     &  pdl%scal0d_oa(ip, v2locsclg(iv) )
        enddo loop_per_re
        endif

            ! scal_ + real or imag
            vnam_oa_r = get_vnam_full( vnamfmt='(a5,a4)'
     &              ,vnams='',verbose=verbose_oa)

        endif if_oa_re_im

        if ( v2locsclg(iv) /= -99 ) then
        loop_for_per_axis : do ip=1,nper_sclg(iv)

                       per0d(ip,v2locsclg(iv)) = 
     &  pdl%per0d_oa( ip, v2locsclg(iv) )

#ifdef OA_TRACES_XIOS             
            MPI_master_only write(stdout,fmt='(a,i3,f16.2)') 
     &      '...OA-Scalogram Period axis ',ip,per0d(ip,v2locsclg(iv))
#endif

#ifdef OA_TRACES
#if defined MILES || defined IGW || defined TANK
                       if(verbose_oa>=verbl) then
                       call traces2_oa_scal0d(ip) 
                       !call traces_oa_scal0d(ip) 
                       endif
#endif
#endif

        enddo loop_for_per_axis
        endif

            pnam_oa = get_pnam_full( pnamfmt='(a5)'
     &              ,pnamdp=pervnam(imm),verbose=verbose_oa )

       inam_oa = get_vnam_full( vnamfmt='(a5,a1)'
     &              ,vnams='i',verbose=verbose_oa)
       jnam_oa = get_vnam_full( vnamfmt='(a5,a1)'
     &              ,vnams='j',verbose=verbose_oa)
       knam_oa = get_vnam_full( vnamfmt='(a5,a1)'
     &              ,vnams='k',verbose=verbose_oa)

            pdl => null()

          endif is_a_scalogram

        enddo the_oa_var_loop

      enddo oa_config_loop

      contains

      subroutine traces_oa1(verbose)        

        implicit none
        integer, intent(in)    :: verbose
        if ( verbose >=6 ) then
            MPI_master_only write(stdout,fmt='(a13,(3(i3.3,1x)))')     
     &     'Cfg-#var-Var ',ic,ivc,iv
        end if
        return
      end subroutine traces_oa1

      subroutine traces_oa_scal0d(ip)        
        implicit none
        integer, intent(in) :: ip

            if (ip<=9) then

                io_node = (10+ip)*1000+mynode
                write (io_node,fmt='(a,i3.3,2(1x,ES22.15E2))') 
     &          'ip, scal0d_oa = '
     &         ,ip
     &         ,REAL(DBLE ( pdl%scal0d_oa(ip,v2locsclg(iv)) ))
     &         ,REAL(DIMAG( pdl%scal0d_oa(ip,v2locsclg(iv)) ))
               
               io_node = (20+ip)*1000+mynode
               write (io_node,fmt='(a,i3.3,2(1x,ES22.15E2))') 
     &         'ip, scal0d_re/im = '
     &         ,ip
     &         ,scal0d_re(ip,v2locsclg(iv))
     &         ,scal0d_im(ip,v2locsclg(iv))

            endif
        return
      end subroutine traces_oa_scal0d

      subroutine traces2_oa_scal0d(ip)        
        implicit none
        integer, intent(in) :: ip

            if (ip<=9) then

                io_node = (10+ip)*1000+mynode
                write (io_node,fmt='(a,f16.2,2(1x,ES22.15E2))') 
     &          'ip, scal0d_oa = '
     &         ,pdl%per0d_oa(             ip,v2locsclg(iv))
     &         ,REAL(DBLE ( pdl%scal0d_oa(ip,v2locsclg(iv)) ))
     &         ,REAL(DIMAG( pdl%scal0d_oa(ip,v2locsclg(iv)) ))
               
               io_node = (20+ip)*1000+mynode
               write (io_node,fmt='(a,f16.2,2(1x,ES22.15E2))') 
     &         'ip, scal0d_re/im = '
     &         ,per0d(    ip,v2locsclg(iv))
     &         ,scal0d_re(ip,v2locsclg(iv))
     &         ,scal0d_im(ip,v2locsclg(iv))

            endif
        return
      end subroutine traces2_oa_scal0d

      function get_vnam_full( vnamfmt, vnams,verbose)

        implicit none
        character(len=*), intent(in) :: vnamfmt, vnams
        integer, intent(in)           :: verbose
        character(len=9)  :: get_vnam_full, tmpvnam !, tmpfmt 
        
        write( tmpvnam, fmt=trim(vnamfmt)) 
     &         tgvnam_oa(tv_oa(iv)),trim(vnams)

#ifdef OA_TRACES       
        if (verbose >=16 ) then
         MPI_master_only write(stdout,fmt='(a9)') 
     &   tmpvnam
         MPI_master_only write(stdout,fmt='(i3.3)') 
     &   tvar_oa(tc_oa(ic),tvc_oa(ic),ivc)
         MPI_master_only write(stdout,fmt='(3(i3.3,1x))') 
     &   tc_oa(ic), tvc_oa(ic), ivc
        end if
#endif

        get_vnam_full = tmpvnam 

        return
      end function get_vnam_full

      function get_pnam_full( pnamfmt
     &                  ,pnamdp,verbose)

        implicit none
        character(len=*), intent(in) :: pnamfmt, pnamdp
        integer, intent(in)           :: verbose
        character(len=5)  :: get_pnam_full, tmppnam !, tmpfmt 

        write( tmppnam, fmt=trim(pnamfmt)) 
     &         trim(pnamdp)

        get_pnam_full = tmppnam 

      end function get_pnam_full

      function get_pnam( pnamfmt
     &                  ,pnamdp,verbose)

        implicit none
        character(len=*), intent(in) :: pnamfmt, pnamdp
        integer, intent(in)           :: verbose
        character(len=22)  :: get_pnam, tmppnam !, tmpfmt 

        write( tmppnam, fmt=trim(pnamfmt)) 
     &         trim(pnamdp),ic,'_',iv

        get_pnam = tmppnam 

      end function get_pnam

!     function get_vnam( vnamfmt='(a5,a5,i3.3,a1,i3.3,a5)'
!    &                  ,vnamdp='3d_r_',vnams='_real', verbose=6)
      function get_vnam( vnamfmt
     &                  ,vnamdp,vnams,verbose)

        implicit none
        character(len=*), intent(in) :: vnamfmt, vnamdp, vnams
        integer, intent(in)           :: verbose
        !character(len=17)  intent(in) :: vnam_oa 
        !character(len=22)  intent(in) :: vnam_oa_r, vnam_oa_i 
        character(len=22)  :: get_vnam, tmpvnam !, tmpfmt 
        

!            write( vnam_oa_r, fmt='(a5,a5,i3.3,a1,i3.3,a5)') 
!& tgvnam_oa(tv_oa(iv)),'3d_r_',ic,'_',iv,'_real'
!            write( vnam_oa_i, fmt='(a5,a5,i3.3,a1,i3.3,a5)') 
!& tgvnam_oa(tv_oa(iv)),'3d_r_',ic,'_',iv,'_imag'
!            write( vnam_oa, fmt='(a5,a5,i3.3,a1,i3.3)') 
!& tgvnam_oa(tv_oa(iv)),'3d_r_',ic,'_',iv

        !tmpfmt=trim(vnamfmt)
        write( tmpvnam, fmt=trim(vnamfmt)) 
     &         tgvnam_oa(tv_oa(iv)),trim(vnamdp),ic,'_',iv,trim(vnams)
       
#ifdef OA_TRACES
        if (verbose >=16 ) then
         MPI_master_only write(stdout,fmt='(a22)') 
     &   tmpvnam
         MPI_master_only write(stdout,fmt='(i3.3)') 
     &   tvar_oa(tc_oa(ic),tvc_oa(ic),ivc)
         MPI_master_only write(stdout,fmt='(3(i3.3,1x))') 
     &   tc_oa(ic), tvc_oa(ic), ivc
        end if
#endif

        get_vnam = tmpvnam 

        return
      end function get_vnam
      
      end subroutine scal0d_oa_out_loc 

#else /* ONLINE_ANALYSIS */

      subroutine scal0d_oa_out_loc_empty
      return
      end subroutine scal0d_oa_out_loc_empty

#endif /* ONLINE_ANALYSIS */
