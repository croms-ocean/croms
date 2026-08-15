#include "cppdefs.h"

module submassbalance

#if defined SUBSTANCE
    USE module_substance
    USE comsubstance
#ifdef MUSTANG
    USE comMUSTANG
#endif
!======================================================================
!                   ***  MODULE  submassbalance  ***
! Ocean dynamics Bio Sed:  Mass budget in one or several zones 
!                          and fluxes threw one or several borders
!
! Adapted from MARS model in 2023
!======================================================================

implicit none

private

public submassbalance_readdomain   ! routine called by substance_read_alloc
public submassbalance_flxcum_kiv   ! routine called by step3d_t
public submassbalance_psource      ! routine called by step3d_t
public submassbalance_flxcum_ws    ! routine called by step3d_t
public submassbalance_main         ! routine called by step

! Shared module variables
integer                                            :: submassbalance_nb_close
integer                                            :: submassbalance_nb_open
integer                                            :: submassbalance_iout
integer, dimension(:,:,:), allocatable             :: submassbalance_mask_bud_in
integer, dimension(:,:,:), allocatable             :: submassbalance_mask_bord_north 
integer, dimension(:,:,:), allocatable             :: submassbalance_mask_bord_east
integer, dimension(:), allocatable                 :: submassbalance_iclose
real(kind = rlg)                                   :: submassbalance_t_out  ! time output budgets
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_bil
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_stok_sed
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_stok_wat
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_stok_wat_fix
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flx_obc_cum
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flx_obczon_cum
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flx_bdl_cum
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flx_bdlzon_cum
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flx_in_cum
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flx_ws_cum
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_stokw_total
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_stoks_total
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_stokwfix_total
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flxobc_total
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flxobczon_total
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flxbdl_total
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flxbdlzon_total
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flxin_total
real(kind = rlg), dimension(:,:), allocatable      :: submassbalance_flxws_total
character(len = lchain), dimension(:), allocatable :: submassbalance_name   

integer :: submassbalance_ncid
integer :: submassbalance_t_dimid, submassbalance_t_varid
integer :: submassbalance_next_record
integer :: submassbalance_border_dimid, submassbalance_border_varid
integer :: submassbalance_budget_dimid, submassbalance_budget_varid
integer :: submassbalance_tracer_dimid, submassbalance_tracer_varid
integer :: submassbalance_tracer_fix_dimid, submassbalance_tracer_fix_varid
integer :: submassbalance_border_flux_varid, submassbalance_budget_tot_varid
integer :: submassbalance_budget_stwat_varid, submassbalance_budget_fix_varid
integer :: submassbalance_budget_flux_obc_varid, submassbalance_budget_flux_in_varid
#ifdef MUSTANG
integer :: submassbalance_budget_stsed_varid, submassbalance_budget_flux_ws_varid
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
integer :: submassbalance_border_flux_bdl_varid, submassbalance_budget_flux_bdl_varid
#endif
#endif

contains

!==========================================================================
subroutine submassbalance_readdomain()
    !----------------------------------------------------------------------
    ! Purpose : read file describing zones and border 
    !           where we want to compute budgets and fluxes
    ! Description : usefull for dissolved, biologic, 
    !           sedimentologic... variables
    ! Called by : substance_read_alloc
    !----------------------------------------------------------------------

    ! Local declarations
    logical :: ex, l_border_close
    integer :: eof, ibord_r, iborder, ind_white, ib, ism1,    &
               nb_seg, is, izon_bud,                          &
               ig_bord, id_bord, jb_bord, jh_bord, i, j, ifile
    integer, dimension(:), allocatable  :: i1, i2, i3
    integer, dimension(1:LLm, 1:MMm)    :: mask_bud_E_tmp, mask_bud_W_tmp
    integer, dimension(1:LLm, 1:MMm)    :: mask_bud_S_tmp, mask_bud_N_tmp
    integer, dimension(:,:,:), allocatable      :: mask_bud_in_tmp
    integer, dimension(:,:,:), allocatable      :: mask_bord_N_tmp
    integer, dimension(:,:,:), allocatable      :: mask_bord_E_tmp
    character(len = lchain)                     :: name_read
    character(len = 5)                          :: comment
    character(len = 1), dimension(:), allocatable :: rose

    ! Executable part
    MPI_master_only  write(stdout, *) ' '
    MPI_master_only  write(stdout, *) ' '
    MPI_master_only  write(stdout, *) '**************************************************'
    MPI_master_only  write(stdout, *) '*************** submassbalance *******************'
    MPI_master_only  write(stdout, *) '**************************************************'
    MPI_master_only  write(stdout, *) ' '
    MPI_master_only  write(stdout, *) 'file defining budget zones and borders for estimate fluxes :'
    MPI_master_only  write(stdout, *) trim(submassbalance_input_file)

    if (trim(submassbalance_input_file) == "" .or. trim(submassbalance_input_file) == "all_domain") then
        ifile = 0  ! file "submassbalance_input_file" is not read
        MPI_master_only  write(stdout, *) 'only one budget zone is taking into account (all the domain)'
        MPI_master_only  write(stdout, *) 'and any fluxes threw open borders '
        submassbalance_nb_close=1
        submassbalance_nb_open=0
    else    
        ! first reading   
        eof = 0
        inquire(file = submassbalance_input_file, exist = ex)
        if (ex) then
        ifile = 1  ! read file "submassbalance_input_file" 
        open(49, file = submassbalance_input_file, form = 'formatted')
        comment = 'debut'
        do while (comment /= '*****')
            read(49,'(a)', iostat = eof) comment
        enddo
        submassbalance_nb_close = 0
        submassbalance_nb_open = 0
        do iborder = 1, submassbalance_nb_border
            read(49,*, iostat = eof) ibord_r ! number of the border
            if (eof == 0) then
            read(49,*) ! name of the border
            read(49,*) l_border_close
            if (l_border_close) then
                submassbalance_nb_close = submassbalance_nb_close + 1
            else
                submassbalance_nb_open = submassbalance_nb_open + 1
            endif
            read(49,*) nb_seg
            do is = 1, nb_seg
                read(49,*)
            enddo
            read(49,*)
            else
            MPI_master_only  write(stdout, *) 'ERROR : the total number of border for the '
            MPI_master_only  write(stdout, *) 'estimation of fluxes and budgets is = ',submassbalance_nb_border
            MPI_master_only  write(stdout, *) '(read in parasubstance file)'
            MPI_master_only  write(stdout, *) 'It is not the good number in the data file = ',trim(submassbalance_input_file)
            stop
            endif
        enddo
        close(49)
        else
            MPI_master_only write(stdout, *) 'ERROR : the data file for budget does not exist'
            MPI_master_only write(stdout, *) 'you have given submassbalance_l=.true. in parasubs.txt'
            MPI_master_only write(stdout, *) 'and the name for the data file =',trim(submassbalance_input_file)
            MPI_master_only write(stdout, *) 'build or rename your data file'
            MPI_master_only write(stdout, *) 'or set submassbalance_l=.false.' 
            MPI_master_only write(stdout, *) 'if you don t want the calculation of budget'
            stop
        endif  ! test on exist file  
    endif   ! test on file or compute on all domain

    ! allocate and initialize masks, stock and fluxes
    allocate(submassbalance_name(submassbalance_nb_border))
    allocate(submassbalance_mask_bord_north(submassbalance_nb_border, GLOBAL_2D_ARRAY))
    allocate(submassbalance_mask_bord_east(submassbalance_nb_border, GLOBAL_2D_ARRAY))
    allocate(mask_bord_N_tmp(submassbalance_nb_border, 1:LLm, 1:MMm))
    allocate(mask_bord_E_tmp(submassbalance_nb_border, 1:LLm, 1:MMm))
    if (submassbalance_nb_open > 0) then
        allocate(submassbalance_flx_obc_cum(submassbalance_nb_open, 1:nv_adv))
        allocate(submassbalance_flxobc_total(submassbalance_nb_open, 1:nv_adv))
#ifdef MUSTANG
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
        allocate(submassbalance_flx_bdl_cum(submassbalance_nb_open, 1:nv_adv))
        allocate(submassbalance_flxbdl_total(submassbalance_nb_open, 1:nv_adv))
#endif
#endif
    endif
    if (submassbalance_nb_close > 0) then
        allocate(submassbalance_iclose(submassbalance_nb_close))
        allocate(submassbalance_mask_bud_in(submassbalance_nb_close, GLOBAL_2D_ARRAY))
        allocate(mask_bud_in_tmp(submassbalance_nb_close, 1:LLm, 1:MMm))
        allocate(submassbalance_flx_in_cum(submassbalance_nb_close, 1:nv_adv))
        allocate(submassbalance_flx_obczon_cum(submassbalance_nb_close, 1:nv_adv))
        allocate(submassbalance_flxin_total(submassbalance_nb_close, 1:nv_adv))
        allocate(submassbalance_flxobczon_total(submassbalance_nb_close, 1:nv_adv))
#ifdef MUSTANG
        allocate(submassbalance_stok_sed(submassbalance_nb_close, 1:nv_adv))
        allocate(submassbalance_flx_ws_cum(submassbalance_nb_close, 1:nv_adv))
        allocate(submassbalance_stoks_total(submassbalance_nb_close, 1:nv_adv))
        allocate(submassbalance_flxws_total(submassbalance_nb_close, 1:nv_adv))
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
        allocate(submassbalance_flx_bdlzon_cum(submassbalance_nb_close, 1:nv_adv))
        allocate(submassbalance_flxbdlzon_total(submassbalance_nb_close, 1:nv_adv))
#endif
#endif
        allocate(submassbalance_stok_wat(submassbalance_nb_close, 1:nv_adv))
        if(nv_fix > 0 ) then
            allocate(submassbalance_stok_wat_fix(submassbalance_nb_close, 1:nv_fix))
            allocate(submassbalance_stokwfix_total(submassbalance_nb_close, 1:nv_fix))
        endif
        allocate(submassbalance_bil(submassbalance_nb_close, 1:nv_adv))
        allocate(submassbalance_stokw_total(submassbalance_nb_close, 1:nv_adv))
    endif

    !initialization
    if (submassbalance_nb_close > 0) then
        submassbalance_mask_bud_in(:,:,:) = 0
        mask_bud_in_tmp(:,:,:) = 0
        submassbalance_flx_in_cum(:,:) = 0.0_rlg
#ifdef MUSTANG
        submassbalance_flx_ws_cum(:,:) = 0.0_rlg
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
        submassbalance_flx_bdlzon_cum(:,:) = 0.0_rlg
#endif
#endif
        submassbalance_flx_obczon_cum(:,:) = 0.0_rlg
    endif
    if (submassbalance_nb_open > 0) then
        submassbalance_flx_obc_cum(:,:) = 0.0_rlg
#ifdef MUSTANG
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
        submassbalance_flx_bdl_cum(:,:) = 0.0_rlg
#endif
#endif
    endif

    submassbalance_mask_bord_north(:,:,:) = 0
    submassbalance_mask_bord_east(:,:,:) = 0
    mask_bord_N_tmp(:,:,:) = 0
    mask_bord_E_tmp(:,:,:) = 0
    submassbalance_dtout = submassbalance_dtout * 3600.0_rlg
    submassbalance_t_out = max(submassbalance_tdeb, start_time)

    ! make mask variables
    if (ifile == 0) then ! no file is reading, 1 budget zone = all domain
        iborder = 1
        submassbalance_name(iborder) = "all_the_domain"
        izon_bud = 1
        submassbalance_iclose(izon_bud) = 1
        mask_bud_in_tmp(1, 1:LLm-1, 1:MMm-1) = 1       ! all cells
        mask_bord_N_tmp(1, 1:LLm, MMm) = -1   ! north border
        mask_bord_N_tmp(1, 1:LLm, 1) = 1      ! south border
        mask_bord_E_tmp(1, 1, 1:MMm) = 1      ! west border
        mask_bord_E_tmp(1, LLm, 1:MMm) = -1   ! east border

    else ! read file 

        open(49,file = submassbalance_input_file, form = 'formatted')
        comment = 'debut'
        do while (comment /= '*****')
            read(49,'(a)',iostat=eof) comment
        enddo
        izon_bud=0
        ! FIRST open borders
        do iborder = 1, submassbalance_nb_open
            read(49, *, iostat = eof) ibord_r
            read(49, '(a)', iostat = eof) name_read
            ind_white = INDEX(name_read, ' ')
            submassbalance_name(iborder) = trim(ADJUSTL(ADJUSTR(name_read(1:ind_white))))
            read(49,*) l_border_close
            if (l_border_close) then
MPI_master_only  write(stdout, *) 'ERROR in the data file for budget :',trim(submassbalance_input_file)
MPI_master_only  write(stdout, *) 'you must give first all the open borders '
MPI_master_only  write(stdout, *) 'before closed borders '
stop   
            endif
            read(49,*) nb_seg
            allocate(i1(nb_seg), i2(nb_seg), i3(nb_seg))
            allocate(rose(nb_seg))
            do is = 1, nb_seg
                read(49,*) i1(is), i2(is), i3(is), rose(is)       
            enddo
            do is = 1, nb_seg
                ism1 = is - 1
                if (rose(is) == 'N') then
                    ! verification
                    if(is > 1) then
                        if( rose(ism1) /= 'N' .AND. (i1(is) .NE. i1(ism1) .OR. i3(is) .NE. i3(ism1))) then
    MPI_master_only  write(stdout, *)'Definition of open border number :',iborder,' name ',trim(submassbalance_name(iborder))
    MPI_master_only  write(stdout, *)'segment NORD number :',is,'i1,i2,i3 ',i1(is),i2(is),i3(is)
    MPI_master_only  write(stdout, *)'is not contiguous to the previous one whom i1,i2,i3 are :', i1(ism1),i2(ism1),i3(ism1)
    MPI_master_only  write(stdout, *)'i1 and i3 of segment',is,'must be = to i1 and i3 respectively of previous segment'
    stop
                        endif
                    endif
                    if (i3(is) >= 1 .and. i3(is) <= MMm) then
                        ig_bord = MIN0(i1(is), i2(is))
                        id_bord = MAX0(i1(is), i2(is))
                        do i = MAX0(1, ig_bord), MIN0(LLm, id_bord-1)
                            mask_bord_N_tmp(iborder, i, i3(is)) = -1   
                        enddo
                    endif
                else if (rose(is) == 'S') then
                    ! verification
                    if(is > 1) then
                        if(rose(ism1) /= 'S' .AND. (i1(is) .NE. i1(ism1) .OR. i3(is) .NE. i3(ism1))) then
    MPI_master_only write(stdout, *)'Definition of open border number :',iborder,' name ',trim(submassbalance_name(iborder))
    MPI_master_only write(stdout, *)'segment SOUTH number :',is,'i1,i2,i3 ',i1(is),i2(is),i3(is)
    MPI_master_only write(stdout, *)'is not contiguous to the previous one whom i1,i2,i3 are :', i1(ism1),i2(ism1),i3(ism1)
    MPI_master_only write(stdout, *)'i1 and i3 of segment',is,'must be = to i1 and i3 respectively of previous segment'
    stop
                        endif
                    endif
                    if (i3(is) >= 1 .and. i3(is) <= MMm) then
                        ig_bord = MIN0(i1(is), i2(is))
                        id_bord = MAX0(i1(is), i2(is))
                        do i = MAX0(1, ig_bord), MIN0(LLm, id_bord-1)
                            mask_bord_N_tmp(iborder, i, i3(is)) = 1
                        enddo
                    endif
                else if (rose(is) == 'E') then
                    ! verification
                    if(is > 1) then
                        if(rose(ism1) /= 'E' .AND. (i1(is) .NE. i2(ism1) .OR. i2(is) .NE. i3(ism1))) then
    MPI_master_only write(stdout, *)'Definition of open border number :',iborder,' name ',trim(submassbalance_name(iborder))
    MPI_master_only write(stdout, *)'segment East number :',is,'i1,i2,i3 ',i1(is),i2(is),i3(is)
    MPI_master_only write(stdout, *)'is not contiguous to the previous one whom i1,i2,i3 are :', i1(ism1),i2(ism1),i3(ism1)
    MPI_master_only write(stdout, *)'i1 and i2 of segment',is,'must be = to i2 and i3 respectively of previous segment'
    stop
                        endif
                    endif
                    if (i1(is) >= 1 .and. i1(is) <= LLm) then
                        jb_bord = MIN0(i3(is), i2(is))
                        jh_bord = MAX0(i3(is), i2(is))
                        do j = MAX0(1, jb_bord), MIN0(MMm, jh_bord-1)
                            mask_bord_E_tmp(iborder, i1(is), j) = -1
                        enddo
                    endif
                else if (rose(is) == 'W') then
                    ! verification
                    if(is > 1) then
                        if(rose(ism1) /= 'W' .AND. (i1(is) .NE. i2(ism1) .OR. i2(is) .NE. i3(ism1))) then
    MPI_master_only write(stdout, *)'Definition of open border number :',iborder,' name ',trim(submassbalance_name(iborder))
    MPI_master_only write(stdout, *)'segment West number :',is,'i1,i2,i3 ',i1(is),i2(is),i3(is)
    MPI_master_only write(stdout, *)'is not contiguous to the previous one whom i1,i2,i3 are :', i1(ism1),i2(ism1),i3(ism1)
    MPI_master_only write(stdout, *)'i1 and i2 of segment',is,'must be = to i2 and i3 respectively of previous segment'
    stop
                        endif
                    endif
                    if (i1(is) >= 1 .and. i1(is) <= LLm) then
                        jb_bord = MIN0(i3(is), i2(is))
                        jh_bord = MAX0(i3(is), i2(is))
                        do j = MAX0(1, jb_bord), MIN0(MMm, jh_bord-1)
                            mask_bord_E_tmp(iborder, i1(is), j) = 1
                        enddo
                    endif
                endif
            enddo
            deallocate(i1, i2, i3, rose)
            read(49,*)
        enddo

        ! SECOND closed border
        do iborder = submassbalance_nb_open + 1, submassbalance_nb_border
            mask_bud_S_tmp(:,:) = 0
            mask_bud_N_tmp(:,:) = 0
            mask_bud_W_tmp(:,:) = 0
            mask_bud_E_tmp(:,:) = 0
            read(49, *, iostat = eof) ibord_r
            read(49,'(a)', iostat = eof) name_read
            ind_white = INDEX(name_read, ' ')
            submassbalance_name(iborder) = trim(ADJUSTL(ADJUSTR(name_read(1:ind_white))))
            read(49,*) l_border_close
            izon_bud = izon_bud+1
            submassbalance_iclose(izon_bud) = iborder
            read(49,*) nb_seg
            allocate(i1(nb_seg), i2(nb_seg), i3(nb_seg))
            allocate(rose(nb_seg))
            do is=1,nb_seg
                read(49,*) i1(is), i2(is), i3(is), rose(is)
            enddo
            do is = 1, nb_seg
                ism1 = is-1
                if (is == 1) ism1 = nb_seg
                if (rose(is) == 'S') then
                ! verification
                if(rose(ism1) /= 'S' .AND. (i1(is) .NE. i1(ism1) .OR. i3(is) .NE. i3(ism1))) then
    MPI_master_only write(stdout, *)'Definition of closed border number :',iborder,' name ',trim(submassbalance_name(iborder))
    MPI_master_only write(stdout, *)'segment SOUTH number :',is,'i1,i2,i3 ',i1(is),i2(is),i3(is)
    MPI_master_only write(stdout, *)'is not contiguous to the previous one whom i1,i2,i3 are :', i1(ism1),i2(ism1),i3(ism1)
    MPI_master_only write(stdout, *)'i1 and i3 of segment',is,'must be = to i1 and i3 respectively of previous segment'
    stop
                endif
                if (i3(is) >= 1 .and.i3(is) <= MMm) then
                    ig_bord = MIN0(i1(is), i2(is))
                    id_bord = MAX0(i1(is), i2(is))
                    do i = MAX0(1, ig_bord), MIN0(LLm, id_bord-1)
                        mask_bord_N_tmp(iborder, i, i3(is)) = 1
                        do j = 1, i3(is)-1
                            mask_bud_S_tmp(i,j) = mask_bud_S_tmp(i,j)-1
                        enddo
                        do j = i3(is), MMm
                            mask_bud_S_tmp(i,j) = mask_bud_S_tmp(i,j)+1
                        enddo
                    enddo
                endif
                else if (rose(is) == 'N') then
                ! verification
                if(rose(ism1) /= 'N' .AND. (i1(is) .NE. i1(ism1) .OR. i3(is) .NE. i3(ism1))) then
    MPI_master_only write(stdout, *)'Definition of closed border number :',iborder,' name ',trim(submassbalance_name(iborder))
    MPI_master_only write(stdout, *)'segment NORD number :',is,'i1,i2,i3 ',i1(is),i2(is),i3(is)
    MPI_master_only write(stdout, *)'is not contiguous to the previous one whom i1,i2,i3 are :', i1(ism1),i2(ism1),i3(ism1)
    MPI_master_only write(stdout, *)'i1 and i3 of segment',is,'must be = to i1 and i3 respectively of previous segment'
    stop
                endif
                if (i3(is) >= 1 .and.i3(is) <= MMm) then
                    ig_bord = MIN0(i1(is), i2(is))
                    id_bord = MAX0(i1(is), i2(is))
                    do i=MAX0(1, ig_bord),MIN0(LLm, id_bord-1)
                        mask_bord_N_tmp(iborder, i, i3(is)) = -1    
                        do j = 1, i3(is)-1
                            mask_bud_N_tmp(i,j) = mask_bud_N_tmp(i,j)+1
                        enddo
                        do j = i3(is), MMm
                            mask_bud_N_tmp(i,j) = mask_bud_N_tmp(i,j)-1
                        enddo
                    enddo
                endif
                else if (rose(is) == 'E') then
                ! verification
                if(rose(ism1) /= 'E' .AND. (i1(is) .NE. i2(ism1) .OR. i2(is) .NE. i3(ism1))) then
    MPI_master_only write(stdout, *)'Definition of closed border number :',iborder,' name ',trim(submassbalance_name(iborder))
    MPI_master_only write(stdout, *)'segment East number :',is,'i1,i2,i3 ',i1(is),i2(is),i3(is)
    MPI_master_only write(stdout, *)'is not contiguous to the previous one whom i1,i2,i3 are :', i1(ism1),i2(ism1),i3(ism1)
    MPI_master_only write(stdout, *)'i1 and i2 of segment',is,'must be = to i2 and i3 respectively of previous segment'
    stop
                endif
                if (i1(is) >= 1 .and.i1(is) <= LLm) then
                    jb_bord = MIN0(i3(is), i2(is))
                    jh_bord = MAX0(i3(is), i2(is))
                    do j = MAX0(1, jb_bord), MIN0(MMm, jh_bord-1)
                        mask_bord_E_tmp(iborder, i1(is), j) = -1
                        do i = 1, i1(is)-1
                            mask_bud_E_tmp(i,j) = mask_bud_E_tmp(i,j)+1
                        enddo
                        do i = i1(is), LLm
                            mask_bud_E_tmp(i,j) = mask_bud_E_tmp(i,j)-1
                        enddo
                    enddo
                endif
                else if (rose(is) == 'W') then
                ! verification
                if(rose(ism1) /= 'W' .AND. (i1(is) .NE. i2(ism1) .OR. i2(is) .NE. i3(ism1))) then
    MPI_master_only write(stdout, *)'Definition of closed border number :',iborder,' name ',trim(submassbalance_name(iborder))
    MPI_master_only write(stdout, *)'segment West number :',is,'i1,i2,i3 ',i1(is),i2(is),i3(is)
    MPI_master_only write(stdout, *)'is not contiguous to the previous one whom i1,i2,i3 are :', i1(ism1),i2(ism1),i3(ism1)
    MPI_master_only write(stdout, *)'i1 and i2 of segment',is,'must be = to i2 and i3 respectively of previous segment'
    stop
                endif
                if (i1(is) >= 1 .and.i1(is) <= LLm) then
                    jb_bord = MIN0(i3(is), i2(is))
                    jh_bord = MAX0(i3(is), i2(is))
                    do j = MAX0(1, jb_bord), MIN0(MMm, jh_bord-1)
                        mask_bord_E_tmp(iborder, i1(is), j) = 1
                        do i = 1, i1(is)-1
                            mask_bud_W_tmp(i,j) = mask_bud_W_tmp(i,j)-1
                        enddo
                        do i = i1(is), LLm
                            mask_bud_W_tmp(i,j) = mask_bud_W_tmp(i,j)+1
                        enddo
                    enddo
                endif
                endif
            enddo

            do i = 1, LLm
                do j = 1, MMm
                if (mask_bud_W_tmp(i,j) + mask_bud_E_tmp(i,j) == 2  .and. &
                    mask_bud_N_tmp(i,j) + mask_bud_S_tmp(i,j) == 2) then
                    mask_bud_in_tmp(izon_bud, i, j) = 1
                else
                    mask_bud_in_tmp(izon_bud, i, j) = 0
                endif
                enddo
            enddo
            deallocate(i1, i2, i3, rose)
            read(49,*)
        enddo
    endif ! test on ifile

#ifdef MPI
    !  one keeps only mask into the zone corresponding to each processor
    submassbalance_mask_bord_east(:, 1:Lmmpi, 1:Mmmpi) = mask_bord_E_tmp(:, iminmpi:imaxmpi, jminmpi:jmaxmpi)
    submassbalance_mask_bord_north(:, 1:Lmmpi, 1:Mmmpi) = mask_bord_N_tmp(:, iminmpi:imaxmpi, jminmpi:jmaxmpi)
    if (submassbalance_nb_close > 0) then
        submassbalance_mask_bud_in(:,1:Lmmpi, 1:Mmmpi) = mask_bud_in_tmp(:, iminmpi:imaxmpi, jminmpi:jmaxmpi)  
    endif
#else
    submassbalance_mask_bord_east(:,1:Lm,1:Mm) = mask_bord_E_tmp(:,:,:)
    submassbalance_mask_bord_north(:,1:Lm,1:Mm) = mask_bord_N_tmp(:,:,:)
    if (submassbalance_nb_close > 0) then
        submassbalance_mask_bud_in(:,1:Lm,1:Mm) = mask_bud_in_tmp(:,:,:)
    endif
#endif

    submassbalance_iout = -1
end subroutine submassbalance_readdomain

!==========================================================================
subroutine submassbalance_flxcum_kiv(Istr, Iend, Jstr, Jend, FX, FE, iv)
    !----------------------------------------------------------------------
    ! Purpose : flux cumulating threw each border line defined by the user 
    !           for dissolved substances
    !           this routine is in k-iv-loop
    ! Called by : step3d_t
    !----------------------------------------------------------------------

    ! Arguments
    integer, intent(in)  :: Istr, Iend, Jstr, Jend, iv
    real(kind=rsh), dimension(PRIVATE_2D_SCRATCH_ARRAY), intent(in) :: FX, FE ! in kg/s for sediment

    ! Local declarations
    integer :: i, j, iz, iborder, ib

    ! Executable part
    if(submassbalance_l .AND. (time .GE. submassbalance_tdeb)) then
        do j = Jstr, Jend
            do i = Istr, Iend
                
                do ib = 1, submassbalance_nb_open
                    submassbalance_flx_obc_cum(ib, iv) = submassbalance_flx_obc_cum(ib, iv)   &
                        + real(submassbalance_mask_bord_north(ib, i, j) * dt * FE(i, j), rlg) &
                        + real(submassbalance_mask_bord_east(ib, i, j) * dt * FX(i, j), rlg)
                enddo

                do iz = 1 , submassbalance_nb_close
                    iborder = submassbalance_iclose(iz)
                    submassbalance_flx_obczon_cum(iz, iv) = submassbalance_flx_obczon_cum(iz, iv)  &
                        + real(submassbalance_mask_bord_north(iborder, i, j) * dt * FE(i, j), rlg) &
                        + real(submassbalance_mask_bord_east(iborder, i, j) * dt * FX(i, j), rlg)
                enddo
            enddo
        enddo
    endif

end subroutine submassbalance_flxcum_kiv

!==========================================================================
subroutine submassbalance_flxcum_ws(Istr, Iend, Jstr, Jend)
    !----------------------------------------------------------------------
    ! Purpose : flux cumulating for water to sediment
    ! Called by : step3d_t
    !----------------------------------------------------------------------
    ! Arguments
    integer, intent(in)  :: Istr, Iend, Jstr, Jend
#ifdef MUSTANG
    ! Local declarations
    integer :: i, j, iz, iv, itrc, ib, iborder
    ! Executable part
    if(submassbalance_l .AND. (time .GE. submassbalance_tdeb)) then
        do j = Jstr, Jend
            do i = Istr, Iend
                do iz = 1 , submassbalance_nb_close
                    do iv = 1, nv_adv
                        itrc = iv + itsubs1 - 1
                        submassbalance_flx_ws_cum(iz, iv) = submassbalance_flx_ws_cum(iz, iv)    &
                            + real(submassbalance_mask_bud_in(iz, i, j) * surf_cell(i,j), rlg)   &
                            * real((dt * flx_s2w_CROCO(i, j, itrc)                &
                            - (1. - typdiss(itrc)) * flx_w2s_sum_CROCO(i,j,itrc)) &
                            , rlg)
                    enddo
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
                    do iv = ibedload1,ibedload2
                        iborder = submassbalance_iclose(iz)
                        submassbalance_flx_bdlzon_cum(iz, iv) = submassbalance_flx_bdlzon_cum(iz, iv) &
                            + real(submassbalance_mask_bord_north(iborder, i, j) * &
                            (min(flx_by(iv,i,j), 0.0_rlg) + max(flx_by(iv,i,j-1), 0.0_rlg)), rlg) &
                            + real(submassbalance_mask_bord_east(iborder, i, j) * &
                            (min(flx_bx(iv,i,j), 0.0_rlg) + max(flx_bx(iv,i-1,j), 0.0_rlg)), rlg)  
                            ! flx_bx and flx_by directly in kg from MUSTANGV2_borne_and_apply_erosion_tot
                            ! i-1 and j-1 possible, sed_exchange done at the end of MUSTANG_update
                    enddo 
#endif
                enddo
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
                do ib = 1 , submassbalance_nb_open
                    do iv = ibedload1,ibedload2
                        submassbalance_flx_bdl_cum(ib, iv) = submassbalance_flx_bdl_cum(ib, iv) &
                            + real(submassbalance_mask_bord_north(ib, i, j) * &
                            (min(flx_by(iv,i,j), 0.0_rlg) + max(flx_by(iv,i,j-1), 0.0_rlg)), rlg) &
                            + real(submassbalance_mask_bord_east(ib, i, j) * &
                            (min(flx_bx(iv,i,j), 0.0_rlg) + max(flx_bx(iv,i-1,j), 0.0_rlg)), rlg)   
                            ! flx_bx and flx_by directly in kg from MUSTANGV2_borne_and_apply_erosion_tot
                            ! i-1 and j-1 possible, sed_exchange done at the end of MUSTANG_update
                    enddo 
                enddo
#endif
            enddo
        enddo
    endif
# endif
end subroutine submassbalance_flxcum_ws

!==========================================================================
subroutine submassbalance_psource(Istr, Iend, Jstr, Jend)
    !----------------------------------------------------------------------
    ! Purpose : flux cumulating for river discharge
    ! Called by : step3d_t
    !----------------------------------------------------------------------
    ! Arguments
    integer, intent(in)  :: Istr, Iend, Jstr, Jend
    ! Local declarations
    integer ::  k, iis, jjs, is, iz, itrc, iv
    real    :: phi, cff2
#  if defined PSOURCE || defined PSOURCE_MASS
    if(submassbalance_l .AND. (time .GE. submassbalance_tdeb)) then
        do is = 1, Nsrc
#   ifdef MPI
            iis = Isrc_mpi(is, mynode)
            jjs = Jsrc_mpi(is, mynode)
#   else
            iis = Isrc(is)
            jjs = Jsrc(is)
#   endif
            if (Istr.le.iis .and. iis.le.Iend+1 .and. Jstr.le.jjs .and. jjs.le.Jend+1) then
                do iv = 1, nv_adv
                    phi = 0.
                    itrc = iv + itsubs1 - 1
                    do k = 1, N 
                        if (Lsrc(is, itrc)) then
                            cff2 = Tsrc(is, k, itrc)
                        else
                            cff2 = 0.
                        endif
                        phi = phi + dt * cff2 * abs(Qsrc(is,k))  ! Q could be <0 in croco even if it is entrance flowrate
                    enddo
                    do iz = 1, submassbalance_nb_close
                        submassbalance_flx_in_cum(iz,iv) = submassbalance_flx_in_cum(iz,iv)   &
                            + real(submassbalance_mask_bud_in(iz, iis, jjs) * phi , rlg)
                    enddo
                enddo
            endif
        enddo ! source is
    endif
#  endif
end subroutine submassbalance_psource

!==========================================================================
subroutine submassbalance_main(tile)
    !----------------------------------------------------------------------
    !Purpose : link between tile and imin/imax/jmin/jmax
    !Called by : step
    !----------------------------------------------------------------------
    ! Arguments
    integer :: tile
# include "compute_tile_bounds.h"

    if(submassbalance_l .and. time .ge. submassbalance_tdeb) then
        if (time .ge. submassbalance_t_out) then
            call submassbalance_comp (Istr, Iend, Jstr, Jend)
        endif
    endif
    
end subroutine submassbalance_main

!==========================================================================
subroutine submassbalance_comp(Istr, Iend, Jstr, Jend)
    !----------------------------------------------------------------------
    !Purpose : mass budget computation and writing on output files
    !Called by : sub_budget_main
    !----------------------------------------------------------------------
# ifdef MPI
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE), blank, ierr
#  ifdef XIOS
#include "mpi_cpl.h"
#  endif
# endif

    ! Arguments
    integer, intent(in)  :: Istr, Iend, Jstr, Jend

    ! Local declarations
    integer             :: iv, i, j, k, ksmin, ksmax, iz, iborder, ib, itrc
    integer             :: ierror
    real(kind=rlg)      :: htot ! total water height
    real(kind=rlg)      :: voltot ! total volume of water cell
#ifdef MUSTANG
    real(kind=rlg)      :: volsed ! volume of sediment cell
#endif

    ! Executable part
    submassbalance_t_out = time + submassbalance_dtout

    ! initialisation to 0
    if (submassbalance_nb_close > 0) then
        submassbalance_stok_wat(:, 1:nv_adv) = 0.0_rlg
        if(nv_fix > 0) submassbalance_stok_wat_fix(:, 1:nv_fix) = 0.0_rsh
#ifdef MUSTANG
        submassbalance_stok_sed(:, 1:nv_adv) = 0.0_rlg
#endif
    endif

    do j = Jstr, Jend
        do i = Istr, Iend

# ifdef MASKING
            htot = (z_w(i,j,N) + h(i,j)) * rmask(i,j)
# else
            htot = z_w(i,j,N) + h(i,j)
# endif
            if ( htot > 0.) then 
                do k = 1, N
                    !surf_cell = dx*dy ; Hz(i,j,k) =dz
                    voltot = real( surf_cell(i,j) * Hz(i,j,k) ,rlg)
                    do iv = 1, nv_adv
                        itrc = iv+itsubs1-1
                        do iz = 1, submassbalance_nb_close
                            if (l_subs2D(iv)) then
                                if(k == 1) then  ! warning, not like in MARS, everything is in cell k=1
                                    submassbalance_stok_wat(iz,iv) = submassbalance_stok_wat(iz,iv) + &
                                        real(submassbalance_mask_bud_in(iz,i,j), rlg) *  &
                                        real(t(i,j,k,nnew,itrc) * voltot, rlg)
                                endif
                            else
                                submassbalance_stok_wat(iz,iv) = submassbalance_stok_wat(iz,iv) + &
                                    real(submassbalance_mask_bud_in(iz,i,j), rlg) *  &
                                    real(t(i,j,k,nnew,itrc) * voltot, rlg)
                            endif


                        enddo     
                    enddo     
                enddo     
                if(nv_fix > 0) then
                    do iv = 1, nv_fix
                        itrc = nv_adv+iv+itsubs1-1
                        do iz = 1, submassbalance_nb_close
                            submassbalance_stok_wat_fix(iz,iv) = submassbalance_stok_wat_fix(iz,iv) + &
                                    real(submassbalance_mask_bud_in(iz,i,j), rlg) *  &
                                    real(t(i,j,1,nnew,itrc) * surf_cell(i,j), rlg)     
                        enddo
                    enddo     
                endif
            endif ! htot

#ifdef MUSTANG
# ifdef MASKING
            if ( rmask(i,j) >0) then
# endif
                ksmin = ksmi(i,j)
                ksmax = ksma(i,j)
                do k = ksmin, ksmax
                    volsed = real(surf_cell(i,j) * dzs(k,i,j), rlg)
                    do iv = 1, nvp
                        do iz=1,submassbalance_nb_close
                            submassbalance_stok_sed(iz,iv) = submassbalance_stok_sed(iz,iv) + &
                                real(submassbalance_mask_bud_in(iz,i,j), rlg) *  &
                                real(cv_sed(iv,k,i,j), rlg) * volsed
                        enddo
                    enddo
                    do iv = nvp+1, nv_adv
                        do iz=1,submassbalance_nb_close
                            submassbalance_stok_sed(iz,iv) = submassbalance_stok_sed(iz,iv) + &
                                real(submassbalance_mask_bud_in(iz,i,j), rlg) *  &
                                real(cv_sed(iv,k,i,j) * poro(k,i,j), rlg) * volsed
                        enddo
                    enddo
                enddo  
# ifdef MASKING
            endif ! rmask
# endif 
#endif
        enddo !i
    enddo !j



#ifdef MPI
    do iv = 1, nv_adv
        do iz = 1, submassbalance_nb_close
        CALL MPI_REDUCE( submassbalance_stok_wat(iz,iv), submassbalance_stokw_total(iz,iv),  & 
            1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_REDUCE( submassbalance_flx_in_cum(iz,iv), submassbalance_flxin_total(iz,iv), & 
            1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_REDUCE( submassbalance_flx_obczon_cum(iz,iv), submassbalance_flxobczon_total(iz,iv), & 
            1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        enddo
        do ib=1,submassbalance_nb_open
        CALL MPI_REDUCE( submassbalance_flx_obc_cum(ib,iv), submassbalance_flxobc_total(ib,iv), & 
            1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        enddo
    enddo
#ifdef MUSTANG
    do iv = 1 , nv_adv
        do iz = 1, submassbalance_nb_close
        CALL MPI_REDUCE( submassbalance_flx_ws_cum(iz,iv), submassbalance_flxws_total(iz,iv), 1, MPI_DOUBLE_PRECISION, & 
            MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        CALL MPI_REDUCE( submassbalance_stok_sed(iz,iv), submassbalance_stoks_total(iz,iv), 1, MPI_DOUBLE_PRECISION, & 
            MPI_SUM, 0, MPI_COMM_WORLD, ierror)
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
        CALL MPI_REDUCE( submassbalance_flx_bdlzon_cum(iz,iv), submassbalance_flxbdlzon_total(iz,iv), 1, MPI_DOUBLE_PRECISION, & 
            MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        enddo
        do ib=1,submassbalance_nb_open
            CALL MPI_REDUCE( submassbalance_flx_bdl_cum(ib,iv), submassbalance_flxbdl_total(ib,iv), 1, MPI_DOUBLE_PRECISION, & 
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
#endif
        enddo
    enddo
#endif
    if (nv_fix > 0) then
    do iv = 1 , nv_fix
        do iz=1,submassbalance_nb_close
        CALL MPI_REDUCE( submassbalance_stok_wat_fix(iz,iv),submassbalance_stokwfix_total(iz,iv), 1, MPI_DOUBLE_PRECISION, & 
            MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        if(submassbalance_stokwfix_total(iz,iv) < 1.e-14) submassbalance_stokwfix_total(iz,iv) = 0.0_rlg
        enddo
    enddo
    endif
#else /* not MPI*/
    do iv = 1, nv_adv
        do iz=1,submassbalance_nb_close
            submassbalance_stokw_total(iz,iv) = submassbalance_stok_wat(iz,iv)
            submassbalance_flxin_total(iz,iv) = submassbalance_flx_in_cum(iz,iv)
            submassbalance_flxobczon_total(iz,iv) = submassbalance_flx_obczon_cum(iz,iv)
        enddo
        do ib=1,submassbalance_nb_open
            submassbalance_flxobc_total(ib,iv) = submassbalance_flx_obc_cum(ib,iv)
        enddo
    enddo
#ifdef MUSTANG
    do iv = 1, nv_adv
        do iz=1,submassbalance_nb_close
            submassbalance_flxws_total(iz,iv) = submassbalance_flx_ws_cum(iz,iv)
            submassbalance_stoks_total(iz,iv) = submassbalance_stok_sed(iz,iv)
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
            submassbalance_flxbdlzon_total(iz,iv) = submassbalance_flx_bdlzon_cum(iz,iv)
        enddo
        do ib=1,submassbalance_nb_open
            submassbalance_flxbdl_total(ib,iv) = submassbalance_flx_bdl_cum(ib,iv)
#endif
        enddo
    enddo
#endif
    if (nv_fix > 0) then
        do iv = 1, nv_fix
        do iz = 1, submassbalance_nb_close
            submassbalance_stokwfix_total(iz,iv) = submassbalance_stok_wat_fix(iz,iv)
            if(submassbalance_stokwfix_total(iz,iv) < 1.e-14) submassbalance_stokwfix_total(iz,iv) = 0.0_rlg
        enddo
        enddo
    endif

#endif

    if (submassbalance_iout == -1 ) then
        submassbalance_iout = 0
#if defined MPI
        if (mynode .gt. 0) then
            call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1, &              
                1, MPI_COMM_WORLD, status, ierr)
        endif
#endif
        CALL submassbalance_def_outnc()
#if defined MPI
        if (mynode .lt. NNODES-1) then
            call MPI_Send (blank, 1, MPI_INTEGER, mynode+1, &        
                1, MPI_COMM_WORLD,  ierr)
        endif
#endif
    endif
    MPI_master_only CALL submassbalance_wrt_outnc()

end subroutine submassbalance_comp

!==========================================================================
subroutine submassbalance_def_outnc()
    !----------------------------------------------------------------------
    !Purpose : write mask in netcdf outputfile
    !Called by : sub_budget_readdomain
    !----------------------------------------------------------------------
    USE netcdf

    integer :: i, j, n, ib, iz, ierr, &
            x_dimid, y_dimid, x_varid, y_varid, lchain_dimid, &
            lonr_varid, latr_varid, xr_varid, yr_varid
    integer :: imin, imax, startx, countx
    integer :: jmin, jmax, starty, county
    real(kind = rlg), dimension(Lm+2) :: x
    real(kind = rlg), dimension(Mm+2) :: y
    integer, dimension(2) :: dimids2

    integer :: name_border_varid
    integer, dimension(3) :: start_border, count_border
    real(kind = rlg), dimension(:), allocatable :: border_var
    character(len = lchain), dimension(:), allocatable :: border_name
    integer ::  mask_border_N_varid, mask_border_E_varid
    integer, dimension(2) :: dimids2_border
    integer, dimension(3) :: dimids3_border, dimids3_border_t

    integer :: name_budget_varid
    integer, dimension(3) :: start_budget, count_budget
    real(kind = rlg), dimension(:), allocatable :: budget_var
    character(len = lchain), dimension(:), allocatable :: budget_name
    integer ::  mask_budget_varid, mask_budget_N_varid, mask_budget_E_varid
    integer, dimension(2) :: dimids2_budget
    integer, dimension(3) :: dimids3_budget, dimids3_budget_t

    integer :: name_tracer_varid
    real(kind = rlg), dimension(:), allocatable :: tracer_var
    character(len = lchain), dimension(:), allocatable :: tracer_name
    integer, dimension(2) :: dimids2_tracer

    integer :: name_tracer_fix_varid
    real(kind = rlg), dimension(:), allocatable :: tracer_fix_var
    character(len = lchain), dimension(:), allocatable :: tracer_fix_name
    integer, dimension(2) :: dimids2_tracer_fix

#ifdef MPI
    if (mynode.eq.0) then
#endif
    call submassbalance_check( nf90_create(submassbalance_output_file, NF90_SHARE, submassbalance_ncid) )

    call submassbalance_check( nf90_def_dim(submassbalance_ncid, 'xi_rho',   xi_rho, x_dimid) )
    call submassbalance_check( nf90_def_dim(submassbalance_ncid, 'eta_rho',  eta_rho, y_dimid) )
    call submassbalance_check( nf90_def_var(submassbalance_ncid, "xi_rho", NF90_DOUBLE, x_dimid, x_varid) )
    call submassbalance_check( nf90_def_var(submassbalance_ncid, "eta_rho", NF90_DOUBLE, y_dimid, y_varid) )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, x_varid, "units", "index x axis") )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, y_varid, "units", "index y axis") )
    dimids2 = (/ x_dimid, y_dimid /)

    call submassbalance_check( nf90_def_dim(submassbalance_ncid, 'lchain',  lchain, lchain_dimid) )
    call submassbalance_check( nf90_def_dim(submassbalance_ncid, 'time', NF90_UNLIMITED, submassbalance_t_dimid) )
    call submassbalance_check( nf90_def_var(submassbalance_ncid, "time", NF90_DOUBLE, &
        submassbalance_t_dimid, submassbalance_t_varid) )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, submassbalance_t_varid, "units", "seconds since 1900-01-01") )

#ifdef SPHERICAL
    call submassbalance_check( nf90_def_var(submassbalance_ncid, 'lon_rho', NF90_DOUBLE, dimids2, lonr_varid) )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, lonr_varid, "long_name", "longitude of RHO-points") )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, lonr_varid, "units", "degree_east") )
    call submassbalance_check( nf90_def_var(submassbalance_ncid, 'lat_rho', NF90_DOUBLE, dimids2, latr_varid) )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, latr_varid, "long_name", "latitude of RHO-points") )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, latr_varid, "units", "degree_north") )
#else
    call submassbalance_check( nf90_def_var(submassbalance_ncid, "x_rho", NF90_DOUBLE, dimids2, xr_varid) )
    call submassbalance_check( nf90_def_var(submassbalance_ncid, "y_rho", NF90_DOUBLE, dimids2, yr_varid) )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, xr_varid, "units", "meter") )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, yr_varid, "units", "meter") )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, xr_varid, "long_name", "x-locations of RHO-points") )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, yr_varid, "long_name", "y-locations of RHO-points") )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, xr_varid, "standard_name", "plane_x_coordinate") )
    call submassbalance_check( nf90_put_att(submassbalance_ncid, yr_varid, "standard_name", "plane_y_coordinate") )
#endif

    if (submassbalance_nb_open .gt. 0) then
        call submassbalance_check( nf90_def_dim(submassbalance_ncid, 'border',  &
            submassbalance_nb_open, submassbalance_border_dimid) )
        call submassbalance_check( nf90_def_var(submassbalance_ncid, "border", NF90_DOUBLE, &
            submassbalance_border_dimid, submassbalance_border_varid) )
        dimids2_border = (/ lchain_dimid, submassbalance_border_dimid/)
        dimids3_border = (/ submassbalance_border_dimid, x_dimid, y_dimid/)
        call submassbalance_check( nf90_def_var(submassbalance_ncid, 'border_name', NF90_CHAR, &
            dimids2_border, name_border_varid))
        call submassbalance_check( nf90_def_var(submassbalance_ncid, 'border_mask_N', NF90_INT, &
            dimids3_border, mask_border_N_varid))
        call submassbalance_check( nf90_def_var(submassbalance_ncid, 'border_mask_E', NF90_INT, &
            dimids3_border, mask_border_E_varid))
    endif
    if (submassbalance_nb_close .gt. 0) then
        call submassbalance_check( nf90_def_dim(submassbalance_ncid, 'budget',  &
            submassbalance_nb_close, submassbalance_budget_dimid) )
        call submassbalance_check( nf90_def_var(submassbalance_ncid, "budget", NF90_DOUBLE, &
            submassbalance_budget_dimid, submassbalance_budget_varid) )
        dimids2_budget = (/ lchain_dimid, submassbalance_budget_dimid/)
        dimids3_budget = (/ submassbalance_budget_dimid, x_dimid, y_dimid/)
        call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_name', NF90_CHAR, &
            dimids2_budget, name_budget_varid))
        call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_mask_N', NF90_INT, &
            dimids3_budget, mask_budget_N_varid))
        call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_mask_E', NF90_INT, &
            dimids3_budget, mask_budget_E_varid))
        call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_mask', NF90_INT, &
            dimids3_budget, mask_budget_varid))
    endif
    if (nv_adv .gt. 0) then
        call submassbalance_check( nf90_def_dim(submassbalance_ncid, 'tracer',  &
            nv_adv, submassbalance_tracer_dimid) )
        call submassbalance_check( nf90_def_var(submassbalance_ncid, "tracer", NF90_DOUBLE, &
            submassbalance_tracer_dimid, submassbalance_tracer_varid) )
        dimids2_tracer = (/ lchain_dimid, submassbalance_tracer_dimid/)
        call submassbalance_check( nf90_def_var(submassbalance_ncid, 'tracer_name', NF90_CHAR, &
            dimids2_tracer, name_tracer_varid))

        if (submassbalance_nb_open .gt. 0) then
            dimids3_border_t = (/ submassbalance_border_dimid, submassbalance_tracer_dimid, submassbalance_t_dimid/)
            call submassbalance_check( nf90_def_var(submassbalance_ncid, 'border_flux', NF90_DOUBLE, &
                dimids3_border_t, submassbalance_border_flux_varid))
            call submassbalance_check( nf90_put_att(submassbalance_ncid, submassbalance_border_flux_varid, &
                "description", "FLux through line") )
#ifdef MUSTANG
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
            call submassbalance_check( nf90_def_var(submassbalance_ncid, 'border_flux_bedload', NF90_DOUBLE, &
                dimids3_border_t, submassbalance_border_flux_bdl_varid))
            call submassbalance_check( nf90_put_att(submassbalance_ncid, submassbalance_border_flux_bdl_varid, &
                "description", "FLux from bedload transport") )
#endif
#endif
        endif
        if (submassbalance_nb_close .gt. 0) then
            dimids3_budget_t = (/ submassbalance_budget_dimid, submassbalance_tracer_dimid, submassbalance_t_dimid/)
            call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_total', NF90_DOUBLE, &
                dimids3_budget_t, submassbalance_budget_tot_varid))
            call submassbalance_check( nf90_put_att(submassbalance_ncid, submassbalance_budget_tot_varid, &
                "description", "Global budget (would must be constant if conservative variable)") )
            call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_stwat', NF90_DOUBLE, &
                dimids3_budget_t, submassbalance_budget_stwat_varid))
            call submassbalance_check( nf90_put_att(submassbalance_ncid, submassbalance_budget_stwat_varid, &
                "description", "Stock in water") )
#ifdef MUSTANG
            call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_stsed', NF90_DOUBLE, &
                dimids3_budget_t, submassbalance_budget_stsed_varid))
            call submassbalance_check( nf90_put_att(submassbalance_ncid, submassbalance_budget_stsed_varid, &
                "description", "Stock in sediment") )
            call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_flux_ws', NF90_DOUBLE, &
                dimids3_budget_t, submassbalance_budget_flux_ws_varid))
            call submassbalance_check( nf90_put_att(submassbalance_ncid, submassbalance_budget_flux_ws_varid, &
                "description", "FLux at interface water-sediment (>0 if from sed to wat)") )
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
            call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_flux_bedload', NF90_DOUBLE, &
                dimids3_budget_t, submassbalance_budget_flux_bdl_varid))
            call submassbalance_check( nf90_put_att(submassbalance_ncid, submassbalance_budget_flux_bdl_varid, &
                "description", "FLux from bedload transport") )
#endif
#endif
            call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_flux_obc', NF90_DOUBLE, &
                dimids3_budget_t, submassbalance_budget_flux_obc_varid))
            call submassbalance_check( nf90_put_att(submassbalance_ncid, submassbalance_budget_flux_obc_varid, &
                "description", "FLux from zone boundaries (>0 if in)") )
            call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_flux_source', NF90_DOUBLE, &
                dimids3_budget_t, submassbalance_budget_flux_in_varid))
            call submassbalance_check( nf90_put_att(submassbalance_ncid, submassbalance_budget_flux_in_varid, &
                "description", "FLux from rivers") )
        endif      
    endif
    if (nv_fix .gt. 0) then
        call submassbalance_check( nf90_def_dim(submassbalance_ncid, 'tracer_fix',  &
            nv_fix, submassbalance_tracer_fix_dimid) )
        call submassbalance_check( nf90_def_var(submassbalance_ncid, "tracer_fix", NF90_DOUBLE, &
            submassbalance_tracer_fix_dimid, submassbalance_tracer_fix_varid) )
        dimids2_tracer_fix = (/ lchain_dimid, submassbalance_tracer_fix_dimid/)
        call submassbalance_check( nf90_def_var(submassbalance_ncid, 'tracer_fix_name', NF90_CHAR, &
            dimids2_tracer_fix, name_tracer_fix_varid))
        if (submassbalance_nb_close .gt. 0) then
            dimids3_budget_t = (/ submassbalance_budget_dimid, submassbalance_tracer_dimid, submassbalance_t_dimid/)
            call submassbalance_check( nf90_def_var(submassbalance_ncid, 'budget_fix_stwat', NF90_DOUBLE, &
                dimids3_budget_t, submassbalance_budget_fix_varid))
            call submassbalance_check( nf90_put_att(submassbalance_ncid, submassbalance_budget_fix_varid, &
                "description", "Stock in water") )
        endif
    endif

    call submassbalance_check( nf90_enddef(submassbalance_ncid) )
    call submassbalance_check( nf90_sync(submassbalance_ncid))

#ifdef MPI
    else ! not MASTER
        call submassbalance_check( nf90_open(submassbalance_output_file, NF90_WRITE, submassbalance_ncid) )
        call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "xi_rho", x_varid) )
        call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "eta_rho", y_varid) )
# ifdef SPHERICAL
        call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "lon_rho", lonr_varid) )
        call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "lat_rho", latr_varid) )
# else
        call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "x_rho", xr_varid) )
        call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "y_rho", yr_varid) )
# endif
        if (submassbalance_nb_open .gt.0) then
            call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "border", submassbalance_border_varid) )
            call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "border_name", name_border_varid) )
            call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "border_mask_N", mask_border_N_varid) )
            call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "border_mask_E", mask_border_E_varid) )
        endif
        if (submassbalance_nb_close .gt.0) then
            call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "budget", submassbalance_budget_varid) )
            call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "budget_name", name_budget_varid) )
            call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "budget_mask", mask_budget_varid) )
            call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "budget_mask_N", mask_budget_N_varid) )
            call submassbalance_check( nf90_inq_varid(submassbalance_ncid, "budget_mask_E", mask_budget_E_varid) )
        endif
    endif
#endif


    startx = 1
    starty = 1
    imin = 0
    jmin = 0

#ifdef MPI
    if (ii.gt.0) then
    startx=1-imin+iminmpi
    imin=1
    endif
    if (ii.eq.NP_XI-1) then
    imax=Lmmpi+1
    else
    imax=Lmmpi
    endif
    if (jj.gt.0) then
        starty=1-jmin+jminmpi
        jmin=1
    endif
    if (jj.eq.NP_ETA-1) then
        jmax=Mmmpi+1
    else
        jmax=Mmmpi
    endif   
#else
    imax=Lm+1
    jmax=Mm+1
#endif
    countx=imax-imin+1
    county=jmax-jmin+1

    do i=1,Lm+2
        x(i) = i + startx - 1
    enddo
    do j=1,Mm+2
        y(j) = j + starty - 1
    enddo

    call submassbalance_check( nf90_put_var(submassbalance_ncid, x_varid, x(:), &
        start=(/startx/), count=(/countx/)) )
    call submassbalance_check( nf90_put_var(submassbalance_ncid, y_varid, y(:), &
        start=(/starty/), count=(/county/)) )

#ifdef SPHERICAL
    call submassbalance_check( nf90_put_var(submassbalance_ncid, lonr_varid, lonr(imin:imax,jmin:jmax), &
    start=(/startx, starty/), count=(/countx, county/)) )
    call submassbalance_check( nf90_put_var(submassbalance_ncid, latr_varid, latr(imin:imax,jmin:jmax), &
    start=(/startx, starty/), count=(/countx, county/)) )
#else
    call submassbalance_check( nf90_put_var(submassbalance_ncid, xr_varid, xr(imin:imax,jmin:jmax), &
    start=(/startx, starty/), count=(/countx, county/)) )
    call submassbalance_check( nf90_put_var(submassbalance_ncid, yr_varid, yr(imin:imax,jmin:jmax), &
    start=(/startx, starty/), count=(/countx, county/)) )
#endif


    if (submassbalance_nb_open .gt.0) then
        allocate(border_var(submassbalance_nb_open))
        allocate(border_name(submassbalance_nb_open))
        do j=1,submassbalance_nb_open
            border_var(j) = j
            border_name(j) = submassbalance_name(j)
        enddo
        start_border = (/ 1, startx, starty/)
        count_border = (/ submassbalance_nb_open, countx, county /)
        call submassbalance_check( nf90_put_var(submassbalance_ncid, submassbalance_border_varid, border_var(:), &
            start=(/1/), count=(/submassbalance_nb_open/)) )
        call submassbalance_check( nf90_put_var(submassbalance_ncid, name_border_varid, &
            border_name, &
            start=(/1, 1/), count=(/lchain, submassbalance_nb_open/)) )
        call submassbalance_check( nf90_put_var(submassbalance_ncid, mask_border_N_varid, &
            submassbalance_mask_bord_north(1:submassbalance_nb_open,imin:imax,jmin:jmax), &
            start=start_border, count=count_border)  )
        call submassbalance_check( nf90_put_var(submassbalance_ncid, mask_border_E_varid, &
            submassbalance_mask_bord_east(1:submassbalance_nb_open,imin:imax,jmin:jmax), &
            start=start_border, count=count_border)  )
    endif
    if (submassbalance_nb_close .gt.0) then
        allocate(budget_var(submassbalance_nb_close))
        allocate(budget_name(submassbalance_nb_close))
        do j=1,submassbalance_nb_close
            budget_var(j) = j
            budget_name(j) = submassbalance_name(submassbalance_nb_open+j)
        enddo
        start_budget = (/ 1, startx, starty/)
        count_budget = (/ submassbalance_nb_close, countx, county/)
        call submassbalance_check( nf90_put_var(submassbalance_ncid, submassbalance_budget_varid, budget_var(:), &
            start=(/1/), count=(/submassbalance_nb_close/)) )
        call submassbalance_check( nf90_put_var(submassbalance_ncid, name_budget_varid, &
            budget_name, &
            start=(/1, 1/), count=(/lchain, submassbalance_nb_close/)) )
        call submassbalance_check( nf90_put_var(submassbalance_ncid, mask_budget_varid, &
            submassbalance_mask_bud_in(:,imin:imax,jmin:jmax), &
            start=start_budget, count=count_budget) )
        call submassbalance_check( nf90_put_var(submassbalance_ncid, mask_budget_N_varid, &
            submassbalance_mask_bord_north(1+submassbalance_nb_open:submassbalance_nb_border,imin:imax,jmin:jmax), &
            start=start_budget, count=count_budget) )
        call submassbalance_check( nf90_put_var(submassbalance_ncid, mask_budget_E_varid, &
            submassbalance_mask_bord_east(1+submassbalance_nb_open:submassbalance_nb_border,imin:imax,jmin:jmax), &
            start=start_budget, count=count_budget) )
    endif
    call submassbalance_check( nf90_sync(submassbalance_ncid))

#ifdef MPI
    if (mynode.eq.0) then
#endif
    if (nv_adv .gt.0) then
        allocate(tracer_var(nv_adv))
        allocate(tracer_name(nv_adv))
        do j=1,nv_adv
            tracer_var(j) = j
            tracer_name(j) = name_var(irk_fil(j))
        enddo
        call submassbalance_check( nf90_put_var(submassbalance_ncid, submassbalance_tracer_varid, tracer_var(:), &
            start=(/1/), count=(/nv_adv/)) )
        call submassbalance_check( nf90_put_var(submassbalance_ncid, name_tracer_varid, &
            tracer_name, &
            start=(/1, 1/), count=(/lchain, nv_adv/)) )
    endif
    if( nv_fix .gt.0) then
        allocate(tracer_fix_var(nv_fix))
        allocate(tracer_fix_name(nv_fix))
        do j=1,nv_fix
            tracer_fix_var(j) = j
            tracer_fix_name(j) = name_var(irk_fil(j+nv_adv))
        enddo
        call submassbalance_check( nf90_put_var(submassbalance_ncid, submassbalance_tracer_fix_varid, tracer_fix_var(:), &
            start=(/1/), count=(/nv_fix/)) )
        call submassbalance_check( nf90_put_var(submassbalance_ncid, name_tracer_fix_varid, &
            tracer_fix_name, &
            start=(/1, 1/), count=(/lchain, nv_fix/)) )
    endif
    call submassbalance_check( nf90_sync(submassbalance_ncid))

    
#ifdef MPI
    else ! MASTER
        call submassbalance_check( nf90_close(submassbalance_ncid))
    endif
#endif

    submassbalance_next_record = 1


end subroutine submassbalance_def_outnc

!==========================================================================
subroutine submassbalance_check(status)
    !----------------------------------------------------------------------
    !Purpose : check netcdf function
    !Called by : sub_budget_writemask
    !----------------------------------------------------------------------
    USE netcdf
    integer, intent(in) :: status

    if (status /= nf90_noerr) then
    write(stdout, *) 'submassbalance_check(): '
    write(stdout, *) trim( nf90_strerror(status) )
    stop
    endif

end subroutine submassbalance_check

!==========================================================================
subroutine submassbalance_wrt_outnc()
    !----------------------------------------------------------------------
    !Purpose : mass budget writing into netcdf output file
    !Called by : sub_budget_out 
    !----------------------------------------------------------------------
    USE netcdf

    integer :: iv, iz
    integer, dimension(3) :: start

    ! write time
    call submassbalance_check( nf90_put_var(submassbalance_ncid, &
        submassbalance_t_varid, time, (/ submassbalance_next_record /)) )

    ! write close line
    if (submassbalance_nb_close > 0) submassbalance_bil(:,:) = 0.0_rlg
    do iz=1,submassbalance_nb_close
        ! compute total
        submassbalance_bil(iz,1:nv_adv) = submassbalance_stokw_total(iz,1:nv_adv) &
#ifdef MUSTANG
                            + submassbalance_stoks_total(iz,1:nv_adv)     &
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
                            - submassbalance_flxbdlzon_total(iz,1:nv_adv) &
#endif
#endif
                            - submassbalance_flxobczon_total(iz,1:nv_adv) &
                            - submassbalance_flxin_total(iz,1:nv_adv)
        do iv=1,nv_adv
            start = (/ iz, iv, submassbalance_next_record /)
            call submassbalance_check( nf90_put_var(submassbalance_ncid, &
                submassbalance_budget_tot_varid, submassbalance_bil(iz,iv), start) )
            call submassbalance_check( nf90_put_var(submassbalance_ncid, &
                submassbalance_budget_stwat_varid, submassbalance_stokw_total(iz,iv), start) )
#ifdef MUSTANG
            call submassbalance_check( nf90_put_var(submassbalance_ncid, &
                submassbalance_budget_stsed_varid, submassbalance_stoks_total(iz,iv), start) )
            call submassbalance_check( nf90_put_var(submassbalance_ncid, &
                submassbalance_budget_flux_ws_varid, submassbalance_flxws_total(iz,iv), start) )
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
            call submassbalance_check( nf90_put_var(submassbalance_ncid, &
                submassbalance_budget_flux_bdl_varid, submassbalance_flxbdlzon_total(iz,iv), start) )
#endif
#endif
            call submassbalance_check( nf90_put_var(submassbalance_ncid, &
                submassbalance_budget_flux_obc_varid, submassbalance_flxobczon_total(iz,iv), start) )
            call submassbalance_check( nf90_put_var(submassbalance_ncid, &
                submassbalance_budget_flux_in_varid, submassbalance_flxin_total(iz,iv), start) )
        enddo
        do iv=1,nv_fix
            call submassbalance_check( nf90_put_var(submassbalance_ncid, &
                submassbalance_budget_fix_varid, submassbalance_stokwfix_total(iz,iv), start) )
        enddo
    enddo

    ! write open line
    do iz=1,submassbalance_nb_open
        do iv=1,nv_adv
            start = (/ iz, iv, submassbalance_next_record /)
            call submassbalance_check( nf90_put_var(submassbalance_ncid, &
                submassbalance_border_flux_varid, submassbalance_flxobc_total(iz,iv), start) )
#ifdef MUSTANG
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
            call submassbalance_check( nf90_put_var(submassbalance_ncid, &
                submassbalance_border_flux_bdl_varid, submassbalance_flxbdl_total(iz,iv), start) )
#endif
#endif
        enddo
    enddo

    call submassbalance_check( nf90_sync(submassbalance_ncid))

    submassbalance_next_record = submassbalance_next_record + 1

end subroutine submassbalance_wrt_outnc

!==========================================================================
#endif /* ifdef SUBSTANCE */
!!======================================================================
end module submassbalance
