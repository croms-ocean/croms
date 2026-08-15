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
#include "cppdefs.h"
#ifdef ABL1D
!======================================================================
      subroutine abl_ini (tile)
      implicit none
# include "param.h"
# include "scalars.h"
# include "private_scratch.h"
      integer tile, trd, omp_get_thread_num
# include "abl1d.h"
# include "compute_tile_bounds.h"
! COLD INITIALIZATION FROM LARGE-SCALE DATA
# ifdef ANA_ABL_LSDATA
      call ana_abl_tile (Istr,Iend,Jstr,Jend)
# else
#  ifdef ONLINE
      call set_bulk_tile_online (Istr,Iend,Jstr,Jend)
#  else
      call set_bulk_tile (Istr,Iend,Jstr,Jend)
#  endif
# endif
      call abl_ini_tile(Istr,Iend,Jstr,Jend)
      return
      end
!======================================================================
      subroutine abl_ini_tile(Istr,Iend,Jstr,Jend)
      implicit none
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "abl1d.h"
# if defined ABL_NUDGING && defined ABL_DYN_RESTORE_EQ
#  include "grid.h"
# endif
# include "netcdf.inc"
      integer Istr,Iend,Jstr,Jend, i,j, k,ierr
      integer imin,imax,jmin,jmax
      real    cff,cff1,cff2, eps
      parameter(eps = 1.e-12)
      integer ncid, varid, lvar, lenstr, lstr
# ifdef ANA_ABL_VGRID
      integer jiter, Ntmp,maxiter
      real    theta_s_abl,z1,ds,tol
      real    x,fxp,fx,alpha
      real    sc_r_abl,sc_w_abl
      real    hmax_abl, hc_abl
      logical impose_z1
      parameter( tol = 1.e-12 )
# endif
# if defined ABL_NUDGING && defined ABL_DYN_RESTORE_EQ
      real    fmax
      parameter( fmax = 2.*7.292116e-05*SIN(0.5*pi) )
# endif
# ifdef MPI
      integer err
# endif
!
# include "compute_extended_bounds.h"
      ierr = 0
!
!======================================================================
!++ Generate the vertical grid online
!======================================================================
# ifdef ANA_ABL_VGRID
      hmax_abl     =   2000.
      theta_s_abl  =      2.
      hc_abl       =    400.
      impose_z1    =  .true.
      z1           =     10.
      Ntmp         = N_abl-1
!======================================================================
      if(impose_z1) THEN
        if(z1.LT.10.) THEN
          WRITE(*,*) " ERROR: z1 < 1st ECMWF level height (~10m)"
          stop
        else
        !! Newton iterations to find the appropriate value of theta_s
          maxiter = 1000
          x       = theta_s_abl
          sc_r_abl= (float(1)-0.5)/float(Ntmp)
          alpha   = (z1 - hc_abl*sc_r_abl) / (hmax_abl - hc_abl)
          !
          do jiter=1,maxiter
            fx   = (sinh(sc_r_abl*x)/sinh(x))-alpha
            fxp  = (sc_r_abl*cosh(sc_r_abl*x)-sinh(sc_r_abl*x)
     &                                 *cosh(x)/sinh(x))/sinh(x)
            if( abs(fx) .lt. tol ) THEN
              exit
            endif
            cff  = fx / fxp
            x    = x - cff
          enddo
          theta_s_abl = x
        endif
        ! print*,'theta_s = ',theta_s_abl
      endif
!======================== ==============================================
      ds =1./float(Ntmp)
      cff=(hmax_abl-hc_abl)/sinh(theta_s_abl)
!
      do k = Ntmp,1,-1
        sc_w_abl    = ds*float(k)
        zw_abl(k+1) = hc_abl*sc_w_abl + cff*sinh(theta_s_abl*sc_w_abl)
        sc_r_abl    = ds*(float(k)-0.5)
        zr_abl(k+1) = hc_abl*sc_r_abl + cff*sinh(theta_s_abl*sc_r_abl)
      enddo
!
      zw_abl (1) = 0.
      Hzr_abl(1) = 0.
      zr_abl (1) = 0.
!
      do k = 2,Ntmp+1
        Hzr_abl(k) = zw_abl(k)-zw_abl(k-1)
      enddo
!
      do k=1,Ntmp
        Hzw_abl(k) = zr_abl(k+1)-zr_abl(k)
      enddo
!
      Hzw_abl(Ntmp+1) = zw_abl(Ntmp+1) - zr_abl(Ntmp+1)
      !
      MPI_master_only write(stdout,'(/1x,A/,/1x,A,10x,A/)')
     &                       'Vertical coordinate System:',
     &                       'level   rho-points    w-points',
     &                       'Hzr     Hzw  '
      do k=1,N_abl
        MPI_master_only write(stdout,'(I6,2x,2F12.3,4x,2F7.3)')
     &                   k, zr_abl(k),zw_abl(k), Hzr_abl(k),Hzw_abl(k)
      enddo
# endif
!
!======================================================================
!++ Generate initial conditions
!======================================================================
      nt_n = 1 + MOD( iic  , 2)
      nt_a = 1 + MOD( iic+1, 2)
!======================================================================

!======================================================================
! Extended range  (needed for subgrid scale closures)
!======================================================================
# ifdef EW_PERIODIC
        imin=Istr-2
        imax=Iend+2
# else
      if (WESTERN_EDGE) then
        imin=Istr-1
      else
        imin=Istr-2
      endif
      if (EASTERN_EDGE) then
        imax=Iend+1
      else
        imax=Iend+2
      endif
# endif
# ifdef NS_PERIODIC
        jmin=Jstr-2
        jmax=Jend+2
# else
      if (SOUTHERN_EDGE) then
        jmin=Jstr-1
      else
        jmin=Jstr-2
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend+1
      else
        jmax=Jend+2
      endif
# endif

! Check the presence of u_abl in croco_ini.nc file in case of
! EXACT_RESTART
# ifdef EXACT_RESTART
      lstr=lenstr(ininame)
      ierr=nf_open(ininame(1:lstr), nf_nowrite, ncid)
      lvar=lenstr(vname_abl1d(1,indxabl_u_abl))
      ierr=nf_inq_varid (ncid, vname_abl1d(1,indxabl_u_abl)(1:lvar), varid)
# else
      ierr=nf_noerr+1
# endif
      if (ierr .ne. nf_noerr) then
        do k=1,N_abl
          do j=jmin,jmax
            do i=imin,imax
              u_abl (i,j,k,nt_n      ) = pu_dta(i,j,k)
              v_abl (i,j,k,nt_n      ) = pv_dta(i,j,k)
              tq_abl(i,j,k,nt_n,jp_ta) = pt_dta(i,j,k) + 273.16
              tq_abl(i,j,k,nt_n,jp_qa) = pq_dta(i,j,k)
              u_abl (i,j,k,nt_a      ) = pu_dta(i,j,k)
              v_abl (i,j,k,nt_a      ) = pv_dta(i,j,k)
              tq_abl(i,j,k,nt_a,jp_ta) = pt_dta(i,j,k) + 273.16
              tq_abl(i,j,k,nt_a,jp_qa) = pq_dta(i,j,k)
            enddo
          enddo
        enddo
!
        do j=jstrR,JendR
          do i=istrR,IendR
            ablh(i,j,nt_n)=zw_abl(3)
            ablh(i,j,nt_a)=zw_abl(3)
# if defined ABL_NUDGING && defined ABL_DYN_RESTORE_EQ
!  ABL_DYN_RESTORE_EQ :: nudging for velocities is function of the Coriolis parmeter
            rest_eq(i,j) = SIN( 0.5*pi*( (f(i,j) - fmax)/fmax ) )**8
# endif
          enddo
        enddo
      endif

# ifdef ABL_NUDGING
!!-------------------------------------------------------------------------------------------
!! Compute parameters to build the vertical profile for the nudging term (used in abl_stp())
!!-------------------------------------------------------------------------------------------
#  if defined ABL_NUDGING_TRA || defined ABL_NUDGING_DYN
      pblh_min = zw_abl(      4) / bmin  !<-- at least 3 grid points at the bottom have value rn_ltra_min
      pblh_max = zw_abl(N_abl-3) / bmax  !<-- at least 3 grid points at the top    have value rn_ltra_max
#  endif
#  if defined ABL_NUDGING_TRA
      cff1     = 1. / ( bmax - bmin )**3
! for active tracers
      alp3_tra = -2. * cff1                 * ( ltra_max - ltra_min )
      alp2_tra =  3. * cff1 * (bmax + bmin) * ( ltra_max - ltra_min )
      alp1_tra = -6. * cff1 *  bmax * bmin  * ( ltra_max - ltra_min )
      alp0_tra = cff1 * ( ltra_max * bmin*bmin * (3.*bmax - bmin)
     &                -   ltra_min * bmax*bmax * (3.*bmin - bmax) )
! Check parameters for tracers
      cff  = alp3_tra * bmin**3 + alp2_tra * bmin**2
     &     + alp1_tra * bmin    + alp0_tra
      cff1 = alp3_tra * bmax**3 + alp2_tra * bmax**2
     &     + alp1_tra * bmax    + alp0_tra
      MPI_master_only write(stdout,'(1pe10.3,2x,A)')
     &        cff*dt ,'ABL Minimum value for tracers restoring'
      MPI_master_only write(stdout,'(1pe10.3,2x,A)')
     &        cff1*dt,'ABL Maximum value for tracers restoring'
      IF( cff > cff1  )  THEN
        ierr = ierr + 1
        MPI_master_only write(stdout,'(A)')
     &                    'Invert ltra_min and ltra_max values'
      ENDIF
      IF( cff*dt > 1.+eps .OR. cff*dt < -eps) THEN
        ierr = ierr + 1
        MPI_master_only write(stdout,'(A)')
     &                    'Wrong minimum value for ABL TRA nudging'
      ENDIF
      IF( cff1*dt > 1.+eps .OR. cff1*dt < -eps) THEN
        ierr = ierr + 1
        MPI_master_only write(stdout,'(A)')
     &                    'Wrong maximum value for ABL TRA nudging'
      ENDIF
#  endif
#  if defined ABL_NUDGING_DYN
      cff1     = 1. / ( bmax - bmin )**3
! for dynamics
      alp3_dyn = -2. * cff1                 * ( ldyn_max - ldyn_min )
      alp2_dyn =  3. * cff1 * (bmax + bmin) * ( ldyn_max - ldyn_min )
      alp1_dyn = -6. * cff1 *  bmax * bmin  * ( ldyn_max - ldyn_min )
      alp0_dyn = cff1 * (  ldyn_max * bmin*bmin * (3.*bmax - bmin)
     &                 -   ldyn_min * bmax*bmax * (3.*bmin - bmax) )
! Check parameters for dynamics
      cff  = alp3_dyn * bmin**3 + alp2_dyn * bmin**2
     &     + alp1_dyn * bmin    + alp0_dyn
      cff1 = alp3_dyn * bmax**3 + alp2_dyn * bmax**2
     &     + alp1_dyn * bmax    + alp0_dyn
      MPI_master_only write(stdout,'(1pe10.3,2x,A)')
     &        cff*dt ,'ABL Minimum value for dynamics restoring'
      MPI_master_only write(stdout,'(1pe10.3,2x,A)')
     &        cff1*dt,'ABL Maximum value for dynamics restoring'
      IF( cff > cff1  )  THEN
        ierr = ierr + 1
        MPI_master_only write(stdout,'(A)')
     &                    'Invert ldyn_min and ldyn_max values'
      ENDIF
      IF( cff*dt > 1.+eps .OR. cff*dt < -eps) THEN
        ierr = ierr + 1
        MPI_master_only write(stdout,'(A)')
     &                    'Wrong minimum value for ABL DYN nudging'
      ENDIF
      IF( cff1*dt > 1.+eps .OR. cff1*dt < -eps) THEN
        ierr = ierr + 1
        MPI_master_only write(stdout,'(A)')
     &                    'Wrong maximum value for ABL DYN nudging'
        print*,'cff1*dt = ',cff1*dt
      ENDIF
#  endif
!======================================================================
!      IF(ierr.ne.0) THEN
!#  ifdef MPI
!         call mpi_abort (MPI_COMM_WORLD, err)
!#  else
!         stop                            !-->  EXIT
!#  endif
!      ENDIF
!======================================================================
# endif
      return
      end
!======================================================================

!======================================================================
      subroutine ana_abl_tile(Istr,Iend,Jstr,Jend)
      implicit none
# include "param.h"
# include "scalars.h"
# include "abl1d.h"
# include "forces.h"
      integer Istr,Iend,Jstr,Jend, i,j, k
      integer imin,imax,jmin,jmax
      real    cff
!
# include "compute_extended_bounds.h"
!
!======================================================================
! Extended range  (needed for subgrid scale closures)
!======================================================================
# ifdef EW_PERIODIC
        imin=Istr-2
        imax=Iend+2
# else
      if (WESTERN_EDGE) then
        imin=Istr-1
      else
        imin=Istr-2
      endif
      if (EASTERN_EDGE) then
        imax=Iend+1
      else
        imax=Iend+2
      endif
# endif
# ifdef NS_PERIODIC
        jmin=Jstr-2
        jmax=Jend+2
# else
      if (SOUTHERN_EDGE) then
        jmin=Jstr-1
      else
        jmin=Jstr-2
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend+1
      else
        jmax=Jend+2
      endif
# endif

!======================================================================
!++ Generate large-scale data analytically
!======================================================================
      cff  = 1.e-4 / (grav*itvref)
      do k=1,N_abl
        do j=jmin,jmax
          do i=imin,imax
            pu_dta (i,j,k) = 15.
            pv_dta (i,j,k) = 0.
            pt_dta (i,j,k) = 15.8 + cff*zr_abl(k)  ! in Celsius
            pq_dta (i,j,k) = 0.
            pgu_dta(i,j,k) = 15.
            pgv_dta(i,j,k) = 0.
            radlw  (i,j)   = 0.
            radsw  (i,j)   = 0.
          enddo
        enddo
      enddo
      return
      end
#else
      subroutine abl_ini_empty
      return
      end
#endif /* !ABL1D */
