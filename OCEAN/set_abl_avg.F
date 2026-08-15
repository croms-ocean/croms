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
#if (defined ABL1D && defined AVERAGES)

      subroutine set_abl_avg(tile)

      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_abl_avg_tile(Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_abl_avg_tile(Istr,Iend,Jstr,Jend)
!
! Compute time-averaged fields within a tile.
! ------- ------------- ------ ------ - -----
! Because of syncronization issues, the delayed mode averaging
! procedure is used. This procedure implies that all fields to be
! averaged are sampled during the next time step, rather than at
! the end of the time step when they were computed.
!
! Although this algorithm results in somewhat awkward controlling
! logic it has the advantage that all fields to be sampled
! correspond to exactly the same time, which is time step "n".
! Particularly, this is done this way because vertical velocity
! corresponding to the newly computed horizontal velocities
! becomes available only during the following time step.
! The same applies to the density field.
!
! The algorithm consists of three logical blocks: (1) initialization
! of the averages arrays: when mod(ilc-1,navg).eq.1 the target arrays
! are set to the first contribution; (2) accumulation of averaged
! data, when mod(ilc-1,navg).gt.1; and (3) adding the last
! contribution and scaling.
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j, itrc, k, ilc
      real    cff,cff1
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "abl1d.h"
!
# include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntstart   ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsablavg
! is a positive number
!
      if (ilc.gt.ntsablavg) then

        if (mod(ilc-ntsablavg,nwrtablavg).eq.1) then
          cff =1.0
          cff1=0.0
          if (ZEROTH_TILE) then
            timeablavg=time
!            MPI_master_only write(*,*) 'started averaging',
!     &                  iic, ntsdia_avg,nwrtdia_avg
          endif
        elseif (mod(ilc-ntsablavg,nwrtablavg).gt.1) then
          cff =1.0
          cff1=1.0
          if (ZEROTH_TILE) then
           timeablavg=timeablavg+time
          endif
        elseif (mod(ilc-ntsablavg,nwrtablavg).eq.0) then
          cff =1.0/float(nwrtablavg)
          cff1=1.0
          if (ZEROTH_TILE) then
            timeablavg=cff*(timeablavg+time)
!            MPI_master_only write(*,*) 'finished averaging',
!     &                  iic,ntsdia_avg,nwrtdia_avg
          endif
        endif
!
        do k=1,N_abl
          do j=JstrR,JendR
            do i=IstrR,IendR
!
              u_abl_avg(i,j,k) = cff*( cff1*u_abl_avg(i,j,k)
     &                                 +u_abl(i,j,k,nt_n) )
              v_abl_avg(i,j,k) = cff*( cff1*v_abl_avg(i,j,k)
     &                                 +v_abl(i,j,k,nt_n) )
              t_abl_avg(i,j,k) = cff*( cff1*t_abl_avg(i,j,k)
     &                                 +tq_abl(i,j,k,nt_n,jp_ta) )
              q_abl_avg(i,j,k) = cff*( cff1*q_abl_avg(i,j,k)
     &                                 +tq_abl(i,j,k,nt_n,jp_qa) )
!
              avm_abl_avg(i,j,k) = cff*( cff1*avm_abl_avg(i,j,k)
     &                                 +avm_abl(i,j,k,nt_n) )
              avt_abl_avg(i,j,k) = cff*( cff1*avt_abl_avg(i,j,k)
     &                                 +avt_abl(i,j,k,nt_n) )
              mxld_abl_avg(i,j,k) = cff*( cff1*mxld_abl_avg(i,j,k)
     &                                 +mxld_abl(i,j,k,nt_n) )
              mxlm_abl_avg(i,j,k) = cff*( cff1*mxlm_abl_avg(i,j,k)
     &                                 +mxlm_abl(i,j,k,nt_n) )
              tke_abl_avg(i,j,k) = cff*( cff1*tke_abl_avg(i,j,k)
     &                                 +tke_abl(i,j,k,nt_n) )
!
              pu_dta_avg(i,j,k) = cff*( cff1*pu_dta_avg(i,j,k)
     &                                 +pu_dta(i,j,k) )
              pv_dta_avg(i,j,k) = cff*( cff1*pv_dta_avg(i,j,k)
     &                                 +pv_dta(i,j,k) )
              pt_dta_avg(i,j,k) = cff*( cff1*pt_dta_avg(i,j,k)
     &                                 +pt_dta(i,j,k) )
              pq_dta_avg(i,j,k) = cff*( cff1*pq_dta_avg(i,j,k)
     &                                 +pq_dta(i,j,k) )
              pgu_dta_avg(i,j,k) = cff*( cff1*pgu_dta_avg(i,j,k)
     &                                 +pgu_dta(i,j,k) )
              pgv_dta_avg(i,j,k) = cff*( cff1*pgv_dta_avg(i,j,k)
     &                                 +pgv_dta(i,j,k) )
!
            enddo
          enddo
        enddo
!
        do j=JstrR,JendR
          do i=IstrR,IendR
            ablh_avg(i,j) = cff*( cff1*ablh_avg(i,j)
     &                               +ablh(i,j,nt_n) )
          enddo
        enddo
!
      endif
      return
      end
#else
      subroutine set_abl_avg_empty
      return
      end
#endif /* defined ABL1D && defined AVERAGES*/
