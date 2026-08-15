!=======================================================================
!  CROCO  --  TGRU_VMIX
!  Transformer + GRU Vertical Mixing Scheme
!-----------------------------------------------------------------------
!  Implemented by :: S. Maishal, 2026
!  Indian Institute of Technology Kharagpur, India
!=======================================================================
!
#include "cppdefs.h"
!
!-----------------------------------------------------------------------
!  // tgru_vmix_init //
!-----------------------------------------------------------------------
      subroutine tgru_vmix_init()
        use tgru_vmix_mod, only: tgru_vmix_init_impl
        implicit none
        call tgru_vmix_init_impl()
      end subroutine tgru_vmix_init

!-----------------------------------------------------------------------
!  // tgru_vmix(tile) //
!-----------------------------------------------------------------------
      subroutine tgru_vmix(tile)
        use tgru_vmix_mod, only: tgru_vmix_tile_exec
        implicit none
        integer tile
#include "param.h"
#include "scalars.h"
#include "grid.h"
#include "ocean3d.h"
#include "mixing.h"
#include "compute_tile_bounds.h"
!---------------------------------------------------
        call tgru_vmix_tile_exec(
     &        Istr, Iend, Jstr, Jend,
     &        N,
     &        nstp,
     &        itemp,
     &        isalt,
     &        dble(time),
     &        Lm, Mm,
     &        Hz,
     &        u,
     &        v,
     &        rho,
     &        z_r,
     &        z_w,
     &        Akv,
     &        Akt)
!---------------------TODO sm-loop--------------------
#ifdef EXCHANGE_MIXING
        call exchange_w3d_tile(Istr,Iend,Jstr,Jend, Akv)
        call exchange_w3d_tile(Istr,Iend,Jstr,Jend,
     &                          Akt(START_2D_ARRAY,0,itemp))
# ifdef SALINITY
        call exchange_w3d_tile(Istr,Iend,Jstr,Jend,
     &                          Akt(START_2D_ARRAY,0,isalt))
# endif
#endif
!*******************************************************
!-------***(Compatibility issues F77/F90)***------------
!-----------------Currently NEED Wrappe-----------------
!*******************************************************
      end subroutine tgru_vmix
