#include "cppdefs.h"
#if defined MPI && defined SOLVE3D

# ifndef MP_3PTS
#  ifndef MP_M3FAST_SEDLAYERS
      subroutine MessPass3D_tile (Istr,Iend,Jstr,Jend, A, nmax)
# else
       subroutine MessPass3D_sl_tile (Istr,Iend,Jstr,Jend, A, 
     &  nmax)
#  endif
# else
#  ifndef MP_M3FAST_SEDLAYERS
      subroutine MessPass3D_3pts_tile (Istr,Iend,Jstr,Jend, A, nmax)
#  else
       subroutine MessPass3D_3pts_sl_tile (Istr,Iend,Jstr,Jend, A,
     &  nmax)
#  endif
# endif
!
!======================================================================
!
! This subroutine is designed for CROCO-MPI code. It exchanges domain
! boundary information, including 2 (or 3) ghost-cells in each
! direction.
!
! Ping Wang 9/15/99.
! Patrick Marchesiello 2012: generic number of ghost-cells Npts
!
!======================================================================
!
!     implicit none
# include "param.h"
# include "scalars.h"
# include "mpi_cpl.h"
      include 'mpif.h'
      logical, parameter :: use_host_device = .TRUE.
# include "mpi_buffers.h"
!
! Nb of boundary points involved in communication
!
      integer Npts,ipts,jpts
# ifndef MP_3PTS
      parameter (Npts=2)
# else
      parameter (Npts=3)
# endif
# undef  NP1      
# define NP1 N+1
#ifdef MP_M3FAST_SEDLAYERS
# undef  NP1
# define NP1 N+N_sl+1
#endif 
      integer nmax
# ifdef  BAND_DEBUG         
      real,dimension(Npts)        :: checksum
# endif
      real A(GLOBAL_2D_ARRAY,nmax)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN) BLOCK_CLAUSE
      integer Istr,Iend,Jstr,Jend, i,j,k, isize,jsize,ksize,
     &        req(8), status(MPI_STATUS_SIZE,8), ierr
      integer iter, mdii, mdjj
!
# include "compute_message_bounds.h"

# ifdef AGRIF
      size_Z=MAXNPTS*MAXNPTS*(NP1)
      sub_X=(Lm+NSUB_X-1)/NSUB_X
      size_X=(NP1)*MAXNPTS*(sub_X+2*MAXNPTS)
      sub_E=(Mm+NSUB_E-1)/NSUB_E
      size_E=(NP1)*MAXNPTS*(sub_E+2*MAXNPTS)

      if(.not. allocated(buf_snd4_3D)) then
      Allocate(buf_snd4_3D(size_Z),  ibuf_sndN_3D(size_X),  
     &         buf_snd2_3D(size_Z),
     &         buf_rev4_3D(size_Z),  ibuf_revN_3D(size_X),  
     &         buf_rev2_3D(size_Z),
     &        jbuf_sndW_3D(size_E),                     jbuf_sndE_3D(size_E),
     &        jbuf_revW_3D(size_E),                     jbuf_revE_3D(size_E),
     &         buf_snd1_3D(size_Z),  ibuf_sndS_3D(size_X),  
     &         buf_snd3_3D(size_Z),
     &         buf_rev1_3D(size_Z),  ibuf_revS_3D(size_X),  
     &         buf_rev3_3D(size_Z))
      endif
# endif

      ksize=Npts*Npts*nmax               ! message sizes for
      isize=Npts*ishft*nmax              ! corner messages and sides
      jsize=Npts*jshft*nmax              ! in XI and ETA directions
!
! Prepare to receive and send: sides....
!
                            !  Message passing split into two stages
                            !  in order to optimize Send-Recv pairing
                            !  in such a way that if one subdomain
      do iter=0,1           !  sends message to, say, its WESTERN
        mdii=mod(ii+iter,2) !  neighbor, that neighbor is preparing
        mdjj=mod(jj+iter,2) !  to receive this message first (i.e.
                            !  message coming from its EASTERN side),
                            !  rather than send his WEST
                            !  bound message, similarly to the first
                            !  subdomain.

!
! Prepare to receive and send: sides....
        if (mdii.eq.0) then
          if ((WEST_INTER2)) then
!$acc parallel loop if(compute_on_device) collapse(3)
            do k=1,nmax
              do j=jmin,jmax
                do ipts=1,Npts
                  jbuf_sndW_3D(k+nmax*(j-jmin+(ipts-1)*jshft))=
     &                  A(ipts,j,k)
                enddo
              enddo
            enddo
!$acc host_data if(compute_on_device)
!$acc&   use_device(jbuf_revW_3D,jbuf_sndW_3D)
            call MPI_Irecv (jbuf_revW_3D, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 2, MPI_COMM_WORLD, req(1), ierr)
            call MPI_Send  (jbuf_sndW_3D, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 1, MPI_COMM_WORLD,         ierr)
!$acc end host_data
          endif
        else
          if (EAST_INTER2) then
!$acc parallel loop if(compute_on_device) collapse(3)
            do k=1,nmax
              do j=jmin,jmax
                do ipts=1,Npts
                  jbuf_sndE_3D(k+nmax*(j-jmin+(ipts-1)*jshft))=
     &                                       A(Lmmpi-Npts+ipts,j,k)
                enddo
              enddo
            enddo
!$acc host_data if(compute_on_device)
!$acc&   use_device(jbuf_revE_3D,jbuf_sndE_3D)
            call MPI_Irecv (jbuf_revE_3D, jsize, MPI_DOUBLE_PRECISION,
     &                        p_E, 1, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Send  (jbuf_sndE_3D, jsize, MPI_DOUBLE_PRECISION,
     &                        p_E, 2, MPI_COMM_WORLD,         ierr)
!$acc end host_data
          endif
        endif

        if (mdjj.eq.0) then
          if (SOUTH_INTER2) then
!$acc parallel loop if(compute_on_device) 
!$acc& default(present) collapse(3)
            do k=1,nmax
              do jpts=1,Npts
                do i=imin,imax
                  ibuf_sndS_3D(k+nmax*(i-imin+(jpts-1)*ishft))=
     &                          A(i,jpts,k)
                enddo
              enddo
            enddo
!$acc update host(ibuf_sndS_3D(1:isize)) if(.not.use_host_device)
!$acc host_data if(compute_on_device.and.use_host_device)
!$acc&   use_device(ibuf_revS_3D,ibuf_sndS_3D)
            call MPI_Irecv (ibuf_revS_3D, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 4, MPI_COMM_WORLD, req(3), ierr)
            call MPI_Send  (ibuf_sndS_3D, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 3, MPI_COMM_WORLD,         ierr)
!$acc end host_data
          endif
        else
          if (NORTH_INTER2) then
!$acc parallel loop if(compute_on_device)
!$acc& default(present) collapse(3)
            do k=1,nmax
              do jpts=1,Npts
                do i=imin,imax
                  ibuf_sndN_3D(k+nmax*(i-imin+(jpts-1)*ishft))=
     &                                         A(i,Mmmpi-Npts+jpts,k)
                enddo
              enddo
            enddo
!$acc update host(ibuf_sndN_3D(1:isize)) if(.not.use_host_device)
!$acc host_data if(compute_on_device.and.use_host_device)
!$acc&   use_device(ibuf_revN_3D,ibuf_sndN_3D)
            call MPI_Irecv (ibuf_revN_3D, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 3, MPI_COMM_WORLD, req(4), ierr)
            call MPI_Send  (ibuf_sndN_3D, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 4, MPI_COMM_WORLD,         ierr)
!$acc end host_data
          endif
        endif
!
! ...corners:
!
        if (mdii.eq.0) then
          if (CORNER_SW) then
!$acc parallel loop if(compute_on_device) collapse(3)
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd1_3D(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                 A(ipts,jpts,k)
                enddo
              enddo
            enddo
!$acc host_data if(compute_on_device)
!$acc&   use_device(buf_rev1_3D,buf_snd1_3D)
            call MPI_Irecv (buf_rev1_3D, ksize, MPI_DOUBLE_PRECISION, p_SW,
     &                                6, MPI_COMM_WORLD, req(5),ierr)
            call MPI_Send  (buf_snd1_3D, ksize, MPI_DOUBLE_PRECISION, p_SW,
     &                                5, MPI_COMM_WORLD,        ierr)
!$acc end host_data
          endif
        else
          if (CORNER_NE) then
!$acc parallel loop if(compute_on_device) collapse(3)
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd2_3D(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                            A(Lmmpi+ipts-Npts,Mmmpi+jpts-Npts,k)
                enddo
              enddo
            enddo
!$acc host_data if(compute_on_device)
!$acc&   use_device(buf_rev2_3D,buf_snd2_3D)
            call MPI_Irecv (buf_rev2_3D, ksize, MPI_DOUBLE_PRECISION, 
     &                  p_NE, 5, MPI_COMM_WORLD, req(6),ierr)
            call MPI_Send  (buf_snd2_3D, ksize, MPI_DOUBLE_PRECISION, 
     &                  p_NE, 6, MPI_COMM_WORLD,        ierr)
!$acc end host_data
          endif
        endif

        if (mdii.eq.1) then
          if (CORNER_SE) then
!$acc parallel loop if(compute_on_device) collapse(3)
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd3_3D(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                                A(Lmmpi+ipts-Npts,jpts,k)
                enddo
              enddo
            enddo
!$acc host_data if(compute_on_device)
!$acc&   use_device(buf_rev3_3D,buf_snd3_3D)
            call MPI_Irecv (buf_rev3_3D, ksize, MPI_DOUBLE_PRECISION, 
     &                  p_SE, 8, MPI_COMM_WORLD, req(7),ierr)
            call MPI_Send  (buf_snd3_3D, ksize, MPI_DOUBLE_PRECISION, 
     &                  p_SE, 7, MPI_COMM_WORLD,        ierr)
!$acc end host_data
          endif
        else
          if (CORNER_NW) then
!$acc parallel loop if(compute_on_device) collapse(3)
            do k=1,nmax
              do jpts=1,Npts
                do ipts=1,Npts
                  buf_snd4_3D(k+nmax*(ipts-1+Npts*(jpts-1)))=
     &                                A(ipts,Mmmpi+jpts-Npts,k)
                enddo
              enddo
            enddo
!$acc host_data if(compute_on_device)
!$acc&   use_device(buf_rev4_3D,buf_snd4_3D)
            call MPI_Irecv (buf_rev4_3D, ksize, MPI_DOUBLE_PRECISION, 
     &                  p_NW, 7, MPI_COMM_WORLD, req(8),ierr)
            call MPI_Send  (buf_snd4_3D, ksize, MPI_DOUBLE_PRECISION, 
     &                  p_NW, 8, MPI_COMM_WORLD,        ierr)
!$acc end host_data
          endif
        endif
      enddo   !<-- iter

!
! Wait for completion of receive and fill ghost points: sides...
!
      if (WEST_INTER2) then
        call MPI_Wait (req(1),status(1,1),ierr)
!$acc parallel loop if(compute_on_device) collapse(3) 
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
# ifdef  BAND_DEBUG         
c            if ( chkbandname .ne. 'none' ) then
            if ( A(ipts-Npts,j,k) == 
     &           jbuf_revW_3D(k+nmax*(j-jmin+(ipts-1)*jshft)) ) 
     &           checksum(ipts) = checksum(ipts)+1
c            endif
# endif
             A(ipts-Npts,j,k)=
     &          jbuf_revW_3D(k+nmax*(j-jmin+(ipts-1)*jshft))
            enddo
          enddo
        enddo
# ifdef  BAND_DEBUG         
c         if ( chkbandname .ne. 'none') then
          do ipts=1,Npts
c             print *,"checksum(ipts)=",checksum(ipts)
c             print *,"jnin,jmax,nmax=",jmin,jmax,nmax,":",ipts
            if( nint(checksum(ipts)) == nmax*(jmax-jmin+1) ) then
               print *,mynode,"ipts=",int(ipts,1),"WEST_INTER2 "//
     &          trim(chkbandname)//':ghostcell  3d'
c                stop
            endif
          enddo
          checksum(:) = 0.
c          endif
# endif        
      endif
      if (WEST_INTER .and. .not. WEST_INTER2) then
# ifdef OPENACC
        print *,'NOT READY FOR OPENACC'
        stop
# endif
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
              A(ipts-Npts,j,k)=A(ipts,j,k)
            enddo
          enddo
        enddo
      endif

      if (EAST_INTER2) then
        call MPI_Wait (req(2),status(1,2),ierr)
# ifdef  BAND_DEBUG         
        checksum(:) = 0.
# endif        
!$acc parallel loop if(compute_on_device) collapse(3) 
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
# ifdef  BAND_DEBUG         
c            if ( chkbandname .ne. 'none' ) then
            if ( A(Lmmpi+ipts,j,k) ==
     &                         jbuf_revE_3D(k+nmax*(j-jmin+(ipts-1)*jshft)) ) 
     &           checksum(ipts) = checksum(ipts)+1
c            endif
# endif
              A(Lmmpi+ipts,j,k)=
     &                         jbuf_revE_3D(k+nmax*(j-jmin+(ipts-1)
     &                         *jshft))
            enddo
          enddo
        enddo
# ifdef  BAND_DEBUG         
c         if ( chkbandname .ne. 'none') then
          do ipts=1,Npts
c             print *,"checksum(ipts)=",checksum(ipts)
c             print *,"jnin,jmax=",jmin,jmax,":",ipts
            if( nint(checksum(ipts)) == nmax*(jmax-jmin+1) ) then
               print *,mynode,"ipts=",int(ipts,1),"EAST_INTER2 "//
     &          trim(chkbandname)//':ghostcell  3d'
c                stop
            endif
          enddo
          checksum(:) = 0.
c          endif
# endif        
      endif
      if (EAST_INTER .and. .not. EAST_INTER2) then
!$acc parallel loop if(compute_on_device) collapse(3) 
        do k=1,nmax
          do j=jmin,jmax
            do ipts=1,Npts
              A(Lmmpi+ipts,j,k)=A(Lmmpi+ipts-Npts,j,k)
            enddo
          enddo
        enddo
      endif


      if (SOUTH_INTER2) then
        call MPI_Wait (req(3),status(1,3),ierr)
!$acc update device(ibuf_revS_3D(1:isize)) if(.not.use_host_device)
!$acc parallel loop if(compute_on_device) collapse(3) 
# ifdef  BAND_DEBUG         
        checksum(:) = 0.
# endif        
        do k=1,nmax
          do jpts=1,Npts
            do i=imin,imax
# ifdef  BAND_DEBUG         
c            if ( chkbandname .ne. 'none' ) then
            if (  A(i,jpts-Npts,k) == 
     &            ibuf_revS_3D(k+nmax*(i-imin+(jpts-1)*ishft)) )
     &          checksum(jpts) = checksum(jpts)+1
c            endif
# endif
              A(i,jpts-Npts,k)=ibuf_revS_3D(k+nmax*(i-imin+(jpts-1)
     &                        *ishft))
            enddo
          enddo
        enddo
      endif
      if (SOUTH_INTER .and. .not. SOUTH_INTER2) then
!$acc parallel loop if(compute_on_device) collapse(3) 
       do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,jpts-Npts,k)=A(i,jpts,k)
            enddo
          enddo
        enddo
# ifdef  BAND_DEBUG         
c         if ( chkbandname .ne. 'none') then
          do jpts=1,Npts
c             print *,"checksum(ipts)=",checksum(ipts)
c             print *,"jnin,jmax=",jmin,jmax,":",ipts
            if( nint(checksum(jpts)) == nmax*(imax-imin+1) ) then
               print *,mynode,"jpts=",int(jpts,1),"SOUTH_INTER2 "//
     &          trim(chkbandname)//':ghostcell  3d'
c                stop
            endif
          enddo
          checksum(:) = 0.
c          endif
# endif        
      endif

      if (NORTH_INTER2) then
        call MPI_Wait (req(4),status(1,4),ierr)
!$acc update device(ibuf_revN_3D(1:isize)) if(.not.use_host_device)
!$acc parallel loop if(compute_on_device) collapse(3)
# ifdef  BAND_DEBUG         
        checksum(:) = 0.
# endif        
        do k=1,nmax
          do jpts=1,Npts
            do i=imin,imax
# ifdef  BAND_DEBUG         
c            if ( chkbandname .ne. 'none' ) then
            if ( A(i,Mmmpi+jpts,k)==
     &                         ibuf_revN_3D(k+nmax*(i-imin+(jpts-1)*ishft)) ) 
     &          checksum(jpts) = checksum(jpts)+1
c            endif
# endif
              A(i,Mmmpi+jpts,k)=
     &          ibuf_revN_3D(k+nmax*(i-imin+(jpts-1)*ishft))
            enddo
          enddo
        enddo
# ifdef  BAND_DEBUG         
c         if ( chkbandname .ne. 'none') then
          do jpts=1,Npts
c             print *,"checksum(ipts)=",checksum(ipts)
c             print *,"jnin,jmax=",jmin,jmax,":",ipts
            if( nint(checksum(jpts)) == nmax*(imax-imin+1) ) then
               print *,mynode,"jpts=",int(jpts,1),"NORTH_INTER2 "//
     &          trim(chkbandname)//':ghostcell  3d'
c                stop
            endif
          enddo
          checksum(:) = 0.
c          endif
# endif        
      endif
      if (NORTH_INTER .and. .not. NORTH_INTER2) then
!$acc parallel loop if(compute_on_device) collapse(3) 
        do k=1,nmax
          do i=imin,imax
            do jpts=1,Npts
              A(i,Mmmpi+jpts,k)=A(i,Mmmpi+jpts-Npts,k)
            enddo
          enddo
        enddo
      endif
!
! ...corners:
!
      if (CORNER_SW) then
        call MPI_Wait (req(5),status(1,5),ierr)
!$acc parallel loop if(compute_on_device) collapse(3) 
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,jpts-Npts,k)=
     &                           buf_rev1_3D(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_SW  .and.
     &   SOUTH_INTER .and.  WEST_INTER ) then
# ifdef OPENACC
        print *,'NOT READY FOR OPENACC'
        stop
# endif
!$acc parallel loop if(compute_on_device) collapse(3) 
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,jpts-Npts,k)=
     &                           A(ipts,jpts,k)
            enddo
          enddo
        enddo
       endif

      if (CORNER_NE) then
        call MPI_Wait (req(6),status(1,6),ierr)
!$acc parallel loop if(compute_on_device) collapse(3) 
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,Mmmpi+jpts,k)=
     &                           buf_rev2_3D(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_NE  .and.
     &   NORTH_INTER .and.  EAST_INTER ) then
!$acc parallel loop if(compute_on_device) collapse(3) 
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,Mmmpi+jpts,k)=
     &              A(Lmmpi+ipts-Npts,Mmmpi+jpts-Npts,k)
            enddo
          enddo
        enddo
       endif


      if (CORNER_SE) then
        call MPI_Wait (req(7),status(1,7),ierr)
!$acc parallel loop if(compute_on_device) collapse(3) 
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(Lmmpi+ipts,jpts-Npts,k)=
     &                           buf_rev3_3D(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_SE .and.
     &   SOUTH_INTER .and.  EAST_INTER ) then
!$acc parallel loop if(compute_on_device) collapse(3)
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
               A(Lmmpi+ipts,jpts-Npts,k)=
     &                           A(Lmmpi+ipts-Npts,jpts,k)
            enddo
          enddo
        enddo
       endif


      if (CORNER_NW) then
        call MPI_Wait (req(8),status(1,8),ierr)
!$acc parallel loop if(compute_on_device) collapse(3) 
        do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,Mmmpi+jpts,k)=
     &                           buf_rev4_3D(k+nmax*(ipts-1+Npts*(jpts-1)))
            enddo
          enddo
        enddo
      endif
      if (.not. CORNER_NW  .and.
     &   NORTH_INTER .and.  WEST_INTER ) then
!$acc parallel loop if(compute_on_device) collapse(3) 
       do k=1,nmax
          do jpts=1,Npts
            do ipts=1,Npts
              A(ipts-Npts,Mmmpi+jpts,k)=A(ipts,Mmmpi+jpts-Npts,k)
            enddo
          enddo
        enddo
       endif


      return
      end

# ifndef MP_3PTS
#  define MP_3PTS
#  include "MessPass3D.F"
#  undef MP_3PTS
# endif

# ifndef MP_3PTS
# ifdef M3FAST_SEDLAYERS      
# ifndef MP_M3FAST_SEDLAYERS
#  define MP_M3FAST_SEDLAYERS
#  include "MessPass3D.F"
#  undef  MP_M3FAST_SEDLAYERS
# endif
# endif      
# endif

# undef  NP1

#else
      subroutine MessPass3D_empty
      return
      end
#endif

