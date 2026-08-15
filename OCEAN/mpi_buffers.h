#if defined MPI
# undef NP1
# define NP1 N+1
#ifdef M3FAST_SEDLAYERS
# undef  NP1
# define NP1 N+N_sl+1
#endif 
#define MAXNPTS 3
      integer sub_X,size_X, sub_E,size_E   
# ifndef AGRIF
      parameter (sub_X=Lm,  size_X=MAXNPTS*(sub_X+2*MAXNPTS)-1,
     &           sub_E=Mm,  size_E=MAXNPTS*(sub_E+2*MAXNPTS)-1)

      real ibuf_sndN(0:size_X), ibuf_revN(0:size_X),
     &     ibuf_sndS(0:size_X), ibuf_revS(0:size_X),
     &     jbuf_sndW(0:size_E), jbuf_sndE(0:size_E),
     &     jbuf_revW(0:size_E), jbuf_revE(0:size_E)
# else
      real, dimension(:), allocatable :: ibuf_sndN, ibuf_revN,
     &     ibuf_sndS, ibuf_revS,
     &     jbuf_sndW, jbuf_sndE,
     &     jbuf_revW, jbuf_revE
# endif     

      real buf_snd4(MAXNPTS*MAXNPTS),     buf_snd2(MAXNPTS*MAXNPTS),
     &     buf_rev4(MAXNPTS*MAXNPTS),     buf_rev2(MAXNPTS*MAXNPTS),
     &     buf_snd1(MAXNPTS*MAXNPTS),     buf_snd3(MAXNPTS*MAXNPTS), 
     &     buf_rev1(MAXNPTS*MAXNPTS),     buf_rev3(MAXNPTS*MAXNPTS)

      common/buf_mpi/buf_snd4,buf_snd2,buf_rev4,buf_rev2,
     &               buf_snd1,buf_snd3,buf_rev1,buf_rev3,
     &     ibuf_sndN, ibuf_revN,
     &     ibuf_sndS, ibuf_revS,
     &     jbuf_sndW, jbuf_sndE,
     &     jbuf_revW, jbuf_revE


# if defined SOLVE3D
      integer size_Z,sub_X_3D,size_X_3D, sub_E_3D,size_E_3D

#  ifndef AGRIF
      parameter (size_Z=MAXNPTS*MAXNPTS*(NP1),
     &       sub_X_3D=(Lm+NSUB_X-1)/NSUB_X,
     &     size_X_3D=(NP1)*MAXNPTS*(sub_X_3D+2*MAXNPTS),
     &     sub_E_3D=(Mm+NSUB_E-1)/NSUB_E,
     &     size_E_3D=(NP1)*MAXNPTS*(sub_E_3D+2*MAXNPTS))

      real ibuf_sndN_3D(size_X_3D), ibuf_revN_3D(size_X_3D),
     &     ibuf_sndS_3D(size_X_3D), ibuf_revS_3D(size_X_3D),
     &     jbuf_sndW_3D(size_E_3D), jbuf_sndE_3D(size_E_3D),
     &     jbuf_revW_3D(size_E_3D), jbuf_revE_3D(size_E_3D)
      real buf_snd4_3D(size_Z), buf_snd2_3D(size_Z),
     &     buf_rev4_3D(size_Z), buf_rev2_3D(size_Z),
     &     buf_snd1_3D(size_Z), buf_snd3_3D(size_Z), 
     &     buf_rev1_3D(size_Z), buf_rev3_3D(size_Z)
#  else
      real, dimension(:), allocatable :: ibuf_sndN_3D, ibuf_revN_3D,
     &     ibuf_sndS_3D, ibuf_revS_3D,
     &     jbuf_sndW_3D, jbuf_sndE_3D,
     &     jbuf_revW_3D, jbuf_revE_3D    
      real, dimension(:), allocatable :: buf_snd4_3D, buf_snd2_3D,
     &     buf_rev4_3D, buf_rev2_3D,
     &     buf_snd1_3D, buf_snd3_3D, 
     &     buf_rev1_3D, buf_rev3_3D
#  endif     

      common/buf_mpi_3D/buf_snd4_3D, buf_snd2_3D,
     &                  buf_rev4_3D, buf_rev2_3D,
     &                  buf_snd1_3D, buf_snd3_3D,
     &                  buf_rev1_3D, buf_rev3_3D,
     &     ibuf_sndN_3D, ibuf_revN_3D,
     &     ibuf_sndS_3D, ibuf_revS_3D,
     &     jbuf_sndW_3D, jbuf_sndE_3D,
     &     jbuf_revW_3D, jbuf_revE_3D
# endif
#undef NP1
#endif
