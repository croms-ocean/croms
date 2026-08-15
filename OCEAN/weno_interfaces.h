# if defined BIO_HADV_WENO5 || \
     defined TS_HADV_WENO5  || defined TS_VADV_WENO5 || \
     defined UV_HADV_WENO5  || defined UV_VADV_WENO5 || \
     defined W_HADV_WENO5   || defined W_VADV_WENO5 || \
     defined BEDLOAD_WENO5
      INTERFACE
      ! flux3_weno
      function flux3_weno( q_im2, q_im1, q_i, q_ip1, ua)
!$acc routine(flux3_weno) seq
      REAL    :: flux3_weno
      REAL    :: q_im2, q_im1, q_i, q_ip1, ua
      end function flux3_weno 

      ! flux5_weno
      function flux5_weno(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
!$acc routine(flux5_weno) seq
      implicit none
      REAL :: flux5_weno
      REAL, INTENT(IN) :: q_im3, q_im2, q_im1,
     &                           q_i, q_ip1, q_ip2, ua
      end function flux5_weno
      END INTERFACE
#endif
