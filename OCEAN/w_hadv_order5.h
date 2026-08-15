!
!===============================================================
!
! Compute wz 5th order horizontal advection
!
!===============================================================
!
!----------------------------------------------------------------------
!  j loop: WFe
!----------------------------------------------------------------------
!
!$acc loop independent
          DO j = max(Jstr,Jmin),min(Jend+1,jmax)  !use full stencil
              DO i = Istr,Iend
                flx5 = Hvom_w(i,j,k)*FLUX5(
     &             wz(i,j-3,k,nrhs), wz(i,j-2,k,nrhs),
     &             wz(i,j-1,k,nrhs), wz(i,j  ,k,nrhs),
     &             wz(i,j+1,k,nrhs), wz(i,j+2,k,nrhs),  Hvom_w(i,j,k) )
#  ifdef MASKING
                flx3 = Hvom_w(i,j,k)*FLUX3(
     &             wz(i,j-2,k,nrhs), wz(i,j-1,k,nrhs),
     &             wz(i,j  ,k,nrhs), wz(i,j+1,k,nrhs),  Hvom_w(i,j,k) )
                flx2 = Hvom_w(i,j,k)*FLUX2(
     &             wz(i,j-1,k,nrhs), wz(i,j,k,nrhs), Hvom_w(i,j,k), cdif)
#   ifdef UP5_MASKING
                mask0=rmask(i,j-1)*rmask(i,j)
                mask2=rmask(i,j-2)*mask0*rmask(i,j+1)
                IF (Hvom_w(i,j,k).gt.0) THEN
                  mask1=rmask(i,j-2)*mask0
                  mask3=rmask(i,j-3)*mask2
                ELSE
                  mask1=rmask(i,j+1)*mask0
                  mask3=rmask(i,j+2)*mask2
                ENDIF
                WFe(i,j)=mask3*flx5+(1-mask3)*mask1*flx3+
     &                             (1-mask3)*(1-mask1)*mask0*flx2
#   else
                mask1=rmask(i,j-2)*rmask(i,j+1)
                mask2=rmask(i,j-3)*rmask(i,j+2)
                mask0=mask1*mask2
                WFe(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                         (1-mask0)*(1-mask1)*flx2
#   endif /* UP5_MASKING */
#  else
                WFe(i,j)=flx5
#  endif /* MASKING */
            ENDDO
          ENDDO

!$acc loop independent
            DO j = Jstr,Jend+1  !j_loop_y_flux_5
              DO i = Istr,Iend
               IF ( j.eq.jmin-2 ) THEN
                 WFe(i,j) = Hvom_w(i,j,k)*FLUX2(
     &             wz(i,j-1,k,nrhs), wz(i,j,k,nrhs), Hvom_w(i,j,k), cdif)
                                                               !
              ELSE IF ( j.eq.jmin-1 .and. jmax.ge.jmin ) THEN  ! 3rd of 4th-order flux 2 in
                                                               ! from south boundary
                flx3 = Hvom_w(i,j,k)*FLUX3(
     &             wz(i,j-2,k,nrhs), wz(i,j-1,k,nrhs),
     &             wz(i,j  ,k,nrhs), wz(i,j+1,k,nrhs),  Hvom_w(i,j,k) )
#  ifdef MASKING
                flx2 = Hvom_w(i,j,k)*FLUX2(
     &             wz(i,j-1,k,nrhs), wz(i,j,k,nrhs), Hvom_w(i,j,k), cdif)
                mask1=rmask(i,j-2)*rmask(i,j+1)
                WFe(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                WFe(i,j)=flx3
#  endif
                                          !
               ELSE IF ( j.eq.jmax+2 ) THEN  ! 2nd order flux next to north
                                          ! boundary
                WFe(i,j) = Hvom_w(i,j,k)*FLUX2(
     &             wz(i,j-1,k,nrhs), wz(i,j,k,nrhs), Hvom_w(i,j,k), cdif)
                                            !
              ELSE IF ( j.eq.jmax+1 ) THEN  ! 3rd or 4th order flux 2 in from
                                            ! north boundary
                flx3 = Hvom_w(i,j,k)*FLUX3(
     &             wz(i,j-2,k,nrhs), wz(i,j-1,k,nrhs),
     &             wz(i,j  ,k,nrhs), wz(i,j+1,k,nrhs),  Hvom_w(i,j,k) )
#  ifdef MASKING
                flx2 = Hvom_w(i,j,k)*FLUX2(
     &             wz(i,j-1,k,nrhs), wz(i,j,k,nrhs), Hvom_w(i,j,k), cdif)
                mask1=rmask(i,j-2)*rmask(i,j+1)
                WFe(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                WFe(i,j)=flx3
#  endif
              ENDIF
            ENDDO
          ENDDO ! j_loop_y_flux_5
!
!----------------------------------------------------------------------
!  i loop: WFx
!----------------------------------------------------------------------
!
          DO j = Jstr,Jend
!$acc loop independent
            DO i = max(Istr,imin),min(Iend+1,imax)  ! i_loop_x_flux_5 use full stencil
              flx5 = Huon_w(i,j,k)*FLUX5(
     &             wz(i-3,j,k,nrhs), wz(i-2,j,k,nrhs),
     &             wz(i-1,j,k,nrhs), wz(i  ,j,k,nrhs),
     &             wz(i+1,j,k,nrhs), wz(i+2,j,k,nrhs),  Huon_w(i,j,k) )
#  ifdef MASKING
              flx3 = Huon_w(i,j,k)*FLUX3(
     &             wz(i-2,j,k,nrhs), wz(i-1,j,k,nrhs),
     &             wz(i  ,j,k,nrhs), wz(i+1,j,k,nrhs),  Huon_w(i,j,k) )
              flx2 = Huon_w(i,j,k)*FLUX2(
     &             wz(i-1,j,k,nrhs), wz(i,j,k,nrhs), Huon_w(i,j,k), cdif)
#   ifdef UP5_MASKING
                mask0=rmask(i-1,j)*rmask(i,j)
                mask2=rmask(i-2,j)*mask0*rmask(i+1,j)
                IF (Huon_w(i,j,k).gt.0) THEN
                  mask1=rmask(i-2,j)*mask0
                  mask3=rmask(i-3,j)*mask2
                ELSE
                  mask1=rmask(i+1,j)*mask0
                  mask3=rmask(i+2,j)*mask2
                ENDIF
                WFx(i,j)=mask3*flx5+(1-mask3)*mask1*flx3+
     &                             (1-mask3)*(1-mask1)*mask0*flx2
#   else
                mask1=rmask(i-2,j)*rmask(i+1,j)
                mask2=rmask(i-3,j)*rmask(i+2,j)
                mask0=mask1*mask2
                WFx(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                         (1-mask0)*(1-mask1)*flx2
#   endif /* UP5_MASKING */
#  else
                WFx(i,j)=flx5
#  endif /* MASKING */
            ENDDO
          ENDDO

          DO j = Jstr,Jend
!$acc loop independent
          DO i = Istr,Iend+1  !i_loop_x_flux_5                                !
            IF ( i.eq.imin-2 ) THEN   ! 2nd order flux next to south
                                           ! boundary
                WFx(i,j) = Huon_w(i,j,k)*FLUX2(
     &             wz(i-1,j,k,nrhs), wz(i,j,k,nrhs), Huon_w(i,j,k), cdif)                                                   !
            ELSE IF ( i.eq.imin-1 .and. imax.ge.imin ) THEN  ! 3rd of 4th order flux 2 in
                                                             ! from south boundary
                flx3 = Huon_w(i,j,k)*FLUX3(
     &             wz(i-2,j,k,nrhs), wz(i-1,j,k,nrhs),
     &             wz(i  ,j,k,nrhs), wz(i+1,j,k,nrhs), Huon_w(i,j,k) )
#  ifdef MASKING
                flx2 = Huon_w(i,j,k)*FLUX2(
     &             wz(i-1,j,k,nrhs), wz(i,j,k,nrhs), Huon_w(i,j,k), cdif)
                mask1=rmask(i-2,j)*rmask(i+1,j)
                WFx(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                WFx(i,j)=flx3
#  endif
                                          !
            ELSE IF ( i.eq.imax+2 ) THEN  ! 2nd order flux next to north
                                          ! boundary
                WFx(i,j) = Huon_w(i,j,k)*FLUX2(
     &             wz(i-1,j,k,nrhs), wz(i,j,k,nrhs), Huon_w(i,j,k), cdif)

                                            !
              ELSE IF ( i.eq.imax+1 ) THEN  ! 3rd or 4th order flux 2 in from
                                            ! north boundary
                flx3 = Huon_w(i,j,k)*FLUX3(
     &             wz(i-2,j,k,nrhs), wz(i-1,j,k,nrhs),
     &             wz(i  ,j,k,nrhs), wz(i+1,j,k,nrhs),  Huon_w(i,j,k) )
#  ifdef MASKING
                flx2 = Huon_w(i,j,k)*FLUX2(
     &             wz(i-1,j,k,nrhs), wz(i,j,k,nrhs), Huon_w(i,j,k), cdif)
                mask1=rmask(i-2,j)*rmask(i+1,j)
                WFx(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                WFx(i,j)=flx3
#  endif
              ENDIF
            ENDDO ! i_loop_x_flux_5
          ENDDO
