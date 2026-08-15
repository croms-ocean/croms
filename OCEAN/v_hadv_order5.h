!
!===============================================================
!
! Compute 5th order horizontal advection
!
!===============================================================
!
#  ifdef NS_PERIODIC
          jmin=0
          jmax=LOCALMM
#  else
#   ifdef MPI
          if (SOUTH_INTER) then
            jmin=0
          else
            jmin=3
          endif
          if (NORTH_INTER) then
            jmax=Mmmpi
          else
            jmax=Mmmpi-2
          endif
#   else
          jmin=3
          jmax=Mm-2
#   endif
#  endif
#  ifdef EW_PERIODIC
          imin=1
          imax=LOCALLM+1
#  else
#   ifdef MPI
          if (WEST_INTER) then
            imin=1
          else
            imin=3
          endif
          if (EAST_INTER) then
            imax=Lmmpi+1
          else
            imax=Lmmpi-1
          endif
#   else
          imin=3
          imax=Lm-1
#   endif
#  endif
!
!----------------------------------------------------------------------
!  j loop: VFe
!----------------------------------------------------------------------
!
!$acc loop independent
          DO j = max(JstrV-1,jmin),min(Jend,jmax) ! j_loop_y_flux_5
                                                  ! use full stencil
            DO i = Istr,Iend
              vel = flux6(Hvom(i,j-2,k),Hvom(i,j-1,k),Hvom(i,j  ,k),
     &                    Hvom(i,j+1,k),Hvom(i,j+2,k),Hvom(i,j+3,k),1.)
              flx5 = vel*FLUX5(
     &             v(i,j-2,k,nrhs), v(i,j-1,k,nrhs), 
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs),
     &             v(i,j+2,k,nrhs), v(i,j+3,k,nrhs), vel )
#  ifdef MASKING 
              vel = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
              flx3 = vel*FLUX3(
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs), vel ) 
              flx2 = vel*FLUX2(
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
#   ifdef UP5_MASKING
              mask2=vmask(i,j-1)*vmask(i,j+2)
              IF (vel.gt.0) THEN
                mask1=vmask(i,j-1)
                mask0=vmask(i,j-2)*mask2          
              ELSE
                mask1=vmask(i,j+2)
                mask0=vmask(i,j+3)*mask2
              ENDIF
#   else
              mask1=vmask(i,j-1)*vmask(i,j+2)
              mask2=vmask(i,j-2)*vmask(i,j+3)
              mask0=mask1*mask2
#   endif
              VFe(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                            (1-mask0)*(1-mask1)*flx2
#  else
              VFe(i,j)=flx5
#  endif /* MASKING */
            ENDDO
          ENDDO ! j_loop_y_flux_5

!$acc loop independent
          DO j = JstrV-1,Jend         ! j_loop_y_flux_5
            DO i = Istr,Iend          !
              IF ( j.eq.jmin-2 ) THEN ! 2nd order flux
                                      ! next to south boundary
                vel = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
                VFe(i,j) = vel*FLUX2(
     &             v(i,j,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
                                                              !
              ELSE IF ( j.eq.jmin-1 .and. jmax.ge.jmin ) THEN ! 3rd of 4th order flux
                                                              ! 2 in from south boundary
                vel = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
                flx3 = vel*FLUX3(
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs), vel )
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
                mask1=vmask(i,j-1)*vmask(i,j+2)
                VFe(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                VFe(i,j)=flx3
#  endif
              ELSE IF ( j.eq.jmax+2 ) THEN  ! 2nd order flux
                                            !  next to north boundary
                vel = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
                VFe(i,j) = vel*FLUX2(
     &             v(i,j,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
                                           !
              ELSE IF ( j.eq.jmax+1 ) THEN ! 3rd or 4th order flux
                                           ! 2 in from north boundary
                vel = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
                flx3 = vel*FLUX3(
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs), vel )
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
                mask1=vmask(i,j-1)*vmask(i,j+2)
                VFe(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                VFe(i,j)=flx3
#  endif
            ENDIF
            ENDDO
          ENDDO ! j_loop_y_flux_5
!
!----------------------------------------------------------------------
!  i loop: VFx
!----------------------------------------------------------------------
!
          DO j = JstrV,Jend
!$acc loop independent
            DO i = max(Istr,imin),min(Iend+1,imax) ! i_loop_x_flux_5
                                                   ! use full stencil
              if ( j.ge.jmin+1 .and. j.le.jmax ) then
                vel = flux6(Huon(i,j-3,k),Huon(i,j-2,k),Huon(i,j-1,k),
     &                      Huon(i,j  ,k),Huon(i,j+1,k),Huon(i,j+2,k),1.)
              else
                vel = 0.5*(Huon(i,j-1,k)+Huon(i,j,k))
              endif
              flx5 = vel*FLUX5(
     &             v(i-3,j,k,nrhs), v(i-2,j,k,nrhs),
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs),
     &             v(i+1,j,k,nrhs), v(i+2,j,k,nrhs), vel )
#  ifdef MASKING
              vel = 0.5*(Huon(i,j-1,k)+Huon(i,j,k))
              flx3 = vel*FLUX3(
     &           v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &           v(i  ,j,k,nrhs), v(i+1,j,k,nrhs), vel )
              flx2 = vel*FLUX2(
     &           v(i-1,j,k,nrhs), v(i  ,j,k,nrhs), vel, cdif)
#   ifdef UP5_MASKING
              mask2=vmask(i-2,j)*vmask(i+1,j)
              IF (vel.gt.0) THEN
                mask1=vmask(i-2,j)
                mask0=vmask(i-3,j)*mask2          
              ELSE
                mask1=vmask(i+1,j)
                mask0=vmask(i+2,j)*mask2
              ENDIF
#   else
              mask1=vmask(i-2,j)*vmask(i+1,j)
              mask2=vmask(i-3,j)*vmask(i+2,j)
              mask0=mask1*mask2
#   endif
              VFx(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                            (1-mask0)*(1-mask1)*flx2
#  else
              VFx(i,j)=flx5
#  endif /* MASKING */
            ENDDO
          ENDDO ! i_loop_x_flux_5

          DO j = JstrV,Jend
!$acc loop independent
            DO i = Istr,Iend+1         ! i_loop_x_flux_5
              IF ( i.eq.imin-2 ) THEN  ! 2nd order flux
                                       ! next to south boundary
                vel = 0.5*(Huon(i,j-1,k)+Huon(i,j,k))
                VFx(i,j) = vel*FLUX2(
     &             v(i-1,j,k,nrhs), v(i,j,k,nrhs), vel, cdif)
                                                              !
              ELSE IF ( i.eq.imin-1 .and. imax.ge.imin ) THEN ! 3rd of 4th order flux
                                                              ! 2 in from south boundary
                vel = 0.5*(Huon(i,j-1,k)+ Huon(i,j,k))
                flx3 = vel*FLUX3(
     &             v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &             v(i  ,j,k,nrhs), v(i+1,j,k,nrhs), vel )
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs), vel, cdif)
                mask1=vmask(i-2,j)*vmask(i+1,j)
                VFx(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                VFx(i,j)=flx3
#  endif
              ELSE IF ( i.eq.imax+2 ) THEN ! 2nd order flux
                                           ! next to north boundary
                vel = 0.5*(Huon(i,j-1,k)+ Huon(i,j,k))
                VFx(i,j) = vel*FLUX2(
     &             v(i-1,j,k,nrhs), v(i,j,k,nrhs), vel, cdif)
                                            !
              ELSE IF ( i.eq.imax+1 ) THEN  ! 3rd or 4th order flux
                                            ! 2 in from north boundary
                vel = 0.5*(Huon(i,j-1,k)+ Huon(i,j,k))
                flx3 = vel*FLUX3(
     &             v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &             v(i  ,j,k,nrhs), v(i+1,j,k,nrhs),  vel )
#  ifdef MASKING
                flx2 = vel*FLUX2(
     &             v(i-1,j,k,nrhs), v(i,j,k,nrhs), vel, cdif)
                mask1=vmask(i-2,j)*vmask(i+1,j)
                VFx(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                VFx(i,j)=flx3
#  endif
              ENDIF
            ENDDO
          ENDDO ! i_loop_x_flux_5

