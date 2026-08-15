!
!====================================================================
!           Wave maker for wave-resolving simulations:
!              monochromatic, bichromatic or JONSWAP
!====================================================================
!
#ifdef WAVE_MAKER_EAST
# define ZBRY zetabry_east
# define UBARBRY ubarbry_east
# define VBARBRY vbarbry_east
# define UBRY ubry_east
# define VBRY vbry_east
# define WBRY wbry_east
# define TBRY tbry_east
# define UNBQBRY unbqbry_east
# define VNBQBRY vnbqbry_east
# define WNBQBRY wnbqbry_east
# define IB0 Lm+1
# define IB1 Lm
#else
# define ZBRY zetabry_west
# define UBARBRY ubarbry_west
# define VBARBRY vbarbry_west
# define UBRY ubry_west
# define VBRY vbry_west
# define WBRY wbry_west
# define TBRY tbry_west
# define UNBQBRY unbqbry_west
# define VNBQBRY vnbqbry_west
# define WNBQBRY wnbqbry_west
# define IB0 0
# define IB1 1
#endif
!
!--------------------------------------------------------------------
!  Configurations
!--------------------------------------------------------------------
!
!  Set configuration parameters or get from croco.in file
!
#ifdef ROGUE_WAVE
# define WAVE_MAKER_SPECTRUM

#elif defined RIP && !defined MRL_WCI
# ifdef WAVE_MAKER_SPECTRUM
#  define WAVE_MAKER_JONSWAP
#  undef  WAVE_MAKER_GAUSSIAN
# else
#  undef  WAVE_MAKER_BICHROMATIC
#  undef  STOKES_WAVES
# endif
        wp=11.            ! period
        wa=0.4            ! amplitude
        wd=-10.           ! incidence angle (deg)
        wds=30.           ! directional spread (deg)
                          !  -> crest length = wl/(2*sin(wds))
#elif defined SWASH
# ifdef SWASH_GLOBEX_B2
#  define WAVE_MAKER_BICHROMATIC
        wf1=2*pi*0.42     ! GLOBEX B2
        wf2=2*pi*0.462
        wa1=0.09
        wa2=0.01
# elif defined SWASH_GLOBEX_B3
#  define WAVE_MAKER_BICHROMATIC
        wf1=2*pi*0.42     ! GLOBEX B3
        wf2=2*pi*0.462
        wa1=0.07
        wa2=0.03
# elif defined SWASH_GLOBEX_A3
#  define WAVE_MAKER_JONSWAP
        wp=2.25           ! period
        wa=0.0354         ! amplitude
        gamma=20.         ! JONSWAP peakedness parameter
# else
        wp=2.
        wa=0.0442
# endif
        wd=0.
        wds=0.
#elif defined SANDBAR && !defined MRL_WCI
# define WAVE_MAKER_JONSWAP
# ifdef SANDBAR_OFFSHORE
        wp=5.             ! period
        wa=0.49           ! amplitude
        gamma=3.3         ! JONSWAP peakedness parameter
# else
        wp=8.             ! period
        wa=0.21           ! amplitude
        gamma=3.3         ! JONSWAP peakedness parameter
# endif
        wd=0.
        wds=0.
#elif defined DUCK3D
# define WAVE_MAKER_JONSWAP
        wp=14.            ! period
        wa=0.5            ! amplitude
        gamma=3.3         ! JONSWAP peakedness parameter
        wd=-10.           ! incidence angle (deg)
        wds=30.           ! directional spread (deg)
                          !  -> crest length = wl/(2*sin(wds))
#else
!
!  get parameters from croco.in
!
        wa=wmaker_amp     ! amplitude
        wp=wmaker_prd     ! period
        wd=wmaker_dir     ! incidence angle
        wds=wmaker_dsp    ! directional spread
        gamma=wmaker_fsp  ! JONSWAP peakedness parameter
#endif
!
        wf=2*pi/wp        ! frequency
!
!  Time & space origins
!
#ifdef ROGUE_WAVE
        x0=14.1
        y0=0.
        time0=64.
#else
        x0=xr(IB0,0)
        y0=0.
        time0=0.
#endif
!
!  Convert angles to rad
!
#ifdef WAVE_MAKER_EAST
        wd=wd+180.0
#endif
        wd =wd *deg2rad  ! incidence angle
        wds=wds*deg2rad  ! directional spread
!
        coswd=cos(wd)
        sinwd=sin(wd)
        coswds=cos(wds)
        sinwds=sin(wds)
!
!  Peak frequency and wavenumber
!
        h0=h(IB0,1)
        wf=2*pi/wp          ! peak frequency
        khd=h0*wf*wf/g      ! peak wavenumber
        wk=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                   khd*(K5+K6*khd)))))) )/h0
!
#if defined WAVE_MAKER_JONSWAP || defined WAVE_MAKER_GAUSSIAN
# define WAVE_MAKER_SPECTRUM
#endif
!
!--------------------------------------------------------------------
!  Initialisation
!--------------------------------------------------------------------
!
#ifdef ROGUE_WAVE
!
!  Read file
!
        if (FIRST_TIME_STEP) then
!         open(117,file='datwaves_CORR1.txt',form='formatted',status='old')
          open(117,file='datwaves.txt',form='formatted',status='old')
          do k=1,Nfrq !--> Nfrq=320 in forces.h
            read(117,*) wa_bry(k), wf_bry(k), wpha_bry(k), wk_bry(k)
            wa_bry(k)=wa_bry(k)*0.154/0.05  ! correct amplitude
            ! wpha_bry(k)=wpha_bry(k) + 1.5*pi
!            khd=h0*wf_bry(k)**2/g      ! recompute wavenumber
!            kh=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
!     &                                       khd*(K5+K6*khd)))))) )
!            wk_bry(k)=kh/h0
          enddo
        endif
        ramp=tanh(dt/2.*float(iic-ntstart))

#elif defined WAVE_MAKER_SPECTRUM
!
!  Build wave spectrum
!
        if (FIRST_TIME_STEP) then

          fmin=0.2*wf  ! frequency spread
          fmax=5.0*wf
          df=(fmax-fmin)/Nfrq
          cff2=0.
          do iw=1,Nfrq
            wf_bry(iw)=fmin+float(iw)*df
            khd=h0*wf_bry(iw)**2/g   ! wavenumber
            kh=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                       khd*(K5+K6*khd)))))) )
            wk_bry(iw)=kh/h0
          enddo
# ifdef WAVE_MAKER_JONSWAP
          sumspec=0.
          do iw=1,Nfrq
            sigma=0.5*( 0.09*(1.+sign(1.,wf_bry(iw)-wf))+
     &                  0.07*(1.-sign(1.,wf_bry(iw)-wf)) )
            cff0=exp(-0.5*((wf_bry(iw)-wf)/(sigma*wf))**2)
            cff1=0.3119*(wf**4)/(wf_bry(iw)**5)
     &                 *exp(-1.25*(wf/wf_bry(iw))**4)*gamma**cff0
            cff2=16.*cff1*df  ! integral must be 1
            wa_bry(iw)=cff2
            sumspec=sumspec+cff2
          enddo
# elif defined WAVE_MAKER_GAUSSIAN
          cff2=0.
          do iw=1,Nfrq
            cff1=exp(-((wf_bry(iw)-wf)/0.05)**2)
            wa_bry(iw)=cff1
            sumspec=sumspec+cff1
          enddo
# endif
# ifdef WAVE_MAKER_DSPREAD
!
! Single-sum wavemaker (Treillou et al. 2024)
!
          jwtt0=minloc(abs(wf_bry-wf),DIM=1)
          jwtt=MOD(jwtt0,Ndir)
          cff1=0.
          do jw=1,Nfrq
            wd_bry(jw)=MOD(float(jw-jwtt),float(Ndir))
            if (wd_bry(jw) .le. 0.) wd_bry(jw)=wd_bry(jw)+float(Ndir)
            wd_bry(jw)=(-1.)**jw * 0.5*pi*( -1. + 
     &                 (floor(wd_bry(jw)-1.))/(float(Ndir)-1.) )
            wd_bry(jw)=wd_bry(jw)+wd
#   ifdef WAVE_MAKER_EAST
            if (wd_bry(jw) .ge.  1.5*pi) wd_bry(jw)= 1.5*pi  ! Wave dir. constrained
            if (wd_bry(jw) .le.  0.5*pi) wd_bry(jw)= 0.5*pi  ! between [-90,90] deg of mean
#   else
            if (wd_bry(jw) .ge.  0.5*pi) wd_bry(jw)= 0.5*pi
            if (wd_bry(jw) .le. -0.5*pi) wd_bry(jw)=-0.5*pi
#   endif
            wa_bry_d(jw)=exp(-((wd_bry(jw)-wd)/
     &                          max(1.5*wds,1.e-12))**2)
            cff1=cff1+wa_bry_d(jw)*wa_bry(jw)
          enddo
          cff2=sumspec/cff1
          do jw=1,Nfrq
            wa_bry_d(jw)=wa_bry_d(jw)*cff2   ! Normalize dir spectrum
            wa_bry_f(jw)=wa_bry(jw)/sumspec  ! Normalize frq spectrum
            wa_bry(jw)=wa*sqrt(wa_bry_f(jw)  ! Finalize amplitude to
     &                        *wa_bry_d(jw)) ! normalize wave energy
          enddo
#  ifdef WAVE_MAKER_DSPREAD_PER
!
! Force periodicity on all wave components such that
! k_i*sin(theta_i)=p*(2*pi/Ly) with p an integer.
!
          wkdom=2.*pi/el       ! domain wavenumber
          cff0=0.
          do jw=1,Nfrq
            cff0=cff0+wd_bry(jw)*wa_bry_d(jw)*wa_bry_f(jw)
            wd_bry_tmp(jw)=wd_bry(jw)
          enddo
          do jw=1,Nfrq         ! Correct each angle of the spectrum
            cff1=wd_bry(jw)    !   angle before correction
            cff2=wk_bry(jw)    !   associated wavenumber
            wky=cff2*sin(cff1)/wkdom
            mindiff=abs(wky - nint(wky))
            if (wk_bry(jw) .gt. wkdom) then
              iw=1
              do while ((iw .lt. 10000) .and. (mindiff .gt. 1.e-6))  
                cff4 = (cff2 * sin(cff1 + 1.e-4))/wkdom
                cff5 = (cff2 * sin(cff1 - 1.e-4))/wkdom
                if (abs(cff4-nint(cff4)) .lt. mindiff) then
                  mindiff = abs(cff4-nint(cff4))
                  cff1 = cff1 + 1.e-4
                else if (abs(cff5-nint(cff5)) .lt. mindiff) then         
                  mindiff = abs(cff5-nint(cff5))
                  cff1 = cff1 - 1.e-4
                endif
                iw=iw+1
              enddo
              wd_bry_tmp(jw)=cff1
            endif 
          enddo
          cff1=0.              ! Test correction on mean angle
          do jw=1,Nfrq
            cff1=cff1+wd_bry_tmp(jw)*wa_bry_d(jw)*wa_bry_f(jw)
          enddo
          cff2=abs(cff0-cff1)*180/pi
          cff3=abs(cff0-cff1)/max(eps,abs(cff0))
          if (cff2 .lt. 0.5 .or. cff3 .lt. 0.3) then  
            do jw=1,Nfrq                ! correction applied        
              wd_bry(jw)=wd_bry_tmp(jw) ! only if change of mean
            enddo                       ! angle <0.5deg or <30%
          endif
#    ifdef WAVE_MAKER_EAST
          cff0=cff0-pi
          cff1=cff1-pi
#    endif
          MPI_master_only write(stdout,'(6x,A,3(f10.5,1x),f10.3/)')
     &            'Mean wave angle correction (Bef. Aft. Diff %):',
     &                        cff0*180/pi,cff1*180/pi,cff2,cff3*100
#  endif /* WAVE_MAKER_DSPREAD_PER  */

          CALL RANDOM_SEED(SIZE=nseed)
          ALLOCATE(seed(nseed))
          seed = 12345  ! Fix seed for reproducibility
          CALL RANDOM_SEED(PUT=seed)
          call RANDOM_NUMBER(wpha_bry)  ! random phase
          do iw=1,Nfrq
             wpha_bry(iw)=wpha_bry(iw)*2.*pi
             wkx_bry(iw)=wk_bry(iw)*cos(wd_bry(iw))
             wky_bry(iw)=wk_bry(iw)*sin(wd_bry(iw))
          enddo

# else
          CALL RANDOM_SEED(SIZE=nseed)
          ALLOCATE(seed(nseed))
          seed = 12345  ! Fix seed for reproducibility
          CALL RANDOM_SEED(PUT=seed)
          call RANDOM_NUMBER(wpha_bry)  ! random phase
          do iw=1,Nfrq
            wa_bry(iw)=wa*sqrt(wa_bry(iw)/sumspec) ! normalize
            wpha_bry(iw)=wpha_bry(iw)*2.*pi
            wkx_bry(iw)=wk_bry(iw)*cos(wd)
            wky_bry(iw)=wk_bry(iw)*sin(wd)
          enddo
# endif /* WAVE_MAKER_DSPREAD */

        endif ! FIRST_TIME_STEP

        ramp=tanh(dt/wp*float(iic-ntstart))

#elif defined WAVE_MAKER_BICHROMATIC
!
!  Bichromatic waves
!
        ramp=tanh(dt/2.*float(iic-ntstart))
        wa1=ramp*wa1
        wa2=ramp*wa2
        khd=h0*wf1*wf1/g   ! compute wavenumber
        wk1=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                    khd*(K5+K6*khd)))))) )/h0
        khd=h0*wf2*wf2/g
        wk2=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                    khd*(K5+K6*khd)))))) )/h0
        wk=0.5*(wk1+wk2);
        wa=0.707*sqrt(wa1**2+wa2**2)
        wf=0.5*(wf1+wf2)
        wp=2*pi/wf
#else
!
!  Monochromatic waves
!
        ramp=tanh(dt/wp*float(iic-ntstart))
        wa=wa*ramp
#endif /* ROGUE_WAVE ... */
!
!--------------------------------------------------------------------
!  Sea level zetabry
!--------------------------------------------------------------------
!
#ifdef Z_FRC_BRY
        do j=JstrR,JendR
          h0=h(IB0,j)
# ifdef WAVE_MAKER_SPECTRUM
#  ifdef DUCK94
          ZBRY(j)=0.7  ! add tidal level
#  else
          ZBRY(j)=0.
#  endif
          do iw=1,Nfrq   ! frequency spread
            theta=(xr(IB0,j)-x0)*wkx_bry(iw)
     &           +(yr(IB0,j)-y0)*wky_bry(iw)
     &              -(time-time0)*wf_bry(iw)
     &                         -wpha_bry(iw)
            ZBRY(j)=ZBRY(j) + ramp*wa_bry(iw)
     &                            *cos(theta)
          enddo

# elif defined WAVE_MAKER_BICHROMATIC /* eg GLOBEX B2/B3 */
          cff1=tanh(wk1*h0)
          cff2=tanh(wk2*h0)
          cff3=wa1*wa1*wk1*(3.-cff1**2)/(4.*cff1**3)
          cff4=wa2*wa2*wk2*(3.-cff2**2)/(4.*cff2**3)
          ZBRY(j)=wa1*cos(wf1*time)
     &           +wa2*cos(wf2*time)
     &          +cff3*cos(2*wf1*time)
     &          +cff4*cos(2*wf2*time)

# else
#  ifdef DUCK94
          ZBRY(j)=0.7
#  else
          ZBRY(j)=0.
#  endif
          time0=0.
          khd=h0*wf*wf/g   ! compute wavenumber
          wk=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                     khd*(K5+K6*khd)))))) )/h0
          sigma=tanh(wk*h0)
          theta=(xr(IB0,j)-x0)*wk*coswd*coswds
     &         +(yr(IB0,j)-y0)*wk*sinwd*coswds
     &         -(time-time0)*wf
          cff_spread=cos(-(xr(IB0,j)-x0)*sinwd*sinwds*wk
     &                   +(yr(IB0,j)-y0)*coswd*sinwds*wk)
          ZBRY(j)=ZBRY(j) + ( wa*cos(theta)
#  ifdef STOKES_WAVES
     &                      +wk*wa*wa*(3.-sigma**2)/
     &                   (4.*sigma**3)*cos(2.*theta)
#  endif
     &                     )*cff_spread
# endif /* ROGUE_WAVE ... */
        enddo  ! j loop
#endif /* Z_FRC_BRY */
!
!--------------------------------------------------------------------
!  XI velocity components ubry and ubarbry
!--------------------------------------------------------------------
!
#ifdef M3_FRC_BRY
        do j=JstrR,JendR
          h0=0.5*(h(IB0,j)+h(IB1,j))
          Du=h0
# ifdef ROGUE_WAVE
          do k=1,N
            UBRY(j,k)=0.
            Zu=Du+0.5*(z_r(IB0,j,k)+z_r(IB1,j,k))
            do iw=1,Nfrq
              UBRY(j,k)=UBRY(j,k) + ramp*wa_bry(iw)*wf_bry(iw)*
     &                                         cosh(wk_bry(iw)*Zu)/
     &                                         sinh(wk_bry(iw)*Du)
     &          *cos((0.5*(xr(IB0,j)+xr(IB1,j))-x0)*wk_bry(iw)
     &                                -(time-time0)*wf_bry(iw)
     &                                           -wpha_bry(iw))
            enddo
          enddo

# elif defined WAVE_MAKER_SPECTRUM
          khd=h0*wf*wf/g     ! compute mean wavenumber
          wk=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                     khd*(K5+K6*khd)))))) )/h0
          do k=1,N
            Zu=Du+0.5*(z_r(IB0,j,k)+z_r(IB1,j,k))
            UBRY(j,k)=ZBRY(j)*wf*cos(wd)*cosh(wk*Zu)
     &                                  /sinh(wk*Du)
          enddo

# elif defined WAVE_MAKER_BICHROMATIC
          cff1=wa1*cos(wf1*time)*wf1/sinh(wk1*Du)
          cff2=wa2*cos(wf2*time)*wf2/sinh(wk2*Du)
          cff3=0.75*wa1*wa1*wk1*wf1/(sinh(wk1*Du))**4
     &                             *cos(2*wf1*time)
          cff4=0.75*wa2*wa2*wk2*wf2/(sinh(wk2*Du))**4
     &                             *cos(2*wf2*time)
          do k=1,N
            Zu=Du+0.5*(z_r(IB0,j,k)+z_r(IB1,j,k))
            UBRY(j,k)=cff1*cosh(wk1*Zu)
     &               +cff2*cosh(wk2*Zu)
     &               +cff3*cosh(2*wk1*Zu)
     &               +cff4*cosh(2*wk2*Zu)
          enddo

# else
          xu=0.5*(xr(IB0,j)+xr(IB1,j))-x0
          yu=0.5*(yr(IB0,j)+yr(IB1,j))-y0
          theta=xu*coswd*coswds*wk
     &         +yu*sinwd*coswds*wk
     &            -(time-time0)*wf
          cff_spread=cos(-xu*sinwd*sinwds*wk
     &                   +yu*coswd*sinwds*wk)
          cff=wa*wf*coswd/sinh(wk*Du)*cos(theta)*cff_spread
#  ifdef STOKES_WAVES
          cff2=cff*cos(2*theta)/cos(theta)
     &            *0.75*wk*wa/(sinh(wk*Du))**3
#  endif
          do k=1,N
            Zu=Du+0.5*(z_r(IB0,j,k)+z_r(IB1,j,k))
            UBRY(j,k)=cff*cosh(wk*Zu)
#  ifdef STOKES_WAVES
     &              +cff2*cosh(2*wk*Zu)
#  endif
          enddo
# endif /* ROGUE_WAVE */

        enddo  ! j loop

        do j=JstrR,JendR                  ! compensation flow
          Du=0.5*(h(IB0,j)+h(IB1,j))
          cff1=0.5*g*wa*wa*wk/(wf*Du)
          do k=1,N
            UBRY(j,k)=UBRY(j,k) - cff1
          enddo
        enddo

#endif /* M3_FRC_BRY */

#ifdef M2_FRC_BRY
        do j=JstrR,JendR
          cff4=0.
          cff5=0.
          do k=1,N
            cff3=Hz(IB0,j,k)+Hz(IB1,j,k)
            cff4=cff4+UBRY(j,k)*cff3
            cff5=cff5+cff3
          enddo
          UBARBRY(j)=cff4/cff5
        enddo
#endif /* M2_FRC_BRY */
!
!--------------------------------------------------------------------
!  ETA velocity components vbry and vbarbry
!--------------------------------------------------------------------
!
#ifdef M3_FRC_BRY
        do j=JstrV,JendR
          h0=0.5*(h(IB0,j)+h(IB0,j-1))
          Dv=h0
# if defined ROGUE_WAVE || defined SWASH || defined SANDBAR
          do k=1,N
            VBRY(j,k)=0.
          enddo

# elif defined WAVE_MAKER_SPECTRUM || defined WAVE_MAKER_BICHROMATIC
          khd=h0*wf*wf/g   ! compute mean wavenumber
          wk=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                     khd*(K5+K6*khd)))))) )/h0
          do k=1,N
            Zv=Dv+0.5*(z_r(IB0,j,k)+z_r(IB0,j-1,k))
            VBRY(j,k)=ZBRY(j)*wf*sin(wd)*cosh(wk*Zv)
     &                                  /sinh(wk*Dv)
          enddo

# else
          xv=0.5*(xr(IB0,j)+xr(IB0,j-1))-x0
          yv=0.5*(yr(IB0,j)+yr(IB0,j-1))-y0
          theta=xv*coswd*coswds*wk
     &         +yv*sinwd*coswds*wk
     &            -(time-time0)*wf
          cff_spread=cos(-xv*sinwd*sinwds*wk
     &                   +yv*coswd*sinwds*wk)
          cff=wa*wf*sinwd/sinh(wk*Dv)*cos(theta)*cff_spread
#  ifdef STOKES_WAVES
          cff2=cff*cos(2*theta)/cos(theta)
     &            *0.75*wk*wa/(sinh(wk*Dv))**3
#  endif
          do k=1,N
            Zv=Dv+0.5*(z_r(IB0,j,k)+z_r(IB0,j-1,k))
            VBRY(j,k)=cff*cosh(wk*Zv)
#  ifdef STOKES_WAVES
     &              +cff2*cosh(2*wk*Zv)
#  endif
          enddo
# endif /* WAVE_MAKER_SPECTRUM */
        enddo  ! j loop
#endif /* M3_FRC_BRY */

#ifdef M2_FRC_BRY
        do j=JstrV,JendR
          cff4=0.
          cff5=0.
          do k=1,N
            cff4=cff4+VBRY(j,k)*(Hz(IB0,j,k)+Hz(IB0,j-1,k))
            cff5=cff5+(Hz(IB0,j,k)+Hz(IB0,j-1,k))
          enddo
          VBARBRY(j)=cff4/cff5
        enddo
#endif /* M2_FRC_BRY */
!
!--------------------------------------------------------------------
!  Z velocity component wbry
!--------------------------------------------------------------------
!
#ifdef W_FRC_BRY
        do j=JstrR,JendR
          Dr=h(IB0,j)
# ifdef WAVE_MAKER_SPECTRUM
          do k=1,N
            WBRY(j,k)=0.
#  ifndef WAVE_MAKER_DSPREAD
            Zr=Dr+z_w(IB0,j,k)
            do iw=1,Nfrq
              WBRY(j,k)=WBRY(j,k) + 
     &                    ramp*wa_bry(iw)*wf_bry(iw)*
     &                               sinh(wk_bry(iw)*Zr)/
     &                               sinh(wk_bry(iw)*Dr)
     &                *sin((xr(IB0,j)-x0)*wk_bry(iw)
     &                      -(time-time0)*wf_bry(iw)
     &                                 -wpha_bry(iw))
            enddo
#  endif
          enddo

# elif defined WAVE_MAKER_BICHROMATIC
          cff1=wa1*cos(wf1*time)*wf1/sinh(wk1*Dr)
          cff2=wa2*cos(wf2*time)*wf2/sinh(wk2*Dr)
          cff3=0.75*wa1*wa1*wk1*wf1/(sinh(wk1*Dr))**4
     &                               *cos(2*wf1*time)
          cff4=0.75*wa1*wa1*wk1*wf1/(sinh(wk2*Dr))**4
     &                               *cos(2*wf2*time)
          do k=0,N
            Zr=Dr+z_w(IB0,j,k)
            WBRY(j,k)=cff1*sinh(wk1*Zr)
     &               +cff2*sinh(wk2*Zr)
     &               +cff3*sinh(2*wk1*Zr)
     &               +cff4*sinh(2*wk2*Zr)
          enddo

# else
          theta=(xr(IB0,j)-x0)*coswd*coswds*wk
     &         +(yr(IB0,j)-y0)*sinwd*coswds*wk - time*wf
          cff=wa*wf/sinh(wk*Dr)*cos(theta)
          cff2=0.75*cff*cos(2*theta)/cos(theta)*wk*wa/(sinh(wk*Dr))**3
          do k=1,N
            Zr=Dr+z_w(IB0,j,k)
            WBRY(j,k)=cff*sinh(wk*Zr)
     &              +cff2*sinh(2*wk*Zr)
          enddo
# endif /* WAVE_MAKER_SPECTRUM ... */

        enddo  ! j loop
#endif /* W_FRC_BRY */
!
!--------------------------------------------------------------------
!  TRACERS tbry
!--------------------------------------------------------------------
!
#ifdef T_FRC_BRY
        if (FIRST_TIME_STEP) then
          do k=1,N
            do j=JstrR,JendR
              do itrc=1,NT
                TBRY(j,k,itrc)=t(IB0,j,k,1,itrc)
              enddo
            enddo
          enddo
        endif
#endif
!
!--------------------------------------------------------------------
!  NBQ variables: unqbry, vnbqbru and wnbqbry
!--------------------------------------------------------------------
!
#ifdef NBQ_FRC_BRY
        do k=1,N
          do j=JstrR,JendR
            UNBQBRY(j,k)=0.5*(Hz(IB0,j,k)+Hz(IB1,j,k))
     &                                          *UBRY(j,k)
          enddo
          do j=JstrV,JendR
            VNBQBRY(j,k)=0.5*(Hz(IB0,j,k)+Hz(IB0,j-1,k))
     &                                          *VBRY(j,k)
          enddo
        enddo
        do k=1,N-1
          do j=JstrR,JendR
            WNBQBRY(j,k)=0.5*(Hz(IB0,j,k)+Hz(IB0,j,k+1))
     &                                          *WBRY(j,k)
          enddo
        enddo
        k=N
        do j=JstrR,JendR
          WNBQBRY(j,k)=0.5*Hz(IB0,j,k)*WBRY(j,k)
        enddo
#endif

#undef ZBRY
#undef UBARBRY
#undef VBARBRY
#undef UBRY
#undef VBRY
#undef WBRY
#undef TBRY
#undef UNBQBRY
#undef VNBQBRY
#undef WNBQBRY
#undef IB0
#undef IB1

