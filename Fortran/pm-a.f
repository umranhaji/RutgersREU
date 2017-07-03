c     --- Given the average proper motion (mu_RA, mu_Dec) in mas/cent,
c         the program converts this proper motion to
c         heliocentric and galactocentric proper motions. In adition,
c         it also calculates the space velocity of the dSph.  The main
c	  purpose of this program is to write a file of the Monte Carlo
c	  Pi, Theta, and Z that can be used to Monte Carlo orbit calculations.
c     --- 16-Jun-2016:  This version uses the LSR and solar motion from
c         Bland-Hawthorn & Gerhard 2016 ARAA.  Note that we only care
c         about the total motion of the Sun in the V direction, not what
c         the circular velocity is at the LSR.  So vlsr includes any
c         possible motion of the LSR wrt V_circ. B-H & G conclude that
c         V_g,sun = 248 +- 3 km/s.
c------------------------------------------------------------------------
c------------------------------------------------------------------------
      IMPLICIT NONE
      REAl*8 alphagp,deltagp,lncp,bncp,gcgan
      REAL*4 asecrad,secyear,kmpc
      REAL*4 pscale1,pscale2,dlsr
      REAL*4 mas100
c     --- asecrad Conversion factor from arc-seconds to radians.
c     --- secyear Conversion factor from year to seconds.
c     --- kmpc    Conversion factor from pc to km.
c     --- yr      Conversion factor from days to years.
c     --- ralsr   RA of the apex of the LSR motion (ie., l=90, b=0)
c     --- declsr  DEC of the apex of the LSR motion (J2000.0)
c     --- Note that rasa and decsa should be consistent with usun,vsun,wsun
c     --- rasa    RA of the apex of the solar peculiar motion (wrt LSR) J2000.0
c     --- decsa   DEC of the apex of the solar peculiar motion (wrt LSR)
c     --- alphagp RA in degrees of the Galactic Pole (J2000.0)
c     --- deltagp DEC in degrees of the Galactic Pole (J2000.0)
c     --- lncp Galactic longitude of the NCP (J2000.0)
c     --- bncp Galactic latitude of the NCP  (J2000.0)
c     --- gcgan Angle along the galactic equator from (l,b)=(0,0) to
c         the galactic ascending node.
c     --- usun U-component of Sun's velocity w.r.t. LSR; Dehnen & Binney
c     --- vsun V-component of Sun's velocity w.r.t. LSR; 1998, MNRAS,
c     --- wsun W-component of Sun's velocity w.r.t. LSR; 298, 387
c     --- pscale1 arcsec per pixel after distortion correction WFPC)
c     --- pscale2 arcsec per pixel after distortion correction STIS)
c     --- vlsr Circular velocity of LSR w.r.t. Galaxy center in km/s
c     --- dlsr Distance of LSR from the Galaxy center in pc.
      REAL*8 ralsr,declsr,rasa,decsa,sind,cosd
      REAL*4 yr,usun,vsun,wsun,vlsr
      PARAMETER(asecrad=206264.8,secyear=3.155693E+7,kmpc=3.08568E+13)
      PARAMETER(mas100=1.0e5)
      PARAMETER(yr=365.2422,ralsr=318.00425d0,declsr=48.329611d0)
      PARAMETER (rasa=251.51,decsa=10.129)
c     --- alphgp, deltagp, gcgan taken from Perryman et al 1997 (Hipparcos
c         catalog description).
      PARAMETER(alphagp=192.85948d0,deltagp=27.12825d0)
      PARAMETER(lncp=122.9319d0,bncp=27.1283d0,gcgan=32.93192d0)
ccc      PARAMETER(alphagp=192.85933,deltagp=27.1282528)
ccc      PARAMETER(lncp=122.9319,bncp=27.1283,gcgan=32.936807)
      PARAMETER(usun=-10.0,vsun=11.0,wsun=7.0,vlsr=237.0,dlsr=8200.0)
c     --- uncertainties: usun: 0.36, vsun: 0.62, wsun: 0.38 km/s
cccc      PARAMETER(pscale1=0.046,pscale2=0.05071,vlsr=220.0,dlsr=8500.0)
      PARAMETER(pscale1=0.046,pscale2=0.05071)
      INTEGER nsim
      PARAMETER (nsim=100000)
c     *** For all of the variables described below, the first element of
c     *** the array is the value and the second element is the uncertainty.
c     --- muxpy  Proper motion in alpha (RA) direction in pixels per year.
c     --- muypy  Proper motion in delta (DEC) direction in pixels per year.
c     --- mux    Proper motion in alpha (RA) direction in arcsec/year.
c     --- muy    Proper motion in delta (DEC) direction in arcsec/year.
c     --- mul    Proper motion in l (galactic) direction in arcsec/year.
c     --- mub    Proper motion in b (galactic) direction in arcsec/year.
c     --- ug     Velcity in Pi direction of the galaxy w.r.t. LSR
c     --- vg     Velcity in Theta direction of the galaxy w.r.t. LSR
c     --- wg     Velcity in Z direction of the galaxy w.r.t. LSR
c     --- upg    Velocity in Pi direction of the galaxy.   (Schweitzer et al 
c     --- vpg    Velocity in Theta direction of the galaxy.     1995, AJ, 
c     --- wpg    Velocity in Z direction of the galaxy.         110,2747)
c     --- uppg   Velocity in the radial direction as seen from the center
c                of the Galaxy, i.e. radial velocity of a galaxy.
c     --- vppg   Velocity of a galaxy in the direction of rotation of the
c                Galactic disk.
c     --- wppg   Velocity of a galaxy up and perpendicular to the radial
c                velocity.
      REAL*4 mux(2),muy(2),mul(2),mub(2)
      REAL*4 ug(2),vg(2),wg(2),upg(2),vpg(2),wpg(2)
      REAL*4 uppg(2),vppg(2),wppg(2)
      REAL*4 pa,d,vrs(2)
      REAL*4 vtx(2),vty(2),vt(2),vtl(2),vtb(2),vtg(2)
      REAL*8 ra,dec,lg,bg
      REAL*4 muxc(2),muyc(2),mulc(2),mubc(2)
ccc      REAL*4 e,pos,rt
      CHARACTER S*1
      REAL*4 mux0(2),muy0(2),normal,x0,y0,z0
      REAL*4 smux,smux2,smuy,smuy2,smul,smul2,smub,smub2,supg,supg2
      REAL*4 svpg,svpg2,swpg,swpg2,suppg,suppg2,svtan,svtan2
      REAL*4 smuxc,smuxc2,smuyc,smuyc2,smulc,smulc2,smubc,smubc2
      REAL*4 vrs0(2),rand,vtan
      REAL*4 vpi0,vtheta0,vz0,vrad0,vtan0
      INTEGER i,iran
ccc 100  format(a32)
      OPEN(UNIT=1,FILE='pm-a.par',STATUS='OLD',form='FORMATTED')
      read(1,'(a1)') S           ! E or G coordinate system for mu
      read(1,*) ra,dec           ! RA & Dec in degrees
      read(1,*) mux0(1),mux0(2)  ! 1st component of proper motion
      read(1,*) muy0(1),muy0(2)  ! 2nd component of proper motion
      read(1,*) d                ! Distance in kpc
      read(1,*) vrs0(1),vrs0(2)  ! Radial velocity + error in km/s
      read (1,*) iran            ! Random seed
      CLOSE(UNIT=1)
      d=d*1000.0                 ! Convert distance to pc
      mux0(1)=mux0(1)/mas100     ! Convert to arcsec/year
      mux0(2)=mux0(2)/mas100
      muy0(1)=muy0(1)/mas100
      muy0(2)=muy0(2)/mas100
      pa=rand(iran)
      do i=1,1000
         pa=normal(0)
      end do
      write(6,101) ra
101   format('Average RA:  ',f9.4,' degrees')
      write(6,102) dec
102   format('Average DEC: ',f9.5,' degrees')
c---------------------------------------------------------------------
c     --- Calculate (l,b) of the galaxy.
      CALL FINDLGBG(ra,dec,alphagp,deltagp,gcgan,LG,BG)
c-----------------------------------------------------------------------      
      write(6,*) 'The galactic coordinates of the galaxy are:'
      write(6,103) lg
103   format('l= ',f9.4,' degrees')
      write(6,104) bg
104   format('b= ',f9.5,' degrees')
      write(6,*)''

c     --- Now loop over the following, drawing observed muxpy and muypy
      smux   =0.0
      smux2  =0.0
      smuy   =0.0
      smuy2  =0.0
c
      smul   =0.0
      smul2  =0.0
      smub   =0.0
      smub2  =0.0
c
      supg   =0.0
      supg2  =0.0
      svpg   =0.0
      svpg2  =0.0
      swpg   =0.0
      swpg2  =0.0
c
      suppg  =0.0
      suppg2 =0.0
      svtan  =0.0
      svtan2 =0.0
      smuxc  =0.0
      smuxc2 =0.0
      smuyc  =0.0
      smuyc2 =0.0
      smulc  =0.0
      smulc2 =0.0
      smubc  =0.0
      smubc2 =0.0
      OPEN(UNIT=20,file='ivel.data',form='formatted',status='unknown')
      do i=1,nsim
         if (i.eq.1) then
            mux(1)=mux0(1)
            mux(2)=mux0(2)
            muy(1)=muy0(1)
            muy(2)=muy0(2)
            else
               mux(1)=mux0(1)+normal(0)*mux0(2)
               mux(2)=mux0(2)
               muy(1)=muy0(1)+normal(0)*muy0(2)
               muy(2)=muy0(2)
         end if
         smux =smux +mux(1)
         smux2=smux2+mux(1)*mux(1)
         smuy =smuy +muy(1)
         smuy2=smuy2+muy(1)*muy(1)
         if (i.eq.1) then
            if(S(1:1).eq.'E') then
               write(6,*) 'RA-dir of proper motion in mas/cent:  ',
     &                     mux(1)*mas100,'+-',mux(2)*mas100
               write(6,*) 'DEC-dir of proper motion in mas/cent: ',
     &                     muy(1)*mas100,'+-',muy(2)*mas100
               write(6,*)
            else if(S(1:1).eq.'G') then
               write(6,*) 'l-dir of proper motion in mas/cent:  ',
     &                     mux(1)*mas100,'+-',mux(2)*mas100
               write(6,*) 'b-dir of proper motion in mas/cent: ',
     &                     muy(1)*mas100,'+-',muy(2)*mas100
               write(6,*)
            end if
         end if
c---------------------------------------------------------------------
c        --- Convert the proper motion from the equatorial coordinate
c             system to the galactic system.
         if (S(1:1).eq.'E') then
            CALL GETMULMUB(ra,dec,lg,bg,alphagp,deltagp,mux,muy,MUL,MUB)
         else if(S(1:1).eq.'G') then
            mul(1)=mux(1)
            mul(2)=mux(2)
            mub(1)=muy(1)
            mub(2)=muy(2)
         end if
         smul=smul+mul(1)
         smul2=smul2+mul(1)*mul(1)
         smub=smub+mub(1)
         smub2=smub2+mub(1)*mub(1)
         if (i.eq.1.and.S(1:1).eq.'E') then
            write(6,*) 'l-dir of proper motion in mas/cent: ',
     &                  mul(1)*mas100,'+-',mul(2)*mas100
            write(6,*) 'b-dir of proper motion in mas/cent: ',
     &                  mub(1)*mas100,'+-',mub(2)*mas100
            write(6,*)
         end if
c
c--------------------------------------------------------------------
c     --- Calculate the initial Galactic X,Y,Z position
      if (i.eq.1) then
         x0=dlsr-d*cosd(bg)*cosd(lg)
         y0=-d*cosd(bg)*sind(lg)
         z0=d*sind(bg)
         write (6,*) 'initial x,y,z (kpc) = ',x0/1000.,y0/1000.,z0/1000.
         write (6,*)
      end if
c
c--------------------------------------------------------------------
c     ---Calculate the heliocentric tangential velocity and its components.
         if(S(1:1).eq.'E'.and.i.eq.1) then
            CALL GETVT(mux,muy,d,asecrad,secyear,kmpc,VTX,VTY,VT)
            write(6,*) 'Tangential velocity in RA dir:  ',
     &                  vtx(1),'+-',vtx(2),' km/s'
            write(6,*) 'Tangential velocity in DEC dir: ',
     &                  vty(1),'+-',vty(2),' km/s'
            write(6,*) 'Tangential velocity magnitude:  ',
     &                  vt(1),'+-',vt(2), ' km/s'
            write(6,*)
         end if
         CALL GETVT(mul,mub,d,asecrad,secyear,kmpc,VTL,VTB,VTG)
         if (i.eq.1) then
            write(6,*) 'Tangential velocity in l dir:  ',
     &                  vtl(1),'+-',vtl(2),' km/s'
            write(6,*) 'Tangential velocity in b dir:  ',
     &                  vtb(1),'+-',vtb(2),' km/s'
            write(6,*) 'Tangential velocity magnitude: ',
     &                  vtg(1),'+-',vtg(2),' km/s'
            write(6,*)
         end if
c--------------------------------------------------------------------     
         if (i.eq.1) then
            vrs(1)=vrs0(1)
            vrs(2)=vrs0(2)
            else
               vrs(1)=vrs0(1)+normal(0)*vrs0(2)
               vrs(2)=vrs0(2)
         end if
c-------------------------------------------------------------------
c        --- Calculate ug,vg,and wg components with respect to LSR
         CALL GETUGVGWG(lg,bg,vtl,vtb,vrs,usun,vsun,wsun,UG,VG,WG)
         if (i.eq.1) then
            write(6,*)'u, v, and w velocities of the galaxy w.r.t. LSR'
            write(6,*) 'U:',ug(1),'+-',ug(2),' km/s'
            write(6,*) 'V:',vg(1),'+-',vg(2),' km/s'
            write(6,*) 'W:',wg(1),'+-',wg(2),' km/s'
            write(6,*)
         end if
c--------------------------------------------------------------------
c     --- Calculate upg,vpg,and wpg components in a coordinate system
c         used by Schweitzer et al. 1995, 110, 2747
         CALL GETUPGVPGWPG(dlsr,vlsr,d,lg,bg,ug,vg,wg,UPG,VPG,WPG)
         supg=supg+upg(1)
         supg2=supg2+upg(1)*upg(1)
         svpg=svpg+vpg(1)
         svpg2=svpg2+vpg(1)*vpg(1)
         swpg=swpg+wpg(1)
         swpg2=swpg2+wpg(1)*wpg(1)
         if (i.eq.1) then
            vpi0=upg(1)
            vtheta0=vpg(1)
            vz0=wpg(1)
            write(6,*) 'PI   :',upg(1),'+-',upg(2),' km/s'
            write(6,*) 'THETA:',vpg(1),'+-',vpg(2),' km/s'
            write(6,*) 'Z    :',wpg(1),'+-',wpg(2),' km/s'
            write(6,*) ''
         end if
c        --- Write out the velocity components for orbit Monte Carlo
         write (20,*) upg(1),upg(2),vpg(1),vpg(2),wpg(1),wpg(2)
c-------------------------------------------------------------------
c        --- Calculate uppg, vppg, and wppg components.
         CALL GETUPPGVPPGWPPG(dlsr,d,lg,bg,upg,vpg,wpg,UPPG,VPPG,WPPG)
         suppg=suppg+uppg(1)
         suppg2=suppg2+uppg(1)*uppg(1)
         vtan=sqrt(vppg(1)*vppg(1)+wppg(1)*wppg(1))
         svtan=svtan+vtan
         svtan2=svtan2+vtan*vtan
         if (i.eq.1) then
            vrad0=uppg(1)
            vtan0=sqrt(vppg(1)*vppg(1)+wppg(1)*wppg(1))
            write(6,*)'In galactocentric coordinate system...'
            write(6,*)''
            write(6,*) 'Radial velocity, u:                  ',
     &                  uppg(1),'+-',uppg(2),' km/s'
            write(6,*) 'Velocity in dir. of disk rotation, v:',
     &                  vppg(1),'+-',vppg(2),' km/s'
            write(6,*) 'Velocity up and perp. to radial      ',
     &                  wppg(1),'+-',wppg(2),' km/s'
            write(6,*) 'Tangential velocity                  ',
     &      sqrt(vppg(1)*vppg(1)+wppg(1)*wppg(1)),'+-',
     &          (1.0/sqrt(vppg(1)*vppg(1)+wppg(1)*wppg(1)))*
     &      sqrt(vppg(1)*vppg(1)*vppg(2)*vppg(2)+
     &           wppg(1)*wppg(1)*wppg(2)*wppg(2)),' km/s'
            write(6,*)
         end if
c
c--------------------------------------------------------------------
c     ---Remove the effect of the Sun's motion (both peculiar and LSR)
c        from the proper motions in the equatorial system
         if(S(1:1).eq.'E') then
            CALL GETVNUVTAU(ra,dec,d,ralsr,declsr,vlsr,
     &                     rasa,decsa,usun,vsun,wsun,mux,muy,MUXC,MUYC)
            smuxc=smuxc+muxc(1)
            smuxc2=smuxc2+muxc(1)*muxc(1)
            smuyc=smuyc+muyc(1)
            smuyc2=smuyc2+muyc(1)*muyc(1)
            if (i.eq.1) then
               write(6,*) 
     &         'Proper motion in mas/cent corrected for solar motion:'  
                write(6,*) 'RA component:  ',
     &                      muxc(1)*mas100,'+-',muxc(2)*mas100
                write(6,*) 'DEC component: ',
     &                      muyc(1)*mas100,'+-',muyc(2)*mas100
                write(6,*)
            end if
c---------------------------------------------------------------------
c           --- Convert the corrected proper motion from the equatorial 
c              coordinate system to the galactic system.
            CALL GETMULMUB(ra,dec,lg,bg,alphagp,deltagp,muxc,muyc,
     &                                                 MULC,MUBC)
            smulc=smulc+mulc(1)
            smulc2=smulc2+mulc(1)*mulc(1)
            smubc=smubc+mubc(1)
            smubc2=smubc2+mubc(1)*mubc(1)
            if (i.eq.1) then
               write(6,*) 
     &        'Proper motion in mas/cent corrected for solar motion:'
               write(6,*) 'l component: ',
     &         mulc(1)*mas100,'+-',mulc(2)*mas100
               write(6,*) 'b component: ',
     &         mubc(1)*mas100,'+-',mubc(2)*mas100
               write(6,*)
            end if
         else if(S(1:1).eq.'G') then
         end if
      end do
      close(UNIT=20)
c
c     --- Calculate the means and rms's
      if(S(1:1).eq.'E') then
         smux=smux/float(nsim)
         smux2=mas100*sqrt((smux2-float(nsim)*smux*smux)/
     &                                        float(nsim-1))
         smux=mas100*smux
         write (6,*) 'mu_RA       (mas/cent):',smux,'+-',smux2
         smuy=smuy/float(nsim)
         smuy2=mas100*sqrt((smuy2-float(nsim)*smuy*smuy)/
     &                                        float(nsim-1))
         smuy=mas100*smuy
         write (6,*) 'mu_Dec      (mas/cent):',smuy,'+-',smuy2
         write(6,*)
      end if
c
      smul=smul/float(nsim)
      smul2=mas100*sqrt((smul2-float(nsim)*smul*smul)/
     &                                         float(nsim-1))
      smul=mas100*smul
      write (6,*) 'mu_l        (mas/cent):',smul,'+-',smul2
      smub=smub/float(nsim)
      smub2=mas100*sqrt((smub2-float(nsim)*smub*smub)/
     &                                         float(nsim-1))
      smub=mas100*smub
      write (6,*) 'mu_b        (mas/cent):',smub,'+-',smub2
      write(6,*)
c
      supg=supg/float(nsim)
      supg2=sqrt((supg2-float(nsim)*supg*supg)/
     &                                         float(nsim-1))
      write (6,*) 'Pi              (km/s):',supg,'+-',supg2,
     &          '  bias =',supg-vpi0
      svpg=svpg/float(nsim)
      svpg2=sqrt((svpg2-float(nsim)*svpg*svpg)/
     &                                         float(nsim-1))
      write (6,*) 'Theta           (km/s):',svpg,'+-',svpg2,
     &          '  bias =',svpg-vtheta0
      swpg=swpg/float(nsim)
      swpg2=sqrt((swpg2-float(nsim)*swpg*swpg)/
     &                                         float(nsim-1))
      write (6,*) 'Z               (km/s):',swpg,'+-',swpg2,
     &          '  bias =',swpg-vz0
      write(6,*)
c
      suppg=suppg/float(nsim)
      suppg2=sqrt((suppg2-float(nsim)*suppg*suppg)/
     &                                         float(nsim-1))
      write (6,*) 'V_R             (km/s):',suppg,'+-',suppg2,
     &          '  bias =',suppg-vrad0
      svtan=svtan/float(nsim)
      svtan2=sqrt((svtan2-float(nsim)*svtan*svtan)/
     &                                         float(nsim-1))
      write (6,*) 'V_T             (km/s):',svtan,'+-',svtan2,
     &          '  bias =',svtan-vtan0
      write(6,*)
c
      if(S(1:1).eq.'E') then
         smuxc=smuxc/float(nsim)
         smuxc2=mas100*sqrt((smuxc2-float(nsim)*smuxc*smuxc)/
     &                                         float(nsim-1))
         smuxc=mas100*smuxc
         write (6,*) 'mu_RA(GRF)  (mas/cent):',smuxc,'+-',smuxc2
         smuyc=smuyc/float(nsim)
         smuyc2=mas100*sqrt((smuyc2-float(nsim)*smuyc*smuyc)/
     &                                         float(nsim-1))
         smuyc=mas100*smuyc
         write (6,*) 'mu_Dec(GRF) (mas/cent):',smuyc,'+-',smuyc2
         write(6,*)
c
         smulc=smulc/float(nsim)
         smulc2=mas100*sqrt((smulc2-float(nsim)*smulc*smulc)/
     &                                         float(nsim-1))
         smulc=mas100*smulc
         write (6,*) 'mu_l(GRF)   (mas/cent):',smulc,'+-',smulc2
         smubc=smubc/float(nsim)
         smubc2=mas100*sqrt((smubc2-float(nsim)*smubc*smubc)/
     &                                         float(nsim-1))
         smubc=mas100*smubc
         write (6,*) 'mu_b(GRF)   (mas/cent):',smubc,'+-',smubc2
      end if
c
      stop
      end
c********************************************************************
C     FUNCTION SUBROUTINE TO RETURN A NORMALLY DISTRIBUTED
C     RANDOM NUMBER WITH UNIT STANDARD DEVIATION.
C
      REAL FUNCTION NORMAL(IRAN)
C
      IMPLICIT NONE
      REAL*4 NORMAL,U1,U2,S,TMP
      INTEGER ISW,IRAN
      SAVE ISW,TMP
C
      IF (ISW.EQ.1) GOTO 2
1     U1=2.00*RAND(IRAN)-1.00
      U2=2.00*RAND(IRAN)-1.00
      S=U1*U1+U2*U2
      IF (S.GE.1.0) GOTO 1
      TMP=SQRT(-2.00*ALOG(S)/S)
      NORMAL=U1*TMP
      TMP=U2*TMP
      ISW=1
      RETURN
C
2     NORMAL=TMP
      ISW=0
      RETURN
      END

c********************************************************************
      SUBROUTINE GETUPPGVPPGWPPG(dlsr,d,lg,bg,upg,vpg,wpg,
     &                                     UPPG,VPPG,WPPG)
      IMPLICIT NONE
      REAL*8 lg,bg,sind,cosd
      REAL*4 dlsr,d,upg(2),vpg(2),wpg(2)
      REAL*4 uppg(2),vppg(2),wppg(2)
c--------------------------------------------------------------------
c                 --- The transformation equations: ---
c       uppg =  upg*cos(beta) + wpg*sin(beta)
c       vppg =  vpg
c       wppg = -upg*sin(beta) + wpg*cos(beta)
c
c     --- beta is an angle between the line-of-sight to a galaxy and
c         line x (see subroutine GETUPGVPGWPG for definition).
      REAL*4 x, cosbeta,sinbeta
c     --- From the law of cosines...
      x=dlsr*dlsr+d*d*cosd(bg)*cosd(bg)-2.0*dlsr*d*cosd(bg)*cosd(lg)
      x=sqrt(x)
      cosbeta=x/sqrt(x*x+d*d*sind(bg)*sind(bg))
      sinbeta=d*sind(bg)/sqrt(x*x+d*d*sind(bg)*sind(bg))
      uppg(1) =  upg(1)*cosbeta + wpg(1)*sinbeta
      vppg(1) =  vpg(1)
      wppg(1) = -upg(1)*sinbeta + wpg(1)*cosbeta
      uppg(2) =  sqrt((upg(2)*cosbeta)**2 + (wpg(2)*sinbeta)**2)
      vppg(2) =  vpg(2)
      wppg(2) =  sqrt((upg(2)*sinbeta)**2 + (wpg(2)*cosbeta)**2)
      return
      end

c*********************************************************************

      SUBROUTINE GETUPGVPGWPG(dlsr,vlsr,d,lg,bg,ug,vg,wg,UPG,VPG,WPG)
      IMPLICIT NONE
      REAL*8 lg,bg,sind,cosd
      REAL*4 dlsr,vlsr,d,ug(2),vg(2),wg(2),upg(2),vpg(2),wpg(2)
c----------------------------------------------------------------------------
c                 --- The transformation equations: ---
c
c       upg =  ug*cos(alpha) + (vg+vlsr)*sin(alpha) ;Pi direction
c       vpg = -ug*sin(alpha) + (vg+vlsr)*cos(alpha) ;Theta direction
c       wpg =  wg                                   ;Z direction
c
c     ---alpha is an angle between dlsr (a line between the Galactic center
c        and the LSR (dlsr=8500 pc) and x (a line between the Galactic center
c        and a point on the Galactic disk corresponding to a vertical 
c        projection of a galaxy onto the disk. 
c        
c     --- The circular velocity of the L.S.R. is removed.
c----------------------------------------------------------------------------
      REAL*4 x,cosalpha,sinalpha
c     --- From the law of cosines...
      x=dlsr*dlsr+d*d*cosd(bg)*cosd(bg)-2.0*dlsr*d*cosd(bg)*cosd(lg)
      x=sqrt(x)
c     --- From the law of sines...
      sinalpha=d*cosd(bg)*sind(lg)/x
c     --- From the law of cosines...
      cosalpha=(x*x+dlsr*dlsr-d*d*cosd(bg)*cosd(bg))/(2.0*x*dlsr)
c     --- Apply the transformation.
      upg(1) =  ug(1)*cosalpha + (vg(1)+vlsr)*sinalpha
      vpg(1) = -ug(1)*sinalpha + (vg(1)+vlsr)*cosalpha
      wpg(1) =  wg(1)
      upg(2) =  sqrt((ug(2)*cosalpha)**2 + (vg(2)*sinalpha)**2)
      vpg(2) =  sqrt((ug(2)*sinalpha)**2 + (vg(2)*cosalpha)**2)
      wpg(2) =  wg(2)
      return
      end
c***************************************************************************
      SUBROUTINE GETVNUVTAU(ra,dec,d,ralsr,declsr,vlsr,
     &     rasa,decsa,usun,vsun,wsun,mux,muy,MUXC,MUYC)
      IMPLICIT NONE
      REAL*8 ra,dec,ralsr,declsr,rasa,decsa,sind,cosd,acosd
      REAL*4 d,vlsr,usun,vsun,wsun
      REAL*4 mux(2),muy(2),MUXC(2),MUYC(2)
c     --- lambda is the angle between the apex and the direction to the galaxy;
c         psi is the angle between the great circles connecting the galaxy
c         to the NCP and the apex
      REAL*8 lambda,sinlambda
      REAL*4 sinpsi,cospsi
      REAL*4 muxlsr,muylsr,muxsa,muysa,velsun
c
c     --- All of the following is based on Binney & Tremaine pp 411-413
c
c     --- First find the correction for the LSR motion
      lambda=acosd(sind(dec)*sind(declsr)+
     &                       cosd(dec)*cosd(declsr)*cosd(ra-ralsr))
      sinlambda=sind(lambda)
      if (sinlambda.eq.0.0) then
c        --- the galaxy is at the apex or the antapax and the correction is 0
         muxlsr=0.0
         muylsr=0.0
         MUXC(1)=mux(1)
         MUYC(1)=muy(1)
      else
c        --- find the correction
         cospsi=(cosd(dec)*sind(declsr) -
     &          sind(dec)*cosd(declsr)*cosd(ra-ralsr))/sinlambda
         sinpsi=cosd(declsr)*sind(ra-ralsr)/sinlambda
         muxlsr=vlsr*sinlambda*sinpsi/(4.7410*d)
         muylsr=-vlsr*sinlambda*cospsi/(4.7410*d)
         MUXC(1)=mux(1)-muxlsr
         MUYC(1)=muy(1)-muylsr
      end if
c
c     --- Similarly, find the correction for the solar peculiar motion
      lambda=acosd(sind(dec)*sind(decsa)+
     &                       cosd(dec)*cosd(decsa)*cosd(ra-rasa))
      sinlambda=sind(lambda)
      if (sinlambda.eq.0.0) then
c        --- the galaxy is at the apex or the antapax and the correction is 0
         muxsa=0.0
         muysa=0.0
         MUXC(1)=MUXC(1)
         MUYC(1)=MUYC(1)
      else
c        --- find the correction
         cospsi=(cosd(dec)*sind(decsa) -
     &          sind(dec)*cosd(decsa)*cosd(ra-rasa))/sinlambda
         sinpsi=cosd(decsa)*sind(ra-rasa)/sinlambda
         velsun=sqrt(usun*usun+vsun*vsun+wsun*wsun)
         muxsa=velsun*sinlambda*sinpsi/(4.7410*d)
         muysa=-velsun*sinlambda*cospsi/(4.7410*d)
         MUXC(1)=MUXC(1)-muxsa
         MUYC(1)=MUYC(1)-muysa
      end if
c
c     -- and, finally, set the uncertainties
      MUXC(2)=mux(2)
      MUYC(2)=muy(2)
c
      return
      end

c*******************************************************************
      SUBROUTINE GETUGVGWG(lg,bg,vtl,vtb,vrs,usun,vsun,wsun,
     &                                               UG,VG,WG)
      IMPLICIT NONE
      REAL*8 lg,bg,sind,cosd
      real*4 vtl(2),vtb(2),vrs(2),usun,vsun,wsun
      real*4 ug(2),vg(2),wg(2)
c     --- The subroutine finds the velocity of the galaxy w.r.t. LSR
c     --- vrs    Heliocentric velocity of the galaxy.
c     --- Find the three components of velocity
c                         Theta direction
      vg(1) =  vrs(1)*cosd(bg)*sind(lg) - vtb(1)*sind(bg)*sind(lg) +
     &                   vtl(1)*cosd(lg) + vsun
c                         PI direction
      ug(1) = -vrs(1)*cosd(bg)*cosd(lg) + vtb(1)*sind(bg)*cosd(lg) +
     &                   vtl(1)*sind(lg) + usun
c                         Z direction
      wg(1) =  vrs(1)*sind(bg) + vtb(1)*cosd(bg) + wsun
c
c     --- and the uncertainties of the three components
      vg(2) = sqrt( (vrs(2)*cosd(bg)*sind(lg))**2 +
     &              (vtb(2)*sind(bg)*sind(lg))**2 +
     &              (vtl(2)*cosd(lg))**2 )
      ug(2) = sqrt( (vrs(2)*cosd(bg)*cosd(lg))**2 +
     &              (vtb(2)*sind(bg)*cosd(lg))**2 +
     &              (vtl(2)*sind(lg))**2 )
      wg(2) = sqrt( (vrs(2)*sind(bg))**2 +
     &              (vtb(2)*cosd(bg))**2 )
      return
      end
c*****************************************************************
      SUBROUTINE GETMULMUB(ra,dec,lg,bg,alphagp,deltagp,mux,muy,
     &                                                  MUL,MUB)
      IMPLICIT NONE
      REAL*8 ra,dec,lg,bg,alphagp,deltagp,sind,cosd,acosd
      REAL*4 mux(2),muy(2),mul(2),mub(2)
c     --- The subroutine converts the proper motion from the equatorial
c         to the galactic coordinate system.
c     --- alphagp  RA of the galactic pole
c     --- deltagp  DEC of the galactic pole
c     --- ngpg     Angle between the NGP and the galaxy
c     --- eta      Angle NCP-G-NGP
      REAL*8 eta
      REAL*4 cosdeta,sindeta
c     --- Calculate eta
      if(alphagp.lt.ra) then
         if ((ra-alphagp).lt.180.d0) then
c           --- eta is positive
            eta=acosd((sind(deltagp)-sind(dec)*sind(bg))/
     &                   (cosd(dec)*cosd(bg)))
         else
c           --- eta is negative
            eta=-acosd((sind(deltagp)-sind(dec)*sind(bg))/
     &                   (cosd(dec)*cosd(bg)))
         end if
      else if (alphagp.eq.ra) then
         if (dec.lt.deltagp) then
            eta=0.0d0
         else
            eta=180.0d0
         end if
      else
         if ((alphagp-ra).gt.180.) then
c           --- eta is positive
            eta=acosd((sind(deltagp)-sind(dec)*sind(bg))/
     &                   (cosd(dec)*cosd(bg)))
         else
c           --- eta is negative
            eta=-acosd((sind(deltagp)-sind(dec)*sind(bg))/
     &                   (cosd(dec)*cosd(bg)))
         end if
      end if
c
c     write(6,*)'The angle eta (NCP-G-NGP):',eta,' degrees'
c     --- Calculate l and b components of the proper motion.
c                         *** Note ***
c         The transformation equations from Spherical Astronomy
c         by W. M. Smart (sixth ed.), page 276 are:
c
c         mul=  mux*cosd(eta)*cosd(dec)/cosd(bg) + muy*sind(eta)/cosd(bg)
c         mub= -mux*sind(eta)*cosd(dec)+muy*cosd(eta)
c
c         However, since our mux and muy already have the cos(dec) factor
c         included, the terms cosd(dec) and cosd(bg) should be then dropped.
c                             ***
      cosdeta=cosd(eta)
      sindeta=sind(eta)
      mul(1) =  mux(1)*cosdeta+muy(1)*sindeta
      mub(1) = -mux(1)*sindeta+muy(1)*cosdeta
      mul(2) = sqrt((mux(2)*cosdeta)**2+(muy(2)*sindeta)**2)
      mub(2) = sqrt((mux(2)*sindeta)**2+(muy(2)*cosdeta)**2)
c     --- Check if the magnitude of the proper motion in both coordinates
c         agree. If not, warn.
      if(abs((mul(1)*mul(1)+mub(1)*mub(1))-
     &                (mux(1)*mux(1)+muy(1)*muy(1))).gt.1.0e-4) then
         write(6,*) 
     &   'WARNING: the magnitude of the proper motion not preserved'
         write(6,*) mul(1)*mul(1)+mub(1)*mub(1),
     &                           mux(1)*mux(1)+muy(1)*muy(1)
      end if
      return
      end
      
c*****************************************************************

      SUBROUTINE GETVT(mux,muy,d,asecrad,secyear,kmpc,vtx,vty,vt)
      IMPLICIT NONE
      REAL*4 mux(2),muy(2),d,asecrad,secyear,kmpc
      REAL*4 vtx(2),vty(2),vt(2)
      vtx(1)=mux(1)*d*kmpc/(secyear*asecrad)
      vty(1)=muy(1)*d*kmpc/(secyear*asecrad)
      vt(1)=vtx(1)*vtx(1)+vty(1)*vty(1)
      vt(1)=sqrt(vt(1))
      vtx(2)=mux(2)*d*kmpc/(secyear*asecrad)
      vty(2)=muy(2)*d*kmpc/(secyear*asecrad)
      vt(2)=sqrt(vtx(2)*vtx(2)+vty(2)*vty(2))
      return
      end

c*****************************************************************

      SUBROUTINE ARCSECPERYEAR(pscale,muxpy,muypy,mux,muy)
      IMPLICIT NONE
      REAL*4 pscale,muxpy(2),muypy(2),mux(2),muy(2)
      mux(1)=pscale*muxpy(1)
      muy(1)=pscale*muypy(1)
      mux(2)=pscale*muxpy(2)
      muy(2)=pscale*muypy(2)
      return
      end

c****************************************************************

      SUBROUTINE PROJECT(pa,muxpy,muypy)
      IMPLICIT NONE
      REAL*8 cosd,sind
      REAL*4 pa,px(2),py(2),muxpy(2),muypy(2)
      REAL*4 cosdpa,sindpa
c     --- Given the position angle, project muxpy and muypy to the 
c         equatorial coordinate system. The positive position angle is 
c         measured from North in the Eastward direction.
      px(1)=muxpy(1)
      px(2)=muxpy(2)
      py(1)=muypy(1)
      py(2)=muypy(2)
      cosdpa=cosd(dble(pa))
      sindpa=sind(dble(pa))
      muxpy(1) = -px(1)*cosdpa+py(1)*sindpa
      muypy(1) =  px(1)*sindpa+py(1)*cosdpa
      muxpy(2) = sqrt((cosdpa*px(2))**2+(sindpa*py(2))**2)
      muypy(2) = sqrt((sindpa*px(2))**2+(cosdpa*py(2))**2)
      return
      end

c*****************************************************************

      SUBROUTINE FINDLGBG(ra,dec,alphagp,deltagp,gcgan,lg,bg)
      IMPLICIT NONE
      REAL*8 ra,dec,alphagp,deltagp,gcgan,lg,bg
c     --- The subroutine finds the galactic coordinates of an object
c         given its RA and DEC.
      REAL*8 x,y,sind,cosd,atan2d,asind
      bg=asind(sind(dec)*cosd(90.0d0-deltagp)-cosd(dec)*
     &     sind(ra-(alphagp+90.0d0))*sind(90.0d0-deltagp))
c
      y=cosd(dec)*sind(ra-(alphagp+90.0d0))*cosd(90.0d0-deltagp)+
     &                           sind(dec)*sind(90.0d0-deltagp)
      x=cosd(dec)*cosd(ra-(alphagp+90.0d0))
      lg=atan2d(y,x)+gcgan
      if(lg.lt.0.0d0) lg=lg+360.0d0
      return
      end

c     ***************************************************
      REAL*8 FUNCTION sind(x)
      IMPLICIT NONE
      REAL*8 sind,x,d2r
      parameter (d2r=1.74532925199d-02)
c
      sind=sin(d2r*x)
c
      return
      end

c     ***************************************************
      REAL*8 FUNCTION cosd(x)
      IMPLICIT NONE
      REAL*8 cosd,x,d2r
      parameter (d2r=1.74532925199d-02)
c
      cosd=cos(d2r*x)
c
      return
      end

c     ***************************************************
      REAL*8 FUNCTION asind(x)
      IMPLICIT NONE
      REAL*8 asind,x,r2d
      parameter (r2d=57.2957795132d0)
c
      asind=r2d*asin(x)
c
      return
      end

c     ***************************************************
      REAL*8 FUNCTION acosd(x)
      IMPLICIT NONE
      REAL*8 acosd,x,r2d
      parameter (r2d=57.2957795132d0)
c
      acosd=r2d*acos(x)
c
      return
      end

c     ***************************************************
      REAL*8 FUNCTION atan2d(y,x)
      IMPLICIT NONE
      REAL*8 atan2d,x,y,r2d
      parameter (r2d=57.2957795132d0)
c
      atan2d=r2d*atan2(y,x)
c
      return
      end
