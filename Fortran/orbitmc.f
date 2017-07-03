c     --- Started writing on October 3, 2002.
c         ORBIT.F integrates an orbit in a given Galactic potential.
c     --- This version uses a file of 1000 input v_pi, v_theta, v_z generated
c         by pm_pmq_mcorb.
c     --- 13-Jun-2016: Implemented the MWPotential2014 from Bovy 2015,
c            ApJS 216:29.  Added an energy calculation for the N potential.
c            Using vcon that is the inverse of the one before, which I think
c            is correct.
c     --- 15-June-2016:  Added a potential from Kafle et al. 2014, ApJ, 794:59,
c            which is similar in form to MWPot2014 but has V_c(R_0) of about
c            240 kpc and an M_vir that is slightly above Bovy's --- though on
c            the low end of the spectrum of derived values.  Also added two
c            potentials (V and X) that have an NFW halo consistent with the
c            one adopted by van der Marel et al. 2012 ApJ 753:8 on the basis
c            of LG timing.  Actual values come from the Bland-Hawthorne and
c            Gerhard 2016 ARAA MW review.  V uses a Hernquist model for the
c            halo (which has finite mass) as described in the van der Marel
c            paper (appendix), while X uses an NFW halo.
c    --- 17-Jun-2016:  Changed dlsr from 8.0 to 8.2 kpc (Bland-Hawthorn and
c            Gerhard 2016).
c
c     --- a      Galactic potential parameter See KJ et al. 2002,AJ,124,127
c     --- b      Galactic potential parameter
c     --- c      Galactic potential parameter
c     --- d      Galactic potential parameter
c     --- g      Gravitational constant in kpc^3/(M_sun Gyr)
c     --- md     Mass of the disk
c     --- ms     Mass of the spherical component
c     --- vh     V-Helo parameter
c     --- fx     X component of the force
c     --- fy     Y component of the force
c     --- fz     Z component of the force
c     --- bg     Galactic latidude of the galaxy.
c     --- dg     Heliocentric distance to the galaxy.
c     --- dlsr   Galactocentric distance to LSR in kpc.
c     --- lg     Galactic longitude of the galaxy.
c     --- x0     Initial X-coordinate of the galaxy.
c     --- y0     Initial Y-coordinate of the galaxy.
c     --- z0     Initial Y-coordinate of the galaxy.
c     --- x      X position of the galaxy at time t
c     --- y      Y position of the galaxy at time t
c     --- z      Z position of the galaxy at time t
c     --- vx0    Initial velocity in X direction of the galaxy, in km/s
c     --- vx0u   Uncertainty in vx0, in km/s
c     --- vy0    Initial velocity in Y direction of the galaxy, in km/s
c     --- vy0u   Uncertainty in vy0, in km/s
c     --- vz0    Initial velocity in Z direction of the galaxy, in km/s
c     --- vz0u   Uncertainty in vz0, in km/s
c     --- vx     Velocity x at time t
c     --- vy     Velocity x at time t 
c     --- vz     Velocity x at time t
c     --- vcon   Conversion factor from kpc/Gyr to km/s
c     --- ra     Apogalacticon
c     --- rp     Perigalacticon
c     --- e      Eccentricity
c     --- theta  Inclination angle
      IMPLICIT NONE
      REAL*8 bg,lg,dg,x0,y0,z0
      REAL*8 vx0,vx0u,vy0,vy0u,vz0,vz0u
      REAL*8 dlsr,vcon
      REAL*8 ratoam,pi,d2r
      INTEGER NN,NS,NC
      PARAMETER (NN=1000,NS=500000,NC=100)
ccc      PARAMETER (dlsr=8.5d0)
      PARAMETER (dlsr=8.2d0)
ccc      PARAMETER (vcon=0.977786932d0)
      PARAMETER (vcon=1.0226903)
      PARAMETER (ratoam=3437.747,PI=3.14159)
      parameter (d2r=1.74532925199D-02)
      INTEGER NAP,idap(NC),NPE,idpe(NC)
      REAL*8 vpi,vpiu,vth,vthu
      REAL*8 fx,fy,fz,x,y,z,vx,vy,vz
      REAL*8 t,dt,ts,rmin,rmax
      REAL*8 rp(NN),ra(NN),e(NN),rt1(NN),rt2(NN),porb(NN),ef(NN),e0(NN)
      REAL*8 theta(NN),omega(NN),ap,bp,cp,dp,lpole(NN),bpole(NN)
      REAL*8 dumb,ave,sd
      REAL*8 mol1,mol2,ldsph,mmw,rap,rpe,ec,fe
      REAL*8 p(3,3)
      REAL*8 time(NS),r2(NS)
      REAL*8 xt(NS),yt(NS),zt(NS)
      REAL*8 cl95(2),rtlim,rpval
      REAL*8 dlam,dbam,lam(NN),bam(NN),lamf
      REAL*4 amx0,amy0,amz0,am0,amx,amy,amz,am
      REAL*4 hrpmax,hramax,hemax,htmax,hphimax,homegamax
      CHARACTER pot*1,fname*32
      INTEGER NSTEP,NMC,i,j,k,NCOUNT,NRP
ccc      INTEGER iseed,iranf
ccc 100  format(1x,f6.3,1x,f7.2,1x,f7.2,1x,f7.2,1x,f7.2,1x,f7.2,1x,f7.2)
      OPEN(1,FILE='orbitmc.par',STATUS='OLD',form='FORMATTED')
      read(1,'(a32)') fname  !name of ivel.data file
      read(1,'(a1)') pot    !Either D,J,N,B,K,V, or X.
      read(1,*) lg ! Read in Galactic longitude [degrees].
      read(1,*) bg ! Read in Galactic latitude [degrees].
      read(1,*) dg ! Read in the heliocentric distance to the galaxy [kpc].
      write (6,*) 'l,b,distance: ',lg,bg,dg
c
      read(1,*) t      ! Total integration time in Gyr
      read(1,*) NSTEP  ! Number of integration steps
      if(NSTEP.gt.NS) then
         write(6,*) 'ERROR 1 in MAIN: Too many steps, increase NS'
         stop
      end if
      dt=t/dble(NSTEP) ! Time-step size
      read(1,*) NMC    ! Read in the number of Monte-Carlo experiments.
      if(NMC.gt.NN) then
         write(6,*) 'Error 1 in the MAIN'
         stop
      end if
c     --- luminosity of the dSph and the min and max M/L
      read(1,*) ldsph
      read(1,*) mol1,mol2
c     --- the limiting r_t; prints the fraction of orbits with smaller values
      read(1,*) rtlim
c     --- a value for perigalacticon, records the fraction of orbits with rp < rpval
      read(1,*) rpval
c     --- read maxima for histogram plots
      read(1,*) hrpmax
      read(1,*) hramax
      read(1,*) hemax
      read(1,*) htmax
      read(1,*) hphimax
      read(1,*) homegamax
      k=1
      NCOUNT=0
      NRP=0
      dlam=0.0d0
      dbam=0.0d0
ccc      open (20, file='ivel.data', form='formatted',status='old')
      write (6,*) 'opening starting velocities file: ',fname
      open (20, file=fname, form='formatted',status='old')
      do j=1,NMC
         rmin=3.2d32
         rmax=0.0d0
c        --- Get velocities from the file and xform to vx, vy, vz
         read (20,*) vpi,vpiu,vth,vthu,vz0,vz0u
         CALL GETINITPOS(dlsr,lg,bg,dg,x0,y0,z0,
     &                         vpi,vpiu,vth,vthu,vx0,vx0u,vy0,vy0u)
c
ccc         write (6,*) j,vpi,vth,vz0
ccc         CALL ENERGY(a,b,c,d,md,ms,vh,g,vcon,x0,y0,z0,vx0,vy0,vz0,e0(k))
c        -- find the initial energy
         if(pot.eq.'J') then
            call JENERGY(vcon,x0,y0,z0,vx0,vy0,vz0,e0(k))
         else if(pot.eq.'D') then
            CALL DENERGY(vcon,x0,y0,z0,vx0,vy0,vz0,e0(k))
         else if(pot.eq.'B') then
            CALL BENERGY(vcon,x0,y0,z0,vx0,vy0,vz0,e0(k))
         else if(pot.eq.'N') then
            CALL NENERGY(vcon,x0,y0,z0,vx0,vy0,vz0,e0(k))
         else if(pot.eq.'K') then
            CALL KENERGY(vcon,x0,y0,z0,vx0,vy0,vz0,e0(k))
         else if(pot.eq.'V') then
            CALL VENERGY(vcon,x0,y0,z0,vx0,vy0,vz0,e0(k))
         else if(pot.eq.'X') then
            CALL XENERGY(vcon,x0,y0,z0,vx0,vy0,vz0,e0(k))
         end if
c
c        --- Calculate the initial angular momentum vector = r x v
c            and the (l,b) that it points
         amx0=y0*vz0-z0*vy0
         amy0=z0*vx0-x0*vz0
         amz0=x0*vy0-y0*vx0
         am0=sqrt(amx0*amx0+amy0*amy0+amz0*amz0)
         bam(k)=asin(amz0/am0)
         lam(k)=atan2(-amy0,-amx0)
         if(lam(k).lt.0.0d0) lam(k)=lam(k)+2.0d0*pi
         if (j.eq.1) then
            write (6,*) 'x0,y0,z0: ',x0,y0,z0
            write (6,*) 'vx0,vy0,vz0: ',vx0,vy0,vz0
            write (6,*) 'amx0,amy0,amz0: ',amx0,amy0,amz0
            write (6,*) 'lam0, bam0: ',lam(k)/d2r,bam(k)/d2r
            write (6,*) ''
ccc        else
ccc            write (6,*) k,lam(i)/d2r,bam(k)/d2r
        end if
c
c        --- integrate the orbit
         x=x0
         y=y0
         z=z0
         vx=vx0
         vy=vy0
         vz=vz0
         ts=0.0d0
         time(1)=ts
         r2(1)=x*x+y*y+z*z
         xt(1)=x
         yt(1)=y
         zt(1)=z
         if(pot.eq.'K') then
            CALL KFORCE(x,y,z,fx,fy,fz)
         else if(pot.eq.'X') then
            CALL XFORCE(x,y,z,fx,fy,fz)
         else if(pot.eq.'V') then
            CALL VFORCE(x,y,z,fx,fy,fz)
         else if(pot.eq.'J') then
            CALL JFORCE(x,y,z,fx,fy,fz)
         else if(pot.eq.'B') then
            CALL BFORCE(x,y,z,fx,fy,fz)
         else if(pot.eq.'D') then
            CALL DFORCE(x,y,z,fx,fy,fz)
         else if(pot.eq.'N') then
            CALL NFORCE(x,y,z,fx,fy,fz)
         end if
         do i=1, NSTEP-1
            vx=vx+0.5d0*dt*fx/vcon
            vy=vy+0.5d0*dt*fy/vcon
            vz=vz+0.5d0*dt*fz/vcon
            x=x+dt*vx*vcon
            y=y+dt*vy*vcon
            z=z+dt*vz*vcon
            ts=ts+dt
            if(pot.eq.'K') then
               CALL KFORCE(x,y,z,fx,fy,fz)
            else if(pot.eq.'X') then
               CALL XFORCE(x,y,z,fx,fy,fz)
            else if(pot.eq.'V') then
               CALL VFORCE(x,y,z,fx,fy,fz)
            else if(pot.eq.'J') then
               CALL JFORCE(x,y,z,fx,fy,fz)
            else if(pot.eq.'B') then
               CALL BFORCE(x,y,z,fx,fy,fz)
            else if(pot.eq.'D') then
               CALL DFORCE(x,y,z,fx,fy,fz)
            else if(pot.eq.'N') then
               CALL NFORCE(x,y,z,fx,fy,fz)
            end if
            vx=vx+0.5d0*dt*fx/vcon
            vy=vy+0.5d0*dt*fy/vcon
            vz=vz+0.5d0*dt*fz/vcon
            time(i+1)=ts
            r2(i+1)=x*x+y*y+z*z
            xt(i+1)=x
            yt(i+1)=y
            zt(i+1)=z
ccc            write (6,*) i,x,y,z,vx,vy,vz
         end do
c        -- find the final energy (the B energy is approximate at small radii)
         if(pot.eq.'J') then
            call JENERGY(vcon,x,y,z,vx,vy,vz,ef(k))
         else if(pot.eq.'B') then
            CALL BENERGY(vcon,x,y,z,vx,vy,vz,ef(k))
         else if(pot.eq.'N') then
            CALL NENERGY(vcon,x,y,z,vx,vy,vz,ef(k))
         else if(pot.eq.'D') then
            CALL DENERGY(vcon,x,y,z,vx,vy,vz,ef(k))
         else if(pot.eq.'K') then
            CALL KENERGY(vcon,x,y,z,vx,vy,vz,ef(k))
         else if(pot.eq.'V') then
            CALL VENERGY(vcon,x,y,z,vx,vy,vz,ef(k))
         else if(pot.eq.'X') then
            CALL XENERGY(vcon,x,y,z,vx,vy,vz,ef(k))
         end if
ccc         CALL ENERGY(a,b,c,d,md,ms,vh,g,vcon,x,y,z,vx,vy,vz,ef(k))
ccc         write (6,*) 'k, e0, (e0-ef)/e0: ',k,e0(k),(e0(k)-ef(k))/e0(k)
c        --- Find the final r x v and compare
         amx=y*vz-z*vy
         amy=z*vx-x*vz
         amz=x*vy-y*vx
         am=sqrt(amx*amx+amy*amy+amz*amz)
         dbam=dbam + (dble(asin(amz/am))-bam(k))
         lamf=atan2(-amy,-amx)
         if(lamf.lt.0.0d0) lamf=lamf+2.0d0*pi
c        --- check to make sure have not crossed 0/2*pi
         if (abs(lamf-lam(k)).gt.pi) then
            if (lam(k).lt.pi/2.0d0) then
               lamf=lamf-2.0d0*pi
            else
               lamf=lamf+2.0d0*pi
            end if
         end if
         dlam=dlam + (lamf-lam(k))
c
         CALL FINDIDORB(NC,NS,NSTEP,r2,nap,idap,npe,idpe)
         if(NAP.lt.2) go to 666
c        --- Find average R_p and store it.
         rp(k)=0.0d0
         do i=1,NPE
            rp(k)=rp(k)+sqrt(r2(idpe(i)))
         end do
         rp(k)=rp(k)/real(NPE)
         if (rp(k).lt.rpval) then
            NRP=NRP+1
         end if
c        --- Find average R_a and store it.
         ra(k)=0.0d0
         do i=1,NAP
            ra(k)=ra(k)+sqrt(r2(idap(i)))
         end do
         ra(k)=ra(k)/real(NAP)
c        --- Find average period from R_a to R_a and store it.
         porb(k)=0.0d0
         do i=2,NAP
            porb(k)=porb(k)+time(idap(i))-time(idap(i-1))
         end do
         porb(k)=porb(k)/real(NAP-1)

c        --- Find average eccentricity and store it.
         e(k)=0.0d0
         do i=1,min(NPE,NAP)
            e(k)=e(k)+(sqrt(r2(idap(i)))-sqrt(r2(idpe(i))))/
     &                (sqrt(r2(idap(i)))+sqrt(r2(idpe(i))))
         end do
         e(k)=e(k)/real(min(NPE,NAP))
c        --- Calculate the tidal radius at pericenter and store it
         do i=1,min(NPE,NAP)
            rap=sqrt(r2(idap(i)))
            rpe=sqrt(r2(idpe(i)))
            ec=(rap-rpe)/(rap+rpe)
            fe=((1.00-ec)**2)/
     &         ((((1.00+ec)**2)/(2.0*ec))*dlog((1.0+ec)/(1.0-ec))+1.00)
            mmw=1.1e10*(rap+rpe)/2.0
            rt1(k)=rt1(k)+((rap+rpe)/2.0)*(fe*mol1*ldsph/mmw)**0.333333
            rt2(k)=rt2(k)+((rap+rpe)/2.0)*(fe*mol2*ldsph/mmw)**0.333333
         end do
         rt1(k)=rt1(k)/real(min(NPE,NAP))
         rt2(k)=rt2(k)/real(min(NPE,NAP))
         rt1(k)=atan(rt1(k)/dg)*ratoam
         rt2(k)=atan(rt2(k)/dg)*ratoam
         theta(k)=0.0d0
         omega(k)=0.0d0
         do i=1,min(NPE,NAP)
            CALL PICK3POINTS(NC,NS,idpe,idap,i,xt,yt,zt,p)
            CALL FINDPLANE(p,ap,bp,cp,dp)
            CALL FINDINC(ap,bp,cp,dumb)
            theta(k)=theta(k)+dumb
            CALL FINDLONGITUDE(ap,bp,dumb)
            omega(k)=omega(k)+dumb
         end do
         theta(k)=theta(k)/real(min(NPE,NAP))
         omega(k)=omega(k)/real(min(NPE,NAP))
         bpole(k)=theta(k)-(pi/2.0d0)
         lpole(k)=omega(k)+0.50d0*pi
         if (lpole(k).gt.2.0d0*pi) then
            lpole(k)=lpole(k)-2.0d0*pi
         end if
c        -- do some checking to avoid splitting omega values across 0/2pi
         if (abs(omega(k)-omega(1)).gt.pi) then
            if (omega(1).lt.pi) then
               omega(k)=omega(k)-2.0d0*pi
            else
               omega(i)=omega(k)+2.0d0*pi
            end if
         end if
c        -- do some checking to avoid splitting lpole values across 0/2pi
         if (abs(lpole(k)-lpole(1)).gt.pi) then
            if (lpole(1).lt.pi) then
               lpole(k)=lpole(k)-2.0d0*pi
            else
               lpole(k)=lpole(k)+2.0d0*pi
            end if
         end if
c        -- and the same for lam
         if (abs(lam(k)-lam(1)).gt.pi) then
            if (lam(1).lt.pi) then
               lam(k)=lam(k)-2.0d0*pi
            else
               lam(k)=lam(k)+2.0d0*pi
            end if
         end if
         k=k+1
         NCOUNT=NCOUNT+1
 666     continue
      end do
      close(unit=20)
ccc      CALL STATISTICS(porb,NN,NCOUNT,porbmean,sdporb)
ccc      CALL MIN1D(porb,NN,NCOUNT,porbmn)
ccc      CALL MAX1D(porb,NN,NCOUNT,porbmx)
ccc      CALL FIND95(porb,NN,NCOUNT,porb95)
      close(unit=1)

 2004 format(1x,'Actual Perigalacticon, in kpc : ',f8.3)
      write (6,2004) rp(1)
 200  format(1x,'Mean Perigalacticon, in kpc   : ',f8.3,1x,'+-',1x,f5.2)
      CALL STATISTICS(rp,NN,NCOUNT,ave,sd)
      write(6,200) ave,sd
 2001 format(1x,'Max. Perigalacticon, in kpc   : ',f8.3)
      CALL MAX1D(rp,NN,NCOUNT,dumb)
      write(6,2001) dumb
 2002 format(1x,'Min. Perigalacticon, in kpc   : ',f8.3)
      CALL MIN1D(rp,NN,NCOUNT,dumb)
      write(6,2002) dumb
 2003 format(1x,'95% Confidence limits, in kpc : ',f8.3,1x,f8.3)
      CALL FIND95(rp,NN,NCOUNT,cl95)
      write(6,2003) cl95(1),cl95(2)
 2005 format(1x,'fraction of orbits with Peri < ',f8.3,' is ',f6.4)
      write (6,2005) rpval,float(NRP)/float(NCOUNT)
      write(6,*)''

 2014 format(1x,'Actual Apogalacticon,  in kpc : ',f8.3)
      write (6,2014) ra(1)
 201  format(1x,'Mean Apogalacticon,  in kpc   : ',f8.3,1x,'+-',1x,f5.2)
      CALL STATISTICS(ra,NN,NCOUNT,ave,sd)
      write(6,201) ave,sd
 2011 format(1x,'Max. Apogalacticon,  in kpc   : ',f8.3)
      CALL MAX1D(ra,NN,NCOUNT,dumb)
      write(6,2011) dumb
 2012 format(1x,'Min. Apogalacticon,  in kpc   : ',f8.3)
      CALL MIN1D(ra,NN,NCOUNT,dumb)
      write(6,2012) dumb
 2013 format(1x,'95% Confidence limits, in kpc : ',f8.3,1x,f8.3)
      CALL FIND95(ra,NN,NCOUNT,cl95)
      write(6,2013) cl95(1),cl95(2)
      write(6,*)''

 2024 format(1x,'Actual eccentricity           : ',f8.3)
      write (6,2024) e(1)
 202  format(1x,'Mean eccentricity             : ',f8.3,1x,'+-',1x,f5.2)
      CALL STATISTICS(e,NN,NCOUNT,ave,sd)
      write(6,202) ave,sd
 2021 format(1x,'Max. eccentricity             : ',f8.3)
      CALL MAX1D(e,NN,NCOUNT,dumb)
      write(6,2021) dumb
 2022 format(1x,'Min. eccentricity             : ',f8.3)
      CALL MIN1D(e,NN,NCOUNT,dumb)
      write(6,2022) dumb
 2023 format(1x,'95% Confidence limits         : ',f8.3,1x,f8.3)
      CALL FIND95(e,NN,NCOUNT,cl95)
      write(6,2023) cl95(1),cl95(2)
      write(6,*)''

 2044 format(1x,'Actual orbital period, in Gyr : ',f8.3)
      write (6,2044) porb(1)
 204  format(1x,'Mean orbital period, in Gyr   : ',f8.3,1x,'+-',1x,f5.2)
      CALL STATISTICS(porb,NN,NCOUNT,ave,sd)
      write(6,204) ave,sd
 2041 format(1x,'Max. orbital period, in Gyr   : ',f8.3)
      CALL MAX1D(porb,NN,NCOUNT,dumb)
      write(6,2041) dumb
 2042 format(1x,'Min. orbital period, in Gyr   : ',f8.3)
      CALL MIN1D(porb,NN,NCOUNT,dumb)
      write(6,2042) dumb
 2043 format(1x,'95% Confidence limits, in Gyr : ',f8.3,1x,f8.3)
      CALL FIND95(porb,NN,NCOUNT,cl95)
      write(6,2043) cl95(1),cl95(2)
      write(6,*)''

 2055 format(1x,'Actual r_t1 [arcmin]          : ',f8.3)
      write (6,2055) rt1(1)
 205  format(1x,'Mean r_t1 [arcmin]            : ',f8.3,1x,'+-',1x,f5.2)
      CALL STATISTICS(rt1,NN,NCOUNT,ave,sd)
      write(6,205) ave,sd
 2051 format(1x,'Max. r_t1 [arcmin]            : ',f8.3)
      CALL MAX1D(rt1,NN,NCOUNT,dumb)
      write(6,2051) dumb
 2052 format(1x,'Min. r_t1 [arcmin]            : ',f8.3)
      CALL MIN1D(rt1,NN,NCOUNT,dumb)
      write(6,2052) dumb
 2053 format(1x,'95% Confidence limits         : ',f8.3,1x,f8.3)
      CALL FIND95(rt1,NN,NCOUNT,cl95)
      write(6,2053) cl95(1),cl95(2)
      CALL FINDSMALLER(rt1,rtlim,NN,NCOUNT,dumb)
      write (6,2054) rtlim,dumb
 2054 format(1x,'fraction smaller than ',f5.2,'   :    ',f6.4)
      write(6,*)''

 2064 format(1x,'Actual r_t2 [arcmin]          : ',f8.3)
      write (6,2064) rt2(1)
 206  format(1x,'Mean r_t2 [arcmin]            : ',f8.3,1x,'+-',1x,f5.2)
      CALL STATISTICS(rt2,NN,NCOUNT,ave,sd)
      write(6,206) ave,sd
 2061 format(1x,'Max. r_t2 [arcmin]            : ',f8.3)
      CALL MAX1D(rt2,NN,NCOUNT,dumb)
      write(6,2061) dumb
 2062 format(1x,'Min. r_t2 [arcmin]            : ',f8.3)
      CALL MIN1D(rt2,NN,NCOUNT,dumb)
      write(6,2062) dumb
 2063 format(1x,'95% Confidence limits         : ',f8.3,1x,f8.3)
      CALL FIND95(rt2,NN,NCOUNT,cl95)
      write(6,2063) cl95(1),cl95(2)
      CALL FINDSMALLER(rt2,rtlim,NN,NCOUNT,dumb)
      write (6,2054) rtlim,dumb
      write(6,*)''

 2074 format(1x,'Actual Inclination [degrees]  : ',f8.3)
      write (6,2074) theta(1)/d2r
 207  format(1x,'Mean Inclination [degrees]    : ',f8.3,1x,'+-',1x,f5.2)
      CALL STATISTICS(theta,NN,NCOUNT,ave,sd)
      write(6,207) ave/d2r,sd/d2r
 2071 format(1x,'Max. Inclination [degrees]    : ',f8.3)
      CALL MAX1D(theta,NN,NCOUNT,dumb)
      write(6,2071) dumb/d2r
 2072 format(1x,'Min. Inclination [degrees]    : ',f8.3)
      CALL MIN1D(theta,NN,NCOUNT,dumb)
      write(6,2072) dumb/d2r
 2073 format(1x,'95% Confidence limits         : ',f8.3,1x,f8.3)
      CALL FIND95(theta,NN,NCOUNT,cl95)
      write(6,2073) cl95(1)/d2r,cl95(2)/d2r
      write(6,*)''

 2084 format(1x,'Actual Longitude [degrees]    : ',f8.3)
      write (6,2084) omega(1)/d2r
 208  format(1x,'Mean Longitude [degrees]      : ',f8.3,1x,'+-',1x,f7.3)
      CALL STATISTICS(omega,NN,NCOUNT,ave,sd)
      write(6,208) ave/d2r,sd/d2r
 2081 format(1x,'Max. Longitude [degrees]      : ',f8.3)
      CALL MAX1D(omega,NN,NCOUNT,dumb)
      write(6,2081) dumb/d2r
 2082 format(1x,'Min. Longitude [degrees]      : ',f8.3)
      CALL MIN1D(omega,NN,NCOUNT,dumb)
      write(6,2082) dumb/d2r
 2083 format(1x,'95% Confidence limits         : ',f8.3,1x,f8.3)
      CALL FIND95(omega,NN,NCOUNT,cl95)
      write(6,2083) cl95(1)/d2r,cl95(2)/d2r
      write(6,*)''

 2094 format(1x,'Actual lpole [degrees]        : ',f8.3)
      write (6,2094) lpole(1)/d2r
 209  format(1x,'Mean lpole [degrees]          : ',f8.3,1x,'+-',1x,f5.2)
      CALL STATISTICS(lpole,NN,NCOUNT,ave,sd)
      write(6,209) ave/d2r,sd/d2r
 2091 format(1x,'Max. lpole [degrees]          : ',f8.3)
      CALL MAX1D(lpole,NN,NCOUNT,dumb)
      write(6,2091) dumb/d2r
 2092 format(1x,'Min. lpole [degrees]          : ',f8.3)
      CALL MIN1D(lpole,NN,NCOUNT,dumb)
      write(6,2092) dumb/d2r
 2093 format(1x,'95% Confidence limits         : ',f8.3,1x,f8.3)
      CALL FIND95(lpole,NN,NCOUNT,cl95)
      write(6,2093) cl95(1)/d2r,cl95(2)/d2r
      write(6,*)''

 2104 format(1x,'Actual bpole [degrees]        : ',f8.3)
      write (6,2104) bpole(1)/d2r
 210  format(1x,'Mean bpole [degrees]          : ',f8.3,1x,'+-',1x,f7.3)
      CALL STATISTICS(bpole,NN,NCOUNT,ave,sd)
      write(6,210) ave/d2r,sd/d2r
 2101 format(1x,'Max. bpole [degrees]          : ',f8.3)
      CALL MAX1D(bpole,NN,NCOUNT,dumb)
      write(6,2101) dumb/d2r
 2102 format(1x,'Min. bpole [degrees]          : ',f8.3)
      CALL MIN1D(bpole,NN,NCOUNT,dumb)
      write(6,2102) dumb/d2r
 2103 format(1x,'95% Confidence limits         : ',f8.3,1x,f8.3)
      CALL FIND95(bpole,NN,NCOUNT,cl95)
      write(6,2103) cl95(1)/d2r,cl95(2)/d2r
      write(6,*)''

 2114 format(1x,'Actual lam [degrees]          : ',f8.3)
      write (6,2114) lam(1)/d2r
 211  format(1x,'Mean lam [degrees]            : ',f8.3,1x,'+-',1x,f5.2)
      CALL STATISTICS(lam,NN,NCOUNT,ave,sd)
      write(6,211) ave/d2r,sd/d2r
 2111 format(1x,'Max. lam [degrees]            : ',f8.3)
      CALL MAX1D(lam,NN,NCOUNT,dumb)
      write(6,2111) dumb/d2r
 2112 format(1x,'Min. lam [degrees]            : ',f8.3)
      CALL MIN1D(lam,NN,NCOUNT,dumb)
      write(6,2112) dumb/d2r
 2113 format(1x,'95% Confidence limits         : ',f8.3,1x,f8.3)
      CALL FIND95(lam,NN,NCOUNT,cl95)
      write(6,2113) cl95(1)/d2r,cl95(2)/d2r
      write(6,*)''

 2124 format(1x,'Actual bam [degrees]          : ',f8.3)
      write (6,2124) bam(1)/d2r
 212  format(1x,'Mean bam [degrees]            : ',f8.3,1x,'+-',1x,f7.3)
      CALL STATISTICS(bam,NN,NCOUNT,ave,sd)
      write(6,212) ave/d2r,sd/d2r
 2121 format(1x,'Max. bam [degrees]            : ',f8.3)
      CALL MAX1D(bam,NN,NCOUNT,dumb)
      write(6,2121) dumb/d2r
 2122 format(1x,'Min. bam [degrees]            : ',f8.3)
      CALL MIN1D(bam,NN,NCOUNT,dumb)
      write(6,2122) dumb/d2r
 2123 format(1x,'95% Confidence limits         : ',f8.3,1x,f8.3)
      CALL FIND95(bam,NN,NCOUNT,cl95)
      write(6,2123) cl95(1)/d2r,cl95(2)/d2r
      write(6,*)''

      write (6,213) dlam/(d2r*float(nmc)),dbam/(d2r*float(nmc))
 213  format(1x,'average change lam,bam [degrees] :',2(1x,f7.4))
      write (6,*) ''

 305  format(1x,'Number of experiments            : ',i5)
      write(6,305) nmc
 306  format(1x,'Number of sucessful experiments  : ',i5)
      write(6,306) ncount
      call pgbegin(0,'orbitmc.ps/ps',1,1)
      call pgsch(1.0)
      call pgslw(2)
      call pgscf(2)
      CALL PLOTENERGY(NN,NCOUNT,e0,ef)
      call pgpage
      CALL PLOTENERGYPERI(NN,NCOUNT,e0,ef,rp)
      call pgpage
      CALL PLOTENERGYPERIOD(NN,NCOUNT,e0,ef,porb)
      call pgpage
      CALL PLOTENERGYECC(NN,NCOUNT,e0,ef,e)
      call pgpage
      CALL PLOTHIST(NN,NCOUNT,rp,ra,e,porb,theta,omega,
     &                 hrpmax,hramax,hemax,htmax,hphimax,homegamax)
      call pgend
      stop
      end

      SUBROUTINE FINDLONGITUDE(a,b,omega)
      IMPLICIT NONE
      REAL*8 a,b,omega,pi
      parameter (pi=3.14159265359d00)
c     -- longitude of the ascending node measured in the direction
c        opposite galactic rotation from x=0
      omega=atan2(b,a)+pi/2.0d0
      if(omega.lt.0.0d0) omega=omega+2.0d0*pi
      return
      end
      
      SUBROUTINE FINDINC(a,b,c,theta)
      IMPLICIT NONE
      REAL*8 a,b,c,theta
c     --- Zero is in the plane and orbiting in the direction of
c         galactic rotation
      theta=acos(-c/sqrt(a*a+b*b+c*c))
      return
      end

      SUBROUTINE FINDPLANE(p,a,b,c,d)
      IMPLICIT NONE
      REAL*8 p(3,3),a,b,c,d
      REAL*8 p1p2(3),p1p3(3)
      p1p2(1)=p(2,1)-p(1,1)
      p1p2(2)=p(2,2)-p(1,2)
      p1p2(3)=p(2,3)-p(1,3)
      p1p3(1)=p(3,1)-p(1,1)
      p1p3(2)=p(3,2)-p(1,2)
      p1p3(3)=p(3,3)-p(1,3)
      a = p1p2(2)*p1p3(3) - p1p3(2)*p1p2(3)
      b =-p1p2(1)*p1p3(3) + p1p3(1)*p1p2(3)
      c = p1p2(1)*p1p3(2) - p1p3(1)*p1p2(2)
      d = a*(-p(1,1)) + b*(-p(1,2)) + c*(-p(1,3))
      return
      end

      SUBROUTINE PICK3POINTS(NC,NS,idpe,idap,index,x,y,z,p)
      IMPLICIT NONE
      INTEGER NS,NC
      REAL*8 x(NS),y(NS),z(NS),p(3,3)
      INTEGER index,idpe(NC),idap(NC)
      if (idpe(index).lt.idap(index)) then
c        --- perigalacticon occurs first
c        --- Perigalacticon first
         p(1,1)=x(idpe(index))
         p(1,2)=y(idpe(index))
         p(1,3)=z(idpe(index))
c        --- Apogalacticon last
         p(3,1)=x(idap(index))
         p(3,2)=y(idap(index))
         p(3,3)=z(idap(index))
c        --- Between is between
         p(2,1)=x((idpe(index)+idap(index))/2)
         p(2,2)=y((idpe(index)+idap(index))/2)
         p(2,3)=z((idpe(index)+idap(index))/2)
      else
c        --- apogalacticon occurs first
c        --- Apogalacticon first
         p(1,1)=x(idap(index))
         p(1,2)=y(idap(index))
         p(1,3)=z(idap(index))
c        --- Perigalacticon last
         p(3,1)=x(idpe(index))
         p(3,2)=y(idpe(index))
         p(3,3)=z(idpe(index))
c        --- Between is between
         p(2,1)=x((idpe(index)+idap(index))/2)
         p(2,2)=y((idpe(index)+idap(index))/2)
         p(2,3)=z((idpe(index)+idap(index))/2)
      end if
      return
      end         

      SUBROUTINE PLOTENERGYPERIOD(NN,NCOUNT,e0,ef,porb)
      IMPLICIT NONE
      INTEGER NN,NCOUNT
      REAL*8 e0(NN),ef(NN),porb(NN)
      REAL xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt
      REAL s(NN),t(NN)
      INTEGER i
      do i=1,NCOUNT
         s(i)=porb(i)
         t(i)=(ef(i)-e0(i))/e0(i)
      end do
      CALL MIN1(s,NN,NCOUNT,xwl)
      CALL MAX1(s,NN,NCOUNT,xwr)
      CALL MIN1(t,NN,NCOUNT,ywb)
      CALL MAX1(t,NN,NCOUNT,ywt)
      xvl=0.2
      xvr=0.8
      yvb=0.2
      yvt=0.8
      call pgvport(xvl,xvr,yvb,yvt)
      call pgwindow(xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),
     &             ywb-0.05*(ywt-ywb),ywt+0.05*(ywt-ywb))
      call pgbox('BCNT',0.0,0,'BCNT',0.0,0)
      call pgmtxt('L',2.4,0.5,0.5,'(E\\df\\u - E\\di\\u)/E\\di\\u')
      call pgmtxt('B',2.4,0.5,0.5,'T [Gyr]')      
      call pgpoint(NCOUNT,s,t,1)
      call pgpage
      do i=1,NCOUNT
         s(i)=e0(i)
         t(i)=porb(i)
      end do
      CALL MIN1(s,NN,NCOUNT,xwl)
      CALL MAX1(s,NN,NCOUNT,xwr)
      CALL MIN1(t,NN,NCOUNT,ywb)
      CALL MAX1(t,NN,NCOUNT,ywt)
      xvl=0.2
      xvr=0.8
      yvb=0.2
      yvt=0.8
      call pgvport(xvl,xvr,yvb,yvt)
      call pgwindow(xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),
     &             ywb-0.05*(ywt-ywb),ywt+0.05*(ywt-ywb))
      call pgbox('BCNT',0.0,0,'BCNT',0.0,0)
      call pgmtxt('B',2.4,0.5,0.5,'E\\di\\u')
      call pgmtxt('L',2.4,0.5,0.5,'T [Gyr]')      
      call pgpoint(NCOUNT,s,t,1)
      return
      end

      SUBROUTINE PLOTENERGYPERI(NN,NCOUNT,e0,ef,rp)
      IMPLICIT NONE
      INTEGER NN,NCOUNT
      REAL*8 e0(NN),ef(NN),rp(NN)
      REAL xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt
      REAL s(NN),t(NN)
      INTEGER i
      do i=1,NCOUNT
         s(i)=rp(i)
         t(i)=(ef(i)-e0(i))/e0(i)
      end do
      CALL MIN1(s,NN,NCOUNT,xwl)
      CALL MAX1(s,NN,NCOUNT,xwr)
      CALL MIN1(t,NN,NCOUNT,ywb)
      CALL MAX1(t,NN,NCOUNT,ywt)
      xvl=0.2
      xvr=0.8
      yvb=0.2
      yvt=0.8
      call pgvport(xvl,xvr,yvb,yvt)
      call pgwindow(xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),
     &             ywb-0.05*(ywt-ywb),ywt+0.05*(ywt-ywb))
      call pgbox('BCNT',0.0,0,'BCNT',0.0,0)
      call pgmtxt('L',2.4,0.5,0.5,'(E\\df\\u - E\\di\\u)/E\\di\\u')
      call pgmtxt('B',2.4,0.5,0.5,'R\\dp\\u [kpc]')      
      call pgpoint(NCOUNT,s,t,1)
      return
      end

      SUBROUTINE PLOTENERGYECC(NN,NCOUNT,e0,ef,e)
      IMPLICIT NONE
      INTEGER NN,NCOUNT
      REAL*8 e0(NN),ef(NN),e(NN)
      REAL xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt
      REAL s(NN),t(NN)
      INTEGER i
      do i=1,NCOUNT
         s(i)=e(i)
         t(i)=(ef(i)-e0(i))/e0(i)
      end do
      CALL MIN1(s,NN,NCOUNT,xwl)
      CALL MAX1(s,NN,NCOUNT,xwr)
      CALL MIN1(t,NN,NCOUNT,ywb)
      CALL MAX1(t,NN,NCOUNT,ywt)
      xvl=0.2
      xvr=0.8
      yvb=0.2
      yvt=0.8
      call pgvport(xvl,xvr,yvb,yvt)
      call pgwindow(xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),
     &             ywb-0.05*(ywt-ywb),ywt+0.05*(ywt-ywb))
      call pgbox('BCNT',0.0,0,'BCNT',0.0,0)
      call pgmtxt('L',2.4,0.5,0.5,'(E\\df\\u - E\\di\\u)/E\\di\\u')
      call pgmtxt('B',2.4,0.5,0.5,'Eccentricity')      
      call pgpoint(NCOUNT,s,t,1)
      return
      end
c******************************************************************
      SUBROUTINE PLOTENERGY(NN,NCOUNT,e0,ef)
      IMPLICIT NONE
      INTEGER NN,NCOUNT
      REAL*8 e0(NN),ef(NN)
      REAL xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt
      REAL s(NN),t(NN)
      INTEGER i
      do i=1,NCOUNT
         s(i)=e0(i)
         t(i)=(ef(i)-e0(i))/e0(i)
      end do
      CALL MIN1(s,NN,NCOUNT,xwl)
      CALL MAX1(s,NN,NCOUNT,xwr)
      CALL MIN1(t,NN,NCOUNT,ywb)
      CALL MAX1(t,NN,NCOUNT,ywt)
      xvl=0.2
      xvr=0.8
      yvb=0.2
      yvt=0.8
      call pgvport(xvl,xvr,yvb,yvt)
      call pgwindow(xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),
     &             ywb-0.05*(ywt-ywb),ywt+0.05*(ywt-ywb))
      call pgbox('BCNT',0.0,0,'BCNT',0.0,0)
      call pgmtxt('L',2.4,0.5,0.5,'(E\\df\\u - E\\di\\u)/E\\di\\u')
      call pgmtxt('B',2.4,0.5,0.5,'E\\di\\u')      
      call pgpoint(NCOUNT,s,t,1)
      return
      end
      
c********************************************************************
      SUBROUTINE PLOTHIST(NN,NCOUNT,rp,ra,e,porb,theta,omega,
     &              hrpmax,hramax,hemax,htmax,hphimax,homegamax)
      IMPLICIT NONE
      INTEGER NN,NCOUNT
      REAL*8 rp(NN),ra(NN),e(NN),porb(NN),theta(NN)
      REAL*8 omega(NN),d2r
      REAL s(NN)
      REAL xwl,xwr,ywb,ywt
      REAL xvl,xvr,yvb,yvt
      REAL hrpmax,hramax,hemax,htmax,hphimax,homegamax
      INTEGER i
      parameter (d2r=1.74532925199D-02)
c     --- Plot histogram of perigalacticons
      do i=1,NCOUNT
        s(i)=rp(i)
      end do
      CALL MIN1(s,NN,NCOUNT,xwl)
      CALL MAX1(s,NN,NCOUNT,xwr)
      ywb=0.0
      ywt=hrpmax
      xvl=0.1
      xvr=0.45
      yvb=0.55
      yvt=0.9
      call pgvport(xvl,xvr,yvb,yvt)
      call pgwindow(xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),ywb,ywt)
      call pgbox('BCNT',0.0,0,'BCNT',0.0,0)
      call pghist(NCOUNT,s,xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),20,1)
      call pgmtxt('L',2.4,0.5,0.5,'No.')
      call pgmtxt('B',2.4,0.5,0.5,'R\\dp\\u [kpc]')
c     --- Plot histogram of apogalacticons
      do i=1,NCOUNT
        s(i)=ra(i)
      end do
      CALL MIN1(s,NN,NCOUNT,xwl)
      CALL MAX1(s,NN,NCOUNT,xwr)
      ywb=0.0
      ywt=hramax
      xvl=0.55
      xvr=0.9
      yvb=0.55
      yvt=0.9
      call pgvport(xvl,xvr,yvb,yvt)
      call pgwindow(xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),ywb,ywt)
      call pgbox('BCNT',0.0,0,'BCNT',0.0,0)
      call pghist(NCOUNT,s,xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),20,1)
      call pgmtxt('L',2.4,0.5,0.5,'No.')
      call pgmtxt('B',2.4,0.5,0.5,'R\\da\\u [kpc]')
c     --- Plot histograms of eccentricities
      do i=1,NCOUNT
        s(i)=e(i)
      end do
      CALL MIN1(s,NN,NCOUNT,xwl)
      CALL MAX1(s,NN,NCOUNT,xwr)
      ywb=0.0
      ywt=hemax
      xvl=0.1
      xvr=0.45
      yvb=0.1
      yvt=0.45
      call pgvport(xvl,xvr,yvb,yvt)
      call pgwindow(xwl,xwr,ywb,ywt)
      call pgbox('BCNT',0.0,0,'BCNT',0.0,0)
      call pghist(NCOUNT,s,xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),20,1)
      call pgmtxt('L',2.4,0.5,0.5,'No.')
      call pgmtxt('B',2.4,0.5,0.5,'e')
c     --- Plot histograms of periods
      do i=1,NCOUNT
        s(i)=porb(i)
      end do
      CALL MIN1(s,NN,NCOUNT,xwl)
      CALL MAX1(s,NN,NCOUNT,xwr)
      ywb=0.0
      ywt=htmax
      xvl=0.55
      xvr=0.9
      yvb=0.1
      yvt=0.45
      call pgvport(xvl,xvr,yvb,yvt)
      call pgwindow(xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),ywb,ywt)
      call pgbox('BCNT',0.0,0,'BCNT',0.0,0)
      call pghist(NCOUNT,s,xwl-0.05*(xwr-xwl),
     &                     xwr+0.05*(xwr-xwl),20,1)      
      call pgmtxt('L',2.4,0.5,0.5,'No.')
      call pgmtxt('B',2.4,0.5,0.5,'T [Gyr]')
      call pgpage
c     --- Plot histogram of Inclinations
      do i=1,NCOUNT
        s(i)=theta(i)/d2r
      end do
      CALL MIN1(s,NN,NCOUNT,xwl)
      CALL MAX1(s,NN,NCOUNT,xwr)
      ywb=0.0
      ywt=hphimax
      xvl=0.1
      xvr=0.45
      yvb=0.55
      yvt=0.9
      call pgvport(xvl,xvr,yvb,yvt)
      call pgwindow(xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),ywb,ywt)
      call pgbox('BCNT',0.0,0,'BCNT',0.0,0)
      call pghist(NCOUNT,s,xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),20,1)
      call pgmtxt('L',2.4,0.5,0.5,'No.')
      call pgmtxt('B',2.4,0.5,0.5,'\\gF [degrees]')
c     --- Plot histogram of longitudes
      do i=1,NCOUNT
        s(i)=omega(i)/d2r
      end do
      CALL MIN1(s,NN,NCOUNT,xwl)
      CALL MAX1(s,NN,NCOUNT,xwr)
      ywb=0.0
      ywt=homegamax
      xvl=0.55
      xvr=0.9
      yvb=0.55
      yvt=0.9
      call pgvport(xvl,xvr,yvb,yvt)
      call pgwindow(xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),ywb,ywt)
      call pgbox('BCNT',0.0,0,'BCNT',0.0,0)
      call pghist(NCOUNT,s,xwl-0.05*(xwr-xwl),xwr+0.05*(xwr-xwl),20,1)
      call pgmtxt('L',2.4,0.5,0.5,'No.')
      call pgmtxt('B',2.4,0.5,0.5,'\\gW [degrees]')
      return
      end
c*********************************************************************
c     I don't remember where this came from, but it has an NFW halo,
c     a Hernquist model bulge, and a Miyamoto & Nagai disk.
      SUBROUTINE NFORCE(x,y,z,fx,fy,fz)
      IMPLICIT NONE
      REAL*8 a,b,c,md,ms,g,x,y,z,fx,fy,fz
      PARAMETER (a=5.0d0,b=0.26d0,c=0.7d0,md=1.0d11)
      PARAMETER (ms=0.8d10,g=4.498933261d-6)
      REAL*8 f_r,fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3
      REAL*8 f_x,f_c,c_n,M_vir,x_n,r_v
      PARAMETER (M_vir=1.0d12,c_n=12.0d0,r_v=258.0d0)
      REAL*8 r,rd
      r  = sqrt(x*x+y*y+z*z)
      rd = sqrt(x*x+y*y)
      x_n = c_n*r/r_v
      f_x = dlog(1.0d0+x_n) - x_n/(1.0d0+x_n)
      f_c = dlog(1.0d0+c_n) - c_n/(1.0d0+c_n)
      f_r=-g*M_vir*f_x/(r*r*f_c)
c     --- Find the X components of the force.
      fx1=-g*md*x/(sqrt(rd*rd+(a+sqrt(z*z+b*b))**2))**3 !F_x disk
      fx2=-g*ms*x/(r*(r+c)*(r+c))                       !F_x bulge
      fx3= x*f_r/r                                      !F_x halo
      fx=fx1+fx2+fx3
c     --- Find the Y components of the force.
      fy1=-g*md*y/(sqrt(rd*rd+(a+sqrt(z*z+b*b))**2))**3 !F_y disk
      fy2=-g*ms*y/(r*(r+c)*(r+c))                       !F_y bulge
      fy3= y*f_r/r                                      !F_y halo
      fy=fy1+fy2+fy3
c     --- Find the Z components of the force.
      fz1=-g*md*z/(sqrt(rd*rd+(a+sqrt(z*z+b*b))**2))**3 
      fz1=fz1*(1.0+a/sqrt(z*z+b*b))                     !F_z disk
      fz2=-g*ms*z/(r*(r+c)*(r+c))                       !F_z bulge
      fz3=z*f_r/r                                       !F_z halo
      fz=fz1+fz2+fz3
      return
      end
c***************************************************************
      SUBROUTINE NENERGY(vcon,x,y,z,vx,vy,vz,e)
      IMPLICIT NONE
      REAL*8 x,y,z,vx,vy,vz,e,vcon,g
      PARAMETER (g=4.498933261d-6)
      REAL*8 a,b,md,rc,ms
c     --- bulge parameters
      PARAMETER (rc=0.7d0,ms=0.8d10)
c     --- disk parameters
      PARAMETER (a=5.0d0,b=0.26d0,md=1.0d11)
      REAL*8 m_vir,c_n,r_v
c     --- halo parameters
      PARAMETER (m_vir=1.0d12,c_n=12.0d0,r_v=258.0d0)
c
      REAL*8 r,rd2,ke,phi1,phi2,phi3
c
c     -- The kinetic energy per unit mass
      ke=vcon*vcon*(vx*vx+vy*vy+vz*vz)/2.0d0
c
c     -- The potential energy per unit mass
      r  = sqrt(x*x+y*y+z*z)
      rd2 = x*x+y*y
c     -- disk
      phi1=-g*md/sqrt(rd2+(a+sqrt(z*z+b*b))**2)
c     -- bulge
      phi2 = -g*ms/(r+rc)
c     -- halo
      phi3 = -g*(m_vir/(log(1.0d0+c_n)-c_n/(1.0d0+c_n))) *
     &                                          log(1.0d0+r*c_n/r_v)/r
c
      e=ke+phi1+phi2+phi3
c      write (6,*) ke,phi1,phi2,phi3,e
c
      return
      end
c**************************************************************
c     Comes from Dehnen, McLaughlin, & Sacania 2006, MNRAS, 269, 1688
c     Has a Dehnen & Binney 1998 disk, a Dehnen & McLaughlin 2005 halo,
c     and a Hernquist model bulge.
      SUBROUTINE DFORCE(x,y,z,fx,fy,fz)
      IMPLICIT NONE
      REAL*8 a,b,c,md,ms,g,x,y,z,fx,fy,fz,Mh,r0
      PARAMETER (a=6.5d0,b=0.26d0,c=0.7d0,md=1.0d11)
      PARAMETER (ms=3.4d10,g=4.498933261d-6)
      PARAMETER (Mh=1.10d13,r0=40.5d0)
      REAL*8 f_r,fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3,p
      REAL*8 r,rd
      r  = sqrt(x*x+y*y+z*z)
      rd = sqrt(x*x+y*y)
      p=4.0d0/9.0d0
      f_r=-(g*Mh/r**2.0d0)*((r**(p))/(r**(p) + r0**(p)))**5.0d0
c     --- Find the X components of the force.
      fx1=-g*md*x/(sqrt(rd*rd+(a+sqrt(z*z+b*b))**2))**3 !F_x disk
      fx2=-g*ms*x/(r*(r+c)*(r+c))                       !F_x bulge
      fx3= x*f_r/r                                      !F_x halo
      fx=fx1+fx2+fx3
c     --- Find the Y components of the force.
      fy1=-g*md*y/(sqrt(rd*rd+(a+sqrt(z*z+b*b))**2))**3 !F_y disk
      fy2=-g*ms*y/(r*(r+c)*(r+c))                       !F_y bulge
      fy3= y*f_r/r                                      !F_y halo
      fy=fy1+fy2+fy3
c     --- Find the Z components of the force.
      fz1=-g*md*z/(sqrt(rd*rd+(a+sqrt(z*z+b*b))**2))**3 
      fz1=fz1*(1.0+a/sqrt(z*z+b*b))                     !F_z disk
      fz2=-g*ms*z/(r*(r+c)*(r+c))                       !F_z bulge
      fz3=z*f_r/r                                       !F_z halo
      fz=fz1+fz2+fz3
      return
      end
c*************************************************************
      SUBROUTINE DENERGY(vcon,x,y,z,vx,vy,vz,e)
      IMPLICIT NONE
      REAL*8 a,b,c,md,ms,vh,g,x,y,z,vx,vy,vz,e,vcon,Mh,r0
      PARAMETER (a=6.5d0,b=0.26d0,c=0.7d0,md=1.0d11)
      PARAMETER (ms=3.4d10,vh=128.0d0,g=4.498933261d-6)
      PARAMETER (Mh=1.10d13,r0=40.5)
      REAL*8 r,rd,ke,phi1,phi2,phi3,p
c
c     -- The kinetic energy per unit mass
      ke=vcon*vcon*(vx*vx+vy*vy+vz*vz)/2.0d0
c
c     -- The potential energy per unit mass
      r  = sqrt(x*x+y*y+z*z)
      rd = sqrt(x*x+y*y)
c     -- disk
      phi1=-g*md/sqrt(rd*rd+(a+sqrt(z*z+b*b))**2)
c     -- spheroid
      phi2=-g*ms/(r+c)
c     -- halo
      p=4.0d0/9.0d0
      phi3=-(g*Mh/r)*((r**(p))/(r**(p) + r0**(p)))**5.0d0
c
      e=ke+phi1+phi2+phi3
      write(6,*) 'got here'
      write (6,*) ke,phi1,phi2,phi3,e
c
      return
      end
c**************************************************************
c     From Johnston, Spergel, and Hernquist 1995, ApJ, 451, 598.
c     Has Miyamoto and Nagai disk, Hernquist model bulge, and a halo
c     with a logarithmic potential. 
      SUBROUTINE JFORCE(x,y,z,fx,fy,fz)
      IMPLICIT NONE
      REAL*8 a,b,c,d,md,ms,vh,g,x,y,z,fx,fy,fz
      PARAMETER (a=6.5d0,b=0.26d0,c=0.7d0,d=12.0d0,md=1.0d11)
      PARAMETER (ms=3.4d10,vh=128.0d0,g=4.498933261d-6)
      REAL*8 fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3
      REAL*8 r,rd
      r  = sqrt(x*x+y*y+z*z)
      rd = sqrt(x*x+y*y)
c     --- Find the X components of the force.
      fx1=-g*md*x/(sqrt(rd*rd+(a+sqrt(z*z+b*b))**2))**3
      fx2=-g*ms*x/(r*(r+c)*(r+c))
      fx3=-2.0d0*vh*vh*x/(r*r+d*d)
      fx=fx1+fx2+fx3
c     --- Find the Y components of the force.
      fy1=-g*md*y/(sqrt(rd*rd+(a+sqrt(z*z+b*b))**2))**3
      fy2=-g*ms*y/(r*(r+c)*(r+c))
      fy3=-2.0d0*vh*vh*y/(r*r+d*d)
      fy=fy1+fy2+fy3
c     --- Find the Z components of the force.
      fz1=-g*md*z/(sqrt(rd*rd+(a+sqrt(z*z+b*b))**2))**3
      fz1=fz1*(1.0+a/sqrt(z*z+b*b))
      fz2=-g*ms*z/(r*(r+c)*(r+c))
      fz3=-2.0d0*vh*vh*z/(r*r+d*d)
      fz=fz1+fz2+fz3
      return
      end
c***********************************************************************
      SUBROUTINE JENERGY(vcon,x,y,z,vx,vy,vz,e)
      IMPLICIT NONE
      REAL*8 a,b,c,d,md,ms,vh,g,x,y,z,vx,vy,vz,e,vcon
      PARAMETER (a=6.5d0,b=0.26d0,c=0.7d0,d=12.0d0,md=1.0d11)
      PARAMETER (ms=3.4d10,vh=128.0d0,g=4.498933261d-6)
      REAL*8 r,rd,ke,phi1,phi2,phi3
c
c     -- The kinetic energy per unit mass
      ke=vcon*vcon*(vx*vx+vy*vy+vz*vz)/2.0d0
c
c     -- The potential energy per unit mass
      r  = sqrt(x*x+y*y+z*z)
      rd = sqrt(x*x+y*y)
c     -- disk
      phi1=-g*md/sqrt(rd*rd+(a+sqrt(z*z+b*b))**2)
c     -- spheroid
      phi2=-g*ms/(r+c)
c     -- halo
      phi3=vh*vh*dlog(r*r+d*d)
c
      e=ke+phi1+phi2+phi3
c      write (6,*) ke,phi1,phi2,phi3,e
c
      return
      end
c**************************************************************
c     Implements the MWPotential2014 from Bovy 2015, ApJS, 216:29.
c     Has a Miyamoto and Nagai disk, a cutoff r^-1.8 bulge, and an NFW halo.
      SUBROUTINE BFORCE(x,y,z,fx,fy,fz)
      IMPLICIT NONE
      REAL*8 x,y,z,fx,fy,fz
      REAL*8 mb,rc,md,a,b,mv,c,ra,g
      PARAMETER (g=4.498933261d-6)
c     --- bulge parameters
      PARAMETER (mb=5.0d9,rc=1.9)
c     --- disk parameters
      PARAMETER (a=3.0d0,b=0.28d0,md=6.8d10)
c     --- halo parameters
      PARAMETER (mv=8.0d11,c=15.3d0,ra=16.0d0)
      REAL*8 r,rd2,gammq,fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3,fb_r,fh_r
      REAL*8 mb_r,mnfw_r
      r  = sqrt(x*x+y*y+z*z)
      rd2 = x*x+y*y
c     --- bulge radial force: -GM_b(r)/r*r
      mb_r=mb*(1.00d0 - 1.0891d0*exp(-r/rc)*(r/rc)**0.2 -
     &                                    gammq(0.2d0,r/rc))
      fb_r = -g*mb_r/(r*r)
c     --- NFW halo radial force: -GM_nfw(r)/r*r
      mnfw_r=mv*(log(1.0d0+r/ra)-r/(ra+r))/(log(1.0d0+c)-c/(1.0d0+c))
      fh_r = -g*mnfw_r/(r*r)
c     --- Find the X components of the force.
      fx1=-g*md*x/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fx2=x*fb_r/r
      fx3=x*fh_r/r
      fx=fx1+fx2+fx3
c     --- Find the Y components of the force.
      fy1=-g*md*y/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fy2=y*fb_r/r
      fy3=y*fh_r/r
      fy=fy1+fy2+fy3
c     --- Find the Z components of the force.
      fz1=-g*md*z/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fz1=fz1*(1.0+a/sqrt(z*z+b*b))
      fz2=z*fb_r/r
      fz3=z*fh_r/r
      fz=fz1+fz2+fz3
      return
      end
c***************************************************************
      SUBROUTINE BENERGY(vcon,x,y,z,vx,vy,vz,e)
      IMPLICIT NONE
      REAL*8 x,y,z,vx,vy,vz,e,vcon,g
      PARAMETER (g=4.498933261d-6)
      REAL*8 mb,rc,md,a,b,mv,c,ra
c     --- bulge parameters
      PARAMETER (mb=5.0d9,rc=1.9)
c     --- disk parameters
      PARAMETER (a=3.0d0,b=0.28d0,md=6.8d10)
c     --- halo parameters
      PARAMETER (mv=8.0d11,c=15.3d0,ra=16.0d0)
c
      REAL*8 r,rd2,ke,phi1,phi2,phi3
c
c     -- The kinetic energy per unit mass
      ke=vcon*vcon*(vx*vx+vy*vy+vz*vz)/2.0d0
c
c     -- The potential energy per unit mass
      r  = sqrt(x*x+y*y+z*z)
      rd2 = x*x+y*y
c     -- disk
      phi1=-g*md/sqrt(rd2+(a+sqrt(z*z+b*b))**2)
c     -- bulge; this isn't right, but is accurate for r >> rc and that
c        is the region of interest for all or almost all of our range of
c        integration
      if (r.gt.rc/2.0d0) then
         phi2 = -g*mb/r
      else
         phi2 = -g*2.0d0*mb/rc
      end if
c     -- halo
      phi3 = -g*(mv/(log(1.0d0+c)-c/(1.0d0+c))) * log(1.0d0+r/ra)/r
c
      e=ke+phi1+phi2+phi3
c      write (6,*) ke,phi1,phi2,phi3,e
c
      return
      end
c**************************************************************
c     Implements the MW model from Kafle et al., ApJ, 794:59.
c     Has a Miyamoto and Nagai disk, a cutoff r^-1.8 bulge, and an NFW halo.
c     I am ignoring the complication that the model was fitted with a
c     flattened bulge -- the bulge doesn't really matter for our integrations.
      SUBROUTINE KFORCE(x,y,z,fx,fy,fz)
      IMPLICIT NONE
      REAL*8 x,y,z,fx,fy,fz
      REAL*8 mb,rc,md,a,b,mv,c,ra,g
      PARAMETER (g=4.498933261d-6)
c     --- bulge parameters
      PARAMETER (mb=9.1d9,rc=2.1d0)
c     --- disk parameters
      PARAMETER (a=4.9d0,b=0.30d0,md=9.5d10)
c     --- halo parameters
      PARAMETER (mv=8.0d11,c=21.1d0,ra=11.33d0)
      REAL*8 r,rd2,gammq,fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3,fb_r,fh_r
      REAL*8 mb_r,mnfw_r
      r  = sqrt(x*x+y*y+z*z)
      rd2 = x*x+y*y
c     --- bulge radial force: -GM_b(r)/r*r
      mb_r=mb*(1.00d0 - 1.0891d0*exp(-r/rc)*(r/rc)**0.2 -
     &                                    gammq(0.2d0,r/rc))
      fb_r = -g*mb_r/(r*r)
c     --- NFW halo radial force: -GM_nfw(r)/r*r
      mnfw_r=mv*(log(1.0d0+r/ra)-r/(ra+r))/(log(1.0d0+c)-c/(1.0d0+c))
      fh_r = -g*mnfw_r/(r*r)
c     --- Find the X components of the force.
      fx1=-g*md*x/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fx2=x*fb_r/r
      fx3=x*fh_r/r
      fx=fx1+fx2+fx3
c     --- Find the Y components of the force.
      fy1=-g*md*y/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fy2=y*fb_r/r
      fy3=y*fh_r/r
      fy=fy1+fy2+fy3
c     --- Find the Z components of the force.
      fz1=-g*md*z/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fz1=fz1*(1.0+a/sqrt(z*z+b*b))
      fz2=z*fb_r/r
      fz3=z*fh_r/r
      fz=fz1+fz2+fz3
      return
      end
c***************************************************************
      SUBROUTINE KENERGY(vcon,x,y,z,vx,vy,vz,e)
      IMPLICIT NONE
      REAL*8 x,y,z,vx,vy,vz,e,vcon,g
      PARAMETER (g=4.498933261d-6)
      REAL*8 mb,rc,md,a,b,mv,c,ra
c     --- bulge parameters
      PARAMETER (mb=9.1d9,rc=2.1d0)
c     --- disk parameters
      PARAMETER (a=4.9d0,b=0.30d0,md=9.5d10)
c     --- halo parameters
      PARAMETER (mv=8.0d11,c=21.1d0,ra=11.33d0)
c
      REAL*8 r,rd2,ke,phi1,phi2,phi3
c
c     -- The kinetic energy per unit mass
      ke=vcon*vcon*(vx*vx+vy*vy+vz*vz)/2.0d0
c
c     -- The potential energy per unit mass
      r  = sqrt(x*x+y*y+z*z)
      rd2 = x*x+y*y
c     -- disk
      phi1=-g*md/sqrt(rd2+(a+sqrt(z*z+b*b))**2)
c     -- bulge; this isn't right, but is accurate for r >> rc and that
c        is the region of interest for all or almost all of our range of
c        integration
      if (r.gt.rc/2.0d0) then
         phi2 = -g*mb/r
      else
         phi2 = -g*2.0d0*mb/rc
      end if
c     -- halo
      phi3 = -g*(mv/(log(1.0d0+c)-c/(1.0d0+c))) * log(1.0d0+r/ra)/r
c
      e=ke+phi1+phi2+phi3
c      write (6,*) ke,phi1,phi2,phi3,e
c
      return
      end
c**************************************************************
c     Implements an approximate MW model using a dark halo from van der
c     Marel 2014 and a disk and bulge from Kafle et al., ApJ, 794:59.
c     Has a Miyamoto and Nagai disk, a cutoff r^-1.8 bulge, and an NFW halo.
c     I am ignoring the complication that the model was fitted with a
c     flattened bulge -- the bulge doesn't really matter for our integrations.
c     Following van der Marel, the halo is represented with a Hernquist
c     model matched to have the same M_vir as an NFW halo with c_vir=10
c     and to have the same density as r->0.  The Bland-Hawthorn & Gerhard
c     2016 ARAA review of the MW suggests r_vir = 282 kpc and M_vir =
c     1.3 (+-0.3) x 10^12 M_sun.  Reducing M_vir to 1.2 x 10^12 to account
c     for the disk and bulge mass, then using the formulae from appendix
c     of van der Marel (2014) gives ra=2.09(r_vir/c)=0.209r_vir = 58.9 kpc
c     and M_h = 1.46 M_vir = 1.75x10^12 M_sun.
      SUBROUTINE VFORCE(x,y,z,fx,fy,fz)
      IMPLICIT NONE
      REAL*8 x,y,z,fx,fy,fz
      REAL*8 mb,rc,md,a,b,mh,ra,g
      PARAMETER (g=4.498933261d-6)
c     --- bulge parameters
      PARAMETER (mb=9.1d9,rc=2.1d0)
c     --- disk parameters
ccc      PARAMETER (a=4.9d0,b=0.30d0,md=9.5d10)
c     --- increase disk mass to get higher v_c at R=R_0.
      PARAMETER (a=4.9d0,b=0.30d0,md=1.3d11)
c     --- halo parameters
      PARAMETER (mh=1.7d12,ra=58.9d0)
      REAL*8 r,rd2,gammq,fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3,fb_r,fh_r
      REAL*8 mb_r
ccc      REAL*8 mh_r
      r  = sqrt(x*x+y*y+z*z)
      rd2 = x*x+y*y
c     --- bulge radial force: -GM_b(r)/r*r
      mb_r=mb*(1.00d0 - 1.0891d0*exp(-r/rc)*(r/rc)**0.2 -
     &                                    gammq(0.2d0,r/rc))
      fb_r = -g*mb_r/(r*r)
c     --- Hernquist halo radial force: -GM_h(r)/r*r
ccc      mh_r=mh*r*r/(ra+r)**2
      fh_r = -g*mh/(ra+r)**2
c     --- Find the X components of the force.
      fx1=-g*md*x/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fx2=x*fb_r/r
      fx3=x*fh_r/r
      fx=fx1+fx2+fx3
c     --- Find the Y components of the force.
      fy1=-g*md*y/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fy2=y*fb_r/r
      fy3=y*fh_r/r
      fy=fy1+fy2+fy3
c     --- Find the Z components of the force.
      fz1=-g*md*z/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fz1=fz1*(1.0+a/sqrt(z*z+b*b))
      fz2=z*fb_r/r
      fz3=z*fh_r/r
      fz=fz1+fz2+fz3
      return
      end
c***************************************************************
      SUBROUTINE VENERGY(vcon,x,y,z,vx,vy,vz,e)
      IMPLICIT NONE
      REAL*8 x,y,z,vx,vy,vz,e,vcon,g
      PARAMETER (g=4.498933261d-6)
      REAL*8 mb,rc,md,a,b,mh,ra
c     --- bulge parameters
      PARAMETER (mb=9.1d9,rc=2.1d0)
c     --- disk parameters
c     --- increase disk mass to get higher v_c at R=R_0.
      PARAMETER (a=4.9d0,b=0.30d0,md=1.3d11)
c     --- halo parameters
      PARAMETER (mh=1.7d12,ra=58.9d0)
c
      REAL*8 r,rd2,ke,phi1,phi2,phi3
c
c     -- The kinetic energy per unit mass
      ke=vcon*vcon*(vx*vx+vy*vy+vz*vz)/2.0d0
c
c     -- The potential energy per unit mass
      r  = sqrt(x*x+y*y+z*z)
      rd2 = x*x+y*y
c     -- disk
      phi1=-g*md/sqrt(rd2+(a+sqrt(z*z+b*b))**2)
c     -- bulge; this isn't right, but is accurate for r >> rc and that
c        is the region of interest for all or almost all of our range of
c        integration
      if (r.gt.rc/2.0d0) then
         phi2 = -g*mb/r
      else
         phi2 = -g*2.0d0*mb/rc
      end if
c     -- halo
      phi3 = -g*mh/(ra+r)
c
      e=ke+phi1+phi2+phi3
c      write (6,*) ke,phi1,phi2,phi3,e
c
      return
      end
c**************************************************************
c     Implements an approximate MW model using a dark halo from van der
c     Marel 2014 and a disk and bulge from Kafle et al., ApJ, 794:59.
c     Has a Miyamoto and Nagai disk, a cutoff r^-1.8 bulge, and an NFW halo.
c     I am ignoring the complication that the model was fitted with a
c     flattened bulge -- the bulge doesn't really matter for our integrations.
c     This version actually uses an NFW model for the halo instead of a
c     Hernquist model (used by van der Marel).  The Bland-Hawthorn & Gerhard
c     2016 ARAA review of the MW suggests r_vir = 282 kpc and M_vir =
c     1.3 (+-0.3) x 10^12 M_sun with c=10 (from van der Marel et al.).
c     Actually reduced M_vir to 1.2 x 10^12 Msun to account for the bulge
c     and disk mass.  Actually had to increase the disk mass to get V_c
c     at R_0 high enough.
      SUBROUTINE XFORCE(x,y,z,fx,fy,fz)
      IMPLICIT NONE
      REAL*8 x,y,z,fx,fy,fz
      REAL*8 mb,rc,md,a,b,mv,ra,c,g
      PARAMETER (g=4.498933261d-6)
c     --- bulge parameters
      PARAMETER (mb=9.1d9,rc=2.1)
c     --- disk parameters
ccc      PARAMETER (a=4.9d0,b=0.30d0,md=9.5d10)
c     --- increase disk mass to get higher v_c at R=R_0.
      PARAMETER (a=4.9d0,b=0.30d0,md=1.35d11)
c     --- halo parameters
      PARAMETER (mv=1.16d12,c=10.0d0,ra=28.2d0)
      REAL*8 r,rd2,gammq,fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3,fb_r,fh_r
      REAL*8 mb_r,mnfw_r
      r  = sqrt(x*x+y*y+z*z)
      rd2 = x*x+y*y
c     --- bulge radial force: -GM_b(r)/r*r
      mb_r=mb*(1.00d0 - 1.0891d0*exp(-r/rc)*(r/rc)**0.2 -
     &                                    gammq(0.2d0,r/rc))
      fb_r = -g*mb_r/(r*r)
c     --- NFW halo radial force: -GM_nfw(r)/r*r
      mnfw_r=mv*(log(1.0d0+r/ra)-r/(ra+r))/(log(1.0d0+c)-c/(1.0d0+c))
      fh_r = -g*mnfw_r/(r*r)
c     --- Find the X components of the force.
      fx1=-g*md*x/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fx2=x*fb_r/r
      fx3=x*fh_r/r
      fx=fx1+fx2+fx3
c     --- Find the Y components of the force.
      fy1=-g*md*y/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fy2=y*fb_r/r
      fy3=y*fh_r/r
      fy=fy1+fy2+fy3
c     --- Find the Z components of the force.
      fz1=-g*md*z/(sqrt(rd2+(a+sqrt(z*z+b*b))**2))**3
      fz1=fz1*(1.0+a/sqrt(z*z+b*b))
      fz2=z*fb_r/r
      fz3=z*fh_r/r
      fz=fz1+fz2+fz3
      return
      end
c***************************************************************
      SUBROUTINE XENERGY(vcon,x,y,z,vx,vy,vz,e)
      IMPLICIT NONE
      REAL*8 x,y,z,vx,vy,vz,e,vcon,g
      PARAMETER (g=4.498933261d-6)
      REAL*8 mb,rc,md,a,b,mv,c,ra
c     --- bulge parameters
      PARAMETER (mb=9.1d9,rc=2.1)
c     --- disk parameters
c     --- increase disk mass to get higher v_c at R=R_0.
      PARAMETER (a=4.9d0,b=0.30d0,md=1.35d11)
c     --- halo parameters
      PARAMETER (mv=1.16d12,c=10.0d0,ra=28.2d0)
c
      REAL*8 r,rd2,ke,phi1,phi2,phi3
c
c     -- The kinetic energy per unit mass
      ke=vcon*vcon*(vx*vx+vy*vy+vz*vz)/2.0d0
c
c     -- The potential energy per unit mass
      r  = sqrt(x*x+y*y+z*z)
      rd2 = x*x+y*y
c     -- disk
      phi1=-g*md/sqrt(rd2+(a+sqrt(z*z+b*b))**2)
c     -- bulge; this isn't right, but is accurate for r >> rc and that
c        is the region of interest for all or almost all of our range of
c        integration
      if (r.gt.rc/2.0d0) then
         phi2 = -g*mb/r
      else
         phi2 = -g*2.0d0*mb/rc
      end if
c     -- halo
      phi3 = -g*(mv/(log(1.0d0+c)-c/(1.0d0+c))) * log(1.0d0+r/ra)/r
c
      e=ke+phi1+phi2+phi3
c      write (6,*) ke,phi1,phi2,phi3,e
c
      return
      end
c***********************************************************************
c
      SUBROUTINE GETINITPOS(dlsr,lg,bg,dg,x0,y0,z0,
     &                         vpi,vpiu,vth,vthu,vx0,vx0u,vy0,vy0u)
c     --- Given the Galactic coordinates of the galaxy, its keliocentric
c         distance, and the Galactocentric distance of the LSR, the
c         subroutine determines the initial position (x0,y0,z0) of the
c         galaxy in a coordinate system centered on the Milky Way.
c         This version also finds vx and vy from vpi and vtheta.
c         The Coordinate System (X,Y,Z)
c         --- X axis is in the Galactic plane and the Sun is at +dlsr
c             along it today.
c         --- Y axis is in the Galactic plane and forms a right-handed
c             coordinate system with X.  The Sun's galactic rotation is
c             towards -Y.
c         --- Z axis is perpendicular to the Galactic plane and points
c             to the N Galactic Pole.
      IMPLICIT NONE
      REAL*8 dlsr,lg,bg,dg,x0,y0,z0
      REAL*8 vpi,vth,vx0,vy0,costh,sinth,vpiu,vthu,vx0u,vy0u
      REAL*8 d2r
      parameter (d2r=1.74532925199D-02)
      x0=dlsr-dg*cos(d2r*bg)*cos(d2r*lg)
      y0=-dg*cos(d2r*bg)*sin(d2r*lg)
      z0=dg*sin(d2r*bg)
c
c     -- determine vx and vy
      costh=x0/sqrt(x0*x0+y0*y0)
      sinth=-y0/sqrt(x0*x0+y0*y0)
      vx0=vpi*costh-vth*sinth
      vx0u=sqrt((vpiu*costh)**2+(vthu*sinth)**2)
      vy0=-vpi*sinth-vth*costh
      vy0u=sqrt((vpiu*sinth)**2+(vthu*costh)**2)
      return
      end
c**********************************************************************
      SUBROUTINE FINDIDORB(NC,NS,NSTEP,r2,nap,idap,npe,idpe)
      IMPLICIT NONE
      INTEGER NC,ns,NSTEP
      REAL*8 r2(NS)
      INTEGER nap,idap(NC),npe,idpe(NC)
      INTEGER i,kp,ka
      NAP=0
      NPE=0
      ka=1
      kp=1
      do i=2,NSTEP-1
         if((r2(i-1).gt.r2(i)).and.(r2(i+1).gt.r2(i))) then
            idpe(kp)=i
            kp=kp+1
            NPE=NPE+1
            if (NPE.gt.NC) then
               write (6,*) 'overflowing idpe'
               stop
            end if
         end if
         if((r2(i-1).lt.r2(i)).and.(r2(i+1).lt.r2(i))) then
            idap(ka)=i
            ka=ka+1
            NAP=NAP+1
            if (NAP.gt.NC) then
               write (6,*) 'overflowing idap'
               stop
            end if
         end if
      end do
      return
      end

      SUBROUTINE STATISTICS(x,n,m,ave,sdev)
      IMPLICIT NONE
      INTEGER N,M,i
      REAL*8 x(n),ave,sdev
c     --- Calculate the mean.
      ave=0.0
      do i=1,m
         ave=ave+x(i)
      end do
      ave=ave/real(m)
c     --- Calculate the standard deviation.
      sdev=0.0
      do i=1,m
         sdev=sdev+(x(i)-ave)*(x(i)-ave)
      end do
      sdev=sdev/real(m-1)
      sdev=sqrt(sdev)
      return
      end

      SUBROUTINE FIND95(x,NN,NCOUNT,x95)
      IMPLICIT NONE
      REAL*8 pnf
      PARAMETER (pnf=0.95d0)
      REAL*8 x(NN),x95(2)
      REAL s(NN)
      INTEGER NN,NCOUNT,i,index(2)
      do i=1,NCOUNT
         s(i)=x(i)
      end do
      call sort(NCOUNT,s)
      index(1)=nint(dble(NCOUNT)*((1.00d0-pnf)/2.0d0))
      index(2)=nint(dble(NCOUNT)*(pnf+((1.00d0-pnf)/2.0d0)))
      x95(1)=dble(s(index(1)))
      x95(2)=dble(s(index(2)))
      return
      end

      SUBROUTINE FINDSMALLER(x,value,NN,NCOUNT,frac)
c     --- Finds the fraction of values in x that are smaller than value
      IMPLICIT NONE
      REAL*8 pnf
      PARAMETER (pnf=0.95d0)
      REAL*8 x(NN),frac,value
      REAL s(NN)
      INTEGER NN,NCOUNT,i
      do i=1,NCOUNT
         s(i)=x(i)
      end do
      call sort(NCOUNT,s)
      i=NCOUNT
      do while ((i.gt.0).and.(s(i).ge.value))
         i=i-1
      end do
      frac=dble(i)/dble(NCOUNT)
      return
      end

      SUBROUTINE MAX1(x,N,NP,xmx)
      IMPLICIT NONE
      INTEGER N,NP,i
      REAL x(N),xmx
      xmx=-9.9d31
      do i=1,NP
         if(x(i).gt.xmx) xmx=x(i)
      end do
      return
      end

      SUBROUTINE MIN1(x,N,NP,xmn)
      IMPLICIT NONE
      INTEGER N,NP,i
      REAL x(N),xmn
      xmn=9.9d31
      do i=1,NP
         if(x(i).lt.xmn) xmn=x(i)
      end do
      return
      end

      SUBROUTINE MAX1D(x,N,NP,xmx)
      IMPLICIT NONE
      INTEGER N,NP,i
      REAL*8 x(N),xmx
      xmx=-9.9d31
      do i=1,NP
         if(x(i).gt.xmx) xmx=x(i)
      end do
      return
      end

      SUBROUTINE MIN1D(x,N,NP,xmn)
      IMPLICIT NONE
      INTEGER N,NP,i
      REAL*8 x(N),xmn
      xmn=9.9d31
      do i=1,NP
         if(x(i).lt.xmn) xmn=x(i)
      end do
      return
      end
c
c     *** Numerical Recipes routines for incomplete gamma function
c
c     --- Returns "upper" incomplete gamma function 1 - Gamma(a,x)
      FUNCTION gammq(a,x)
      IMPLICIT NONE
      REAL*8 a,gammq,x
c     --- uses gcf and gser
      REAL*8 gammcf,gamser,gln
      if ((x.lt.0.d0).or.(a.le.0.)) then
         write (6,*) 'bad arguments in gammq: a,x=',a,x
         stop
      end if
      if (x.lt.a+1.0d0) then
c        --- Use the series representation
         call gser(gamser,a,x,gln)
c        --- and take its complement
         gammq=1.0d0-gamser
      else
c        --- use the continued fraction representation
         call gcf(gammcf,a,x,gln)
         gammq=gammcf
      end if
      return
      end
c
      SUBROUTINE GSER(GAMSER,A,X,GLN)
C
C     RETURNS THE INCOMPLETE GAMMA FUNCTION P(A,X) EVALUTATED BY ITS SERIES
C     REPRESENTATION AS GAMSER.  ALSO RETURNS LN(GAMMA(A)) AS GLN.
c     Uses gammln.
C
      IMPLICIT NONE
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
      INTEGER n
      REAL*8 ap,del,sum,gammln
c
      GLN=GAMMLN(A)
      IF (X.LE.0.d0) THEN
        IF (X.LT.0.d0) pause 'x < 0 in gser'
        GAMSER=0.d0
        RETURN
      END IF
      AP=A
      SUM=1.0d0/A
      DEL=SUM
      DO N=1,ITMAX
        AP=AP+1.0d0
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF (ABS(DEL).LT.ABS(SUM)*EPS) GOTO 1
      end do
      PAUSE 'A TOO LARGE, ITMAX TOO SMALL'
1     continue
      GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
C
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
C
C     RETURNS THE INCOMPLETE GAMMA FUNCTION Z(A,X) EVALUATED BY ITS CONTINUED
C     FRACTION REPRESENTATION AS GAMMCF.  ALSO RETURNS LN(GAMMA(A)) AS GLN.
c     Parameters: ITMAX is the maximum allowed number of iterations; EPS is the
c     relative accuracy; FPMIN is a number near the smallest representable floating
c     point number.
c     Calls gammln
C
      IMPLICIT NONE
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100, EPS=3.d-7, FPMIN=1.0d-30)
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
c
      GLN=GAMMLN(A)
c
c     --- Set up for evaluating the continued fraction by modified Lentz's
c         method with b_0=0.
      b=x+1.0d0-a
      c=1.0d0/FPMIN
      d=1.0d0/b
      h=d
      do i=1,ITMAX
         an=-i*(i-a)
         b=b+2.0d0
         d=an*d+b
         if (abs(d).lt.FPMIN) d=FPMIN
         c=b+an/c
         if (abs(c).lt.FPMIN) c=FPMIN
         d=1.0d0/d
         del=d*c
         h=h*del
         if (abs(del-1.0d0).lt.EPS) goto 1
      end do
      pause 'a too large, ITMAX too small in gcf'
1     continue
C      -- put factors in front
      gammcf=exp(-x+a*log(x)-gln)*h
      return
      end
c
c     --- Returns the value of ln(Gamma(xx)) for xx > 0.
c
      FUNCTION gammln(xx)
      IMPLICIT NONE
      REAL*8 gammln,xx
      INTEGER j
      REAL*8 ser,stp,tmp,x,y,cof(6)
c     --- Internal arithmetic will be done in double precision, a nicety
c         that you can omit if five-figure accuracy is good enough.
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     &   24.01409824083091d0,-1.231739572450155d0,0.1208650973866179d-2,
     &   -0.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
         y=y+1.0d0
         ser=ser+cof(j)/y
      end do
      gammln=tmp+log(stp*ser/x)
      return
      end
