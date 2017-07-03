c     --- Started writing on October 3, 2002.
c     --- 06-29-2007 Changed so that a user can choose a type of
c         potential for the Milky Way.
c     --- 06-07-2014 Made dlsr = 8.0 kpc
c     --- 06-13-2016 Added the Bovy MWPotential2014 (type B).
c                    Changed vcon to the inverse of the value that was
c                    being used, which I think is correct.
c     --- 17-Jun-2016: Added other potentials (see the discussions in
c             the subroutines.  Also changed dlsr to 8.2 kpc (Bland-Hawthorn
c             & Gebhard 2016).
c         ORBIT.F integrates an orbit in a given Galactic potential.
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
c     --- vx0    Initial velocity in X direction of the galaxy in km/s
c     --- vy0    Initial velocity in Y direction of the galaxy in km/s
c     --- vz0    Initial velocity in Z direction of the galaxy in km/s
c     --- vx     Velocity x at time t
c     --- vy     Velocity x at time t 
c     --- vz     Velocity x at time t
c     --- vcon   Conversion factor from kpc/Gyr to km/s
c     --- Supportet potentials
c     J -- From Johnston et al. 1995 (JFORCE,JENERGY)
c     D -- From Dehnen et al. 2006 (DFORCE,DENERGY)
c     N -- Is an NFW model from somewhere.
c     B -- Is the MWPotential2014 model from Bovy 2015 ApJS, 216:29
      IMPLICIT NONE
      REAL*8 bg,lg,dg,x0,y0,z0
      REAL*8 vx0,vy0,vz0,vpi,vth
      REAL*8 dlsr,vcon
ccc      PARAMETER (dlsr=8.5d0)
      PARAMETER (dlsr=8.2d0)
      PARAMETER (vcon=1.0226903)
ccc      PARAMETER (vcon=0.977786932d0)
      REAL*8 fx,fy,fz,x,y,z,vx,vy,vz,e0,ef
      REAL*8 t,dt,ts
      REAL R_E,R,E_T,dR,R_CE,v_co,v_toto,v_c2
      REAL R_C,V_C,F_RX,F_RY,F_RZ,factor,dV,C_o,J_o,J_c
      CHARACTER tdir*1,pot*1
      INTEGER NSTEP,i
      OPEN(1,FILE='orbit.par',STATUS='OLD',form='FORMATTED')
      read(1,'(a1)') pot ! Read a type of potential to use [J,D,N,B,K,V,X]
      read(1,*) lg ! Read in Galactic longitude [degrees].
      read(1,*) bg ! Read in Galactic latitude [degrees].
      read(1,*) dg ! Read in the heliocentric distance to the galaxy [kpc].
      read(1,'(a1)') tdir  !F=forward in time; B=backward in time
      read(1,*) vpi ! Read in initial velocity in pi direction
      read(1,*) vth ! Read in initial velocity in theta direction
      read(1,*) vz0 ! Read in initial velocity in Z direction
      read(1,*) t               ! Total integration time in Gyr
      read(1,*) NSTEP ! Number of integration steps
      read(1,*) R_e   ! Read in radius at which to calculate energy
      close(unit=1)
      if(tdir(1:1).ne.'F'.and.tdir(1:1).ne.'B') then
         write(6,*) 'Error: Specify direction of time with F or B'
         stop
      end if
      if(tdir(1:1).eq.'B') then
         vpi = -vpi
         vth = -vth
         vz0 = -vz0
         write(6,*)'Integration backward in time; time[Gyr] =',sngl(t)
      else if(tdir(1:1).eq.'F') then
         write(6,*)'Integration forward in time; time[Gyr] =',sngl(t)
      end if
      CALL GETINITPOS(dlsr,lg,bg,dg,x0,y0,z0,vpi,vth,vx0,vy0)
      x=x0
      y=y0
      z=z0
      vx=vx0
      vy=vy0
      vz=vz0
c
c     -- find the initial energy
      if(pot.eq.'J') then
         call JENERGY(vcon,x,y,z,vx,vy,vz,e0)
      else if(pot.eq.'B') then
         CALL BENERGY(vcon,x,y,z,vx,vy,vz,e0)
      else if(pot.eq.'N') then
         CALL NENERGY(vcon,x,y,z,vx,vy,vz,e0)
      else if(pot.eq.'D') then
         CALL DENERGY(vcon,x,y,z,vx,vy,vz,e0)
      else if(pot.eq.'K') then
         CALL KENERGY(vcon,x,y,z,vx,vy,vz,e0)
      else if(pot.eq.'V') then
         CALL VENERGY(vcon,x,y,z,vx,vy,vz,e0)
      else if(pot.eq.'X') then
         CALL XENERGY(vcon,x,y,z,vx,vy,vz,e0)
      end if
c
      dt=t/dble(NSTEP) ! Time-step size
      ts=0.0d0
      if(pot.eq.'K') then
         CALL KFORCE(x,y,z,fx,fy,fz)
      else if(pot.eq.'X') then
         CALL XFORCE(x,y,z,fx,fy,fz)
      else if(pot.eq.'V') then
         CALL VFORCE(x,y,z,fx,fy,fz)
      else if(pot.eq.'J') then
         CALL JFORCE(x,y,z,fx,fy,fz)
      else if(pot.eq.'D') then
         CALL DFORCE(x,y,z,fx,fy,fz)
      else if(pot.eq.'N') then
         CALL NFORCE(x,y,z,fx,fy,fz)
      else if(pot.eq.'B') then
         CALL BFORCE(x,y,z,fx,fy,fz)
      end if
c
c     -- open the output file
      if(pot.eq.'J') then
         OPEN(2,FILE='orbit_J.res',STATUS='unknown',form='FORMATTED')
      else if(pot.eq.'D') then
         OPEN(2,FILE='orbit_D.res',STATUS='unknown',form='FORMATTED')
      else if(pot.eq.'N') then
         OPEN(2,FILE='orbit_N.res',STATUS='unknown',form='FORMATTED')
      else if(pot.eq.'B') then
         OPEN(2,FILE='orbit_B.res',STATUS='unknown',form='FORMATTED')
      else if(pot.eq.'K') then
         OPEN(2,FILE='orbit_K.res',STATUS='unknown',form='FORMATTED')
      else if(pot.eq.'V') then
         OPEN(2,FILE='orbit_V.res',STATUS='unknown',form='FORMATTED')
      else if(pot.eq.'X') then
         OPEN(2,FILE='orbit_X.res',STATUS='unknown',form='FORMATTED')
      end if
c
      write(2,*) ts,x,y,z,vx,vy,vz
      dR=5 !kpc
      dV=100.0 ! km/s
      do i=1, NSTEP-1
         vx=vx+0.50d0*dt*fx/vcon
         vy=vy+0.50d0*dt*fy/vcon
         vz=vz+0.50d0*dt*fz/vcon
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
         else if(pot.eq.'D') then
            CALL DFORCE(x,y,z,fx,fy,fz)
         else if(pot.eq.'N') then
            CALL NFORCE(x,y,z,fx,fy,fz)
         else if(pot.eq.'B') then
            CALL BFORCE(x,y,z,fx,fy,fz)
         end if
         vx=vx+0.50d0*dt*fx/vcon
         vy=vy+0.50d0*dt*fy/vcon
         vz=vz+0.50d0*dt*fz/vcon
c     --- Find Energy
         R=sqrt(x*x + y*y + z*z)
         if(abs(R-R_E).lt.dR) then
            dR=abs(R-R_E)
            E_T = (vx*vx + vy*vy + vz*vz)/2.0
            R_CE=R
         end if
c     ---------------------------------------
c     --- Find circularity
         factor=sngl((fx*x + fy*y + fz*z))/(R*R)
         F_RX=factor*x
         F_RY=factor*y
         F_RZ=factor*z
         V_C2=sqrt(R*sqrt(FX**2 + FY**2 + FZ**2))/vcon
         V_C=sqrt(R*sqrt(F_RX**2 + F_RY**2 + F_RZ**2))/vcon
         if(abs(sqrt(vx**2 + vy**2 + vz**2) - V_C).lt.dV) then
            dV=abs(sqrt(vx**2 + vy**2 + vz**2) - V_C)
            J_o=(y*vz-z*vy)**2 + (z*vx - x*vz)**2 + (x*vy - y*vx)**2
            J_o=sqrt(J_o)
            J_c=R*sqrt(vx**2 + vy**2 + vz**2)
            C_o=J_o/J_c
            R_c=R
            V_co=V_C
            V_toto=sqrt(vx**2 + vy**2 + vz**2)
         end if
         write(2,*) ts,x,y,z,vx,vy,vz
c         write (6,*) r,v_c,x,y,z,fx,fy,fz,factor,f_rx
      end do
      close(unit=2)
      write(6,*)'Total Energy and Radii',E_T,R_E,R_CE
      write(6,*)
      write(6,*)'Circularity of the orbit:',C_o,J_o,J_c,dV
      write(6,*) 'R_c, V_c, V_tot: ',r_c,v_co,v_toto
c
c     -- find the final energy and the change
      if(pot.eq.'J') then
         call JENERGY(vcon,x,y,z,vx,vy,vz,ef)
      else if(pot.eq.'B') then
         call BENERGY(vcon,x,y,z,vx,vy,vz,ef)
      else if(pot.eq.'N') then
         call NENERGY(vcon,x,y,z,vx,vy,vz,ef)
      else if(pot.eq.'D') then
         call DENERGY(vcon,x,y,z,vx,vy,vz,ef)
      else if(pot.eq.'K') then
         CALL KENERGY(vcon,x,y,z,vx,vy,vz,ef)
      else if(pot.eq.'V') then
         CALL VENERGY(vcon,x,y,z,vx,vy,vz,ef)
      else if(pot.eq.'X') then
         CALL XENERGY(vcon,x,y,z,vx,vy,vz,ef)
      end if
      write (6,*) 'e0,ef,(e0-ef)/e0=',e0,ef,(e0-ef)/e0
      stop
      end
c**************************************************************
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
c***************************************************************
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
      PARAMETER (mv=8.0d11,c=15.3,ra=16.0)
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
      PARAMETER (mv=8.0d11,c=15.3,ra=16.0)
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
      SUBROUTINE GETINITPOS(dlsr,lg,bg,dg,x0,y0,z0,vpi,vth,vx0,vy0)
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
      REAL*8 dlsr,lg,bg,dg,x0,y0,z0
      REAL*8 vpi,vth,vx0,vy0,costh,sinth
      REAL*8 d2r
      parameter (d2r=1.74532925199D-02)
      x0=dlsr-dg*cos(d2r*bg)*cos(d2r*lg)
      y0=-dg*cos(d2r*bg)*sin(d2r*lg)
      z0=dg*sin(d2r*bg)
      write (6,*) lg,bg,x0,y0,z0
c
c     -- determine vx and vy
      costh=x0/sqrt(x0*x0+y0*y0)
      sinth=-y0/sqrt(x0*x0+y0*y0)
      vx0=vpi*costh-vth*sinth
      vy0=-vpi*sinth-vth*costh
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
