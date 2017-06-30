import numpy as np
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

#MC Parameters
NMC = 1000 #Number of MC experiments 
t = 100 #Gyr
NSTEP = 100000
pot = 'X'

"""
name = "Leo_II"
ra = 168.3700 #degrees
dec = 22.1517 #degrees
distance = 233. #kpc
mu_a = -6.9 #mas/cent
mu_d = -8.7 #mas/cent
rv = 78.0 #km/s
#Uncertainty in the proper motion components
#These are taken as the std. devs. for the Gaussian distributions for the MC
mu_a_std = 3.7  
mu_d_std = 3.9
rv_std = 0.5
#luminosity of the dSph and the min and max M/L
L = .74e6 #1.e3
#r_h = 176 #pc
"""
"""
def linear_to_angular_size(linear_size, distance):
    #Returns angular size in arcminutes
    linear_size=float(linear_size)
    theta = np.arctan(linear_size/distance) #Radians
    return theta * 3437.75 #Arcmin
"""

#rtlim=2.6
#rtlim = 30

rpval = 120. #kpc
mol1 = 1.
mol2 = 3.

name = str(raw_input("Name?:"))
ra = float(raw_input("RA (degrees)?:"))
dec = float(raw_input("DEC (degrees)?:"))
distance = float(raw_input("Distance (kpc)?:"))
mu_a = float(raw_input("mu_a (mas/cent)?:"))
mu_a_std = float(raw_input("mu_a uncertainty?:"))
mu_d = float(raw_input("mu_d (mas/cent)?:"))
mu_d_std = float(raw_input("mu_d uncertainty?:"))
rv = float(raw_input("Radial Velocity (km/s)?:"))
rv_std = float(raw_input("Radial Velocity uncertainty?:"))
L = float(raw_input("Luminosity (solar units)?:"))
rtlim = float(raw_input("r_h (arcmin)?:"))


print "\n Your values are:"
print "Name:", name
print "RA, DEC:", ra, dec
print "Distance:", distance
print "mu_a:", mu_a, "+-", mu_a_std
print "mu_d:", mu_d, "+-", mu_d_std
print "Radial Velocity:", rv, "+-", rv_std
print "Luminosity:", L
print "mol1, mol2:", mol1, mol2
print "rtlim:", rtlim
print "rpval:", rpval


from astropy import units as u
from astropy.coordinates import SkyCoord, Galactocentric, ICRS, Galactic
import gala.coordinates as gc
import sys, time, mpmath
import scipy.stats as st

#Constants
G = 4.498933261e-6 #Gravitational constant in kpc^3/(M_sun Gyr)
vcon = 1.0226903 #Conversion factor from kpc/Gyr to km/s (i.e. multiply by km/s to get kpc/Gyr)
asec_to_rad = 206264.8 #Converts arcseconds to radians
ratoam = 3437.747 #Converts radians to arcmin
year_to_sec = 3.155693e7 #Converts from years to seconds
pc_to_km = 3.08568e13 #Converts pc to km
year_to_day = 365.2422 #Converts years to days
usun = -10.0 #U-component of Sun's velocity w.r.t. LSR; Dehnen & Binney                  
vsun = 11.0 #V-component of Sun's velocity w.r.t. LSR; 1998, MNRAS, 
wsun = 7.0 #W-component of Sun's velocity w.r.t. LSR; 298, 387
vlsr = 237.0 #Circular velocity of LSR w.r.t. Galaxy center in km/s
dlsr = 8.2 #Distance of LSR from the Galaxy center in pc

#Functions to get Positions and Velocities

def sin_deg(x):
    #Sine for degree-valued args
    return np.sin(np.deg2rad(x))

def cos_deg(x):
    #Cosine for degree-valued args
    return np.cos(np.deg2rad(x))

def radec_to_lb(ra, dec):
    #Returns galactic coords (l, b) in degrees, given an object's RA/DEC in degrees
    coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    galactic = coords.galactic
    return galactic.l.deg, galactic.b.deg

def lb_to_galactic_xyz(l, b, d, dlsr=dlsr):
    #Takes l and b in degrees, d (distance to object) in kpc, dlsr in kpc
    x = (dlsr-d*cos_deg(b)*cos_deg(l))
    y = -d*cos_deg(b)*sin_deg(l)
    z = d*sin_deg(b)
    return x, y, z #kpc

def pm_radec_to_lb(ra, dec, d, mu_a, mu_d, rv):
    #Takes ra/dec in degrees, distance in pc, mu_a/mu_d in mas/cent, rv in km/s. Returns mu_l/mu_b in mas/cent
    c = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, distance=d*u.pc)
    pm = [mu_a, mu_d] * u.mas/(100*u.yr) 
    result = gc.pm_icrs_to_gal(c, pm).to(u.mas/(u.yr))*u.yr/u.mas
    mu_l = float(result[0]*100)
    mu_b = float(result[1]*100)
    
    return mu_l, mu_b #mas/cent

def pm_lb_to_galactic_vt(mu_a,mu_d,d,asec_to_rad=asec_to_rad,year_to_sec=year_to_sec,pc_to_km=pc_to_km):
    #mul, mub in mas/cent, distance d in kpc
    #Converts l,b pm into tangential velocity in l-direction and b-direction,
    #and gets magnitude vt of tangential velocity
    
    #Get from mas/cent into arcsec/year
    mu_a = (mu_a*(1e-3)/100)
    mu_d = (mu_d*(1e-3)/100)

    vtl=mu_a*d*1000*pc_to_km/(year_to_sec*asec_to_rad)
    vtb=mu_d*d*1000*pc_to_km/(year_to_sec*asec_to_rad)
    vt=np.sqrt(vtl*vtl+vtb*vtb)
    return vtl, vtb, vt  

def galactic_vt_to_v_wrt_lsr(l, b, vtl, vtb, v_r_sun=rv, usun=usun, vsun=vsun, wsun=wsun):
    #Finds velocity of galaxy wrt LSR. velocity inputs are in km/s
    #Theta direction
    vg = v_r_sun*cos_deg(b)*sin_deg(l) - vtb*sin_deg(b)*sin_deg(l)+vtl*cos_deg(l) + vsun
    #Pi direction
    ug = -v_r_sun*cos_deg(b)*cos_deg(l) + vtb*sin_deg(b)*cos_deg(l) + vtl*sin_deg(l) + usun
    #Z direction
    wg =  v_r_sun*sin_deg(b) + vtb*cos_deg(b) + wsun
    
    return ug, vg, wg #km/s

def v_uvg_to_cyl(l,b,d,ug,vg,wg, dlsr=dlsr, vlsr=vlsr):
    #input velocities in km/s, distances in kpc
    x=np.sqrt(dlsr*dlsr+d*d*cos_deg(b)*cos_deg(b)-2.0*dlsr*d*cos_deg(b)*cos_deg(l))
    sinalpha=d*cos_deg(b)*sin_deg(l)/x
    cosalpha=(x*x+dlsr*dlsr-d*d*cos_deg(b)*cos_deg(b))/(2.0*x*dlsr)

    vpi =  ug*cosalpha + (vg+vlsr)*sinalpha
    vtheta = -ug*sinalpha + (vg+vlsr)*cosalpha
    vz =  wg
    return vpi, vtheta, vz #km/s

def v_cyl_to_xyz(l, b, d, vpi, vth, vz, dlsr=dlsr):
    #distances in kpc, vels in km/s

    x=dlsr-d*cos_deg(b)*cos_deg(l)
    y=-d*cos_deg(b)*sin_deg(l)
    z=d*sin_deg(b)

    costh=x/np.sqrt(x**2+y**2)
    sinth=-y/np.sqrt(x**2+y**2)
    vx=vpi*costh-vth*sinth
    vy=-vpi*sinth-vth*costh
    return vx, vy, vz

def get_v_spherical(l, b, d, vpi, vtheta, vz, dlsr=dlsr):

    x=np.sqrt(dlsr*dlsr+d*d*cos_deg(b)*cos_deg(b)-2.0*dlsr*d*cos_deg(b)*cos_deg(l))
    cosbeta=x/np.sqrt(x*x+d*d*sin_deg(b)*sin_deg(b))
    sinbeta=d*sin_deg(b)/np.sqrt(x*x+d*d*sin_deg(b)*sin_deg(b))
    vr =  vpi*cosbeta + vz*sinbeta
    vrot = vtheta
    v_perp_radial = -vpi*sinbeta + vz*cosbeta
    vt = np.sqrt(vrot**2+v_perp_radial**2)
    
    return vr, vrot, v_perp_radial, vt

def get_vxvyvz(ra, dec, distance, mu_a, mu_d, rv):
    l, b = radec_to_lb(ra, dec)
    x, y, z = lb_to_galactic_xyz(l, b, distance)
    
    mu_l, mu_b = pm_radec_to_lb(ra, dec, distance, mu_a, mu_d, rv)
    vtl, vtb, vt = pm_lb_to_galactic_vt(mu_l, mu_b, distance)    
    ug, vg, wg = galactic_vt_to_v_wrt_lsr(l, b, vtl, vtb, rv)
    vpi, vtheta, vz = v_uvg_to_cyl(l, b, distance, ug, vg, wg)
    vx, vy, vz = v_cyl_to_xyz(l, b, distance, vpi, vtheta, vz)
    
    return vx, vy, vz

def get_xyz(ra, dec, distance):
    l, b = radec_to_lb(ra, dec)
    x, y, z = lb_to_galactic_xyz(l, b, distance)
    return x, y, z

#~~~~~~~~~~~~~~~~~~~~Energy and Force functions for X potential~~~~~~~~~~~~~~~~~~~~

def gammq(a,x):
    #a-1 is power inside the integral
    #x = lower bound of integration
    #mpmath.gammainc defines them differently,
    #so we make the switch when calling it herein
    return float(mpmath.gammainc(z=a, a=x, regularized=True))

def XENERGY(x, y, z, vx, vy, vz, vcon=vcon, G=G):
    #bulge parameters
    mb=9.1e9
    rc=2.1
    #disk parameters
    #increase disk mass to get higher v_c at R=R_0.
    a=4.9e0
    b=0.30e0
    md=1.35e11
    #halo parameters
    mv=1.16e12
    c=10.0
    ra=28.2

    #The kinetic energy per unit mass
    ke=vcon*vcon*(vx*vx+vy*vy+vz*vz)/2.0
    
    #The potential energy per unit mass
    r  = np.sqrt(x*x+y*y+z*z)
    rd2 = x*x+y*y
    
    #disk
    phi1=-G*md/np.sqrt(rd2+(a+np.sqrt(z*z+b*b))**2)

    #bulge; this isn't right, but is accurate for r >> rc and that
    #is the region of interest for all or almost all of our range of
    #integration
    if r > rc/2.0:
        phi2 = -G*mb/r
    else:
        phi2 = -G*2.0*mb/rc
        
    #halo
    phi3 = -G*(mv/(np.log(1.0+c)-c/(1.0+c))) * np.log(1.0+r/ra)/r

    e=ke+phi1+phi2+phi3    
    return e

def XFORCE(x, y, z, G=G):
    #bulge parameters
    mb=9.1e9
    rc=2.1
    #disk parameters
    #increase disk mass to get higher v_c at R=R_0.
    a=4.9
    b=0.30
    md=1.35e11
    
    #halo parameters
    mv=1.16e12
    c=10.0
    ra=28.2
    
    r  = np.sqrt(x*x+y*y+z*z)
    rd2 = x*x+y*y
    
    #bulge radial force: -GM_b(r)/r*r
    mb_r = mb*(1.-1.0891*np.exp(-r/rc)*(r/rc)**.2 - gammq(.2, r/rc))
    #mb_r=mb*(1.00 - 1.0891*np.exp(-r/rc)*(r/rc)**.2 - gammq(.2,r/rc))
    fb_r = -G*mb_r/(r*r)

    #NFW halo radial force: -GM_nfw(r)/r*r
    mnfw_r=mv*(np.log(1.0+r/ra)-r/(ra+r))/(np.log(1.0+c)-c/(1.0+c))
    fh_r = -G*mnfw_r/(r*r)

    #Find the X components of the force.
    fx1=-G*md*x/(np.sqrt(rd2+(a+np.sqrt(z*z+b*b))**2))**3
    fx2=x*fb_r/r
    fx3=x*fh_r/r
    fx=fx1+fx2+fx3

    #Find the Y components of the force.
    fy1=-G*md*y/(np.sqrt(rd2+(a+np.sqrt(z*z+b*b))**2))**3
    fy2=y*fb_r/r
    fy3=y*fh_r/r
    fy=fy1+fy2+fy3
    
    #Find the Z components of the force.
    fz1=-G*md*z/(np.sqrt(rd2+(a+np.sqrt(z*z+b*b))**2))**3
    fz1=fz1*(1.0+a/np.sqrt(z*z+b*b))
    fz2=z*fb_r/r
    fz3=z*fh_r/r
    fz=fz1+fz2+fz3

    return fx, fy, fz

#Create dictionaries which store references to force funtions and energy functions. 
#This gives easy way to get fx, fy, fz, energy for any potential by changing a single parameter. 
#If we add later add force and energy functions for different potentials
#we can simply add their names to the dictionaries below and pass the
#'pot' argument in integrate_orbit() directly to these dictionaries.

forcefunctions = {'X':XFORCE} 
energyfunctions = {'X':XENERGY}

def integrate_orbit(pot, x0, y0, z0, vx0, vy0, vz0, vcon=vcon, t=t, NSTEP=NSTEP):
    
    dt = float(t)/float(NSTEP)
    
    #Initial angular momenta
    amx0=y0*vz0-z0*vy0
    amy0=z0*vx0-x0*vz0
    amz0=x0*vy0-y0*vx0
    am0=np.sqrt(amx0*amx0+amy0*amy0+amz0*amz0)
    bam0=np.arcsin(amz0/am0)
    lam0=np.arctan2(-amy0,-amx0)
    if lam0 < 0.:
        lam0 += 2.*np.pi
    if lam0 > 2.*np.pi:
        lam0 -= 2.*np.pi

    #Initial energy
    initialenergy = energyfunctions['X'](x0, y0, z0, vx0, vy0, vz0)
    
    dt = float(t)/float(NSTEP)
    x=x0
    y=y0
    z=z0
    r2=x**2+y**2+z**2
    vx=vx0
    vy=vy0
    vz=vz0
    
    times, positions = [[] for i in np.arange(2)] 
    timestep = 0
    times.append(timestep)
    positions.append((x,y,z,r2))

    fx, fy, fz = forcefunctions['X'](x, y, z)
    for i in np.arange(0, NSTEP): #-1)?
        timestep = timestep + dt
        vx=vx+0.5*dt*fx/vcon
        vy=vy+0.5*dt*fy/vcon
        vz=vz+0.5*dt*fz/vcon
        x=x+dt*vx*vcon
        y=y+dt*vy*vcon
        z=z+dt*vz*vcon
        r2 = x**2 + y**2 + z**2
        
        fx, fy, fz = forcefunctions['X'](x, y, z)
        vx=vx+0.5*dt*fx/vcon
        vy=vy+0.5*dt*fy/vcon
        vz=vz+0.5*dt*fz/vcon
        
        times.append(timestep)
        positions.append((x,y,z,r2))
        sys.stdout.flush()
        sys.stdout.write('\r'+'Done {0} steps.' .format(i))
    
    sys.stdout.write('\r' + ' '*50)
    #sys.stdout.flush()
    
    #Final angular momenta
    amx=y*vz-z*vy
    amy=z*vx-x*vz
    amz=x*vy-y*vx
    am=np.sqrt(amx*amx+amy*amy+amz*amz)
    bamf = np.arcsin(amz/am)
    dbam = bamf - bam0

    lamf=np.arctan2(-amy,-amx)
    if lamf < 0.:
        lamf += 2.*np.pi
    if lamf > 2.*np.pi:
        lamf -= 2.*np.pi
    dlam=lamf-lam0

    #Final energy
    finalenergy = energyfunctions['X'](x, y, z, vx, vy, vz)
        
    return {'positions':np.array(positions),
            'initialenergy':np.array(initialenergy),
            'finalenergy':np.array(finalenergy), 
            'times':np.array(times), 
            'lam0':np.rad2deg(lam0),
            'bam0':np.rad2deg(bam0),
            'dlam':np.rad2deg(dlam),
            'dbam':np.rad2deg(dbam)}

#Functions to get Orbital Parameters
def find_apo_peri(r2):
    id_pe = [i for i in np.arange(1, len(r2)-1) if r2[i-1] > r2[i] and r2[i+1] > r2[i]]
    id_ap = [i for i in np.arange(1, len(r2)-1) if r2[i-1] < r2[i] and r2[i+1] < r2[i]]
       
    apos = [r2[i] for i in id_ap]
    peris = [r2[i] for i in id_pe]

    if len(apos) >= 1:
        apo_mean = np.sqrt(np.mean(apos))
    else:
        apo_mean = None

    if len(peris) >= 1:
        peri_mean = np.sqrt(np.mean(peris))
    else:
        peri_mean = None

    return id_pe, id_ap, peri_mean, apo_mean

def find_ecc(r2, id_pe, id_ap):
    N_pe = len(id_pe)
    N_ap = len(id_ap)
    eccentricities = []
    for i in np.arange(0, min(N_pe,N_ap)):   
        e = (np.sqrt(r2[id_ap[i]]) - np.sqrt(r2[id_pe[i]]))/(np.sqrt(r2[id_ap[i]])+np.sqrt(r2[id_pe[i]]))
        eccentricities.append(e)
    return np.mean(eccentricities)

def find_period(times, id_ap): 
    #Find average period from R_a to R_a 
    periods = []
    for i in np.arange(1, len(id_ap)):
        period = times[id_ap[i]]-times[id_ap[i-1]]
        periods.append(period)
    return np.mean(periods) #Gyr

def pick3points(idpe,idap,index,x,y,z):
    #x, y, z are arrays, i.e. time series of x, y, z values
    if idpe[index] < idap[index]:
        #perigalacticon occurs first
        #Peri first
        point1 = (x[idpe[index]], y[idpe[index]], z[idpe[index]])
        
        #Apogalacticon last
        point3 = (x[idap[index]], y[idap[index]], z[idap[index]])
        
    else:
        #apogalacticon occurs first
        #Apogalacticon first
        point1 = (x[idap[index]], y[idap[index]], z[idap[index]])
        
        #Perigalacticon last
        point3 = (x[idpe[index]], y[idpe[index]], z[idpe[index]])

    #Point between the two
    indexbetween = (idpe[index]+idap[index])/2
    point2 = (x[indexbetween], y[indexbetween], z[indexbetween])
    return [point1, point2, point3]

def FINDPLANE(points): #points is list of three 3-tuples
    point1, point2, point3 = points
    
    #Vector p1p2 goes from p1 to p2, vector p1p3 goes from p1 to p3
    #The points are ordered in time, so galaxy reaches p1, then p2, then p3
    #Cross product p1p2 X p1p3 thus yields normal vector to 
    #orbital plane, this vector being parallel to the angular velocity vector, 
    #by right hand rule
    
    p1p2 = (point2[0]-point1[0], point2[1]-point1[1], point2[2]-point1[2]) 
    p1p3 = (point3[0]-point1[0], point3[1]-point1[1], point3[2]-point1[2])
    normalvector = np.cross(p1p2, p1p3)
    return normalvector

def FINDINC(normalvector):
    #Normal vector is assumed to be parallel to angular velocity vector
    a, b, c = normalvector 
    
    #Inclination angle of orbit plane to galactic plane.
    #Zero is in galactic plane and orbiting in the direction of galactic rotation
    theta=np.arccos(-c/np.sqrt(a*a+b*b+c*c))
    return theta #radians!

def FINDLONGITUDE(normalvector):     
    #input is normal vector to orbital plane, parallel to angular velocity vector
    #longitude of ascending node measured in direction opposite galactic rotation from x=0
    a = normalvector[0]
    b = normalvector[1]
    omega=np.arctan2(b,a)+np.pi/2.
    if omega < 0:
        omega += 2.*np.pi
    if omega > 2.*np.pi:
        omega -= 2.*np.pi
    return omega #radians

def find_lpole_bpole(theta, omega):
    bpole=theta-(np.pi/2.)
    lpole=omega+0.5*np.pi
    lpole = np.rad2deg(lpole)
    bpole = np.rad2deg(bpole)
    if lpole >= 0:
        lpole = lpole % 360
    else:
        lpole = lpole % -360
    if bpole >= 0:
        bpole = bpole % 360
    else:
        bpole = bpole % -360
    return lpole, bpole

def get_rt(r2, idpe, idap, L=L, mol1=mol1, mol2=mol2, d=distance, ratoam=ratoam):
    #Calculate the tidal radius at pericenter and store it
    N_pe = len(idpe)
    N_ap = len(idap)
    rt1 = []
    rt2 = []
    for i in np.arange(0, min(N_pe, N_ap)): #CHANGED FROM np.arange(1, min...)!!!!
        r_ap = np.sqrt(r2[idap[i]])    
        r_pe=np.sqrt(r2[idpe[i]])
        ec=(r_ap-r_pe)/(r_ap+r_pe)
        fe=((1.-ec)**2)/((((1.+ec)**2)/(2.*ec))*np.log((1.+ec)/(1.-ec))+1.)
        mmw=1.1e10*(r_ap+r_pe)/2.
        rt1.append(((r_ap+r_pe)/2.)*(fe*mol1*L/mmw)**0.333333)
        rt2.append(((r_ap+r_pe)/2.)*(fe*mol2*L/mmw)**0.333333)
        
    rt1=np.mean(rt1)
    rt2=np.mean(rt2)
    rt1=np.arctan(rt1/d)*ratoam
    rt2=np.arctan(rt2/d)*ratoam
    
    return rt1, rt2

#~~~~~~~~~~~~~~~~~~~~Main Code~~~~~~~~~~~~~~~~~~~~

print "Constructing sample of {0} proper motions from Gaussian distributions of mu_a, mu_d, radial velocity..." \
.format(NMC)
mu_a_sample = np.random.normal(loc=mu_a, scale=mu_a_std, size=NMC-1)
mu_d_sample = np.random.normal(loc=mu_d, scale=mu_d_std, size=NMC-1)
rv_sample = np.random.normal(loc=rv, scale=rv_std, size=NMC-1)

pm_sample = np.column_stack([mu_a_sample, mu_d_sample, rv_sample])
pm_sample = np.insert(pm_sample, 0, (mu_a, mu_d, rv), axis=0) 
#Now the first value is a tuple with the measured or "actual" mu_a, mu_d, rv

print "Converting to sample of xyz velocities..."
velsample = [get_vxvyvz(ra, dec, distance, pm[0], pm[1], pm[2]) for pm in pm_sample]

#Make bunch of empty lists for our orbital parameters
apos, peris, eccs, periods, thetas, omegas, lpoles, bpoles, \
lam0s, bam0s, lamfs, bamfs, dlams, dbams, rt1s, rt2s, deltaEs = [[] for i in np.arange(17)]

Nfailed = 0
experiment = 0

#Get initial position
x0, y0, z0 = get_xyz(ra, dec, distance)

print "Commencing MC experiments..."
for velocityvector in velsample:
    experiment += 1
    
    sys.stdout.write('\r'+'Doing experiment {0}...' .format(experiment))
    sys.stdout.flush()
    time.sleep(.5)
    sys.stdout.write('\r'+' '*50)

    vx0, vy0, vz0 = velocityvector
    result = integrate_orbit(pot, x0, y0, z0, vx0, vy0, vz0, t=t, NSTEP=NSTEP)
    
    positions = result['positions']
    times = result['times']
    x = [pos[0] for pos in positions]
    y = [pos[1] for pos in positions]
    z = [pos[2] for pos in positions]
    r2 = [pos[3] for pos in positions]

    #Get apo and peri
    id_pe, id_ap, peri_mean, apo_mean = find_apo_peri(r2)
    N_pe = len(id_pe)
    N_ap = len(id_ap)
    if N_ap < 2:
        #If not at least one complete orbit from apo to apo, can't get parameters
        Nfailed += 1
        print "\r Experiment {0} failed." .format(experiment)
        continue
    apos.append(apo_mean)
    peris.append(peri_mean)
     
    #Get eccentr
    ecc = find_ecc(r2, id_pe, id_ap)
    eccs.append(ecc)
    
    #Get period
    avg_period = find_period(times, id_ap)
    periods.append(avg_period)
    
    #Get Inclination and Longitude
    setsofpoints = [pick3points(id_pe, id_ap, i, x, y, z) for i in np.arange(0,min(N_pe,N_ap))] #Changed from np.range(1, min...)!!!!
    inclinations = [FINDINC(FINDPLANE(pointset)) for pointset in setsofpoints]
    longitudes = [FINDLONGITUDE(FINDPLANE(pointset)) for pointset in setsofpoints]
    avg_theta = np.mean(inclinations)
    avg_omega = np.mean(longitudes)
    thetas.append(np.rad2deg(avg_theta))
    omegas.append(np.rad2deg(avg_omega))
    
    #SAME THING AS ABOVE
    #for i in np.arange(1,min(N_pe,N_ap)):
    #points = pick3points(id_pe, id_ap, i, test_x, test_y, test_z)
    #normalvector = FINDPLANE(points)       
    #theta = FINDINC(normalvector)
    #omega = FINDLONGITUDE(normalvector)
    #print normalvector, theta, omega
    
    #Get lpole, bpole
    lpole, bpole = find_lpole_bpole(avg_theta, avg_omega)
    lpoles.append(lpole)
    bpoles.append(bpole)
    
    #Get initial angular momentum
    lam0 = float(result['lam0'])
    bam0 = result['bam0']
    lam0s.append(lam0)
    bam0s.append(bam0)
    
    #Get rt1, rt2
    rt1, rt2 = get_rt(r2, id_pe, id_ap) 
    rt1s.append(rt1)
    rt2s.append(rt2)
    #print "rt1, 2 =", rt1, rt2
    
    #dlam, dbam
    dlams.append(result['dlam'])
    dbams.append(result['dbam'])

    #Get initial energy E0, final energy Ef
    E0 = result['initialenergy']
    Ef = result['finalenergy']
    #deltaEs.append(np.abs((E0-Ef)/E0))
    deltaEs.append((E0-Ef)/E0)

if Nfailed == NMC:
    print "All experiments failed to get a complete orbit."
    quit()

#Initial galactic coordinates
l, b = radec_to_lb(ra, dec)

#Initial velocity
vx0, vy0, vz0 = velsample[0]

#Initial angular momentum
amx0=y0*vz0-z0*vy0
amy0=z0*vx0-x0*vz0
amz0=x0*vy0-y0*vx0

print

print "Input Galaxy Parameters:"
print "Name \t\t\t\t\t:", name
print "RA, DEC (degrees) \t\t\t:", ra, dec 
print "Heliocentric Distance (kpc) \t\t:", distance
print "mu_alpha (mas/century) \t\t\t:", mu_a, "+-", mu_a_std
print "mu_delta (mas/century) \t\t\t:", mu_d, "+-", mu_d_std
print "Heliocentric Radial Velocity (km/s) \t:", rv, "+-", rv_std
print "Luminosity (solar luminosities) \t:", L
print "Min and Max M/L (solar units) \t\t:", mol1, mol2, '\n'
 
print "l, b (degrees) \t\t\t:", l, b
print "x0, y0, z0 (kpc) \t\t:", x0, y0, z0   
print "vx0, vy0, vz0 (km/s) \t\t:", vx0, vy0, vz0
print "amx0, amy0, amz0 \t\t:", amx0, amy0, amz0   
print "lam0, bam0 (degrees) \t\t:", lam0s[0], bam0s[0], "\n"  

print "Actual Perigalacticon (kpc) \t:", peris[0]
print "Mean Perigalacticon (kpc) \t:", np.mean(peris), "+-", np.std(peris)
print "Max Perigalacticon (kpc) \t:", np.amax(peris)
print "Min Perigalacticon (kpc) \t:", np.amin(peris)
lower, upper = st.t.interval(0.95, len(peris)-1, loc=np.mean(peris), scale=st.sem(peris))
print "95% Confidence Limits (kpc) \t:", lower, upper
print "Fraction of orbits with peri <", rpval, "is", np.count_nonzero(np.array(peris) < rpval)/float(len(peris)), "\n"

print "Actual Apogalacticon (kpc) \t:", apos[0] 
print "Mean Apogalacticon (kpc) \t:", np.mean(apos), "+-", np.std(apos)
print "Max Apogalacticon (kpc) \t:", np.amax(apos)
print "Min Apogalacticon (kpc) \t:", np.amin(apos)
lower, upper = st.t.interval(0.95, len(apos)-1, loc=np.mean(apos), scale=st.sem(apos))
print "95% Confidence Limits (kpc) \t:", lower, upper, '\n'

print "Actual Eccentricity \t\t:", eccs[0] 
print "Mean Eccentricity \t\t:", np.mean(eccs), "+-", np.std(eccs)
print "Max Eccentricity \t\t:", np.amax(eccs)
print "Min Eccentricity \t\t:", np.amin(eccs)
lower, upper = st.t.interval(0.95, len(eccs)-1, loc=np.mean(eccs), scale=st.sem(eccs))
print "95% Confidence Limits \t\t:", lower, upper, '\n'

print "Actual Orbital Period (Gyr) \t:", periods[0] 
print "Mean Orbital Period (Gyr) \t:", np.mean(periods), "+-", np.std(periods)
print "Max Orbital Period (Gyr) \t:", np.amax(periods)
print "Min Orbital Period (Gyr) \t:", np.amin(periods)
lower, upper = st.t.interval(0.95, len(periods)-1, loc=np.mean(periods), scale=st.sem(periods))
print "95% Confidence Limits (Gyr) \t:", lower, upper, '\n'

print "Actual r_t1 (arcmin) \t\t:", rt1s[0] 
print "Mean r_t1 (arcmin) \t\t:", np.mean(rt1s), "+-", np.std(rt1s)
print "Max r_t1 (arcmin) \t\t:", np.amax(rt1s)
print "Min r_t1 (arcmin) \t\t:", np.amin(rt1s)
lower, upper = st.t.interval(0.95, len(rt1s)-1, loc=np.mean(rt1s), scale=st.sem(rt1s))
print "95% Confidence Limits (arcmin) \t:", lower, upper
print "Fraction smaller than", rtlim, "\t\t:", np.count_nonzero(np.array(rt1s) < rtlim)/float(len(rt1s)), "\n"

print "Actual r_t2 (arcmin) \t\t:", rt2s[0] 
print "Mean r_t2 (arcmin) \t\t:", np.mean(rt2s), "+-", np.std(rt2s)
print "Max r_t2 (arcmin) \t\t:", np.amax(rt2s)
print "Min r_t2 (arcmin) \t\t:", np.amin(rt2s)
lower, upper = st.t.interval(0.95, len(rt2s)-1, loc=np.mean(rt2s), scale=st.sem(rt2s))
print "95% Confidence Limits (arcmin) \t:", lower, upper
print "Fraction smaller than", rtlim, "\t\t:", np.count_nonzero(np.array(rt2s) < rtlim)/float(len(rt2s)), "\n"

print "Actual Inclination (degrees) \t:", thetas[0] 
print "Mean Inclination (degrees) \t:", np.mean(thetas), "+-", np.std(thetas)
print "Max Inclination (degrees) \t:", np.amax(thetas)
print "Min Inclination (degrees) \t:", np.amin(thetas)
lower, upper = st.t.interval(0.95, len(thetas)-1, loc=np.mean(thetas), scale=st.sem(thetas))
print "95% Confidence Limits \t\t:", lower, upper, '\n'

print "Actual Longitude (degrees) \t:", omegas[0] 
print "Mean Longitude (degrees) \t:", np.mean(omegas), "+-", np.std(omegas)
print "Max Longitude (degrees) \t:", np.amax(omegas)
print "Min Longitude (degrees) \t:", np.amin(omegas)
lower, upper = st.t.interval(0.95, len(omegas)-1, loc=np.mean(omegas), scale=st.sem(omegas))
print "95% Confidence Limits \t\t:", lower, upper, '\n'

print "Actual lpole (degrees) \t\t:", lpoles[0] 
print "Mean lpole (degrees) \t\t:", np.mean(lpoles), "+-", np.std(lpoles)
print "Max lpole (degrees) \t\t:", np.amax(lpoles)
print "Min lpole (degrees) \t\t:", np.amin(lpoles)
lower, upper = st.t.interval(0.95, len(lpoles)-1, loc=np.mean(lpoles), scale=st.sem(lpoles))
print "95% Confidence Limits \t\t:", lower, upper, '\n'

print "Actual bpole (degrees) \t\t:", bpoles[0] 
print "Mean bpole (degrees) \t\t:", np.mean(bpoles), "+-", np.std(bpoles)
print "Max bpole (degrees) \t\t:", np.amax(bpoles)
print "Min bpole (degrees) \t\t:", np.amin(bpoles)
lower, upper = st.t.interval(0.95, len(bpoles)-1, loc=np.mean(bpoles), scale=st.sem(bpoles))
print "95% Confidence Limits \t\t:", lower, upper, '\n'

print "Actual lam (degrees) \t\t:", lam0s[0] 
print "Mean lam (degrees) \t\t:", np.mean(lam0s), "+-", np.std(lam0s)
print "Max lam (degrees) \t\t:", np.amax(lam0s)
print "Min lam (degrees) \t\t:", np.amin(lam0s)
lower, upper = st.t.interval(0.95, len(lam0s)-1, loc=np.mean(lam0s), scale=st.sem(lam0s))
print "95% Confidence Limits \t\t:", lower, upper, '\n'

print "Actual bam (degrees) \t\t:", bam0s[0] 
print "Mean bam (degrees) \t\t:", np.mean(bam0s), "+-", np.std(bam0s)
print "Max bam (degrees) \t\t:", np.amax(bam0s)
print "Min bam (degrees) \t\t:", np.amin(bam0s)
lower, upper = st.t.interval(0.95, len(bam0s)-1, loc=np.mean(bam0s), scale=st.sem(bam0s))
print "95% Confidence Limits \t\t:", lower, upper, '\n'

print "Actual delE/E \t\t\t:", deltaEs[0]
print "Mean delE/E \t\t\t:", np.mean(deltaEs), "+-", np.std(deltaEs)
print "Max delE/E \t\t\t:", np.amax(deltaEs)
print "Min delE/E \t\t\t:", np.amin(deltaEs)
lower, upper = st.t.interval(0.95, len(deltaEs)-1, loc=np.mean(deltaEs), scale=st.sem(deltaEs))
print "95% Confidence Limits \t\t:", lower, upper, '\n'

print "Average change lam, bam (deg) \t:", np.mean(dlams), np.mean(dbams), '\n'

print "No. of MC Experiments \t\t:", NMC
print "No. of Successful Experiments \t:", NMC - Nfailed

with open('galacticorbitmc_{0}_{1}.res' .format(pot, name), 'w') as f:
    f.write("Input Galaxy Parameters:\n")
    f.write("Name \t\t\t\t\t: {0} \n" .format(name))
    f.write("RA, DEC (degrees) \t\t\t: {0} {1} \n" .format(ra, dec))
    f.write("Heliocentric Distance (kpc) \t\t: {0} \n" .format(distance))
    f.write("mu_alpha (mas/century) \t\t\t: {0} +- {1} \n" .format(mu_a, mu_a_std))
    f.write("mu_delta (mas/century) \t\t\t: {0} +- {1} \n" .format(mu_d, mu_d_std))
    f.write("Heliocentric Radial Velocity (km/s) \t: {0} +- {1} \n" .format(rv, rv_std))
    f.write("Luminosity (solar luminosities) \t: {0} \n" .format(L))
    f.write("Min and Max M/L (solar units) \t\t: {0} {1} \n" .format(mol1, mol2))
    f.write("Limiting Radius (arcmin) \t\t: {0} \n\n" .format(rtlim))

    f.write("Integration and MC Parameters:\n")
    f.write("Time (Gyr) \t\t\t\t: {0} \n" .format(t))
    f.write("Number of Steps \t\t\t: {0} \n" .format(NSTEP)) 
    f.write("Number of MC Experiments \t\t: {0} \n" .format(NMC))
    f.write("Number of Successful Experiments \t: {0} \n\n" .format(NMC - Nfailed))

    f.write("l, b (degrees) \t\t\t: {0} {1} \n" .format(l, b))
    f.write("x0, y0, z0 (kpc) \t\t: {0} {1} {2} \n" .format(x0, y0, z0))
    f.write("vx0, vy0, vz0 (km/s) \t\t: {0} {1} {2} \n" .format(vx0, vy0, vz0))
    f.write("amx0, amy0, amz0 \t\t: {0} {1} {2} \n" .format(amx0, amy0, amz0))
    f.write("lam0, bam0 (degrees) \t\t: {0} {1} \n\n" .format(lam0s[0], bam0s[0]))

    f.write("Actual Perigalacticon (kpc) \t: {0} \n" .format(peris[0]))
    f.write("Mean Perigalacticon (kpc) \t: {0} +- {1} \n" .format(np.mean(peris), np.std(peris)))
    f.write("Max Perigalacticon (kpc) \t: {0} \n" .format(np.amax(peris)))
    f.write("Min Perigalacticon (kpc) \t: {0} \n" .format(np.amin(peris)))
    lower, upper = st.t.interval(0.95, len(peris)-1, loc=np.mean(peris), scale=st.sem(peris))
    f.write("95% Confidence Limits (kpc) \t: {0} {1} \n" .format(lower, upper))
    f.write("Fraction of orbits with peri < {0} is {1} \n\n" .format(rpval, np.count_nonzero(np.array(peris) < rpval)/float(len(peris)))) 

    f.write("Actual Apogalacticon (kpc) \t: {0} \n" .format(apos[0]))
    f.write("Mean Apogalacticon (kpc) \t: {0} +- {1} \n" .format(np.mean(apos), np.std(apos)))
    f.write("Max Apogalacticon (kpc) \t: {0} \n" .format(np.amax(apos)))
    f.write("Min Apogalacticon (kpc) \t: {0} \n" .format(np.amin(apos)))
    lower, upper = st.t.interval(0.95, len(apos)-1, loc=np.mean(apos), scale=st.sem(apos))
    f.write("95% Confidence Limits (kpc) \t: {0} {1} \n\n" .format(lower, upper))

    f.write("Actual Eccentricity \t\t: {0} \n" .format(eccs[0]))
    f.write("Mean Eccentricity \t\t: {0} +- {1} \n" .format(np.mean(eccs), np.std(eccs)))
    f.write("Max Eccentricity \t\t: {0} \n" .format(np.amax(eccs)))
    f.write("Min Eccentricity \t\t: {0} \n" .format(np.amin(eccs)))
    lower, upper = st.t.interval(0.95, len(eccs)-1, loc=np.mean(eccs), scale=st.sem(eccs))
    f.write("95% Confidence Limits \t\t: {0} {1} \n\n" .format(lower, upper))

    f.write("Actual Orbital Period (Gyr) \t: {0} \n" .format(periods[0]))
    f.write("Mean Orbital Period (Gyr) \t: {0} +- {1} \n" .format(np.mean(periods), np.std(periods)))
    f.write("Max Orbital Period (Gyr) \t: {0} \n" .format(np.amax(periods)))
    f.write("Min Orbital Period (Gyr) \t: {0} \n" .format(np.amin(periods)))
    lower, upper = st.t.interval(0.95, len(periods)-1, loc=np.mean(periods), scale=st.sem(periods))
    f.write("95% Confidence Limits (Gyr) \t: {0} {1} \n\n" .format(lower, upper))

    f.write("Actual r_t1 (arcmin) \t\t: {0} \n" .format(rt1s[0]))
    f.write("Mean r_t1 (arcmin) \t\t: {0} +- {1} \n" .format(np.mean(rt1s), np.std(rt1s)))
    f.write("Max r_t1 (arcmin) \t\t: {0} \n" .format(np.amax(rt1s)))
    f.write("Min r_t1 (arcmin) \t\t: {0} \n" .format(np.amin(rt1s)))
    lower, upper = st.t.interval(0.95, len(rt1s)-1, loc=np.mean(rt1s), scale=st.sem(rt1s))
    f.write("95% Confidence Limits (arcmin) \t: {0} {1} \n" .format(lower, upper))
    f.write("Fraction smaller than {0} \t: {1} \n\n" .format(rtlim, np.count_nonzero(np.array(rt1s) < rtlim)/float(len(rt1s))))

    f.write("Actual r_t2 (arcmin) \t\t: {0} \n" .format(rt2s[0]))
    f.write("Mean r_t2 (arcmin) \t\t: {0} +- {1} \n" .format(np.mean(rt2s), np.std(rt2s)))
    f.write("Max r_t2 (arcmin) \t\t: {0} \n" .format(np.amax(rt2s)))
    f.write("Min r_t2 (arcmin) \t\t: {0} \n" .format(np.amin(rt2s)))
    lower, upper = st.t.interval(0.95, len(rt2s)-1, loc=np.mean(rt2s), scale=st.sem(rt2s))
    f.write("95% Confidence Limits (arcmin) \t: {0} {1} \n" .format(lower, upper))
    f.write("Fraction smaller than {0} \t: {1} \n\n" .format(rtlim, np.count_nonzero(np.array(rt2s) < rtlim)/float(len(rt2s))))

    f.write("Actual Inclination (degrees) \t: {0} \n" .format(thetas[0]))
    f.write("Mean Inclination (degrees) \t: {0} +- {1} \n" .format(np.mean(thetas), np.std(thetas)))
    f.write("Max Inclination (degrees) \t: {0} \n" .format(np.amax(thetas)))
    f.write("Min Inclination (degrees) \t: {0} \n" .format(np.amin(thetas)))
    lower, upper = st.t.interval(0.95, len(thetas)-1, loc=np.mean(thetas), scale=st.sem(thetas))
    f.write("95% Confidence Limits \t\t: {0} {1} \n\n" .format(lower, upper))

    f.write("Actual Longitude (degrees) \t: {0} \n" .format(omegas[0]))
    f.write("Mean Longitude (degrees) \t: {0} +- {1} \n" .format(np.mean(omegas), np.std(omegas)))
    f.write("Max Longitude (degrees) \t: {0} \n" .format(np.amax(omegas)))
    f.write("Min Longitude (degrees) \t: {0} \n" .format(np.amin(omegas)))
    lower, upper = st.t.interval(0.95, len(omegas)-1, loc=np.mean(omegas), scale=st.sem(omegas))
    f.write("95% Confidence Limits \t\t: {0} {1} \n\n" .format(lower, upper))

    f.write("Actual lpole (degrees) \t\t: {0} \n" .format(lpoles[0]))
    f.write("Mean lpole (degrees) \t\t: {0} +- {1} \n" .format(np.mean(lpoles), np.std(lpoles)))
    f.write("Max lpole (degrees) \t\t: {0} \n" .format(np.amax(lpoles)))
    f.write("Min lpole (degrees) \t\t: {0} \n" .format(np.amin(lpoles)))
    lower, upper = st.t.interval(0.95, len(lpoles)-1, loc=np.mean(lpoles), scale=st.sem(lpoles))
    f.write("95% Confidence Limits \t\t: {0} {1} \n\n" .format(lower, upper))

    f.write("Actual bpole (degrees) \t\t: {0} \n" .format(bpoles[0]))
    f.write("Mean bpole (degrees) \t\t: {0} +- {1} \n" .format(np.mean(bpoles), np.std(bpoles)))
    f.write("Max bpole (degrees) \t\t: {0} \n" .format(np.amax(bpoles)))
    f.write("Min bpole (degrees) \t\t: {0} \n" .format(np.amin(bpoles)))
    lower, upper = st.t.interval(0.95, len(bpoles)-1, loc=np.mean(bpoles), scale=st.sem(bpoles))
    f.write("95% Confidence Limits \t\t: {0} {1} \n\n" .format(lower, upper))

    f.write("Actual lam (degrees) \t\t: {0} \n" .format(lam0s[0]))
    f.write("Mean lam (degrees) \t\t: {0} +- {1} \n" .format(np.mean(lam0s), np.std(lam0s)))
    f.write("Max lam (degrees) \t\t: {0} \n" .format(np.amax(lam0s)))
    f.write("Min lam (degrees) \t\t: {0} \n" .format(np.amin(lam0s)))
    lower, upper = st.t.interval(0.95, len(lam0s)-1, loc=np.mean(lam0s), scale=st.sem(lam0s))
    f.write("95% Confidence Limits \t\t: {0} {1} \n\n" .format(lower, upper))

    f.write("Actual bam (degrees) \t\t: {0} \n" .format(bam0s[0]))
    f.write("Mean bam (degrees) \t\t: {0} +- {1} \n" .format(np.mean(bam0s), np.std(bam0s)))
    f.write("Max bam (degrees) \t\t: {0} \n" .format(np.amax(bam0s)))
    f.write("Min bam (degrees) \t\t: {0} \n" .format(np.amin(bam0s)))
    lower, upper = st.t.interval(0.95, len(bam0s)-1, loc=np.mean(bam0s), scale=st.sem(bam0s))
    f.write("95% Confidence Limits \t\t: {0} {1} \n\n" .format(lower, upper))

    f.write("Average change lam, bam (deg) \t: {0} {1} \n\n" .format(np.nanmean(dlams), np.nanmean(dbams)))

    f.write("Some statistics on the percent change in energy from start to finish, \n")
    f.write("to serve as a check on the acceptability of the integration precision:\n\n")
   # f.write("integration precision: \n\n")

    f.write("Mean deltaE/E_0 : {0} +- {1} \n" .format(np.mean(deltaEs), np.std(deltaEs)))
    f.write("Max deltaE/E_0 \t: {0} \n" .format(np.amax(deltaEs)))
    f.write("Min deltaE/E_0 \t: {0} \n" .format(np.amin(deltaEs)))
