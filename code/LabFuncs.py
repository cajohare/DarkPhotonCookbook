#================================LabFuncs.py===================================#
# Created by Ciaran O'Hare 2021

# Description:
# Contains an assortment of functions that are all related to the 'Lab' somehow
# e.g. the nuclear form factor, lab velocity etc.


# Contains:

##### Lab velocity
# LabVelocity: Full lab velocity in (N,W,Z) with Earth rotation
# LabVelocitySimple: Simplified Lab velocity in galactic coordinates
# JulianDay: JulianDay at dd-mm-yyyy hh:hh
# EarthVelocity: Earth velocity to second order in eccentricity
# EarthVector: Earth radius vector to second order in eccentricity
# v_infinity: transforms velocity to the value outside the Sun's potential
# v_infinity_alt: same as v_inficity but with a different velocity discretisation
#####

##### Solar direction:
# EarthSunDistance: Distance between Earth and Sun as a function of time
# SolarDirection: Direction of the sun at a given time
#####

##### Co-ordinate transformations
# eqt2lab: Equatorial system to laboratory system
# gal2eqt: Galactic system to equatorial system
# gal2lab: Galactic system to lab system
#####
#==============================================================================#

import numpy as np
from numpy import cos, sin, pi, floor, exp, sqrt, size, zeros, shape, arccos
from numpy import array, trapz
import Params



#==============================Form Factors====================================#
def FormFactorHelm(E_r,A):
    q = sqrt(2*A*931.5*1000*E_r)*1.0e-12/1.97e-7
    c1 = 1.23*A**(1.0/3.0)-0.6
    s = 0.9
    R_1 = sqrt(c1**2 + (7.0/3.0)*pi**2.0*(0.52**2.0) - 5*s**2.0)
    F = (3*(sin(q*R_1) - q*R_1*cos(q*R_1))*exp(-q*q*s*s/2.0)/(q*R_1)**3)
    return F
#------------------------------------------------------------------------------#





#=======================Apply Angular Resolution===============================#
def Smear(x,dR,sig_gamma):
    # x = cartesian vectors
    # dR = value of rate at directions in x
    # sig_gamma = Gaussian width to smear dR by
    npix = size(dR)
    dR_smeared = zeros(shape=shape(dR))
    for i in range(0,npix):
        x0 = x[i,:]
        gamma = x0[0]*x[:,0] + x0[1]*x[:,1] + x0[2]*x[:,2]
        gamma[i] = 1.0
        gamma = arccos(gamma)
        dR_smeared[i] = sum(dR*exp(-gamma**2.0/(2*sig_gamma**2.0)))

    # Make sure it's normalised to what it was before the smearing
    dR_smeared = dR_smeared*sum(dR)/sum(dR_smeared)
    return dR_smeared
#------------------------------------------------------------------------------#

#===========================Apply Energy Res===================================#
def SmearE(E,dR,sig_E):
    # E = energies
    # dR = value of rate at energies in E
    # sig_E = Gaussian width to smear dR by
    nE = size(dR)
    dR_smeared = zeros(shape=shape(dR))
    if size(sig_E)==1:
        sig_E *= ones(shape=shape(dR))
    for i in range(0,nE):
        Ediff = abs(E-E[i])
        dR_smeared[i] = trapz(dR*exp(-Ediff**2.0/(2*sig_E**2.0)),E)

    # Make sure it's normalised to what it was before the smearing
    dR_smeared = dR_smeared*trapz(dR,E)/trapz(dR_smeared,E)
    return dR_smeared
#------------------------------------------------------------------------------#



#==============================Lab Velocity====================================#
# Peculiar velocity
v_pec = array([11.1,12.2,7.3])

# Earth orbital params
vv_earthrev = 29.79
eccentricity = 0.016722
eccentricity_deg = 0.9574
orb_long_ecliptic = 13.0+1.0
lat_ecl_gal = np.array([-5.5303,59.575,29.812])
long_ecl_gal = np.array([266.141,-13.3485,179.3212])
e1 = array([0.9941,0.1088,0.0042])
e2 = array([-0.0504,0.4946,-0.8677])
w_p = 2*pi/365 # orbital freq.
t1 = 79
ve = 29.79 # Earth's revolution
vrot = 0.47 # Earth's rotation

# Other constants
AstronomicalUnit = 1.49597892e11 # Astronomical Unit
EarthRadius = 6371.01*1000.0 # Earth Radius
Msun = 2.0e30 # Solar mass (kg)
bigG = 6.67e-11*(1.0e3)**(-3)
Jan1 = 2458849.5 # Julian date of January 1 2019

#------------------------------------------------------------------------------#
def LabVelocity(day, Loc=Params.Boulby, v_LSR=233.0):
    JD = day+Jan1
    lat = Loc.Latitude
    lon = Loc.Longitude

    # Convert day into phase of Earth rotation t_lab
    UT = 24*(JD+0.5-floor(JD+0.5)) #Universal time
    MJD = JD - 2400000.5 #Modified Julian Day
    T_0 = (floor(MJD)-55197.5)/36525.0
    t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
    t_lab = t_GAST + lon/15
    t_lab = 15*t_lab #Lab time in degrees


    # Galactic (LSR) Rotation
    vtemp = np.array([0.0,v_LSR,0.0])
    v_galrot = gal2lab(vtemp,t_lab, lat) #transform to lab co-ords

    # Peculiar solar Motion
    vtemp1 = v_pec
    v_solar = gal2lab(vtemp1,t_lab, lat) # transform to lab co-ords

    #Earth's revolution (first calculate in galactic frame then transform)
    e = eccentricity
    lambda_0 = orb_long_ecliptic
    L = 281.0298 + 36000.77*T_0 + 0.04107*UT
    g = 357.9258 + 35999.05*T_0 + 0.04107*UT
    lambda_sun = L + (1.915 - 0.0048*T_0)*sin(g*pi/180.0)\
         + 0.020*sin(2*g*pi/180.0)
    beta = lat_ecl_gal
    lambda_i = long_ecl_gal
    v_earthrev1 = vv_earthrev*(1-e*sin(pi/180.0*(lambda_sun-lambda_0)))*\
         (cos(beta*pi/180.0)*sin(pi/180.0*(lambda_sun-lambda_i)))
    v_earthrev = gal2lab(v_earthrev1,t_lab, lat) #transform to lab co-ords

    # Earth's rotation (already in lab co-ords)
    v_earthrot = 0.465102*cos(lat*pi/180)*np.array([0.0,-1.0,0.0])

    # Add them all together (delete as needed)
    v_lab = np.array([0.,0.,0.])
    v_lab += v_earthrot
    v_lab += v_earthrev
    v_lab += v_solar
    v_lab += v_galrot

    return v_lab


def JulianDay(month, day, year, hour): # Calculates time in JD for a given date
    year_r = year+4800-floor((14-month)/12.0)
    month_r = month+12*floor((14-month)/12.0)-3
    JulianDay = day + floor((153*month_r+2)/5.0) + 365*year_r\
                + floor(year_r/4.0) - floor(year_r/100.0)\
                + floor(year_r/400.0) - 32045 + (hour-12.0)/24.0
    return JulianDay

def LabVelocitySimple(day,v_LSR=233.0):
    # day measured from Jan1
    vsun = array([0.0,v_LSR,0.0])+v_pec
    v_lab = vsun + EarthVelocity(day)
    return v_lab

def EarthVelocity(day):
    # Second order in eccentricity
    # day measured from Jan1
    lambda_p = 102.93*pi/180.0
    th = w_p*(day-t1)
    v_E = cos(th)*(e1-2*eccentricity*sin(lambda_p)*e2) \
          +sin(th)*(e2+2*eccentricity*sin(lambda_p)*e1) \
          -eccentricity*(cos(2*th)*(cos(lambda_p)*e1-sin(lambda_p)*e2) \
          +sin(2*th)*(sin(lambda_p)*e1+cos(lambda_p)*e2))
    return vv_earthrev*v_E


def EarthVector(day):
    # Earth's orbital radius vectors
    # day measured from Jan1
    # Second order in Earth's eccentricity
    a_earth = AstronomicalUnit/1.0e3
    tp = 3
    lamb_p = 102*pi/180
    g = w_p*(day-tp)
    nu = g + 2.*eccentricity*sin(g)*(5.0/4.0)+eccentricity**2.0*sin(2*g)
    r = a_earth*(1-eccentricity**2.0)/(1+eccentricity*cos(nu))
    r_earth = r*(-sin(lamb_p+nu)*e1 + cos(lamb_p+nu)*e2)
    return r_earth

def v_infinity(v,costh,phi,day):
    # v_infinity used for Gravitational focusing, this version uses a set of
    # angles costh and phi
    # day measured from Jan1
    x_earth = EarthVector(day)
    r_earth = sqrt(sum(x_earth**2.0))
    x_earth /= r_earth # unit vector towards Earth
    v_earth = EarthVelocity(day)
    uu_esc = 2*bigG*Msun/r_earth # escape speed

    vx = v*sqrt(1.0-costh**2.0)*cos(phi)+v_earth[0]
    vy = v*sqrt(1.0-costh**2.0)*sin(phi)+v_earth[1]
    vz = v*costh+v_earth[2]

    vv_inf = (vx**2.0+vy**2.0+vz**2.0)-uu_esc
    vv_inf = (vv_inf+abs(vv_inf))/2.0

    vdotr = vx*x_earth[0]+vy*x_earth[1]+vz*x_earth[2]

    v_inf = sqrt(vv_inf)

    denom = vv_inf + 0.5*uu_esc - v_inf*vdotr
    v_infx = (vv_inf*vx + 0.5*v_inf*uu_esc*x_earth[0] - v_inf*vx*vdotr)/denom
    v_infy = (vv_inf*vy + 0.5*v_inf*uu_esc*x_earth[1] - v_inf*vy*vdotr)/denom
    v_infz = (vv_inf*vz + 0.5*v_inf*uu_esc*x_earth[2] - v_inf*vz*vdotr)/denom


    return v_infx,v_infy,v_infz

def v_infinity_alt(v3,day):
    # v_infinity used for Gravitational focusing, this version uses a set of
    # angles cartesian velocity vectors in v3 which defines a Healpix
    # discretisation. Tends to be a bit faster and more accurate.
    # day measured from Jan1

    x_earth = EarthVector(day)
    r_earth = sqrt(sum(x_earth**2.0))
    x_earth /= r_earth # unit vector towards Earth
    v_earth = EarthVelocity(day)
    uu_esc = 2*bigG*Msun/r_earth # escape speed

    vx = v3[:,0]+v_earth[0] # galactic x-component
    vy = v3[:,1]+v_earth[1] # galactic y-component
    vz = v3[:,2]+v_earth[2] # galactic z-component

    vv_inf = (vx**2.0+vy**2.0+vz**2.0)-uu_esc
    vv_inf = (vv_inf+abs(vv_inf))/2.0
    #vv_inf[vv_inf<0.0] = 0.0

    vdotr = vx*x_earth[0]+vy*x_earth[1]+vz*x_earth[2] # (v.x_earth)

    v_inf = sqrt(vv_inf)

    denom = vv_inf + 0.5*uu_esc - v_inf*vdotr
    v_infx = (vv_inf*vx + 0.5*v_inf*uu_esc*x_earth[0] - v_inf*vx*vdotr)/denom
    v_infy = (vv_inf*vy + 0.5*v_inf*uu_esc*x_earth[1] - v_inf*vy*vdotr)/denom
    v_infz = (vv_inf*vz + 0.5*v_inf*uu_esc*x_earth[2] - v_inf*vz*vdotr)/denom
    return v_infx,v_infy,v_infz

#==========================Solar direction=====================================#
def EarthSunDistance(JD): # Earth-sun distance at Julian Day (JD)
    D = JD-2451545.0
    g = 357.529 + 0.98560028*D
    g = g*pi/180.0
    r_es = 1.00014 - 0.01671*cos(g) - 0.00014*cos(2*g)
    r_es = r_es*AstronomicalUnit
    return r_es

#------------------------------------------------------------------------------#
def SolarDirection(JD,Loc): # Solar direction in lab coords at Julian Day (JD)

    lat = Loc.Latitude
    lon = Loc.Longitude

    # Compute RA and dec of Sun
    #JD = day+Jan1
    n = JD - 2451545.0
    Omega = 2.1429-0.0010394594*n
    L = 4.8950630 + 0.017202791698*n
    g = 6.2400600 + 0.0172019699*n
    ll = L+0.03341607*sin(g) + 0.00034894*sin(2*g)\
        - 0.0001134 - 0.0000203*sin(Omega)
    ep = 0.4090928 - 6.214e-9*n + 0.0000396*cos(Omega)
    ra = np.arctan2((cos(ep)*sin(ll)),cos(ll)) # Right ascension of Sun
    dec = np.arcsin(sin(ep)*sin(ll)) # Declination of sun

    # Solar vector
    x_sun1 = np.array([0.,0.,0.])
    x_sun1[0] = cos(dec)*cos(ra)
    x_sun1[1] = cos(dec)*sin(ra)
    x_sun1[2] = sin(dec)

    # Lab time conversion
    UT = 24*(JD+0.5-floor(JD+0.5))
    MJD = JD - 2400000.5
    T_0 = (floor(MJD)-55197.5)/36525.0
    t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
    t_lab = t_GAST + lon/15.0
    t_lab = 15*t_lab # DEGREES

    # Convert vector from equatorial system into lab system
    x_sun = eqt2lab(x_sun1,t_lab,lat)
    return x_sun

def EarthSunDistanceMod(JD):
    # Solar neutrinos:
    # Flux is scaled by 1/EarthSunDistance^2 but since Flux is already averaged
    # We need to also divide by Integral(1/R^2) over one year
    # Integral_inv_EarthSun_sq is defined in params.f95
    Integral_inv_EarthSun_sq = 4.468864372000642e-23 # integral(1/R^2) over 1 year
    f = (1.0/Integral_inv_EarthSun_sq)*(1.0/EarthSunDistance(JD)**2.0)
    return f

#------------------------------------------------------------------------------#


#==============================================================================#
#---------------------------Coordinate trans.----------------------------------#
def eqt2lab(vp,t_lab,lat): # Equatorial (x_e,y_e,z_e) to Laboratory (N,W,Z)
    t = t_lab*pi/180.0
    latr = lat*pi/180.0
    v = vp*0.0
    v[0] = -cos(t)*sin(latr)*vp[0] - sin(t)*sin(latr)*vp[1] + cos(latr)*vp[2]
    v[1] = sin(t)*vp[0] - cos(t)*vp[1]
    v[2] = cos(t)*cos(latr)*vp[0] + cos(latr)*sin(t)*vp[1] + sin(latr)*vp[2]
    return v

def gal2eqt(vp): # Galactic (x_g,y_g,z_g) to Equatorial (x_e,y_e,z_e)
    v = 0.0*vp
    v[0] = -0.06699*vp[0] + 0.4927*vp[1] - 0.8676*vp[2]
    v[1] = -0.8728*vp[0]  - 0.4503*vp[1] - 0.1884*vp[2]
    v[2] = -0.4835*vp[0]  + 0.7446*vp[1] + 0.4602*vp[2]
    return v

def lab2gal(v,t_lab, lat): # Laboratory (N,W,Z) to Galactic (x_g,y_g,z_g)
    vp = lab2eqt(v,t_lab,lat)
    return eqt2gal(vp)

def gal2lab(v,t_lab, lat): # Galactic (x_g,y_g,z_g) to Laboratory (N,W,Z)
    vp = gal2eqt(v)
    return eqt2lab(vp, t_lab, lat)
#==============================================================================#




# def BinEvents(Expt,dRfunc,*Args):
#     # Expt = Detector class
#     # dRfunc = differential recoil rate that is being binned
#     # Args = anything else needed by dRfunc
#
#     # Energy and time binning:
#     E_bins = Expt.Energies
#     t_bins = Expt.Times
#     Efficiency = Expt.Efficiency
#     ne = size(E_bins)
#     nt = size(t_bins)
#
#     # DIRECTIONAL LIMITS
#     if Expt.Directional:
#         q = Expt.Directions
#         sig_gamma = Expt.AngularResolution
#         HeadTail = Expt.HeadTailEfficiency
#
#         npix = size(q)/3
#         E_r = zeros(shape=(ne*nt*npix))
#         E = zeros(shape=(ne*nt*npix,3))
#         eff = zeros(shape=(ne*nt*npix))
#         t = zeros(shape=(ne*nt*npix))
#         eff_HT = zeros(shape=(ne*nt*npix))
#         ii = 0
#         for i in range(0,nt):
#             for j in range(0,ne):
#                 for k in range(0,npix):
#                     E_r[ii] = E_bins[j]
#                     E[ii,:] = E_bins[j]*q[k,:]
#                     t[ii] = t_bins[i]
#                     eff[ii] = Efficiency[j]
#                     eff_HT[ii] = HeadTail[j]
#                     ii += 1
#
#         # Correct for Head-Tail
#         if HeadTail[0]>0.99:
#             dR = dRfunc(E,t,Expt,Args[0],Args[1])
#         else:
#             dR = (1.0-eff_HT)*dRfunc(E,t,Expt,Args[0],Args[1])\
#                     +eff_HT*dRfunc(-1.0*E,t,Expt,Args[0],Args[1])
#
#         dR = dR*4*pi/(1.0*npix*nt)
#         # Correct for Efficiency
#         if Efficiency[0]<0.99:
#             dR = dR*eff
#
#         # Correct for Angular resolution
#         if sig_gamma[0]>0.01:
#             i1 = 0
#             dR_smear = zeros(shape=shape(dR))
#             for i in range(0,nt):
#                 for j in range(0,ne):
#                     i2 = i1 + npix - 1
#                     if sum(dR[i1:i2+1])>0.0:
#                         dR_smear[i1:i2+1] = Smear(q,dR[i1:i2+1],sig_gamma[j])
#                     i1 = i2+1
#             dR = dR_smear
#
#         # Bin events
#         i1 = 0
#         RD = zeros(shape=(ne-1)*nt*npix)
#         for i in range(0,nt):
#             for j in range(0,ne-1):
#                 i2 = i1 + npix - 1
#                 dR1 = dR[(t==t_bins[i])&(E_r==E_bins[j])]
#                 dR2 = dR[(t==t_bins[i])&(E_r==E_bins[j+1])]
#                 RD[i1:i2+1] = 0.5*(E_bins[j+1] - E_bins[j])*(dR1+dR2)
#                 i1 = i2+1
#         # Last step: turn energy off if needed
#         if Expt.EnergyOff:
#             i1 = 0
#             RD_reduced = zeros(shape=(ne-1)*nt)
#             it1 = 0
#             for i in range(0,nt):
#                 it2 = it1 + npix -1
#                 for j in range(0,ne-1):
#                     i2 = i1 + npix - 1
#                     RD_reduced[it1:it2+1] += RD[i1:i2]
#                     i1 = i2 + 1
#                 it1 = it2+1
#             RD = RD_reduced
#
#     # Non-directional limits
#     else:
#         E = zeros(shape=(ne*nt))
#         t = zeros(shape=(ne*nt))
#         eff = zeros(shape=(ne*nt))
#         ii = 0
#         for i in range(0,nt):
#             for j in range(0,ne):
#                 E[ii] = E_bins[j]
#                 t[ii] = t_bins[i]
#                 eff[ii] = Efficiency[j]
#                 ii += 1
#
#         dR = dRfunc(E,t,Expt,Args[0],Args[1])
#         dR = dR/(1.0*nt)
#         # Correct for Efficiency
#         if Efficiency[0]<0.99:
#             dR = dR*eff
#
#         # Bin events
#         i1 = 0
#         RD = zeros(shape=(ne-1)*nt)
#         for i in range(0,nt):
#             i2 = i1 + ne - 2
#             dR1 = dR[(t==t_bins[i])]
#             RD[i1:i2+1] = 0.5*(E_bins[1:] - E_bins[0:-1])*(dR1[1:]+dR1[0:-1])
#             i1 = i2 + 1
#
#     RD *= Expt.Exposure
#     return RD
# #------------------------------------------------------------------------------#
#
#
def GravFocusAngles(vv,costh,phi,day,sig=164.75,v_LSR=233.0,v_shift=array([0.0,0.0,0.0])):
    v1 = vv*sqrt(1.0-costh**2.0)*cos(phi)
    v2 = vv*sqrt(1.0-costh**2.0)*sin(phi)
    v3 = vv*costh

    v0 = sqrt(2.0)*sig

    vsun = array([0.0,v_LSR,0.0])+v_pec-v_shift
    vearth = EarthVelocity(day)
    rearth = EarthVector(day)
    r = sqrt(sum(rearth**2.0))
    xearth = rearth/r

    vsum1 = v1+vearth[0]
    vsum2 = v2+vearth[1]
    vsum3 = v3+vearth[2]
    normvsum = sqrt(vsum1**2 + vsum2**2 + vsum3**2)
    xsum1 = vsum1/normvsum
    xsum2 = vsum2/normvsum
    xsum3 = vsum3/normvsum
    rx = 1.0-(xearth[0]*xsum1 + xearth[1]*xsum2 + xearth[2]*xsum3)
    vrx = (v1+vearth[0]+vsun[0])*(xearth[0]-xsum1)\
            +(v2+vearth[1]+vsun[1])*(xearth[1]-xsum2)\
            +(v3+vearth[2]+vsun[2])*(xearth[2]-xsum3)

    J = (-2*bigG*Msun/(r*v0**2.0*vv))*(vrx/rx)
    return J
