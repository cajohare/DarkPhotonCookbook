#================================LabFuncs.py===================================#
# Created by Ciaran O'Hare 2021
#==============================================================================#

from numpy import *

#==============================================================================#
'''
Cos(th)^2 factors for DP polarisation axes defined by
Input:
t = time in days
costh_X,phi_X = DP polarisation axis (arrays of size N)
lat = latitude of expt. in degrees
OUTPUT: cos(theta(t))^2 for each (th_X,phi_X)
'''
# North-pointing experiments
@jit(nopython=True)
def costh2N(t,costh_X,phi_X,lat):
    lat *= pi/180
    wt = 2*pi*t
    th_X = arccos(costh_X)
    return (sin(th_X)*cos(phi_X)*sin(lat)*cos(wt) \
                - sin(th_X)*sin(phi_X)*sin(lat)*sin(wt) \
                + costh_X*cos(lat))**2

# West-pointing experiments
@jit(nopython=True)
def costh2W(t,costh_X,phi_X,lat):
    lat *= pi/180
    wt = 2*pi*t
    th_X = arccos(costh_X)
    return (sin(th_X)*cos(phi_X)*sin(wt) - sin(th_X)*sin(phi_X)*cos(wt))**2

# Zenith-pointing experiments
@jit(nopython=True)
def costh2Z(t,costh_X,phi_X,lat):
    lat *= pi/180
    wt = 2*pi*t
    th_X = arccos(costh_X)
    return (sin(th_X)*cos(phi_X)*cos(lat)*cos(wt)\
            + sin(th_X)*sin(phi_X)*cos(lat)*sin(wt)
            + costh_X*sin(lat))**2

# North-facing experiments
@jit(nopython=True)
def costh2ZW(t,costh_X,phi_X,lat):
    return 1 - costh2N(t,costh_X,phi_X,lat)

# West-facing experiments
@jit(nopython=True)
def costh2ZN(t,costh_X,phi_X,lat):
    return 1 - costh2W(t,costh_X,phi_X,lat)

# Zenith-facing experiments
@jit(nopython=True)
def costh2NW(t,costh_X,phi_X,lat):
    return 1 - costh2Z(t,costh_X,phi_X,lat)

@jit(nopython=True)
def costh_arb(t,costh_X,phi_X,lat,n,w,z):
    # For arbitrary vector (nN + wW + zZ)
    lat *= pi/180
    wt = 2*pi*t
    th_X = arccos(costh_X)
    N = sin(th_X)*cos(phi_X)*sin(lat)*cos(wt) - sin(th_X)*sin(phi_X)*sin(lat)*sin(wt) + costh_X*cos(lat)
    W = sin(th_X)*cos(phi_X)*sin(wt) - sin(th_X)*sin(phi_X)*cos(wt)
    Z = sin(th_X)*cos(phi_X)*cos(lat)*cos(wt) + sin(th_X)*sin(phi_X)*cos(lat)*sin(wt) + costh_X*sin(lat)
    return (n*N + w*W + z*Z)**2
#==============================================================================#


#==============================================================================#
# Get the longitude and latitude of the ecliptic plane
def Ecliptic(npoints):
    # ecliptic
    lon = linspace(0,359,100)*u.deg
    lat = zeros(shape=100)*u.deg
    c_eclp = coord.SkyCoord(lon=lon,lat=lat,frame='geocentrictrueecliptic').transform_to('galactic')
    b_eclp = c_eclp.b.to(u.deg).value
    l_eclp = c_eclp.l.to(u.deg).value
    l_eclp[l_eclp>180] = l_eclp[l_eclp>180]-360
    b_eclp = b_eclp[argsort(l_eclp)]
    l_eclp = l_eclp[argsort(l_eclp)]
    return l_eclp,b_eclp
#==============================================================================#

#==============================================================================#
#------------------------------Coordinate trans.-------------------------------#
@jit(nopython=True)
def lab2eqt(vp,t_lab,lat):
    t = t_lab*pi/180.0
    latr = lat*pi/180.0
    v = zeros(shape=shape(vp))
    v[:,0] = -cos(t)*sin(latr)*vp[:,0] + sin(t)*vp[:,1] + cos(latr)*cos(t)*vp[:,2]
    v[:,1] = -sin(latr)*sin(t)*vp[:,0] + -cos(t)*vp[:,1] + cos(latr)*sin(t)*vp[:,2]
    v[:,2] = cos(latr)*vp[:,0] + 0*vp[:,1] + sin(latr)*vp[:,2]
    return v

@jit(nopython=True)
def eqt2gal(vp):
    v = zeros(shape=shape(vp))
    v[:,0] = -0.066945*vp[:,0] - 0.872755*vp[:,1] - 0.483505*vp[:,2]
    v[:,1] = +0.492754*vp[:,0] - 0.450313*vp[:,1] + 0.744620*vp[:,2]
    v[:,2] = -0.867607*vp[:,0] - 0.188340*vp[:,1] + 0.460194*vp[:,2]
    return v

@jit(nopython=True)
def lab2gal(v,JD, lat, lon):
    # Convert day into phase of Earth rotation t_lab
    UT = 24*(JD+0.5-floor(JD+0.5)) #Universal time
    MJD = JD - 2400000.5 #Modified Julian Day
    T_0 = (floor(MJD)-55197.5)/36525.0
    t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
    t_lab = t_GAST + lon/15
    t_lab = 15*t_lab #Lab time in degrees
    vp = lab2eqt(v,t_lab,lat)
    return eqt2gal(vp)
#==============================================================================#



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
