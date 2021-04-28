#================================HaloFuncs.py==================================#
# Created by Ciaran O'Hare 2019

# Description:
# This file has every function needed to calculate all the various speed and
# velocity distributions for the SHM, Sausage, DarkShards, etc. These are then
# used to make most of the plots in the later sections of the paper

# All speeds are in km/s and all distributions are in units of inverse speed to
# whichever power the dimensionality of the distribution is

# At the moment everything gets calculated in Galactic coordinates centered on
# the earth. Mostly this is expressed in cylindrical coordinates

#==============================================================================#



from numpy import pi, sqrt, exp, zeros, size, shape, linspace, meshgrid, cos, sin
from numpy import trapz, arange, array, flipud, interp, inf, isnan, vstack
from scipy.integrate import cumtrapz
from numpy.linalg import norm
from scipy.special import erf, erfi
from scipy.stats import norm
import LabFuncs
import Params
import healpy as hp



#=============================Discretisations==================================#
# discretisations: These functions require a lot of integrals over angles, I
# doing all of by brute force by summing over trapezoids on a sphere. This Tends
# to be fast enough for the accuracy that is actually needed.

# Healpix discretisation of a sphere (more efficient for these calculations)
nside = 64
npix = 12*nside**2
dpix = 4*pi/(npix*1.0)
x_pix = zeros(shape=(npix,3))
for i in range(0,npix):
    x_pix[i,:] = hp.pix2vec(nside, i)

# Square discretisation of a sphere (less accurate for these calculations)
n = 201
dth = 2.0/(n-1.0)
dph = 2*pi/(2*n*1.0)
cvals = arange(-1.0,1.0,dth)
pvals = arange(0,2*pi-dph,dph)
C,P = meshgrid(cvals,pvals)
#==============================================================================#


#======================== Weighting ============================================#
# Calculates how much each shard needs to be weighted and then groups the Shards
# into the categories listed in the paper
def ShardsWeights(names,pops,Psun):
    # names = Array of Shards names
    # pops = Array of populations
    # Psun = Distance to Sun (in sigma)

    # first use height of gaussian distribution at the distance to the Sun
    weights = norm.pdf(Psun,loc=0.0,scale=1.0)/norm.pdf(0.0,loc=0.0,scale=1.0)
    n = size(pops)
    w = zeros(shape=n)
    w = pops*weights

    # Then loop over each group and lump each individual Shard together
    # Groups: S1, S2, Retrograde, Prograde, Low energy
    for idi in ['S1','S2','R','Ca','N']:
        mask1 = zeros(shape=n)==1
        for i in range(0,n):
            mask1[i] = names[i].startswith(idi)
        w[mask1] = w[mask1]/sum(w[mask1])

    # Then divide by 5 since there are 5 groups
    w /= 5.0

    return w
#==============================================================================#


#==============================================================================#
# Normalisation constants for the SHM's isotropic Gaussian distribution as well
# as the Sausage's triaxial distribution (also used for the Streams)
def Nesc_Isotropic(sig,v_esc):
    # sig = 1d dispersion
    # v_esc = escape velocity
    return erf(v_esc/(sqrt(2)*sig)) - \
    sqrt(2.0/pi)*(v_esc/sig)*exp(-v_esc**2.0/(2.0*sig**2.0))

def Nesc_Triaxial(sigr,sigphi,beta,v_esc):
    # sigr, sigphi = radial and azimuthal dispersions
    # beta = anisotropy parameter
    # v_esc = escape velocity
    N_esc = erf(v_esc/(sqrt(2.0)*sigr)) - sqrt((1.0-beta)/beta)\
            *exp(-v_esc**2.0/(2.0*sigphi**2.0))\
            *erfi(v_esc/(sqrt(2)*sigr)*sqrt(beta/(1-beta)))

    # Make sure the N_esc doesn't blow up (which it does for some streams)
    if (abs(N_esc)==inf) or isnan(N_esc):
        N_esc = 1.0
    return N_esc
#==============================================================================#






#======================= Velocity distributions ===============================#
def VelocityDist_Isotropic(v,day,v_LSR=233.0,sig=164.75,v_esc=528.0,\
                        v_shift=array([0.0,0.0,0.0]),GravFocus=False,\
                        EscapeSpeed=True):
    ## Isotropic velocity distribution
    # v = velocities (in galactic cylindrical coordinates)
    # v_LSR = Local standard of rest
    # sig = 1d dispersion
    # v_esc = escape speed
    # v_shift = any shift to v_lab needed (usually the stream velocity)
    # GravFocus = whether to calculate gravfocus or not
    # EscpaeSpeed = whether to implement v_esc or not (i.e. v_esc=inf)

    # lab velocity
    v_lab = LabFuncs.LabVelocitySimple(day,v_LSR=v_LSR)-v_shift

    v0 = sig*sqrt(2.0) # denominator in exponent
    N_esc = Nesc_Isotropic(sig,v_esc) # normalisation

    vr = v[:,0] # radial comp.
    vphi = v[:,1] # azimuthal comp.
    vz = v[:,2] # vertical comp.

    V = sqrt((vr+v_lab[0])**2.0+(vphi+v_lab[1])**2.0+(vz+v_lab[2])**2.0) # speed

    # Now compute value of distribution
    fv3  = (1.0/(N_esc*sqrt(2*pi)*sig**3.0)*\
          exp(-((vr+v_lab[0])**2.0\
                +(vz+v_lab[2])**2.0\
                +(vphi+v_lab[1])**2.0)/(2*sig**2.0)))*(V<v_esc)
    return fv3


def VelocityDist_Triaxial(v,day,sig3,v_LSR=233.0,v_esc=528.0,\
                        v_shift=array([0.0,0.0,0.0]),GravFocus=False):
    ## Triaxial velocity distribution
    # v = velocities (in galactic cylindrical coordinates)
    # v_LSR = Local standard of rest
    # sig3 = 3d dispersion (array)
    # v_esc = escape speed
    # v_shift = any shift to v_lab needed (usually the stream velocity)
    # GravFocus = whether to calculate gravfocus or not

    v_lab = LabFuncs.LabVelocitySimple(day,v_LSR=v_LSR)-v_shift

    # Dispersions
    sigr = sig3[0]
    sigphi = sig3[1]
    sigz = sig3[2]
    beta = 1.0-(sigphi**2.0+sigz**2.0)/(2*sigr**2.0)

    # Normalisation
    if beta>0.0:
        N_esc = Nesc_Triaxial(sigr,sigphi,beta,v_esc)
    elif beta==0.0:
        N_esc = Nesc_Isotropic(sigr,v_esc)
    else:
        N_esc = 1.0
    N = 1.0/(N_esc*(2*pi)**(1.5)*sigr*sigphi*sigz)

    vr = v[:,0] # radial comp.
    vphi = v[:,1] # azimuthal comp.
    vz = v[:,2] # vertical comp.

    V = sqrt((vr+v_lab[0])**2.0+(vphi+v_lab[1])**2.0+(vz+v_lab[2])**2.0) # Speed

    # Now compute value of distribution
    fv3  = N*exp(-((vr+v_lab[0])**2.0/(2*sigr**2.0))\
               -((vz+v_lab[2])**2.0/(2*sigz**2.0))\
               -((vphi+v_lab[1])**2.0/(2*sigphi**2.0)))*(V<v_esc)
    return fv3
#==============================================================================#




#====================== Speed distributions====================================#
def SpeedDist_Isotropic(v,day,v_LSR=233.0,sig=164.75,v_esc=528.0,\
                        v_shift=array([0.0,0.0,0.0]),GravFocus=False,\
                        EscapeSpeed=True):
    ## Isotropic speed distribution
    # v = speeds (in km/s)
    # v_LSR = Local standard of rest
    # sig = 1d dispersion
    # v_esc = escape speed
    # v_shift = any shift to v_lab needed (usually the stream velocity)
    # GravFocus = whether to calculate gravfocus or not
    # EscpaeSpeed = whether to implement v_esc or not (i.e. v_esc=inf)


    # Analytic form for the speed distribution:
    v_lab = LabFuncs.LabVelocitySimple(day,v_LSR=v_LSR)-v_shift
    v_e = sqrt(sum(v_lab**2.0))
    v0 = sig*sqrt(2.0)
    N_esc = Nesc_Isotropic(sig,v_esc)
    fv1 = (1.0/(N_esc*sqrt(2*pi)))*(v/(v_e*sig))\
        *(exp(-(v**2.0+v_e**2.0-2.0*v*v_e)/(2*sig**2.0))\
        -exp(-(v**2.0+v_e**2.0+2.0*v*v_e)/(2*sig**2.0)))\
        *((v)<(v_esc+v_e))
    f1 = (1.0/(N_esc*sqrt(2*pi)))*(v/(v_e*sig))\
        *(exp(-(v**2.0+v_e**2.0-2.0*v*v_e)/(2*sig**2.0))\
        -exp(-(v**2.0+v_e**2.0+2.0*v*v_e)/(2*sig**2.0)))
    f2 = (1.0/(N_esc*sqrt(2*pi)))*(v/(v_e*sig))\
                    *(exp(-(v**2.0+v_e**2.0-2.0*v*v_e)/(2*sig**2.0))\
                  -(2/sqrt(pi))*exp(-v_esc**2.0/(2*sig**2.0)))
    fv1[v<(v_esc-v_e)] = f1[v<(v_esc-v_e)]
    fv1[v>(v_esc-v_e)] = f2[v>(v_esc-v_e)]
    fv1 /= trapz(fv1,v)

    # Turn off escape if needed
    if not EscapeSpeed:
        N_esc = 1.0
        v_esc = 1000.0

    # Can do analytic grav focusing calculation which uses a pertubative result
    # currently not used in the paper since every example uses the function
    # following this one.
    if GravFocus:
        nvals = size(v)
        fvJ = zeros(shape=nvals)
        # Loop over speeds
        for i in range(0,nvals):
            vv = v[i]
            voff2 = (vv*sqrt(1-C**2)*cos(P)+v_lab[0])**2\
                    + (vv*sqrt(1-C**2.0)*sin(P)+v_lab[1])**2.0\
                    + (vv*C+v_lab[2])**2 # speed
            fv3 = vv**2.0*(1.0/(sqrt(pi)*pi))*(1.0/v0**3)*exp(-voff2/v0**2)
            fJ = fv3*LabFuncs.GravFocusAngles(vv,C,P,day,sig=sig)
            fvJ[i] = dth*dph*sum(sum(fJ))
        fv1 += fvJ

    return fv1



def SpeedDist_Triaxial(v,day,sig3,v_LSR=233.0,v_esc=528.0,\
                        v_shift=array([0.0,0.0,0.0]),GravFocus=False,\
                        GalFrame=False,\
                      EscapeSpeed=True,SmoothCutoff=False):
    ## Triaxial speed distribution
    # v = velocities (in galactic cylindrical coordinates)
    # v_LSR = Local standard of rest
    # sig3 = 3d dispersion (array)
    # v_esc = escape speed
    # v_shift = any shift to v_lab needed (usually the stream velocity)
    # GravFocus = whether to calculate gravfocus or not
    # GalFrame = whether to compute in Galactic Frame or not
    # EscapeSpeed = whether to have a finite escape speed or not
    # SmoothCutoff = whether to have a sharp or smooth escape speed cutoff

    # Dispersions
    sigr = sig3[0]
    sigphi = sig3[1]
    sigz = sig3[2]

    # Normalisation
    beta = 1.0-(sigphi**2.0+sigz**2.0)/(2*sigr**2.0)
    if beta>0.0:
        N_esc = Nesc_Triaxial(sigr,sigphi,beta,v_esc)
    elif beta==0.0:
        N_esc = Nesc_Isotropic(sigr,v_esc)
    else:
        N_esc = 1.0
    if not EscapeSpeed:
        N_esc = 1.0
        v_esc = 1000.0
    N = 1.0/(N_esc*(2*pi)**(1.5)*sigr*sigphi*sigz)



    n = size(v)
    fv1 = zeros(shape=n)

    if GravFocus==False:
        if GalFrame:
            v_off = -v_shift # Shift back by -v_shift to undo shift
            v_max = v_esc
        else:
            v_e = LabFuncs.LabVelocitySimple(day,v_LSR=v_LSR)
            v_max = v_esc+sqrt(sum(v_e**2.0))
            v_off = v_e-v_shift

        # Need a correction Fcorr if a smooth escape speed cutoff is needed
        if SmoothCutoff:
            vr = (v_max)*sqrt(1-C**2.0)*cos(P)+v_off[0]
            vphi = (v_max)*sqrt(1-C**2.0)*sin(P)+v_off[1]
            vz = (v_max)*C+v_off[2]
            V = sqrt(vr**2.0+vphi**2.0+vz**2.0)
            Fcorr = N*exp(-(vr**2.0/(2*sigr**2.0))\
                          -(vz**2.0/(2*sigz**2.0))\
                          -(vphi**2.0/(2*sigphi**2.0)))
        else:
            Fcorr = 0.0

        # Loop over speeds each step doing an angular integral
        for i in range(0,n):
            v1 = v[i]
            vr = v1*sqrt(1-C**2.0)*cos(P)+v_off[0]
            vphi = v1*sqrt(1-C**2.0)*sin(P)+v_off[1]
            vz = v1*C+v_off[2]
            V = sqrt(vr**2.0+vphi**2.0+vz**2.0)

            F  = N*exp(-(vr**2.0/(2*sigr**2.0))\
                       -(vz**2.0/(2*sigz**2.0))\
                       -(vphi**2.0/(2*sigphi**2.0)))*(V<v_esc)

            F = F-Fcorr # take off smooth cutoff correction
            fv1[i] = (v1**2.0)*dth*dph*sum(sum(F)) # sum over angles
        fv1[v>v_max] = 0.0
        fv1 /= trapz(fv1,v) # renormalise for security
    else:

        ##### Gravitation focusing calculation #####
        v_off = LabFuncs.v_pec+array([0.0,v_LSR,0.0])-v_shift
        vv_e = sqrt(sum(v_off**2.0))
        v_max = v_esc+vv_e

        if SmoothCutoff:
            vr = (v_max)*sqrt(1-C**2.0)*cos(P)
            vphi = (v_max)*sqrt(1-C**2.0)*sin(P)
            vz = (v_max)*C
            V = sqrt(vr**2.0+vphi**2.0+vz**2.0)
            Fcorr = N*exp(-(vr**2.0/(2*sigr**2.0))\
                          -(vz**2.0/(2*sigz**2.0))\
                          -(vphi**2.0/(2*sigphi**2.0)))
        else:
            Fcorr = 0.0

        # Loop over speeds each step doing an angular integral
        for i in range(0,n):
            v1 = v[i]
            # Transform velocities to v_infinity (needed to implement Grav focus)
            vr,vphi,vz = LabFuncs.v_infinity(v1,C,P,day)

            vr += v_off[0]
            vphi += v_off[1]
            vz += v_off[2]
            V = sqrt(vr**2.0+vphi**2.0+vz**2.0)
            F  = N*exp(-(vr**2.0/(2*sigr**2.0))\
                       -(vz**2.0/(2*sigz**2.0))\
                       -(vphi**2.0/(2*sigphi**2.0)))*(V<(v_esc))
            F = F-Fcorr
            fv1[i] = (v1**2.0)*dth*dph*sum(sum(F))
        ####################

    return fv1


def SpeedDist_Triaxial_alt(v,day,sig3,v_LSR=233.0,v_esc=528.0,\
                        v_shift=array([0.0,0.0,0.0]),GravFocus=False):
    # Exactly the same as previous function but uses a faster discretisation
    sigr = sig3[0]
    sigphi = sig3[1]
    sigz = sig3[2]

    beta = 1.0-(sigphi**2.0+sigz**2.0)/(2*sigr**2.0)
    if beta>0.0:
        N_esc = Nesc_Triaxial(sigr,sigphi,beta,v_esc)
    elif beta==0.0:
        N_esc = Nesc_Isotropic(sigr,v_esc)
    else:
        N_esc = 1.0

    v_e = LabFuncs.LabVelocitySimple(day,v_LSR=v_LSR)

    N = 1.0/(N_esc*(2*pi)**(1.5)*sigr*sigphi*sigz)
    n = size(v)
    fv1 = zeros(shape=n)

    if GravFocus==False:
        v_off = v_e-v_shift
        vv_e = sqrt(sum(v_off**2.0))

        for i in range(0,n):
            v1 = v[i]
            vr = v1*x_pix[:,0]+v_off[0]
            vphi = v1*x_pix[:,1]+v_off[1]
            vz = v1*x_pix[:,2]+v_off[2]
            V = sqrt(vr**2.0+vphi**2.0+vz**2.0)

            F  = N*exp(-(vr**2.0/(2*sigr**2.0))\
                       -(vz**2.0/(2*sigz**2.0))\
                       -(vphi**2.0/(2*sigphi**2.0)))*(V<(v_esc))
            fv1[i] = (v1**2.0)*sum(F)*dpix
    else:

        v_off = LabFuncs.v_pec+array([0.0,v_LSR,0.0])-v_shift
        vv_e = sqrt(sum(v_off**2.0))

        for i in range(0,n):
            v1 = v[i]
            vr,vphi,vz = LabFuncs.v_infinity_alt(v1*x_pix,day)

            vr += v_off[0]
            vphi += v_off[1]
            vz += v_off[2]
            V = sqrt(vr**2.0+vphi**2.0+vz**2.0)
            F  = N*exp(-(vr**2.0/(2*sigr**2.0))\
                       -(vz**2.0/(2*sigz**2.0))\
                       -(vphi**2.0/(2*sigphi**2.0)))*(V<(v_esc))
            fv1[i] = (v1**2.0)*sum(F)*dpix

    return fv1



def VelocityDist1D_Triaxial(v,day,sig3,v_LSR=233.0,v_esc=528.0,\
                        v_shift=array([0.0,0.0,0.0]),GalFrame=False,\
                      EscapeSpeed=True,SmoothCutoff=False):

    # Same as previous function but finds the 1-dimensional speed distributions
    # along each galactic cylindrical coordinate

    sigr = sig3[0]
    sigphi = sig3[1]
    sigz = sig3[2]

    beta = 1.0-(sigphi**2.0+sigz**2.0)/(2*sigr**2.0)
    N_esc = 1.0

    if not EscapeSpeed:
        N_esc = 1.0
        v_esc = 1000.0

    N = 1.0/(N_esc*(2*pi)**(1.5)*sigr*sigphi*sigz)
    n = size(v)
    fvr = zeros(shape=n)
    fvphi = zeros(shape=n)
    fvz = zeros(shape=n)

    if GalFrame:
        v_off = -v_shift
        v_max = v_esc
    else:
        v_e = LabFuncs.LabVelocitySimple(day,v_LSR=v_LSR)
        v_max = v_esc+sqrt(sum(v_e**2.0))
        v_off = v_e-v_shift

    # discretisation of speeds in each dimension
    nfine = 51
    vfine =linspace(-v_max,v_max,nfine)
    V1,V2 = meshgrid(vfine,vfine)
    dv = vfine[1]-vfine[0]

    if SmoothCutoff:
        vr = (v_max)*sqrt(1-C**2.0)*cos(P)+v_off[0]
        vphi = (v_max)*sqrt(1-C**2.0)*sin(P)+v_off[1]
        vz = (v_max)*C+v_off[2]
        V = sqrt(vr**2.0+vphi**2.0+vz**2.0)
        Fcorr = N*exp(-(vr**2.0/(2*sigr**2.0))\
                      -(vz**2.0/(2*sigz**2.0))\
                      -(vphi**2.0/(2*sigphi**2.0)))
    else:
        Fcorr = 0.0

    # Loop over speeds
    for i in range(0,n):

        # For each component create a matrix of velocities (V1,V2)
        # in the other two components and sum over those.

        # R-component
        vr = v[i]+v_off[0]
        vphi = V1+v_off[1]
        vz = V2+v_off[2]
        V = sqrt(vr**2.0+vphi**2.0+vz**2.0)
        F  = N*exp(-(vr**2.0/(2*sigr**2.0))\
                   -(vz**2.0/(2*sigz**2.0))\
                   -(vphi**2.0/(2*sigphi**2.0)))*(V<v_esc)-Fcorr
        fvr[i] = dv*dv*sum(sum(F))

        # Phi-component
        vr = V1+v_off[0]
        vphi = v[i]+v_off[1]
        vz = V2+v_off[2]
        V = sqrt(vr**2.0+vphi**2.0+vz**2.0)
        F  = N*exp(-(vr**2.0/(2*sigr**2.0))\
                   -(vz**2.0/(2*sigz**2.0))\
                   -(vphi**2.0/(2*sigphi**2.0)))*(V<v_esc)-Fcorr
        fvphi[i] = dv*dv*sum(sum(F))

        # Z-component
        vr = V1+v_off[0]
        vphi = V2+v_off[1]
        vz = v[i]+v_off[2]
        V = sqrt(vr**2.0+vphi**2.0+vz**2.0)
        F  = N*exp(-(vr**2.0/(2*sigr**2.0))\
                   -(vz**2.0/(2*sigz**2.0))\
                   -(vphi**2.0/(2*sigphi**2.0)))*(V<v_esc)-Fcorr
        fvz[i] = dv*dv*sum(sum(F))


    # Implement cutoff
    fvr[v>v_max] = 0.0
    fvphi[v>v_max] = 0.0
    fvz[v>v_max] = 0.0

    # normalise each one
    fvr /= trapz(fvr,v)
    fvphi /= trapz(fvphi,v)
    fvz /= trapz(fvz,v)

    # Stack together for output
    fv3 = vstack((fvr.T,fvphi.T,fvz.T))

    return fv3



#======================== Halo integrals ======================================#
def gvmin_Isotropic(v_min,day,v_LSR=233.0,sig=164.75,v_esc=528.0,\
                    v_shift=array([0.0,0.0,0.0]),GravFocus=False,v_exponent=-1.0):
    # MeanInverse Speed g(vmin) for an isotropic gaussian distribution
    # v_min = minimum wimp speed
    # day = day of the year, where Jan 1 = 0
    # v_exponent = the exponent to raise v_min to in the integral (-1 gives g(vmin))

    if (GravFocus) or (v_exponent!=-1.0):
        # If gravitational focusing is turned on then the integral needs to be
        # done numerically. I use flipud and then cumtrapz to do the integral
        # backwards over the range v_min_fine, which is then interpolated to
        # the input v_min

        v_min_fine = linspace(0.0001,800.0,300)
        fv = flipud((v_min_fine**v_exponent)*\
                    SpeedDist_Isotropic(v_min_fine,day,v_LSR=v_LSR,sig=sig,\
                                 v_esc=v_esc,v_shift=v_shift,GravFocus=GravFocus))
        gvmin_fine = zeros(shape=size(v_min_fine))
        gvmin_fine[0:-1] = flipud(cumtrapz(fv,v_min_fine))
        gvmin_fine[-1] = 0.0
        gvmin = interp(v_min,v_min_fine,gvmin_fine) # interp to input v_min

    else:
        # If Grav focus being ignored and exponent=-1, can use analytic result
        v_lab = LabFuncs.LabVelocitySimple(day,v_LSR=v_LSR)-v_shift
        N_esc = Nesc_Isotropic(sig,v_esc)
        v_e = sqrt(sum(v_lab**2.0))
        v0 = sig*sqrt(2.0)
        x = v_min/v0
        z = v_esc/v0
        y = v_e/v0
        gvmin = zeros(shape=shape(v_min))
        g1 = (1.0/(v0*y))
        g2 = (1.0/(2.0*N_esc*v0*y))*(erf(x+y)-erf(x-y)-(4.0/sqrt(pi))*y*exp(-z**2))
        g3 = (1.0/(2.0*N_esc*v0*y))*(erf(z)-erf(x-y)-(2.0/sqrt(pi))*(y+z-x)*exp(-z**2))
        gvmin[(x<abs(y-z))&(z<y)] = g1
        gvmin[(x<abs(y-z))&(z>y)] = g2[(x<abs(y-z))&(z>y)]
        gvmin[(abs(y-z)<x)&(x<(y+z))] = g3[(abs(y-z)<x)&(x<(y+z))]
        gvmin[(y+z)<x] = 0.0

    return gvmin


def gvmin_Triaxial(v_min,day,sig,v_LSR=233.0,v_esc=528.0,\
                  v_shift=array([0.0,0.0,0.0]),GravFocus=False,v_exponent=-1.0):
    v_min_fine = linspace(0.0001,800.0,300)
    # Same as previous but for the Triaxial Gaussian
    fv = flipud((v_min_fine**v_exponent)*\
                SpeedDist_Triaxial_alt(v_min_fine,day,sig,v_LSR=v_LSR,v_esc=v_esc,\
                                              v_shift=v_shift,GravFocus=GravFocus))
    gvmin_fine = zeros(shape=size(v_min_fine))
    gvmin_fine[0:-1] = flipud(cumtrapz(fv,v_min_fine))
    gvmin_fine[-1] = 0.0
    gvmin = interp(v_min,v_min_fine,gvmin_fine)
    return gvmin

def fhat_Isotropic(v_min,x,day,v_LSR=233.0,sig=164.75,v_esc=528.0,\
                    v_shift=array([0.0,0.0,0.0])):
    # Radon transform: the directional halo integral for isotropic Gaussian
    # v_min = min speed again
    # x = array recoil vectors
    # Important: RECOIL VECTOR x JUST NEEDS TO BE IN SAME COORDINATES AS v_lab

    v_lab = LabFuncs.LabVelocitySimple(day,v_LSR=v_LSR)-v_shift

    v0 = sig*sqrt(2.0)
    N_esc = Nesc_Isotropic(sig,v_esc)

    # dot product of v_lab and recoil vectors
    vlabdotq = (x[:,0]*v_lab[0]+x[:,1]*v_lab[1]+x[:,2]*v_lab[2]) # recoil projection

    fhat = zeros(shape=size(vlabdotq))
    fhat[((v_min+vlabdotq)<(v_esc))] = (1/(N_esc*sqrt(2*pi*sig**2.0)))\
                                        *(exp(-(v_min+vlabdotq[((v_min+vlabdotq)<(v_esc))])\
                                        **2.0/(2*sig**2.0))\
                                        -exp(-v_esc**2.0/(2*sig**2.0)))

    # needs to be divided by 2*pi, then the WIMP rate formula is independent of
    # which halo integral is used
    fhat /= 2*pi
    return fhat


def fhat_Triaxial(v_min,x,day,sig3,v_LSR=233.0,v_esc=528.0,\
                    v_shift=array([0.0,0.0,0.0])):
    # Same as previous but for a Triaxial gaussian

    v_lab = LabFuncs.LabVelocitySimple(day,v_LSR=v_LSR)-v_shift

    sigr = sig3[0]
    sigphi = sig3[1]
    sigz = sig3[2]

    # Normalisation
    beta = 1.0-(sigphi**2.0+sigz**2.0)/(2*sigr**2.0)
    if beta>0.0:
        N_esc = Nesc_Triaxial(sigr,sigphi,beta,v_esc)
    else:
        N_esc = 1.0

    vlabdotq = (x[:,0]*v_lab[0]+x[:,1]*v_lab[1]+x[:,2]*v_lab[2]) # recoil projection

    qsq = (x[:,0]*sigr)**2.0 + (x[:,1]*sigphi)**2.0 + (x[:,2]*sigz)**2.0

    fhat = zeros(shape=size(vlabdotq))
    mask = ((v_min+vlabdotq)<(v_esc))
    fhat[mask] = (1.0/(N_esc*sqrt(2*pi*qsq[mask])))*\
              (exp(-(v_min+vlabdotq[mask])**2.0/(2*qsq[mask]))-\
              exp(-v_esc**2.0/(2*qsq[mask])))

    fhat /= 2*pi
    return fhat


#==============================================================================#
