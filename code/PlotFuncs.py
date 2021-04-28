#================================PlotFuncs.py==================================#
# Created by Ciaran O'Hare 2020

# Description:
# This file has many functions which are used throughout the project, but are
# all focused around the bullshit that goes into making the plots

#==============================================================================#

from numpy import *
from numpy.random import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib import colors
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from scipy.stats import zscore,chi2,multivariate_normal
from scipy.special import erfinv
from scipy.stats import gaussian_kde
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import norm
from Params import *

pltdir = '../plots/'
pltdir_png = pltdir+'plots_png/'
#==============================================================================#
# My custom figure save
def MySaveFig(fig,pltname,pngsave=True):
    fig.savefig(pltdir+pltname+'.pdf',bbox_inches='tight')
    if pngsave:
        fig.savefig(pltdir_png+pltname+'.png',bbox_inches='tight')
#==============================================================================#


#==============================================================================#
# My preferred style of colorbar
def cbar(mappable,extend='neither',\
                label='',lfs=35,labelpad=40,labelrotation=-90):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = fig.colorbar(mappable, cax=cax,extend=extend)
    cbar.ax.tick_params(which='minor',length=8,width=2)
    cbar.ax.tick_params(which='major',length=10,width=3)
    cbar.set_label(label,fontsize=lfs,rotation=labelrotation,labelpad=labelpad)
    return cbar
#==============================================================================#



#==============================================================================#
def MySquarePlot(xlab='',ylab='',\
                 lw=2.5,lfs=45,tfs=25,size_x=13,size_y=12,Grid=False):
    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{mathpazo}']

    fig = plt.figure(figsize=(size_x,size_y))
    ax = fig.add_subplot(111)

    ax.set_xlabel(xlab,fontsize=lfs)
    ax.set_ylabel(ylab,fontsize=lfs)

    ax.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)
    if Grid:
        ax.grid()
    return fig,ax

def MyDoublePlot(xlab1='',ylab1='',xlab2='',ylab2='',\
                 wspace=0.25,lw=2.5,lfs=45,tfs=25,size_x=20,size_y=11,Grid=False):
    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{mathpazo}']


    fig, axarr = plt.subplots(1, 2,figsize=(size_x,size_y))
    gs = gridspec.GridSpec(1, 2)
    gs.update(wspace=wspace)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax1.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax1.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)
    ax2.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax2.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax1.set_xlabel(xlab1,fontsize=lfs)
    ax1.set_ylabel(ylab1,fontsize=lfs)

    ax2.set_xlabel(xlab2,fontsize=lfs)
    ax2.set_ylabel(ylab2,fontsize=lfs)

    if Grid:
        ax1.grid()
        ax2.grid()
    return fig,ax1,ax2

def MyDoublePlot_Vertical(xlab1='',ylab1='',xlab2='',ylab2='',\
                 hspace=0.25,lw=2.5,lfs=45,tfs=35,size_x=20,size_y=10,Grid=False):
    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{mathpazo}']


    fig, axarr = plt.subplots(2,1,figsize=(size_x,size_y))
    gs = gridspec.GridSpec(2, 1)
    gs.update(hspace=hspace)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    ax1.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax1.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax2.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax2.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax1.set_xlabel(xlab1,fontsize=lfs)
    ax1.set_ylabel(ylab1,fontsize=lfs)

    ax2.set_xlabel(xlab2,fontsize=lfs)
    ax2.set_ylabel(ylab2,fontsize=lfs)


    if Grid:
        ax1.grid()
        ax2.grid()
    return fig,ax1,ax2

def MyTriplePlot(xlab1='',ylab1='',xlab2='',ylab2='',xlab3='',ylab3='',\
                 wspace=0.25,lw=2.5,lfs=45,tfs=25,size_x=20,size_y=7,Grid=False,width_ratios=[1,1,1]):

    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{mathpazo}']

    fig, axarr = plt.subplots(1, 3,figsize=(size_x,size_y))
    gs = gridspec.GridSpec(1, 3,width_ratios=width_ratios)
    gs.update(wspace=wspace)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])

    ax1.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax1.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax2.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax2.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax3.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax3.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax1.set_xlabel(xlab1,fontsize=lfs)
    ax1.set_ylabel(ylab1,fontsize=lfs)

    ax2.set_xlabel(xlab2,fontsize=lfs)
    ax2.set_ylabel(ylab2,fontsize=lfs)

    ax3.set_xlabel(xlab3,fontsize=lfs)
    ax3.set_ylabel(ylab3,fontsize=lfs)

    if Grid:
        ax1.grid()
        ax2.grid()
        ax3.grid()
    return fig,ax1,ax2,ax3
#==============================================================================#



def MyTriplePlot_Vertical(xlab1='',ylab1='',xlab2='',ylab2='',xlab3='',ylab3='',\
                 hspace=0.25,lw=2.5,lfs=45,tfs=35,size_x=20,size_y=15,Grid=False):
    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{mathpazo}']


    fig, axarr = plt.subplots(3,1,figsize=(size_x,size_y))
    gs = gridspec.GridSpec(3, 1)
    gs.update(hspace=hspace)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])

    ax1.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax1.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax2.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax2.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax3.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax3.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    ax1.set_xlabel(xlab1,fontsize=lfs)
    ax1.set_ylabel(ylab1,fontsize=lfs)

    ax2.set_xlabel(xlab2,fontsize=lfs)
    ax2.set_ylabel(ylab2,fontsize=lfs)

    ax3.set_xlabel(xlab3,fontsize=lfs)
    ax3.set_ylabel(ylab3,fontsize=lfs)

    if Grid:
        ax1.grid()
        ax2.grid()
        ax3.grid()
    return fig,ax1,ax2,ax3
#==============================================================================#







#==============================================================================#
# Colormap stuff:

from copy import copy
def cmap_setunderwhite(cmap):
    cmap_c = copy(plt.get_cmap(cmap))
    cmap_c.set_under('white', 1.0)
    return cmap_c

def reverse_colourmap(cmap, name = 'my_cmap_r'):
    reverse = []
    k = []

    for key in cmap._segmentdata:
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:
            data.append((1-t[0],t[2],t[1]))
        reverse.append(sorted(data))

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL)
    return my_cmap_r


def col_alpha(col,alpha=0.1):
    rgb = colors.colorConverter.to_rgb(col)
    bg_rgb = [1,1,1]
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(rgb, bg_rgb)]
#==============================================================================#


#==============================================================================#
# Set up axion plot as desired:
def AxionPlotSetup(CAST_text_on=True,Haloscopes_text_on=True,QCD_text_on=True,\
                        StellarBounds_text_on=True):
    fig,ax = AxionPhoton.FigSetup(Shape='Square',\
                    ylab='$|g_{a\gamma}|$ [GeV$^{-1}$]',mathpazo=True,\
                    g_min=1e-12,g_max=1e-9,m_min=1e-6,m_max=1e0,tfs=35)

    AxionPhoton.QCDAxion(ax,text_on=False,thick_lines=True)
    AxionPhoton.StellarBounds(ax,text_on=False)
    AxionPhoton.Helioscopes(ax,projection=False,text_on=False,line_alpha=1)

    m_min=1e-6
    m_max=1e0
    plt.plot([m_min,m_max],[2e-10*0.1*m_min,2e-10*0.1*m_max],'-',\
                        color='brown',lw=1.5,zorder=0,alpha=0.1)
    plt.plot([m_min,m_max],[2e-10*10*m_min,2e-10*10*m_max],'-',\
                        color='brown',lw=1.5,zorder=0,alpha=0.1)

    fs = 25
    AxionPhoton.NeutronStars(ax,col=[0.0, 0.66, 0.42],fs=15,\
                            RescaleByMass=False,text_on=False)
    AxionPhoton.ADMX(ax,col='darkred',fs=fs,text_on=False)
    AxionPhoton.RBF_UF(ax,col='darkred',fs=fs-2,text_on=False)
    AxionPhoton.HAYSTAC(ax,col='darkred',text_on=False)
    AxionPhoton.CAPP(ax,col='darkred',fs=fs-4,text_on=False)
    AxionPhoton.ORGAN(ax,col='darkred',text_on=False)


    if CAST_text_on:
        plt.gcf().text(0.3,0.75,r'{\bf CAST}',color='k',fontsize=50)
        plt.gcf().text(0.295,0.75,r'{\bf CAST}',color='w',fontsize=50)

    if StellarBounds_text_on:
        plt.gcf().text(0.89,0.59,r'{\bf Stellar bounds}',color='k',\
                            fontsize=30,horizontalalignment='right')
        plt.gcf().text(0.888,0.59,r'{\bf Stellar bounds}',color='w',\
                            fontsize=30,horizontalalignment='right')

    if Haloscopes_text_on:
        plt.gcf().text(0.188,0.135,r'{\bf Haloscopes}',color='k',\
                        rotation=90,fontsize=32,rotation_mode='anchor',zorder=0)
        plt.gcf().text(0.185,0.1351,r'{\bf Haloscopes}',color='w',\
                        rotation=90,fontsize=32,rotation_mode='anchor',zorder=0)

    if QCD_text_on:
        cols = cm.get_cmap('YlOrBr')
        plt.gcf().text(0.561,0.13,r'{\bf KSVZ}',rotation=63,\
                                fontsize=32,rotation_mode='anchor')
        plt.gcf().text(0.611,0.13,r'{\bf DFSZ II}',rotation=63,\
                                fontsize=32,rotation_mode='anchor')
        plt.gcf().text(0.56,0.13,r'{\bf KSVZ}',rotation=63,fontsize=32,\
                                    rotation_mode='anchor',color=cols(0.9))
        plt.gcf().text(0.61,0.13,r'{\bf DFSZ II}',rotation=63,fontsize=32,\
                                    rotation_mode='anchor',color=cols(0.9))


    plt.title('space',color='w',fontsize=15)

    return fig,ax


#==============================================================================#
class AxionPhoton():
    def FigSetup(xlab=r'$m_a$ [eV]',ylab='',\
                     g_min = 1.0e-19,g_max = 1.0e-6,\
                     m_min = 1.0e-12,m_max = 1.0e7,\
                     lw=2.5,lfs=45,tfs=25,tickdir='out',\
                     Grid=False,Shape='Rectangular',mathpazo=False,TopAndRightTicks=False):

            plt.rcParams['axes.linewidth'] = lw
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif',size=tfs)

            if mathpazo:
                mpl.rcParams['text.latex.preamble'] = [r'\usepackage{mathpazo}']

            if Shape=='Wide':
                fig = plt.figure(figsize=(16.5,5))
            elif Shape=='Rectangular':
                fig = plt.figure(figsize=(16.5,11))
            elif Shape=='Square':
                fig = plt.figure(figsize=(16,16))

            ax = fig.add_subplot(111)

            ax.set_xlabel(xlab,fontsize=lfs)
            ax.set_ylabel(ylab,fontsize=lfs)

            ax.tick_params(which='major',direction=tickdir,width=2.5,length=13,right=TopAndRightTicks,top=TopAndRightTicks,pad=7)
            ax.tick_params(which='minor',direction=tickdir,width=1,length=10,right=TopAndRightTicks,top=TopAndRightTicks)


            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlim([m_min,m_max])
            ax.set_ylim([g_min,g_max])

            locmaj = mpl.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=50)
            locmin = mpl.ticker.LogLocator(base=10.0, subs=arange(2, 10)*.1,numticks=100)
            ax.xaxis.set_major_locator(locmaj)
            ax.xaxis.set_minor_locator(locmin)
            ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

            locmaj = mpl.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
            locmin = mpl.ticker.LogLocator(base=10.0, subs=arange(2, 10)*.1,numticks=100)
            ax.yaxis.set_major_locator(locmaj)
            ax.yaxis.set_minor_locator(locmin)
            ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

            if Shape=='Rectangular':
                plt.xticks(rotation=20)

            if Grid:
                ax.grid(zorder=0)
            return fig,ax

    def QCDAxion(ax,coupling='Photon',
                      C_logwidth=10,KSVZ_on=True,DFSZ_on=True,
                      cmap='YlOrBr',fs=18,RescaleByMass=False,text_on=True,thick_lines=False):
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        ## QCD Axion band:
        g_min,g_max = ax.get_ylim()
        m_min,m_max = ax.get_xlim()

        # Mass-coupling relation
        def g_x(C_ag,m_a):
            return 2e-10*C_ag*m_a
        KSVZ = 1.92
        DFSZ = 0.75

        if rs1==0:
            # Plot Band
            n = 200
            g = logspace(log10(g_min),log10(g_max),n)
            m = logspace(log10(m_min),log10(m_max),n)
            QCD = zeros(shape=(n,n))
            for i in range(0,n):
                QCD[:,i] = norm.pdf(log10(g)-log10(g_x(1.0,m[i])),0.0,0.8)
            cols = cm.get_cmap(cmap)

            cols.set_under('w') # Set lowest color to white
            vmin = amax(QCD)/(C_logwidth/4.6)
            plt.contourf(m, g, QCD, 50,cmap=cols,vmin=vmin,vmax=0.9,zorder=0)
            plt.contourf(m, g, QCD, 50,cmap=cols,vmin=vmin,vmax=0.9,zorder=0)
            plt.contourf(m, g, QCD, 50,cmap=cols,vmin=vmin,vmax=0.9,zorder=0)

            # QCD Axion models
            rot = 45.0
            trans_angle = plt.gca().transData.transform_angles(array((rot,)),array([[0, 0]]))[0]
            m2 = array([1e-9,5e-8])
            if KSVZ_on:
                if thick_lines:
                    plt.plot(m,g_x(KSVZ,m),'-',linewidth=5,color='k',zorder=0)
                    plt.plot(m,g_x(KSVZ,m),'-',linewidth=3,color=cols(0.7),zorder=0)
                else:
                    plt.plot(m,g_x(KSVZ,m),'-',linewidth=2,color=cols(1.0),zorder=0)
                if text_on:
                    plt.text(1e-8,g_x(KSVZ,1e-8)*1.05,r'{\bf KSVZ}',fontsize=fs,rotation=trans_angle,color=cols(1.0),ha='left',va='bottom',rotation_mode='anchor')
            if DFSZ_on:
                if thick_lines:
                    plt.plot(m,g_x(DFSZ,m),'-',linewidth=5,color='k',zorder=0)
                    plt.plot(m,g_x(DFSZ,m),'-',linewidth=3,color=cols(0.7),zorder=0)
                else:
                    plt.plot(m,g_x(DFSZ,m),'-',linewidth=2,color=cols(1.0),zorder=0)
                if text_on:
                    plt.text(5e-8,g_x(DFSZ,5e-8)/1.5,r'{\bf DFSZ II}',fontsize=fs,rotation=trans_angle,color=cols(1.0),ha='left',va='top',rotation_mode='anchor')
        else:
            C_min,C_max = ax.get_ylim()
            n = 200
            C = logspace(log10(C_min),log10(C_max),n)
            m = logspace(log10(m_min),log10(m_max),n)
            QCD = zeros(shape=(n,n))
            for i in range(0,n):
                QCD[:,i] = norm.pdf(log10(C),0.0,0.8)
            cols = cm.get_cmap(cmap)
            cols.set_under('w') # Set lowest color to white
            vmin = amax(QCD)/(C_logwidth/2)
            plt.contourf(m, C, QCD, 50,cmap=cols,vmin=vmin,vmax=0.9,zorder=0)
            plt.contourf(m, C, QCD, 50,cmap=cols,vmin=vmin,vmax=0.9,zorder=0)
            plt.contourf(m, C, QCD, 50,cmap=cols,vmin=vmin,vmax=0.9,zorder=0)
            if thick_lines:
                plt.plot([1e-9,1e0],[0.75,0.75],'-',lw=5,color='k')
                plt.plot([1e-9,1e0],[0.75,0.75],'-',lw=3,color='k')
            else:
                plt.plot([1e-9,1e0],[0.75,0.75],'-',lw=2,color='k')
            if text_on:
                plt.text(1e-2,0.75/3,r'{\bf DFSZ II}',fontsize=fs,color='k')
        return


    def ADMX(ax,col=[0.8, 0.0, 0.0],projection=False,fs=15,RescaleByMass=False,text_on=True):
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        # 2018: arXiv[1804.05750]
        # 2019: arXiv[1910.08638]
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/ADMX.txt")
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=0.1)
        dat = loadtxt(limit_dir+"AxionPhoton/ADMX2018.txt")
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=0.1)
        dat = loadtxt(limit_dir+"AxionPhoton/ADMX2019_1.txt")
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=0.1)
        dat = loadtxt(limit_dir+"/AxionPhoton/ADMX2019_2.txt")
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=0.1)
        dat = loadtxt(limit_dir+"/AxionPhoton/ADMX_Sidecar.txt")
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=0.1)


        if projection:
            # ADMX arXiv[1804.05750]
            dat = loadtxt(limit_dir+"/AxionPhoton/Projections/ADMX_Projected.txt")
            plt.plot(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),'-',linewidth=1.5,color=col,zorder=0)
            plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=0,alpha=0.1)
            if rs1==0:
                if text_on:
                    plt.text(4e-5,9e-16,r'{\bf ADMX}',fontsize=20,color=col,rotation=0,ha='left',va='top')
                plt.plot([4e-5,3e-5],[9e-16,2.1e-15],'k-',lw=1.5)
            else:
                plt.text(0.9e-6,0.15,r'{\bf ADMX}',fontsize=fs,color=col,rotation=0,ha='left',va='top')
        else:
            if rs1==0:
                if text_on:
                    plt.text(0.7e-6,1e-13,r'{\bf ADMX}',fontsize=fs,color=col,rotation=90,ha='left',va='top')

        return


    def NeutronStars(ax,col=[0.1, 0.5, 0.2],fs=15,RescaleByMass=False,text_on=True):
        # Neutron stars arXiv:[2004.00011]
        y2 = ax.get_ylim()[1]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        dat = loadtxt(limit_dir+'/AxionPhoton/NeutronStars.txt')
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=0.11)
        if rs1==0:
            if text_on:
                plt.text(5e-6,1e-12,r'{\bf Neutron stars}',fontsize=fs,color='w',ha='left')
        else:
            if text_on:
                plt.text(1e-7,4e3,r'{\bf Neutron}',fontsize=fs,color=col,ha='center')
                plt.text(1e-7,1e3,r'{\bf stars}',fontsize=fs,color=col,ha='center')
            plt.plot([3.5e-7,7e-6],[6e3,2e4],lw=1.5,color=col)

    def RBF_UF(ax,col ='darkred',fs=13,RescaleByMass=False,text_on=True):
        # UF: Phys. Rev. D42, 1297 (1990).
        # RBF: Phys. Rev. Lett. 59, 839 (1987).
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/RBF_UF_Haloscopes.txt")
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=0.1)
        if rs1==0:
            if text_on:
                plt.text(0.4e-5,3e-11,r'{\bf RBF+UF}',fontsize=fs,color='w',rotation=-90,ha='left',va='top')
        else:
            if text_on:
                plt.text(0.7e-5,4e3,r'{\bf RBF}',fontsize=fs,color='w',rotation=0,ha='center',va='top')
                plt.text(0.7e-5,1e3,r'{\bf UF}',fontsize=fs,color='w',rotation=0,ha='center',va='top')

        return

    def HAYSTAC(ax,col=[0.88, 0.07, 0.37],fs=13,RescaleByMass=False,projection=True,text_on=True):
        # HAYSTAC arXiv:[1803.03690]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
            zo = 3
        else:
            rs1 = 0.0
            rs2 = 1.0
            zo = 0
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/HAYSTAC.txt")
        if rs1==0:
            plt.plot([dat[0,0],dat[0,0]],[dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),y2/(rs1*2e-10*dat[0,0]+rs2)],color=col,zorder=zo,lw=3)
            if text_on:
                if projection:
                    plt.text(2.4e-5,4e-12,r'{\bf HAYSTAC}',fontsize=fs,color=col,rotation=-90,ha='left',va='top')
                else:
                    plt.text(2.4e-5,5e-13,r'{\bf HAYSTAC}',fontsize=fs,color=col,rotation=-90,ha='left',va='top')
        else:
            plt.plot([dat[0,0],dat[0,0]],[dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),y2/(rs1*2e-10*dat[0,0]+rs2)],color='k',zorder=zo,lw=4)
            plt.plot([dat[0,0],dat[0,0]],[dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),y2/(rs1*2e-10*dat[0,0]+rs2)],color=col,zorder=zo,lw=3)
            if text_on:
                plt.text(dat[0,0],y2*1.2,r'{\bf HAYSTAC}',fontsize=fs,color=col,rotation=40,ha='left',rotation_mode='anchor')
            plt.plot(dat[0,0],dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),'.',markersize=15,color=col,markeredgecolor='k',zorder=zo)
        return

    def CAPP(ax,col=[1, 0.1, 0.37],fs=15,RescaleByMass=False,text_on=True):
        # CAPP arXiv:[2001.05102]
        y2 = ax.get_ylim()[1]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
            zo = 3
        else:
            rs1 = 0.0
            rs2 = 1.0
            zo = 0
        dat = loadtxt(limit_dir+"/AxionPhoton/CAPP-8TB.txt")
        if rs1==0:
            plt.plot([dat[0,0],dat[0,0]],[dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),y2/(rs1*2e-10*dat[0,0]+rs2)],color=col,zorder=zo,lw=3)
            if text_on:
                plt.text(1e-5,1e-13,r'{\bf CAPP}',fontsize=fs,color=col,rotation=-90,ha='center',va='top')
        else:
            plt.plot([dat[0,0],dat[0,0]],[dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),y2/(rs1*2e-10*dat[0,0]+rs2)],color='k',zorder=zo,lw=4)
            plt.plot([dat[0,0],dat[0,0]],[dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),y2/(rs1*2e-10*dat[0,0]+rs2)],color=col,zorder=zo,lw=3)
            if text_on:
                plt.text(dat[0,0],y2*1.8,r'{\bf CAPP}',fontsize=fs,color=col,rotation=40,ha='left',va='top',rotation_mode='anchor')
            plt.plot(dat[0,0],dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),'.',markersize=15,color=col,markeredgecolor='k',zorder=zo)
        return

    def QUAX(ax,col='crimson',fs=15,RescaleByMass=False,text_on=True):
        # QUAX arXiv:[1903.06547]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
            zo = 3
        else:
            rs1 = 0.0
            rs2 = 1.0
            zo = 0
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/QUAX.txt")

        if rs1==0:
            plt.plot([dat[0,0],dat[0,0]],[dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),y2/(rs1*2e-10*dat[0,0]+rs2)],color=col,lw=2,zorder=zo)
            if text_on:
                plt.text(5.2e-5,4e-11,r'{\bf QUAX}',fontsize=fs,color=col,rotation=-90,ha='center',va='top')
        else:
            plt.plot([dat[0,0],dat[0,0]],[dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),y2/(rs1*2e-10*dat[0,0]+rs2)],color='k',lw=4,zorder=zo)
            plt.plot([dat[0,0],dat[0,0]],[dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),y2/(rs1*2e-10*dat[0,0]+rs2)],color=col,lw=3,zorder=zo)
            if text_on:
                plt.text(dat[0,0]*1.3,y2*1.2,r'{\bf QUAX}',fontsize=fs,color=col,rotation=40,ha='left',rotation_mode='anchor')
            plt.plot(dat[0,0],dat[0,1]/(rs1*2e-10*dat[0,0]+rs2),'.',markersize=15,color=col,markeredgecolor='k',zorder=zo)
        return

    def ABRACADABRA(ax,col=[0.83, 0.07, 0.37],fs=15,projection=False,RescaleByMass=False,text_on=True):
        # ABRACADABRA arXiv:[1810.12257]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/ABRACADABRA.txt")
        n = shape(dat)[0]
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=2)
        x = dat[arange(0,n,20),0]
        y = dat[arange(0,n,20),1]
        y[-1] = y2
        plt.plot(x,y/(rs1*2e-10*x+rs2),'k-',lw=1,zorder=10,alpha=0.5)
        if text_on:
            if rs1==0:
                plt.text(1.5e-9,3e-8,r'{\bf ABRA}',fontsize=fs,color='w',rotation=0,ha='center',va='top',zorder=10)
                plt.text(1.5e-9,1e-8,r'10 cm',fontsize=fs,color='w',rotation=0,ha='center',va='top',zorder=10)

        if projection:
            dat = loadtxt(limit_dir+"/AxionPhoton/Projections/ABRACADABRA.txt")
            plt.plot(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),'-',linewidth=1.5,color=col,zorder=0)
            plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=0,alpha=0.1)
            if rs1==0:
                if text_on:
                    plt.text(1e-12,2.5e-18,r'{\bf ABRACADABRA}',fontsize=fs-1,color=col,rotation=13,ha='left',va='top')
            else:
                if text_on:
                    plt.text(1.3e-9,1.0e2,r'{\bf ABRACADABRA}',fontsize=fs-1,color=col,rotation=0,ha='left',va='top')
                plt.plot([dat[-1,0],dat[-1,0]],[dat[-1,1]/(rs1*2e-10*dat[-1,0]+rs2),1e6],lw=1.5,color=col,zorder=0)
        return

    def ORGAN(ax,col='crimson',projection=False,fs=15,RescaleByMass=False,text_on=True):
        # ORGAN arXiv[1706.00209]
        y2 = ax.get_ylim()[1]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        col = [0.8, 0.0, 0.0]
        dat = loadtxt(limit_dir+"/AxionPhoton/ORGAN.txt")
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=col,facecolor=col,zorder=0.1,lw=2)

        if projection:
            dat = loadtxt(limit_dir+"/AxionPhoton/Projections/ORGAN_Projected.txt")
            plt.plot(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),'-',linewidth=1.5,color=col,zorder=0)
            plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=0,alpha=0.1)
            if rs1==0:
                if text_on:
                    plt.text(5e-4,1.15e-14,r'{\bf ORGAN}',fontsize=18,color=col,rotation=0,ha='left',va='top')
                plt.plot([5e-4,1.5e-4],[1.3e-14,6e-13],'k-',lw=1.5)
            else:
                plt.text(1.2e-4,1e3,r'{\bf ORGAN}',fontsize=18,color=col,rotation=-90,ha='left',va='top')

        else:
            if rs1==0:
                if text_on:
                    plt.text(110e-6,6e-11,r'{\bf ORGAN}',fontsize=fs,color=col,rotation=-90,ha='left',va='top')
        return


    def MADMAX(ax,col=[0.6, 0.1, 0.1],fs=18,RescaleByMass=False,text_on=True):
        # MADMAX arXiv[2003.10894]
        y2 = ax.get_ylim()[1]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        dat = loadtxt(limit_dir+"/AxionPhoton/Projections/MADMAX.txt")
        plt.plot(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),'-',linewidth=1.5,color=col,zorder=0)
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=0,alpha=0.2)
        if rs1==0:
            plt.text(1.5e-4,3.5e-15,r'{\bf MADMAX}',fontsize=18,color=col,rotation=0,ha='left',va='top')
            plt.plot([3e-4,1.3e-4],[4.5e-15,2.6e-14],'k-',lw=1.5)
        else:
            if text_on:
                plt.text(4e-5,3.5e-1,r'{\bf MADMAX}',fontsize=fs,color=col,rotation=0,ha='left',va='top')
        return

    def KLASH(ax,col=[0.6, 0.1, 0.2],fs=15,RescaleByMass=False,text_on=True):
        # KLASH arXiv:[1707.06010]
        y2 = ax.get_ylim()[1]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        dat = loadtxt(limit_dir+"/AxionPhoton/Projections/KLASH.txt")
        plt.plot(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),'-',linewidth=1.5,color=col,zorder=0)
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,facecolor=col,zorder=0,alpha=0.3)
        if text_on:
            if rs1==0:
                plt.text(1e-7,1e-12,r'{\bf KLASH}',rotation=90,fontsize=fs,color=col,ha='left',va='top')
            else:
                plt.text(2.5e-7,1.3e0,r'{\bf KLASH}',rotation=90,fontsize=fs,color=col,ha='left',va='top',rotation_mode='anchor')

    def BRASS(ax,col=[0.5, 0.1, 0.2],fs=15,RescaleByMass=False,text_on=True):
        # BRASS http://www.iexp.uni-hamburg.de/groups/astroparticle/brass/brassweb.htm
        y2 = ax.get_ylim()[1]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        dat = loadtxt(limit_dir+"/AxionPhoton/Projections/BRASS.txt")
        plt.plot(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),'-',linewidth=1.5,color=col,zorder=0)
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,facecolor=col,zorder=0,alpha=0.1)
        if text_on:
            if rs1==0:
                plt.text(2.3e-3,0.6e-10,r'{\bf BRASS}',rotation=56,fontsize=fs,color=col,ha='left',va='top')
            else:
                plt.text(1e-3,0.12e3,r'{\bf BRASS}',rotation=15,fontsize=fs,color=col,ha='left',va='top')


    def TOORAD(ax,col=[0.8, 0.1, 0.2],fs=15,RescaleByMass=False,text_on=True):
        # TOORAD arXiv[1807.08810]
        y2 = ax.get_ylim()[1]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        dat = loadtxt(limit_dir+"/AxionPhoton/Projections/TOORAD.txt")
        plt.plot(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),'-',linewidth=1.5,color=col,zorder=0)
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,facecolor=col,zorder=0,alpha=0.1)
        if rs1==0:
            if text_on:
                plt.text(0.5e-2,5e-14,r'{\bf TOORAD}',rotation=0,fontsize=18,color=col,ha='left',va='top')
            plt.plot([0.5e-2,0.21e-2],[5e-14,1e-13],'k-',lw=1.5)
        else:
            if text_on:
                plt.text(0.6e-3,4e-1,r'{\bf TOORAD}',rotation=-25,fontsize=fs,color=col,ha='left',va='top')

    def LAMPOST(ax,col=[0.8, 0.1, 0.2],fs=15,RescaleByMass=False,text_on=True):
        # LAMPOST arXiv[1803.11455]
        y2 = ax.get_ylim()[1]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        dat = loadtxt(limit_dir+"/AxionPhoton/Projections/LAMPOST.txt",delimiter=',')
        plt.plot(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),'-',linewidth=1.5,color=col,zorder=0)
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,facecolor=col,zorder=0,alpha=0.1)

        if rs1==0:
            if text_on:
                plt.text(0.8e-1,5e-12,r'{\bf LAMPOST}',rotation=-90,fontsize=fs,color=col,ha='left',va='top')
        else:
            if text_on:
                plt.text(0.9e-1,1.9e-1,r'{\bf LAMPOST}',rotation=0,fontsize=fs,color=col,ha='left',va='top')


    # Low mass ALP haloscopes
    def DANCE(ax,col=[0.8, 0.1, 0.2],fs=15,text_on=True):
        # DANCE arXiv[1911.05196]
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/Projections/DANCE.txt")
        plt.plot(dat[:,0],dat[:,1],'-',linewidth=1.5,color=col,zorder=0)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0,alpha=0.1)
        if text_on:
            plt.text(1.7e-12,2e-13,r'{\bf DANCE}',rotation=50,fontsize=fs,color=col,ha='left',va='top')

    def aLIGO(ax,col=[0.8, 0.1, 0.2],fs=15,text_on=True):
        # aLIGO arXiv[1903.02017]
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/Projections/aLIGO.txt")
        plt.plot(dat[:,0],dat[:,1],'-',linewidth=1.5,color=col,zorder=0)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0,alpha=0.1)
        if text_on:
            plt.text(0.2e-9,0.35e-13,r'{\bf aLIGO}',rotation=0,fontsize=fs,color=col,ha='left',va='top')

    def ADBC(ax,col=[0.8, 0.1, 0.2],fs=15,text_on=True):
        # ADBC arXiv[1809.01656]
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/Projections/ADBC.txt")
        plt.plot(dat[:,0],dat[:,1],'-',linewidth=1.5,color=col,zorder=0)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0,alpha=0.1)
        if text_on:
            plt.text(1e-11,1.3e-12,r'{\bf ADBC}',rotation=26,fontsize=fs,color=col,ha='left',va='top')

    def SHAFT(ax,col='red',fs=16,text_on=True):
        # SHAFT arXiv:[2003.03348]
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/SHAFT.txt")
        n = shape(dat)[0]
        x = dat[arange(0,n,2),0]
        y = dat[arange(0,n,2),1]
        y[-1] = y2
        plt.plot(x,y,'k-',lw=1,zorder=10,alpha=0.5)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.8)
        if text_on:
            plt.text(0.8e-10,3e-10,r'{\bf SHAFT}',fontsize=fs,color='w',rotation=0,ha='center',va='top',zorder=9)
        return

    def ALPS(ax,projection=True,col=[0.8, 0.25, 0.33],fs=15,RescaleByMass=False,text_on=True):
        # ALPS-I arXiv:[1004.1313]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0

        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/ALPS.txt")
        plt.plot(dat[:,0],dat[:,1],'k-',lw=2.5,zorder=1.53,alpha=0.5)
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor=None,facecolor=col,zorder=1.53,lw=0.01)
        if text_on:
            plt.text(1e-5,7e-8,r'{\bf ALPS-I}',fontsize=20,color='w')
        if projection:
            dat = loadtxt(limit_dir+"/AxionPhoton/Projections/ALPS-II.txt")
            plt.plot(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),'-',lw=1.5,zorder=1.5,color='k',alpha=0.5)
            if text_on:
                if RescaleByMass:
                    plt.text(9e-4,2.5e3,r'{\bf ALPS-II}',fontsize=20,color='k',rotation=20,alpha=0.5)
                else:
                    plt.text(1.5e-3,3e-9,r'{\bf ALPS-II}',rotation=60,fontsize=18,color='w',zorder=10)
        return

    def OSQAR(ax,col=[0.6, 0.2, 0.25],fs=15,text_on=True):
        # OSQAR arXiv:[]
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/OSQAR.txt")
        plt.plot(dat[:,0],dat[:,1],'k-',lw=2.5,zorder=1.52,alpha=0.5)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.52,lw=0.01)
        if text_on:
            plt.text(1e-5,1.5e-8,r'{\bf OSQAR}',fontsize=17,color='w')
        return

    def PVLAS(ax,col=[0.4, 0.2, 0.2],fs=15,text_on=True):
        # PVLAS arXiv:[]
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/PVLAS.txt")
        plt.plot(dat[:,0],dat[:,1],'k-',lw=2.5,zorder=1.51,alpha=0.4)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.51,lw=0.01)
        if text_on:
            plt.text(2e-3,9e-8,r'{\bf PVLAS}',fontsize=17,color='w',rotation=45)
        return


    def CROWS(ax,col=[0.7, 0.2, 0.2],fs=15,text_on=True):
        # CROWS arXiv:[1310.8098]
        y2 = ax.get_ylim()[1]
        dat = loadtxt(limit_dir+"/AxionPhoton/CROWS.txt")
        plt.plot(dat[:,0],dat[:,1],'k-',lw=2.5,zorder=1.54,alpha=0.4)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.54,lw=0.01)
        if text_on:
            plt.text(1e-7,1.5e-7,r'{\bf CROWS}',fontsize=17,color='w',rotation=0)
        return


    ####################################################
    def Helioscopes(ax,col=[0.5, 0.0, 0.13],fs=25,projection=True,RescaleByMass=False,text_on=True,line_alpha=0.5):
        # CAST arXiv:[1705.02290]
        y2 = ax.get_ylim()[1]
        if RescaleByMass:
            rs1 = 1.0
            rs2 = 0.0
        else:
            rs1 = 0.0
            rs2 = 1.0
        dat = loadtxt(limit_dir+"/AxionPhoton/CAST_highm.txt")
        plt.plot(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),'k-',lw=2,zorder=1.49,alpha=line_alpha)
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor='k',facecolor=col,zorder=1.49,lw=0.1)
        dat = loadtxt(limit_dir+"/AxionPhoton/CAST.txt")
        plt.plot(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),'k-',lw=2,zorder=1.5,alpha=line_alpha)
        plt.fill_between(dat[:,0],dat[:,1]/(rs1*2e-10*dat[:,0]+rs2),y2=y2,edgecolor='k',facecolor=col,zorder=1.5,lw=0.1)
        if text_on:
            if rs1==0:
                plt.text(1e-1,1.5e-9,r'{\bf CAST}',fontsize=fs+4,color='w',rotation=0,ha='center',va='top')
            else:
                plt.text(4e-2,5e3,r'{\bf CAST}',fontsize=fs+4,color='w',rotation=0,ha='center',va='top')

        if projection:
            # IAXO arXiv[1212.4633]
            IAXO_col = 'purple'
            IAXO = loadtxt(limit_dir+"/AxionPhoton/Projections/IAXO.txt")
            plt.plot(IAXO[:,0],IAXO[:,1]/(rs1*2e-10*IAXO[:,0]+rs2),'--',linewidth=2.5,color=IAXO_col,zorder=0.5)
            plt.fill_between(IAXO[:,0],IAXO[:,1]/(rs1*2e-10*IAXO[:,0]+rs2),y2=y2,edgecolor=None,facecolor=IAXO_col,zorder=0,alpha=0.3)
            if text_on:
                if rs1==0:
                    plt.text(0.35e-1,0.2e-11,r'{\bf IAXO}',fontsize=fs,color=IAXO_col,rotation=45)
                else:
                    plt.text(0.7e-2,0.15e1,r'{\bf IAXO}',fontsize=fs,color=IAXO_col,rotation=-20)
        return


    def Haloscopes(ax,projection=True,fs=20,text_on=True):
        AxionPhoton.ADMX(ax,projection=projection,fs=fs,text_on=text_on)
        AxionPhoton.RBF_UF(ax,fs=fs-2,text_on=text_on)
        AxionPhoton.HAYSTAC(ax,projection=projection,text_on=text_on)
        AxionPhoton.ABRACADABRA(ax,fs=fs,projection=projection,text_on=text_on)
        AxionPhoton.SHAFT(ax,text_on=text_on)
        AxionPhoton.CAPP(ax,fs=fs-4,text_on=text_on)
        AxionPhoton.ORGAN(ax,projection=projection,text_on=text_on)

        if projection:
            AxionPhoton.MADMAX(ax,text_on=text_on)
            AxionPhoton.KLASH(ax,text_on=text_on)
            AxionPhoton.TOORAD(ax,text_on=text_on)
            AxionPhoton.BRASS(ax,text_on=text_on)
            AxionPhoton.ADBC(ax,text_on=text_on)
            AxionPhoton.DANCE(ax,text_on=text_on)
            AxionPhoton.aLIGO(ax,text_on=text_on)
        else:
            AxionPhoton.QUAX(ax,text_on=text_on)
        return

    def LSW(ax,projection=True,text_on=True):
        AxionPhoton.PVLAS(ax,text_on=text_on)
        AxionPhoton.ALPS(ax,projection=projection,text_on=text_on)
        AxionPhoton.OSQAR(ax,text_on=text_on)
        AxionPhoton.CROWS(ax,text_on=text_on)
        return


    def AstroBounds(ax,projection=True,fs=15,text_on=True):
        y2 = ax.get_ylim()[1]
        ### Astrophysical constraints

        # SN-gamma rays arXiv:[1410.3747]
        SNgamma_col = [0.05, 0.5, 0.06]
        SNgamma = loadtxt(limit_dir+"/AxionPhoton/SN-gamma.txt")
        plt.plot(SNgamma[:,0],SNgamma[:,1],'k-',alpha=0.6,zorder=0.21,lw=2)
        plt.fill_between(SNgamma[:,0],SNgamma[:,1],y2=y2,edgecolor=None,facecolor=SNgamma_col,zorder=0.21)

        # M87 Limits from arXiv:[1703.07354]
        M87_col = [0.0, 0.66, 0.42]
        M87 = loadtxt(limit_dir+"/AxionPhoton/M87.txt")
        plt.plot(M87[:,0],M87[:,1],'k-',lw=2,alpha=0.8,zorder=0.2)
        plt.fill_between(M87[:,0],M87[:,1],y2=y2,edgecolor=None,facecolor=M87_col,zorder=0.2)

        # HYDRA-A arXiv:[1304.0989]
        HYDRA_col = [0.24, 0.71, 0.54]
        HYDRA = loadtxt(limit_dir+"/AxionPhoton/HYDRA_A.txt")
        plt.plot(HYDRA[:,0],HYDRA[:,1],'k-',alpha=0.6,zorder=0.23,lw=2)
        plt.fill_between(HYDRA[:,0],HYDRA[:,1],y2=y2,edgecolor=None,facecolor=HYDRA_col,zorder=0.23)

        # HESS arXiv:[1304.0700]
        HESS_col = [0.0, 0.55, 0.3]
        HESS = loadtxt(limit_dir+"/AxionPhoton/HESS.txt")
        plt.plot(HESS[:,0],HESS[:,1],'k-',alpha=0.6,zorder=0.2,lw=2)
        plt.fill_between(HESS[:,0],HESS[:,1],y2=y2,edgecolor=None,facecolor=HESS_col,zorder=0.2)

        # Fermi NGC1275 arXiv:[1603.06978]
        Fermi_col = [0.0, 0.42, 0.24]
        Fermi1 = loadtxt(limit_dir+"/AxionPhoton/Fermi1.txt")
        Fermi2 = loadtxt(limit_dir+"/AxionPhoton/Fermi2.txt")
        plt.fill_between(Fermi1[:,0],Fermi1[:,1],y2=y2,edgecolor=Fermi_col,facecolor=Fermi_col,zorder=0.24,lw=3)
        plt.fill(Fermi2[:,0],1.01*Fermi2[:,1],edgecolor=Fermi_col,facecolor=Fermi_col,lw=3,zorder=0.24)
        Fermi1 = loadtxt(limit_dir+"/AxionPhoton/Fermi_bound.txt")
        Fermi2 = loadtxt(limit_dir+"/AxionPhoton/Fermi_hole.txt")
        plt.plot(Fermi1[:,0],Fermi1[:,1],'k-',alpha=0.5,lw=1.5,zorder=0.24)
        plt.plot(Fermi2[:,0],Fermi2[:,1],'k-',alpha=0.5,lw=1.5,zorder=0.24)

        # Optical telescope [astro-ph/0611502]
        Telescopes_col = [0.09, 0.45, 0.27]
        Telescopes = loadtxt(limit_dir+"/AxionPhoton/Telescopes.txt")
        plt.fill_between(Telescopes[:,0],Telescopes[:,1],y2=y2,edgecolor=None,facecolor=Telescopes_col,zorder=0.2)
        if text_on:
            plt.text(3.3,4e-12,r'{\bf Telescopes}',fontsize=fs,color=Telescopes_col,rotation=-90,ha='left',va='top')
            plt.text(4.8e-10,1.2e-11,r'{\bf Fermi}',fontsize=fs,color='w',ha='left',va='top')
            plt.text(2e-8,1.6e-11,r'{\bf HESS}',fontsize=fs+1,color=HESS_col,ha='left',va='top')
            plt.text(1.5e-12,4e-11,r'{\bf Hydra}',fontsize=fs-2,color='w',ha='left',va='top')
            plt.text(3e-12,2e-11,r'\quad {\bf A}',fontsize=fs-2,color='w',ha='left',va='top')
            plt.text(1.4e-12,4e-12,r'\quad {\bf M87}',fontsize=fs,color='w',ha='left',va='top')
            plt.text(3e-11,2e-11,r'{\bf SN}-$\gamma$',fontsize=fs,color='w',ha='left',va='top')


        # Chandra arXiv:[1907.05475]
        Chandra_col = [0.0, 0.3, 0.24]
        Chandra = loadtxt(limit_dir+'/AxionPhoton/Chandra.txt')
        plt.plot(Chandra[:,0],Chandra[:,1],'k-',alpha=0.8,lw=2,zorder=0.1)
        plt.fill_between(Chandra[:,0],Chandra[:,1],y2=y2,edgecolor=None,facecolor=Chandra_col,zorder=0.1)
        if text_on:
            if projection==False:
                plt.text(1.1e-11,2e-12,r'{\bf Chandra}',fontsize=fs,color=Chandra_col,rotation=0,ha='left',va='top')

        if projection==True:
            # Fermi nearby SN prospects arXiv:[1609.02350]
            FermiSN = loadtxt(limit_dir+"/AxionPhoton/Projections/FermiSN.txt")
            plt.fill_between(FermiSN[:,0],FermiSN[:,1],y2=y2,edgecolor=Fermi_col,linewidth=1.5,facecolor=Fermi_col,zorder=0.1,alpha=0.2)
            if text_on:
                plt.text(1e-9,4e-12,r'{\bf Fermi SN}',fontsize=fs,color=Fermi_col,rotation=43,ha='left',va='top')


        return

    def Cosmology(ax,fs=30,text_on=True):
        y2 = ax.get_ylim()[1]
        ## Cosmology constraints see arXiv:[1210.3196] for summary
        # Xray Background
        XRAY_col = [0.03, 0.57, 0.82]
        XRAY = loadtxt(limit_dir+"/AxionPhoton/XRAY.txt")
        plt.plot(XRAY[:,0],XRAY[:,1],color='k',alpha=0.5,zorder=0.3,lw=2)
        plt.fill_between(XRAY[:,0],XRAY[:,1],y2=1e-11,edgecolor=None,facecolor=XRAY_col,zorder=0.3)

        # Extragalactic background light
        EBL_col =  [0.0, 0.2, 0.6]
        EBL = loadtxt(limit_dir+"/AxionPhoton/EBL.txt")
        EBL2 = loadtxt(limit_dir+"/AxionPhoton/EBL2.txt")
        plt.plot(EBL[:,0],EBL[:,1],'k',lw=2.5,zorder=0.4,alpha=0.8)
        plt.fill_between(EBL[:,0],EBL[:,1],y2=y2,edgecolor=None,facecolor=EBL_col,zorder=0.5)
        plt.fill_between(EBL2[:,0],EBL2[:,1],y2=y2,edgecolor=None,facecolor=EBL_col,zorder=0.5)

        # Ionisation fraction
        x_ion_col = [0.27, 0.51, 0.71]
        x_ion = loadtxt(limit_dir+"/AxionPhoton/x_ion.txt")
        plt.plot(x_ion[:,0],x_ion[:,1],'k',lw=2.5,zorder=0.4,alpha=0.8)
        plt.fill_between(x_ion[:,0],x_ion[:,1],y2=y2,edgecolor=None,facecolor=x_ion_col,zorder=0.5)

        # BBN+N_eff arXiv:[2002.08370]
        BBN_col = [0.27, 0.51, 0.71]
        BBN = loadtxt(limit_dir+"/AxionPhoton/BBN_Neff.txt")
        plt.plot(BBN[:,0],BBN[:,1],'k',lw=2,zorder=0.4,alpha=0.5)
        plt.fill_between(BBN[:,0],BBN[:,1],y2=y2,edgecolor=None,facecolor=BBN_col,zorder=0.4)

        if text_on:
            plt.text(3e3,0.8e-16,r'{\bf X-rays}',fontsize=fs,color='w',rotation=-50,ha='left',va='top')
            plt.text(1e4,5e-14,r'{\bf EBL}',fontsize=fs+5,color='w',rotation=-55,ha='left',va='top')
            plt.text(100.5744,5.1720e-11,r'{\bf Ionisation}',fontsize=fs-6,color='w',rotation=-90,ha='left',va='top')
            plt.text(40,4.1720e-11,r'{\bf fraction}',fontsize=fs-6,color='w',rotation=-90,ha='left',va='top')
            plt.text(0.05e6,8e-10,r'{\bf BBN}+$N_{\rm eff}$',fontsize=fs,color='w',rotation=-55,ha='left',va='top')


    def StellarBounds(ax,fs=30,text_on=True):
        y2 = ax.get_ylim()[1]
        # Stellar physics constraints

        # Globular clusters arXiv:[1406.6053]
        HB_col = [0.0, 0.66, 0.42]
        HB = loadtxt(limit_dir+"/AxionPhoton/HorizontalBranch.txt")
        plt.plot(HB[:,0],HB[:,1],color='k',alpha=0.5,zorder=1,lw=2)
        plt.fill_between(HB[:,0],HB[:,1],y2=y2,edgecolor=None,facecolor=HB_col,zorder=1)

        # Solar neutrino B8 bound arXiv:[1501.01639]
        SolarNu_col = [0.01, 0.75, 0.24]
        SolarNu = loadtxt(limit_dir+"/AxionPhoton/SolarNu.txt")
        plt.plot(SolarNu[:,0],SolarNu[:,1],color='k',lw=2,alpha=0.5,zorder=1)
        plt.fill_between(SolarNu[:,0],SolarNu[:,1],y2=y2,edgecolor=None,facecolor=SolarNu_col,zorder=1)

        # SN1987A-neutrinos updated arXiv:[1808.10136]
        SN = loadtxt(limit_dir+"/AxionPhoton/SN1987A_2019.txt")
        plt.fill_between(SN[:,0],SN[:,1],y2=y2,edgecolor=None,facecolor='ForestGreen',zorder=0.1)
        # SN1987A-decay arXiv:[1702.02964]
        SN = loadtxt(limit_dir+"/AxionPhoton/SN1987A_decay.txt")
        plt.fill_between(SN[:,0],SN[:,1],y2=y2,edgecolor=None,facecolor='ForestGreen',zorder=0.1)

        if text_on:
            plt.text(0.4e6,6e-7,r'{\bf SN1987A}',fontsize=fs-9,color='w',rotation=-60,ha='left',va='top')
            plt.text(1e1,3e-9,r'{\bf Solar} $\nu$',fontsize=fs+3,color='w')
            plt.text(1.4e0,1.5e-10,r'{\bf Horizontal branch}',fontsize=fs-7,color='w')

#==============================================================================#



from matplotlib import patches
from matplotlib import text as mtext
import numpy as np
import math

class CurvedText(mtext.Text):
    """
    A text object that follows an arbitrary curve.
    """
    def __init__(self, x, y, text, axes, **kwargs):
        super(CurvedText, self).__init__(x[0],y[0],' ', **kwargs)

        axes.add_artist(self)

        ##saving the curve:
        self.__x = x
        self.__y = y
        self.__zorder = self.get_zorder()

        ##creating the text objects
        self.__Characters = []
        for c in text:
            if c == ' ':
                ##make this an invisible 'a':
                t = mtext.Text(0,0,'a')
                t.set_alpha(0.0)
            else:
                t = mtext.Text(0,0,c, **kwargs)

            #resetting unnecessary arguments
            t.set_ha('center')
            t.set_rotation(0)
            t.set_zorder(self.__zorder +1)

            self.__Characters.append((c,t))
            axes.add_artist(t)


    ##overloading some member functions, to assure correct functionality
    ##on update
    def set_zorder(self, zorder):
        super(CurvedText, self).set_zorder(zorder)
        self.__zorder = self.get_zorder()
        for c,t in self.__Characters:
            t.set_zorder(self.__zorder+1)

    def draw(self, renderer, *args, **kwargs):
        """
        Overload of the Text.draw() function. Do not do
        do any drawing, but update the positions and rotation
        angles of self.__Characters.
        """
        self.update_positions(renderer)

    def update_positions(self,renderer):
        """
        Update positions and rotations of the individual text elements.
        """

        #preparations

        ##determining the aspect ratio:
        ##from https://stackoverflow.com/a/42014041/2454357

        ##data limits
        xlim = self.axes.get_xlim()
        ylim = self.axes.get_ylim()
        ## Axis size on figure
        figW, figH = self.axes.get_figure().get_size_inches()
        ## Ratio of display units
        _, _, w, h = self.axes.get_position().bounds
        ##final aspect ratio
        aspect = ((figW * w)/(figH * h))*(ylim[1]-ylim[0])/(xlim[1]-xlim[0])

        #points of the curve in figure coordinates:
        x_fig,y_fig = (
            np.array(l) for l in zip(*self.axes.transData.transform([
            (i,j) for i,j in zip(self.__x,self.__y)
            ]))
        )

        #point distances in figure coordinates
        x_fig_dist = (x_fig[1:]-x_fig[:-1])
        y_fig_dist = (y_fig[1:]-y_fig[:-1])
        r_fig_dist = np.sqrt(x_fig_dist**2+y_fig_dist**2)

        #arc length in figure coordinates
        l_fig = np.insert(np.cumsum(r_fig_dist),0,0)

        #angles in figure coordinates
        rads = np.arctan2((y_fig[1:] - y_fig[:-1]),(x_fig[1:] - x_fig[:-1]))
        degs = np.rad2deg(rads)


        rel_pos = 10
        for c,t in self.__Characters:
            #finding the width of c:
            t.set_rotation(0)
            t.set_va('center')
            bbox1  = t.get_window_extent(renderer=renderer)
            w = bbox1.width
            h = bbox1.height

            #ignore all letters that don't fit:
            if rel_pos+w/2 > l_fig[-1]:
                t.set_alpha(0.0)
                rel_pos += w
                continue

            elif c != ' ':
                t.set_alpha(1.0)

            #finding the two data points between which the horizontal
            #center point of the character will be situated
            #left and right indices:
            il = np.where(rel_pos+w/2 >= l_fig)[0][-1]
            ir = np.where(rel_pos+w/2 <= l_fig)[0][0]

            #if we exactly hit a data point:
            if ir == il:
                ir += 1

            #how much of the letter width was needed to find il:
            used = l_fig[il]-rel_pos
            rel_pos = l_fig[il]

            #relative distance between il and ir where the center
            #of the character will be
            fraction = (w/2-used)/r_fig_dist[il]

            ##setting the character position in data coordinates:
            ##interpolate between the two points:
            x = self.__x[il]+fraction*(self.__x[ir]-self.__x[il])
            y = self.__y[il]+fraction*(self.__y[ir]-self.__y[il])

            #getting the offset when setting correct vertical alignment
            #in data coordinates
            t.set_va(self.get_va())
            bbox2  = t.get_window_extent(renderer=renderer)

            bbox1d = self.axes.transData.inverted().transform(bbox1)
            bbox2d = self.axes.transData.inverted().transform(bbox2)
            dr = np.array(bbox2d[0]-bbox1d[0])

            #the rotation/stretch matrix
            rad = rads[il]
            rot_mat = np.array([
                [math.cos(rad), math.sin(rad)*aspect],
                [-math.sin(rad)/aspect, math.cos(rad)]
            ])

            ##computing the offset vector of the rotated character
            drp = np.dot(dr,rot_mat)

            #setting final position and rotation:
            t.set_position(np.array([x,y])+drp)
            t.set_rotation(degs[il])

            t.set_va('center')
            t.set_ha('center')

            #updating rel_pos to right edge of character
            rel_pos += w-used
