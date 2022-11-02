import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
import bell_checker.tools as tools

def func(x, A, B, C, D):
    """
    This is a sinusoidal model that is used to obtain a fit .
            
        Args:
            x:              A generic variable for this model.
            A, B, C, D:     The parameters for this model.
            
        Output:
            A function func with variable x and parameters A, B, C and D.
    """
    return A * np.sin(B*x + D) + C


def g1D(S=None, t=None, v=None, title="1D plot"):
    """
    Creates a one dimensional plot of any pair of lists, adding the classical limit S.
            
        Args:
            S:     The classical limit in Bell's inequality.
            t:     List containing angle values.
            v:     List containing values to plot
            title: The plot title
            
        Output:
            Plot t vs v.
    """
    if v==None:
        print("There is not any list to plot")
        return
    
    r = [0,np.pi/2,np.pi,3*np.pi/2,2*np.pi]
    
    # Interpolation part
    plt.figure(figsize=(14, 8))
    plt.plot(t,v,"b.",ms=20)
    interp_func = CubicSpline(t,v)
    
    # Fitting part
    plt.plot(np.linspace(0,2*np.pi,100),interp_func(np.linspace(0,2*np.pi,100)),"b",label="Interpolation")
    params, values = curve_fit(func, t, v)
    plt.plot(np.linspace(0,2*np.pi,100),func(np.linspace(0,2*np.pi,100),params[0],params[1],params[2],params[3]),"k--",label="Fit" )
    

    if S!=None:
        plt.axhline(y=S, color="r")
        plt.axhline(y=-S, color="r")

    plt.xticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],fontsize=20);
    plt.yticks(fontsize=20);
    plt.legend(fontsize=20)
    plt.title(title,fontsize=20)
    plt.xlabel(r'$\theta$',fontsize=20)
    plt.grid("on")
    
def g2D(S,x=None,y=None,v=None):
    
    """
    Creates a contour plot and a density plot to observe how Bell's inequality is broken.
            
        Args:
            S:     The classical limit in Bell's inequality.
            x:     List containing horizontal values.
            y:     List containing vertical values.
            v:     Matrix containing the grid values to plot (i.e v = [[f(i,j) for i in x] for j in y]).
            
        Output:
            Contour and density plot in different figures.
    """
    
    if v==None:
        print("There is not any matrix to plot")
        return
    
    data = gaussian_filter(v, 0.5);
    
    levels = [-S,-S/2,-S/8,S/8,S/2,S]
    
    plt.figure(figsize=(24,8))
    
    plt.subplot(1, 2, 1)
    
    r = [0,np.pi/2,np.pi,3*np.pi/2,2*np.pi,3*np.pi,4*np.pi]
    
    plt.xticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$', r'$3\pi$', r'$4\pi$'],fontsize=25);
    plt.yticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$', r'$3\pi$', r'$4\pi$'],fontsize=25);
    CS = plt.contour(np.linspace(x[0],x[-1],np.array(data).shape[1]),np.linspace(y[0],y[-1],np.array(data).shape[0]),data, levels=levels, cmap="seismic");
    plt.clabel(CS, inline=3, fontsize=12)
    plt.title('Contour plot',fontsize=20);
    plt.xlabel('x',fontsize=20)
    plt.ylabel('y',fontsize=20)

    plt.subplot(1, 2, 2)
    plt.imshow(v, cmap="seismic",interpolation='bicubic',extent=[x[0],x[-1], y[0], y[-1]],origin='lower',alpha=0.8, aspect='auto');
    plt.colorbar().ax.tick_params(labelsize=20);
    plt.xticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$', r'$3\pi$', r'$4\pi$'],fontsize=25);
    plt.yticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$', r'$3\pi$', r'$4\pi$'],fontsize=25);
    plt.xlim(x[0],x[-1])
    plt.ylim(y[0],y[-1])
    plt.title('Density plot',fontsize=20);
    plt.xlabel('x',fontsize=20)
    plt.ylabel('y',fontsize=20)
    
    plt.show()
    
def g2D_2(S,x=None,y=None,v=None,title='Contour and density plot'):
    
    """
    Creates a contour plot and a density plot in the same figure to observe how Bell's inequality is broken.
            
        Args:
            S:     The classical limit in Bell's inequality.
            x:     List containing horizontal values.
            y:     List containing vertical values.
            v:     Matrix containing the grid values to plot (i.e v = [[f(i,j) for i in x] for j in y]).
            
        Output:
            Contour and density plot in the same figure.
    """
    if v==None:
        print("There is not any matrix to plot")
        return
    
    data = gaussian_filter(v, 0.5);
    levels = [-S,-S/2,-S/8,S/8,S/2,S]
    r = [0,np.pi/2,np.pi,3*np.pi/2,2*np.pi,3*np.pi,4*np.pi]
    
    plt.figure(figsize=(10,8))
    plt.imshow(v, cmap='RdGy',interpolation='bicubic',extent=[x[0],x[-1], y[0], y[-1]],origin='lower',alpha=0.8, aspect='auto');
    plt.colorbar().ax.tick_params(labelsize=20);
    plt.xticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$', r'$3\pi$', r'$4\pi$'],fontsize=25);
    plt.yticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$', r'$3\pi$', r'$4\pi$'],fontsize=25);
    CS=plt.contour(np.linspace(x[0],x[-1],np.array(data).shape[1]),np.linspace(y[0],y[-1],np.array(data).shape[0]),data, levels=levels, colors='k', alpha=0.8);
    plt.clabel(CS, inline=1, fontsize=12)
    plt.xlim(x[0],x[-1])
    plt.ylim(y[0],y[-1])
    plt.title(title,fontsize=20);
    plt.xlabel('x',fontsize=20)
    plt.ylabel('y',fontsize=20)