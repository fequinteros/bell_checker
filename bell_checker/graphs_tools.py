import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.interpolate import Rbf,Akima1DInterpolator,CubicSpline
from tqdm import tqdm
from scipy.optimize import curve_fit


def func(x, A, B, C, D):
    """Modelo para nuestros datos."""
    return A * np.sin(B*x + D) + C

def g1D(obj=None, Sabxy=None, S=None, t_min=0, t_max=2*np.pi, t=None, v=None):
    """
    The general form of each circuit, evaluate different cases and make the correspond circuit.
            
        Args:
            init_state:  The state to evaluate in QuantumCircuit, array or list form.
            obs:         Observables to measure for each qubit.
            obs_conf:    
            theta_matr:  The array or list of angles that rotate each observable, the combination of angles are given in sub-arrays (i.e [[pi]] --> rotates Alice obs stay fix and Bob obs is rotated in pi ).
            index:       The index of the observables to use in the QuantumCircuit
            
        Output:
            Each QuantumCircuit needed to evaluate the bell inequality.
    """
    if v==None:
        t = np.linspace(t_min,t_max,7)
        v = [obj.witness(init_state=init_state,theta_vec=[i], Sabxy=Sabxy, S=S) for i in t]
    
    r = [0,np.pi/2,np.pi,3*np.pi/2,2*np.pi]
    
    plt.figure(figsize=(14, 8))
    plt.plot(t,v,"b.",ms=20)
    interp_func = CubicSpline(t,v)
    plt.plot(np.linspace(0,2*np.pi,100),interp_func(np.linspace(0,2*np.pi,100)),"b",label="Interpolation")
    params, values = curve_fit(func, t, v)
    plt.plot(np.linspace(0,2*np.pi,100),func(np.linspace(0,2*np.pi,100),params[0],params[1],params[2],params[3]),"k--",label="Fit" )
    

    if S!=None:
        plt.plot(t,[S for i in t],"r")
        plt.plot(t,[-S for i in t],"r")

    plt.xticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],fontsize=20);
    plt.yticks(fontsize=20);
    plt.legend(fontsize=20)
    plt.title("1D plot",fontsize=20)
    plt.grid("on")
    
def g2D(S,Sabxy=None,obj=None,x_min=0,x_max=2*np.pi,y_min=-np.pi/4,y_max=2*np.pi,x=None,y=None,v=None):
    if v==None:
        x = np.linspace(x_min,x_max,20)
        y = np.linspace(y_min,y_max,20)
        v = [[obj.witness(init_state,[i,j], Sabxy, S) for j in x] for i in tqdm(y,ncols=80)]
    
    data = gaussian_filter(v, 0.5);
    
    levels = [-S,-S/2,-S/8,S/8,S/2,S]
    
    plt.figure(figsize=(24,8))
    
    plt.subplot(1, 2, 1)
    
    r = [0,np.pi/2,np.pi,3*np.pi/2,2*np.pi]
    
    plt.xticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],fontsize=25);
    plt.yticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],fontsize=25);
    CS = plt.contour(np.linspace(x[0],x[-1],np.array(data).shape[1]),np.linspace(y[0],y[-1],np.array(data).shape[0]),data, levels=levels, cmap="seismic");
    plt.clabel(CS, inline=3, fontsize=12)
    plt.title('Gráfica contour',fontsize=20);

    plt.subplot(1, 2, 2)
    plt.imshow(v, cmap="seismic",interpolation='bicubic',extent=[x[0],x[-1], y[0], y[-1]],origin='lower',alpha=0.8);
    plt.colorbar().ax.tick_params(labelsize=20);
    plt.xticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],fontsize=25);
    plt.yticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],fontsize=25);
    plt.xlim(x[0],x[-1])
    plt.ylim(y[0],y[-1])
    plt.title('Gráfica contour',fontsize=20);
    
    plt.show()
    
def g2D_2(S,Sabxy=None,obj=None,x_min=0,x_max=2*np.pi,y_min=-np.pi/4,y_max=2*np.pi,x=None,y=None,v=None):
    if v==None:
        if obj==None:
            print("The object argument is required")
            return 0
        else:
            x = np.linspace(x_min,x_max,20)
            y = np.linspace(y_min,y_max,20)
            v = [[obj.witness(init_state,[i,j], Sabxy, S) for j in x] for i in tqdm(y,ncols=80)]
    
    data = gaussian_filter(v, 0.5);
    levels = [-S,-S/2,-S/8,S/8,S/2,S]
    r = [0,np.pi/2,np.pi,3*np.pi/2,2*np.pi]
    
    plt.figure(figsize=(10,8))
    plt.imshow(v, cmap='RdGy',interpolation='bicubic',extent=[x[0],x[-1], y[0], y[-1]],origin='lower',alpha=0.8);
    plt.colorbar().ax.tick_params(labelsize=20);
    plt.xticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],fontsize=25);
    plt.yticks(r,[r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],fontsize=25);
    CS=plt.contour(np.linspace(x[0],x[-1],np.array(data).shape[1]),np.linspace(y[0],y[-1],np.array(data).shape[0]),data, levels=levels, colors='k', alpha=0.8);#, linestyles='-.'
    plt.clabel(CS, inline=1, fontsize=12)
    plt.xlim(x[0],x[-1])
    plt.ylim(y[0],y[-1])
    plt.title('Gráfica contour',fontsize=20);