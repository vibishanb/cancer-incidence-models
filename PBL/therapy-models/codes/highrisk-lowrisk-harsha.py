import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.stats import truncnorm

def lotvol(t,x,r,A_min,A_max,mu,sig):
    #Defining a random Environment variable bounded by 0 and 1
    a, b = (0 - mu) / sig, (1 - mu) / sig
    #Truncated normal distribution sampled at every time step
    E=truncnorm.rvs(a,b,loc=mu,scale=sig)
    #Effective alpha at that time
    A=A_max-E*(A_max-A_min)
    # Returns the array with dy_i/dt at that time
    return np.array(r*x*(1-A@x))

y0=np.array([0.5,0.5]) #initial population level
r=np.array([2,1]) #growth rate
A_min=np.array([[1,1],[1,1]]) #alpha values when E=1 (best environment)
A_max=np.array([[1.5,1.1],[1.5,1.1]]) #alpha values when E=0 (worst environment)
mu,sig = (0.5,0.5)
t_max,dt=(100,0.1)
#Timeseries arrays
t=np.arange(0,t_max,dt)
sol = solve_ivp(lotvol, [0, t_max], y0, args=(r,A_min,A_max,mu,sig),t_eval=t,dense_output=True)
plt.plot(t,sol.y.T)
plt.xlabel("Time (days)")
plt.ylabel("Density")
