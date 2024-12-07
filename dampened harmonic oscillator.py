import numpy as np 
from numpy import sin, cos
from matplotlib import pyplot as plt 
import math 
import pandas as pd

g = 9.81
l = g/(np.pi)**2
b = 0.1
m = 1
delta = 0.02

gamma = b/m

w0 = np.sqrt(g/l)
w = (w0)**2 - (gamma)**2/4


def theta_double_dot(theta):
    theta2prime = -(g/l) * sin(theta)
    return theta2prime

delta = 0.02

def runge_kutta_2_damped(delta,theta,theta_dot):
    
    k1_theta = delta*theta_dot
    
    k1_theta_dot = delta * ((-b/m)*theta_dot  - (g/l) * sin(theta ))
    
    k2_theta = delta * (theta_dot + k1_theta_dot/2)
    
    k2_theta_dot = delta * (((-b/m) * theta_dot + k1_theta_dot/2) - (g/l)*sin(theta + k1_theta/2))
    
    theta_new = theta + k2_theta
    
    theta_dot_new = theta_dot + k2_theta_dot
    
    return theta_new, theta_dot_new

## initial theta value
theta_0 = 0.1

## initial theta_dot vvalue
theta_dot_0 = 0

## length of time
t_n = 20

## num time steps
N = int(t_n/delta)

times = np.linspace(0, t_n + delta, N+1)

thetas_damped = [theta_0]

theta_dots_damped = [theta_dot_0]

for i in range(0,N):
    
    t_curr, td_curr = runge_kutta_2_damped(delta,thetas_damped[-1],theta_dots_damped[-1])
    
    thetas_damped.append(t_curr)
    theta_dots_damped.append(td_curr)
    
def theta(t,theta0, theta_dot0, omega):
    return np.exp(-gamma*t/2)*cos(w*t)
    
plt.figure(figsize = (10,6))

plt.plot(times, thetas_damped, label = "theta_damped")

plt.plot(times, theta(times, theta_0, theta_dot_0, np.sqrt(g/l)), label = 'damped analytical')

plt.xlabel("time", fontsize = 20)
plt.ylabel("theta", fontsize = 20)

plt.legend(fontsize = 15)

plt.show()

data = {"theta damped": thetas_damped , "time":times , "thetadots damped": theta_dots_damped}

df= pd.DataFrame(data)

print(df)

