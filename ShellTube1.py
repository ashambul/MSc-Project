import numpy as np
import matplotlib.pyplot as plt


L = 60 #Pipe length (m)
r1 = 0.1 #Pipe radius (m)
r2 = 0.15 #Outer pipe radius (m)
n = 100 # nodes used

m1 = 3 #mass flow rate (kg/s)
Cp1 = 4180 #heat capacity of fluid water (J/kg*K)
rho1 = 1000 #density of fluid water (kg/m^3)

m2 = 5 #mass flow rate (kg/s)
Cp2 = 4180 #heat capacity of fluid water (J/kg*K)
rho2 = 1000 #density of fluid water (kg/m^3)

pi = 3.14159
Ac1 = pi*r1**2 # cross sectional area of inner pipe (m^2)
Ac2 = pi*(r2**2-r1**2) # cross-sectional area of outer pipe (m^2)

T1i = 400 # inlet temperature (K)
T2i = 800 # pipe inner surface temperature (K)
T0 = 300 # initial temperature of fluid throughout the pipe (K)

U = 340 # overall heat transfer coefficient (W/m^2*K)

dx = L/n # node width (m)

t_final = 1000 # simulation time (s)
dt = 1 # time step (s)


x = np.linspace (dx/2, L-dx/2, n)

T1 = np.ones(n)*T0
T2 = np.ones(n)*T0

dT1dt = np.zeros(n)
dT2dt = np.zeros(n)

t = np.arange(0, t_final, dt)

for j in range(1,len(t)):
    plt.clf()
    
    dT1dt[1:n] = (m1*Cp1*(T1[0:n-1]-T1[1:n])+U*2*pi*r1*dx*(T2[1:n]-T1[1:n]))/(rho1*Cp1*dx*Ac1)
    dT1dt[0] = (m1*Cp1*(T1i-T1[0])+U*2*pi*r1*dx*(T2[0]-T1[0]))/(rho1*Cp1*dx*Ac1)
    
    dT2dt[1:n] = (m2*Cp2*(T2[0:n-1]-T2[1:n])-U*2*pi*r1*dx*(T2[1:n]-T1[1:n]))/(rho2*Cp2*dx*Ac2)
    dT2dt[0] = (m2*Cp2*(T2i-T2[0])-U*2*pi*r1*dx*(T2[0]-T1[0]))/(rho2*Cp2*dx*Ac2)
    
    T1 = T1 + dT1dt*dt
    T2 = T2 + dT2dt*dt
    
    plt.figure(1)
    plt.plot(x,T1, color = 'green', label = 'Inside')
    plt.plot(x,T2, color = 'red', label = 'Outside')
    plt.axis([0, L, 298, 820])
    plt.xlabel('Length (m)')
    plt.ylabel('Temperature(K)')
    plt.legend(loc = 'upper right')
    plt.show()
    plt.pause(0.005)
    
    
    
