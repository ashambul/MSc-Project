import numpy as np
import matplotlib.pyplot as plt

# BOUNDARY CONDITIONS

m1_total = 3*np.ones(200) #mass flow rate (kg/s) 
T1i = 400 # inlet temperature (K)


m2_total = 5*np.ones(200) #mass flow rate (kg/s)
T2i = 600 # pipe inner surface temperature (K)

T0 = 300 # initial temperature of fluid throughout the pipe (K)
# PHYSICAL PROPERTIES
# tube fluid properties
T1 = np.ones(n)*T0
T2 = np.ones(n)*T0

V_shell = np.ones(n)
V_tube = np.ones(n)

h_tube = np.ones(n)
h_shell = np.ones(n)
Cp_tube = np.ones(n)
Cp_shell =np.ones(n)
rho_tube = np.ones(n)
rho_shell = np.ones(n)
Re_tube = np.ones(n)
Re_shell = np.ones(n)
Pr_tube= np.ones(n)
Pr_shell = np.ones(n)
Nu_tube = np.ones(n)
Nu_shell = np.ones(n)

k_tube = (4.2365*10**(-9))*T1**3-(1.1440*10**-5)*T1**2+7.1959*10**(-3)*T1-0.63262  #(Zografos 1987)
k_shell = (4.2365*10**(-9))*T2**3-(1.1440*10**-5)*T2**2+7.1959*10**(-3)*T2-0.63262 # (Zografos 1987)


k_wall = 36 # W/mK thermal conductivity of the wall (Carbon steel @400K)


k_foul = 0.329 #W/mK #fouling layer thermal conductivity (Diaz Bejarano thesis p.41)



# GEOMETRY

L = 6.1 #Pipe length (m)
r_tube_i = 0.009525 # Inner radius of the pipe (m)
r_tube_e = 0.127 #Outer pipe radius (m)
r_shell =  0.0254 #Radius of the shell (m)
r_flow = 0.008525 # Flow radius with fouling inside the tube
n = 100 # nodes used
N_tubes = 200 #number of tubes


#CALCULATIONS

pi = 3.14159
Ac1 = pi*r_tube_i**2 # cross sectional area of inner tube (m^2)
Ac2 = pi*r_shell**2-N_tubes*r_tube_e**2# cross-sectional area of shell without tubes  (m^2)



#Aeff = 2*r_tube_i*pi*L
#t_wall = r_tube_e-r_tube_i 
#t_foul = r_tube_i-r_flow
#A_wall = 2*pi*r_shell+N_tubes*(pi*r_tube_e**2-r_tube_i**2)




m1 = m1_total[1]


m2 = m2_total[1]

#Dittus Boelter

V_shell = m1/(rho_shell*Ac2) #Velocity on shell side m/s based on Ac2
Pr_shell = Cp_shell*mu_shell/k_shell
Re_shell = rho_shell*V_shell*2*r_shell/mu_shell
Nu_shell = 0.023*Re_shell**0.8*Pr_shell**0.4
h_shell = Nu_shell*k_shell/L

V_tube = m2/(rho_tube*Ac1)#Velocity on tube side m/s based on the flow area
Pr_tube = Cp_tube*mu_tube/k_tube
Re_tube = rho_tube*V_tube*2*r_flow/mu_tube
Nu_tube = 0.023*Re_tube**0.8*Pr_tube**0.4
h_tube = Nu_tube*k_tube/L





dx = L/n    #node width (m)

t_final = 1000 # simulation time (s)
dt = 1 # time step (s)


x = np.linspace (dx/2, L-dx/2, n)





dT1dt = np.zeros(n)
dT2dt = np.zeros(n)


t = np.arange(0, t_final, dt)


UA = 1/((1/h_shell*pi*r_tube_e*dx)+((np.log(r_tube_e/r_tube_i))/2*pi*k_wall*dx)+((np.log(r_tube_i/r_flow))/2*pi*k_foul*dx)+(1/(h_tube*2*pi*r_flow)) # overall heat transfer coefficient (W/m^2*K)



for j in range(1,len(t)):
    
    plt.clf()
    

    #tube fluid properties
    
    Cp_tube = (T1**3*1.785*10**-7)-(1.9149*10**-4*T1**2)+(6.7953*10**-2*T1-3.7559) #heat capacity of fluid water (J/kg*K)
    rho_tube = -3.0115*10**(-6)*(T1-252.33)**(-1) #density of fluid water (kg/m^3)
    mu_tube = 3.8208*10**(-2)*(T1-252.33)**-1 #Viscoity offluid in tube side (kg/ms)

    #shell fluid properties

    Cp_shell = T2**(3)*1.785*10**(-7)-1.9149*10**(-4)*(T2**2)+6.7953*10**(-2)*T2-3.7559 #heat capacity of fluid water (J/kg*K)
    rho_shell = -3.0115*10**(-6)*(T2-252.33)**(-1) #density of fluid water (kg/m^3)
    mu_shell = 3.8208*10**(-2)*(T2-252.33)**-1 #Viscosity of fluid in shell side (kg/ms)
    
    dT1dt[1:n] = (m1*Cp_tube*(T1[0:n-1]-T1[1:n])+UA*dx*(T2[1:n]-T1[1:n]))/(rho_tube*Cp_tube*dx*A_tube)
    dT1dt[0] = (m1*Cp_tube*(T1i-T1[0])+U*2*pi*r_tube*dx*(T2[0]-T1[0]))/(rho_tube*Cp_tube*dx*A_tube)
    
    dT2dt[1:n] = (m2*Cp_shell*(T2[0:n-1]-T2[1:n])-UA*dx*(T2[1:n]-T1[1:n]))/(rho_shell*Cp_shell*dx*Ac2)
    dT2dt[0] = (m2*Cp_shell*(T2i-T2[0])-U*2*pi*r1*dx*(T2[0]-T1[0]))/(rho_shell*Cp_shell*dx*Ac2)
    
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
    
    
    
