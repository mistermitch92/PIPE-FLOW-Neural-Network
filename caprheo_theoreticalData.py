#Mitch Delemeester
#5/21/2024
#Capillary Rheometry Theoretical Data Generation


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm


#GIVEN VALUES WE WILL EVENTUALLY TRY TO PREDICT

yield_stress = np.linspace(0, 2.5E6, num=5) #yield stress (pascals)
slip_yield_stress = np.linspace(0, 5.0E5, num=5) #slip yield stress (pascals)
Beta = np.linspace(5E-11, 5E-9, num=5) #slip velocity coefficient
thinning_index = np.linspace(0.6, 0.7, num=5) #shear thinning index
Viscosity_consistency = np.linspace(800, 1300, num=5) #shear thinning consistency (pa.s^n)


#initial conditions:

Diameter = np.array([0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.5, 1.5, 1.5])*10**-3 #diameter (meters)
Radius = Diameter/2 #radius (meters)

aspectRatio = np.array([ 10, 20, 30, 10, 20, 30, 10, 15, 20]) #aspect ratio L/D
Length = Diameter*aspectRatio #length (meters)

G = np.linspace(1.5E10, 1.14E11, 5) #pressure gradient array (pascals/meter)


iterations = len(yield_stress)*len(slip_yield_stress)*len(Beta)*len(thinning_index)*len(Viscosity_consistency)*len(Radius)*len(G)
print("Total iterations: ", iterations)




#empty array to store the data
stored_theory_values = np.array([])
stored_given_values = np.array([])


#calculate the apparent wall shear rate (gamma_dot_apparent) for each case
for ty in tqdm(yield_stress):
    for tys in slip_yield_stress:
        for beta in Beta:
            for n in thinning_index:
                s = 1/n #shear thinning index reciprocal
                for k in Viscosity_consistency:
                    for ii in range(len(Radius)):
                        for g_i in G:
                            #pressure drop for each Radius/Aspect Ratio pair and the given pressure gradient (Pa)
                            dP = g_i*Length[ii]
                        
                            #wall shear stress (Pa)
                            tw = g_i*Radius[ii]/2

                            #Radius of plug (meters)
                            ri = Radius[ii]*ty/tw

                            #thickness of slip layer (meters)
                            delta = Radius[ii]*tys/tw

                            #slip velocity (m/s)
                            if tw > tys:
                                vs0 = beta*tw
                            else:
                                vs0 = 0

                            #plug velocity (m/s)
                            if tw > ty:
                                v0 = (g_i*(Radius[ii]-ri)**(s+1))/(2*k*(s+1))
                            else:
                                v0 = 0

                            #volumetric flow rate from slip (m^3/s)
                            Qslip = np.pi*Radius[ii]**2*vs0*(1-delta/Radius[ii])

                            #volumetric flow rate from plug (m^3/s)
                            Qplug = np.pi*v0*(Radius[ii]**2-2*(Radius[ii]-ri)**2/(s+3))

                            #total volumetric flow rate (m^3/s)
                            Qtot = Qslip + Qplug

                            #apparent wall shear rate (1/s)
                            gamma_dot_apparent = 4*Qtot/(np.pi*Radius[ii]**3)

                            temp_theorydata = np.array([Radius[ii], aspectRatio[ii], dP, g_i, gamma_dot_apparent, Qtot])

                            temp_givendata = np.array([ty, tys, beta, n, k])

                            #store the theoretical data from each iteration in a matrix (turn into a pandas dataframe later)
                            stored_theory_values = np.concatenate((stored_theory_values, temp_theorydata), axis = 0)

                            #store the given data from each iteration in a matrix (turn into a pandas dataframe later)
                            stored_given_values = np.concatenate((stored_given_values, temp_givendata), axis = 0)




#reshape the stored data into a matrix
stored_theory_values = np.reshape(stored_theory_values, (iterations, 6))

#convert the stored data into a pandas dataframe
Theory_Data = pd.DataFrame(stored_theory_values, columns = ['Radius (m)', 'Aspect Ratio', 'Pressure Drop (Pa)', 'Pressure Gradient(Pa m^-1)', 'Apparent Wall Shear Rate (s^-1)', 'Volumetric Flow Rate (m^3 s^-1)'])

#save the data to a csv file
Theory_Data.to_csv('Theory_Data.csv', index = False)

#reshape the given data into a matrix
stored_given_values = np.reshape(stored_given_values, (iterations, 5))

#convert the stored data into a pandas dataframe
Given_Data = pd.DataFrame(stored_given_values, columns = ['Yield Stress (Pa)', 'Slip Yield Stress (Pa)', 'Slip Velocity Coefficient', 'Shear Thinning Index', 'Shear Thinning Consistency (Pa s^n)'])

#save the data to a csv file
Given_Data.to_csv('Given_Data.csv', index = False)