#Mitch Delemeester
#5/21/2024
#Capillary Rheometry Theoretical Data Generation


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm


#TARGETS WE WILL EVENTUALLY TRY TO PREDICT

yield_stress = np.array([5E4]) #yield stress (pascals)
slip_yield_stress = np.array([5E3]) #slip yield stress (pascals)
Beta = np.array([5E-10]) #slip velocity coefficient
thinning_index = np.array([0.67]) #shear thinning index
Viscosity_consistency = np.array([1100]) #shear thinning consistency (pa.s^n)


#initial conditions:

Diameter = np.linspace(0.5, 1.5, 50)*10**-3 #diameter (meters)
Radius = Diameter/2 #radius (meters)

aspectRatio = np.linspace(10,100, 50) #aspect ratio L/D
Length = Diameter*aspectRatio #length (meters)

G = np.linspace(3.33E8, 6.67E9, 50) #pressure gradient array (pascals/meter)


iterations = len(yield_stress)*len(slip_yield_stress)*len(Beta)*len(thinning_index)*len(Viscosity_consistency)*len(Radius)*len(G)
print("Total iterations: ", iterations)




#empty array to store the data
stored_theory_values = np.array([])
stored_target_values = np.array([])
stored_target_values_2 = np.array([])


#calculate the apparent wall shear rate (gamma_dot_apparent) for each case
for g_i in tqdm(G):
     for ii in range(len(Radius)):

        #pressure drop for each Radius/Aspect Ratio pair and the pressure gradient (Pa)
        dP = g_i*Length[ii]

        #wall shear stress (Pa)
        tw = g_i*Radius[ii]/2


        for ty in yield_stress:
            #Radius of plug (meters)
            ri = Radius[ii]*ty/tw

            for tys in slip_yield_stress:
                #thickness of slip layer (meters)
                delta = Radius[ii]*tys/tw

                for beta in Beta:
                    #slip velocity (m/s)
                    if tw > tys:
                        vs0 = beta*(tw-tys)
                    else:
                        vs0 = 0

                    for n in thinning_index:
                        #shear thinning index reciprocal
                        s = 1/n 

                        for k in Viscosity_consistency:
                            
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

                            temp_theorydata = np.array([Radius[ii], L[ii], aspectRatio[ii], dP, g_i, gamma_dot_apparent, Qtot])

                            temp_target_data = np.array([Qslip, Qplug])

                            temp_target_data_2 = np.array([ty, tys, beta, n, k])

                            #store the theoretical data from each iteration in a matrix (turn into a pandas dataframe later)
                            stored_theory_values = np.concatenate((stored_theory_values, temp_theorydata), axis = 0)

                            #store the target values from each iteration in a matrix (turn into a pandas dataframe later)
                            stored_target_values = np.concatenate((stored_target_values, temp_target_data), axis = 0)

                            stored_target_values_2 = np.concatenate((stored_target_values_2, temp_target_data_2), axis = 0)




#reshape the stored data into a matrix
stored_theory_values = np.reshape(stored_theory_values, (iterations, len(temp_theorydata)))

#convert the stored data into a pandas dataframe
Theory_Data = pd.DataFrame(stored_theory_values, columns = ['Radius (m)', 'Length (m)','Aspect Ratio', 'Pressure Drop (Pa)', 'Pressure Gradient(Pa m^-1)', 'Apparent Wall Shear Rate (s^-1)', 'Volumetric Flow Rate (m^3 s^-1)'])

#save the data to a csv file
Theory_Data.to_csv('Theory_Data.csv', index = False)




#reshape the target data into a matrix
stored_target_values = np.reshape(stored_target_values, (iterations, len(temp_target_data)))

#convert the stored data into a pandas dataframe
Target_Data = pd.DataFrame(stored_target_values, columns = ['Slip Flow Rate (m^3 s^-1)', 'Plug Flow Rate (m^3 s^-1)'])

#save the data to a csv file
Target_Data.to_csv('Target_Data.csv', index = False)




#reshape the target data into a matrix
stored_target_values_2 = np.reshape(stored_target_values_2, (iterations, len(temp_target_data_2)))

#convert the stored data into a pandas dataframe
Target_Data_2 = pd.DataFrame(stored_target_values_2, columns = ['Yield Stress (Pa)', 'Slip Yield Stress (Pa)', 'Beta', 'Shear Thinning Index', 'Viscosity Consistency (Pa s^n)'])

#save the data to a csv file
Target_Data_2.to_csv('Target_Data_2.csv', index = False)