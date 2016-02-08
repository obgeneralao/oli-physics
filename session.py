# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 17:44:31 2016

@author: oli
"""



from nbarrier import n_barrier
import scipy as sc



barrier = 2		       	#number of potential barriers
barrier_width = 10e-9		#barrier width (m), 10 Angstrom = 1e-9m = 1nm
well_width    = 4e-9		#well width (m)
V0 = 0.240			#barrier energy (eV),height of barrier
                                  	# 1 eV = 1.60217733e-19 J
npoints = 1000			#number of points in energy plot



initial_value = n_barrier(barrier,barrier_width,well_width,V0,npoints)
#initial_value.potential_barrier_count()

Trans, E = initial_value.transmission()


Emax = V0 + 0.1

import matplotlib.pyplot as plt  



Emax    =  V0 + 0.1	
#fig = plt.figure(figsize=(10, 6))
#ax = fig.add_subplot(1,1,1)   
plt.axis([0,Emax,-20,1])
plt.plot(E,sc.log(Trans),'b.')

### major ticks every 0.2, minor ticks every 0.05                                     
##x_major_ticks = numpy.arange(0, Emax, 0.1)                                              
##x_minor_ticks = numpy.arange(0, Emax, 0.05)   
##y_major_ticks = numpy.arange(0, 5, 1)                                              
##y_minor_ticks = numpy.arange(0, 5, 0.5)       
#
##ax.set_xticks(x_major_ticks)                                                       
##ax.set_xticks(x_minor_ticks, minor=True)                                           
##ax.set_yticks(y_major_ticks)                                                       
##ax.set_yticks(y_minor_ticks, minor=True)                                           
#
### and a corresponding grid                                                       
#
##ax.grid(which='both')                                                            
#
### or if you want differnet settings for the grids:                               
##ax.grid(which='minor', alpha=0.2)                                                
##ax.grid(which='major', alpha=0.5)   
##plt.text(0.12, 1.09, r'$109.33 meV$')
##plt.legend() 
##plt.axvline(x=V0, ymin=-30, ymax = 1, linewidth=2, color='k')
#
#
plt.show()


	
## First maximum of Transmission coeff between 0 eV and 0.1 eV
#for i = 0:len(Trans)
#    if (0 < E(i)) && (E(i)< 0.1)
#         En1(i) = E(i)
#         Tran1(i) = Trans(i)
#    end
#end
#[Transmax1,ind1] = max(Tran1)
#Emax1 = En1(ind1) # First maximum of Transmission coeff
#disp(['First maximum of Transmission coeff has the Energy  ', num2str(Emax1),' eV'])
#disp(' ')

## Second maximum of Transmission coeff between 0.1 eV and 0.24 eV
#for i = 0:len(Trans)
#    if (0.1 < E(i)) && (E(i)< V0)
#         En2(i) = E(i)
#         Tran2(i) = Trans(i)
#    end
#end



#import matplotlib.pyplot as plt  

#plt.plot(E,np.log(Trans),'b.')
#plt.axis([0,Emax,-30,0])
#plt.annotate('test',xy=( 0.15,-5),color='g')



#plt.show()





