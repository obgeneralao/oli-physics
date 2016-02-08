#  Resonant tunneling through 2 barriers with Vectorized
#  code and Using array not matrix

import numpy as np  
import math  
import scipy.constants as sc_const
import time
import scipy
import cmath



barrier = 2;		       	#number of potential barriers
barrier_width = 10e-9;		#barrier width (m), 10 Angstrom = 1e-9m = 1nm
well_width    = 4e-9;		#well width (m)
V0 = 0.240;			#barrier energy (eV),height of barrier
                            	# 1 eV = 1.60217733e-19 J

                    
N=(2*barrier)+1            #number of samples of potential
indj1 = 0
indj2 = 1
dL = np.empty(N,dtype=np.float64)
V = np.empty(N,dtype=np.float64)	#set up position and potential array



while indj1 <= N:
	dL[indj1]=well_width
	V[indj1]=0
	while indj2 < N:
		dL[indj2]=barrier_width
		V[indj2]=V0
		indj2 = indj2+2
	indj1 = indj1+2
	
 #  dL = [well_width,barrier_width,well_width,...]   
 #  V  = [     0    ,       V0    ,     0    ,...]                            
  


Emin    = math.pi *1e-5			#add (pi*1.0e-5) to energy to avoid divide by zero
Emax    = 0.3				#maximum particle energy (eV)
npoints = 10000			#number of points in energy plot
# #  1,000,000 npoints takes about 441.0666 seconds/ 7.35 min to finish
# # 10,000,000  npoints takes about 
# # The Equation on how it takes for the program to finish 
# # is  clock(seconds) = 0.0004*npoints
dE      = V0/npoints		#energy increment (eV)
eye     = 0 + 1j	#square root of -1
#sc_const.m_e is  bare electron mass (9.109382e-31 kg)
m    = 0.06 * sc_const.m_e  #effective electron mass / m0
# Vectorization of codes, necessary for defining E,Trans, etc in the while loop below
E = np.empty(npoints)
Trans = np.empty(npoints)
k = np.empty(N,dtype=np.complex128)
M = np.empty([2,2],dtype=np.complex128)
#index initialization 
k_ind = 0
i = 0
n = 0
ind = 0

tic = time.clock()

while  k_ind < npoints:
	
	E[k_ind]= dE*(k_ind+1)+Emin
	Ms = np.array([[1,0],[0,1]],dtype=np. complex128)	#default value of propagation matrix for the whole system
	if i < N:
		k[i]=scipy.sqrt(2*m*sc_const.e*(E[k_ind]-V[i]))/sc_const.hbar  #wave vector at energy E
		# (sc_const.hbar)Planck's constant (J s), sc_const.e is electron charge (1.6021764e-19 C)
		
		#print V
		#print k
		#print i
		i +=1
		
	else:
		while n <(N-1):			#multiply out propagation matrix
			M[0,0]=0.5*(1+k[n+1]/k[n])*cmath.exp(-eye*k[n]*dL[n])
			M[0,1]=0.5*(1-k[n+1]/k[n])*cmath.exp(-eye*k[n]*dL[n])
			M[1,0]=0.5*(1-k[n+1]/k[n])*cmath.exp(eye*k[n]*dL[n])
			M[1,1]=0.5*(1+k[n+1]/k[n])*cmath.exp(eye*k[n]*dL[n])
			Ms = np.dot(Ms,M)
			n = n+1
		
		Trans[k_ind]=(abs(1/Ms[0,0]))**2	#transmission coefficient,normalized
	#print E
		#print i
		#print k_ind
		i = 0  # reset the value of i
		n = 0  #reset the value of n for E
		k_ind += 1

toc = time.clock()
print " The time to finish the while loop  is %f seconds with %d points" %(toc-tic,npoints)
	
#This creates a file Transmission.out that contains the values of the array of Trans in 10 characters wide and 10 digit precision in exponent form
#np.savetxt('Transmission.out',Trans,'%10.10f','  ',)


import matplotlib.pyplot as plt  


#fig = plt.figure(figsize=(10, 6))
#ax = fig.add_subplot(1,1,1)   
plt.axis([0,Emax,-30,1])
plt.plot(E,np.log(Trans),'b.')

## major ticks every 0.2, minor ticks every 0.05                                     
#x_major_ticks = numpy.arange(0, Emax, 0.1)                                              
#x_minor_ticks = numpy.arange(0, Emax, 0.05)   
#y_major_ticks = numpy.arange(0, 5, 1)                                              
#y_minor_ticks = numpy.arange(0, 5, 0.5)       

#ax.set_xticks(x_major_ticks)                                                       
#ax.set_xticks(x_minor_ticks, minor=True)                                           
#ax.set_yticks(y_major_ticks)                                                       
#ax.set_yticks(y_minor_ticks, minor=True)                                           

## and a corresponding grid                                                       

#ax.grid(which='both')                                                            

## or if you want differnet settings for the grids:                               
#ax.grid(which='minor', alpha=0.2)                                                
#ax.grid(which='major', alpha=0.5)   
#plt.text(0.12, 1.09, r'$109.33 meV$')
#plt.legend() 
#plt.axvline(x=V0, ymin=-30, ymax = 1, linewidth=2, color='k')


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





