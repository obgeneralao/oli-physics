#  Resonant tunneling through 2 barriers with Vectorized
#  code and Using array not matrix

import numpy as np  
import math  
import scipy.constants as sc_const
import scipy 
import cmath
from numba import jit


class n_barrier:
    
    def __init__(self,barrier,barrier_width,well_width,V0,npoints,Emax):
        self.barrier = barrier		       	#number of potential barriers
        self.barrier_width = barrier_width		#barrier width (m), 10 Angstrom = 1e-9m = 1nm
        self.well_width    = well_width		#well width (m)
        self.V0 = V0			#barrier energy (eV),height of barrier
                                    	# 1 eV = 1.60217733e-19 J
        self.npoints = npoints			#number of points in energy plot
        self.N = (2*barrier)+1
        self.Emax = Emax
        
        
#    def potential_barrier_count(self):
#        barrier_width,well_width, V0, N =self.barrier_width,\
#                    self.well_width, self.V0, self.N
#                    #number of samples of potential
#        indj1 = 0
#        indj2 = 1
#        dL = np.empty(N,dtype=np.float64)
#        V = np.empty(N,dtype=np.float64)	#set up position and potential array
#        while indj1 <= N:
#        	dL[indj1]=well_width
#        	V[indj1]=0
#        	while indj2 < N:
#        		dL[indj2]=barrier_width
#        		V[indj2]=V0
#        		indj2 = indj2+2
#        	indj1 = indj1+2
#         
#        return V
	
 #  dL = [well_width,barrier_width,well_width,...]   
 #  V  = [     0    ,       V0    ,     0    ,...]           

               
    #@jit  
    def transmission(self):
        
        barrier_width,well_width, V0, N, npoints,Emax =self.barrier_width,\
                    self.well_width, self.V0, self.N, self.npoints, self.Emax
                    #number of samples of potential
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
         
              

        Emin    = math.pi *1e-5			#add (pi*1.0e-5) to energy to avoid divide by zero
        			#maximum particle energy (eV)
        dE      = Emax/npoints		#energy increment (eV)
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
          
        return Trans, E
        
              	
        #This creates a file Transmission.out that contains the values of the array of Trans in 10 characters wide and 10 digit precision in exponent form
        #np.savetxt('Transmission.out',Trans,'%10.10f','  ',)
               

#    if __name__ == '__main__':
#        potential_barrier_count()
#        transmission()
    