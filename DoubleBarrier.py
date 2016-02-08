#  Resonant tunneling through 2 barriers with Vectorized
#  code and Using array not matrix

import numpy as np  
import math  
import scipy.constants as sc_const
import scipy 
import cmath


class nbarrier:
    
    def __init__(self,barrier,barrier_width,well_width,V0,n_points):
        self.barrier = 2		       	#number of potential barriers
        self.barrier_width = 10e-9		#barrier width (m), 10 Angstrom = 1e-9m = 1nm
        self.well_width    = 4e-9		#well width (m)
        self.V0 = 0.240			#barrier energy (eV),height of barrier
                                    	# 1 eV = 1.60217733e-19 J
        self.npoints = 10000			#number of points in energy plot
        self.N = (2*barrier)+1
        
        
        
    def potential_barrier_count(self):
        barrier_width,well_width, V0, N =self.barrier_width,\
                    self.well_width, self.V0, self.N
                    #number of samples of potential
        indj1 = 0
        indj2 = 1
        dL = np.empty(N,dtype=np.float64)
        V_barrier = np.empty(N,dtype=np.float64)	#set up position and potential array
        while indj1 <= N:
        	dL[indj1]=well_width
        	V_barrier[indj1]=0
        	while indj2 < N:
        		dL[indj2]=barrier_width
        		V_barrier[indj2]=V0
        		indj2 = indj2+2
        	indj1 = indj1+2
         
        return V_barrier
	
 #  dL = [well_width,barrier_width,well_width,...]   
 #  V  = [     0    ,       V0    ,     0    ,...]           

#%%%%%%%%%%%%                 
  
    def transmission(self):
        
        npoints, V0, N = self.npoints, self.V0, self.N
        V = nbarrier.potential_barrier_count()

        Emin    = math.pi *1e-5			#add (pi*1.0e-5) to energy to avoid divide by zero
        			#maximum particle energy (eV)
        dL = np.empty(N,dtype=np.float64)
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
          
          return Trans
        
              	
        #This creates a file Transmission.out that contains the values of the array of Trans in 10 characters wide and 10 digit precision in exponent form
        #np.savetxt('Transmission.out',Trans,'%10.10f','  ',)
