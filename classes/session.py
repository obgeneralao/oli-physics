# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 17:44:31 2016

@author: oli
"""



from nbarrier import n_barrier
from ROOT import TCanvas, TGraph

from ROOT import gROOT
import numpy as np

gROOT.Reset()





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





Emax    =  V0 + 0.1	


c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )

c1.SetFillColor( 42 )
c1.SetGrid()

gr = TGraph( npoints, E, np.log(Trans) )
gr.SetLineColor( 2 )
gr.SetLineWidth( 4 )
#gr.SetMarkerColor( 4 )
#gr.SetMarkerStyle( 21 )
gr.SetTitle( 'Transmission of at n-barrier' )
gr.GetXaxis().SetTitle( 'Energy' )
gr.GetYaxis().SetTitle( 'ln of Transmission' )
gr.GetYaxis().SetRangeUser(-30.0,0.0 )
gr.Draw( 'ACP' )

# TCanvas.Update() draws the frame, after which one can change it
c1.Update()
c1.GetFrame().SetFillColor( 21 )
c1.GetFrame().SetBorderSize( 12 )
c1.Modified()
c1.Update()

	
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





