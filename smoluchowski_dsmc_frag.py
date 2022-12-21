#Direct Simulation Monte Carlo code to simulate the Smoluchowski equation with fragmentation.
#
#Reference 1: "Fragmentation and entanglement limit vimentin intermediate filament assembly".
#Authors: Quang D. Tran, Valerio Sorichetti, Gerard Pehau-Arnaudet, Martin Lenz, Cécile Leduc.
#doi: https://doi.org/10.1101/2022.03.19.484978 (biorXiv).
#
#Author of this code: Valerio Sorichetti (valeriosorichetti@gmail.com)
#
####################################################################################

import os, sys, time, re
from numpy import *

#reaction rate
def k_rate(mi,mj):
	return 1./mi+1./mj

k_max0=2 				#Initial max value of K_ij -- updated during the simulation
Alpha=1 				#This parameter must lay between 0 and 1.
nevery_hist=10000 		#Print histogram every this many time steps
time_max=1e4 			#Stop the simulation when the waiting time is larger than this value.
ntot=10000	 			#Number of monomers
K_d=1.00e-3  			#Equilibrium dissociation constant
density=0.1 			#Monomer number density
f_max0=0.5*K_d*k_max0   #Initial max value of F_ij -- updated during the simulation
myseed=123 				#Seed for random number generation

if(K_d<0 or K_d>1):
	print("INPUT ERROR: Parameter K_d must be in interval [0,1].")
	sys.exit()

#Seed random generator
random.seed(myseed)

#Save initial values because we may want to check if they changed
k_max=k_max0
f_max=f_max0

os.system("mkdir histograms_ntot%d_K_d%.2e_density%.2e"%(ntot,K_d,density))

#The nth element of the array is the mass of molecule n; initially, only monomers are present
masses=zeros(ntot,dtype="int32")+1 

with open("mav_vs_t_ntot%d_K_d%.2e_density%.2e_tmax%.2e.dat"%(ntot,K_d,density,time_max),"w") as fout:
	waiting_time=0.
	step=0
	hist_count=0

	n_poly=len(masses[masses>1]) #Number of masses of size >1; mass=1 can't be fragmented
	n_mol=len(masses[masses!=0]) #Number of molecules

	while(waiting_time<time_max):
		step+=1
		
		#Compute coagulation probability
		p_coag=1./(1+(ntot*n_poly*f_max)/(n_mol*(n_mol-1)*density*k_max))

		if(random.rand()<p_coag):
			#Attempt coagulation

			#Randomly select 2 molecules
			#masses[i] must be not 0 since 0 is for absence of molecule
			#also i must be different from j

			i=random.choice(where(masses>0)[0])
			j=int(random.rand()*ntot)

			while(masses[j]==0 or j==i):
				j=int(random.rand()*ntot)

			mi=masses[i]
			mj=masses[j]

			#Calculate bonding rate
			k_ij=k_rate(mi,mj)

			if(k_ij>k_max): 
				#Update max. rate estimate
				k_max=k_ij

			else:
				if(random.rand()<k_ij/k_max):
					#Coagulation is performed

					#Increment time
					waiting_time+=2*Alpha*ntot/(n_mol*(n_mol-1)*density*k_ij)

					#Update masses array
					masses[j]=mi+mj
					masses[i]=0

					#Update relevant quantities
					n_mol-=1

					if(mi==1 and mj==1):
						n_poly+=1
					elif(mi>1 and mj>1):
						n_poly-=1

					#Compute mean and max mass
					mean_mass=mean(masses[masses!=0])
					reac_extent=1-float(n_mol)/ntot

					#Print relevant quantities
					fout.write("%.4e %.4e %.4e\n"%(waiting_time,mean_mass,reac_extent))

		else:
			if(n_poly>0):
				#Attempt fragmentation

				#Randomly select 1 molecule with mass>1 (mass=1 cannot break)
				k=random.choice(where(masses>1)[0])
				mk=masses[k]

				#Randomly select fragments size:
				m1=1+int(random.rand()*(mk-1))
				m2=mk-m1

				#Calculate fragmentation rate respecting detailed balance:
				f_m1m2=0.5*K_d*k_rate(m1,m2)

				if((mk-1)*f_m1m2>f_max): 
					#Update max. rate estimate
					f_max=(mk-1)*f_m1m2

				else:
					if(random.rand()<(mk-1)*f_m1m2/f_max):	
						#Fragmentation is performed

						#Increment time
						waiting_time+=2*(1-Alpha)/(n_poly*(mk-1)*f_m1m2)

						#Update masses array
						masses[k]=m1
						masses[random.choice(where(masses==0)[0])]=m2

						#Update relevant quantities
						n_mol+=1
						
						if(m1==1 and m2==1):
							n_poly-=1
						elif(m1>1 and m2>1):
							n_poly+=1

						#Compute mean and max mass
						mean_mass=mean(masses[masses!=0])
						reac_extent=1-float(n_mol)/ntot

						#Print relevant quantities
						fout.write("%.4e %.4e %.4e\n"%(waiting_time,mean_mass,reac_extent))

		#Save and print histogram
		if(hist_count%nevery_hist==0):
			masses_values=sort(array(list(set(masses[masses!=0]))))
			masses_occurrences=array([list(masses).count(x) for x in masses_values])
			hist=column_stack((masses_values,masses_occurrences))
			savetxt('histograms_ntot%d_K_d%.2e_density%.2e/n%.10d_t%.3e'%(ntot,K_d,density,hist_count,waiting_time),hist,fmt='%d %d')			
		hist_count+=1

#Save last array so that we can resume simulation if needed
savetxt('masses_final_ntot%d_K_d%.2e_density%.2e_t%.3e.dat'%(ntot,K_d,density,waiting_time),masses,fmt='%d')
