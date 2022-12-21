#Direct Simulation Monte Carlo code to simulate the irreversible Smoluchowski equation with ring formation.
#
#Reference 1: "Runaway Transition in Irreversible Polymer Condensation with Cyclisation"
#Authors: Maria Panoukidou, Simon Weir, Valerio Sorichetti, Yair Gutierrez Fosado, Martin Lenz, Davide Michieletto 
#doi: https://doi.org/10.48550/arXiv.2210.14010
#
#Author of this code: Valerio Sorichetti (valeriosorichetti@gmail.com)
#
####################################################################################

import os, sys, time, re
from numpy import *

#Global_variables
glob_pref=float(input("Enter prefactor for ring formation rate:"))

#Coagulation rate K_ij
def k_rate(mi,mj):
	return (sqrt(mi)+sqrt(mj))*(1./mi+1./mj) #De Gennes kernel for Gaussian chains

#Cyclization rate R_i
def r_rate(mi):
	return glob_pref*mi**(-1.5)

k_max0=4 				#Initial max value of K_ij (coagulaton rate) -- updated during the simulation
r_max0=glob_pref*0.5 	#Initial max value of R_i (cyclization rate) -- updated during the simulation
Alpha=1 				#A parameter that must be between 0 and 1. 1 is the value suggested.
nevery_hist=10000 		#Print histogram every this many time steps
time_max=1e5 			#Stop the simulation when the waiting time is larger than this value.

ntot=int(input("Enter n. of particles:"))
density=float(input("Enter reduced monomer density:"))
myseed=int(input("Enter a seed for the random number generator:"))

#Seed random generator
random.seed(myseed)

#Save initial values because we may want to check if they changed
k_max=k_max0
r_max=r_max0

#The nth element of the array is the mass of molecule n; initially, only monomers are present
masses=zeros(ntot,dtype="int32")+1 
rings=[]
with open("mav_vs_t_ntot%d_density%.2e_pref%.2f.dat"%(ntot,density,glob_pref),"w") as fout_mav:
	with open("nchain_rings_ntot%d_density%.2e_pref%.2f.dat"%(ntot,density,glob_pref),"w") as fout_nchain:
		waiting_time=0.
		step=0
		hist_count=0

		n_rings=len(rings)
		n_chains=len(masses[masses!=0]) #Number of molecules

		while(waiting_time<time_max and n_chains>1):
			step+=1
			
			#Compute coagulation probability
			p_coag=1./(1+(2*ntot*r_max)/((n_chains-1)*density*k_max))

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
						waiting_time+=2*Alpha*ntot/(n_chains*(n_chains-1)*density*k_ij)

						#Update masses array
						masses[j]=mi+mj
						masses[i]=0

						#Update relevant quantities
						n_chains-=1

						#Compute mean chain mass
						mean_mass=float((sum(masses[masses!=0])+sum(rings)))/(len(masses[masses!=0])+len(rings))
						mean_chain_mass=mean(masses[masses!=0])

						#Cumpute mean ring mass
						if(n_rings>0):
							mean_ring_mass=mean(rings)
						else:
							mean_ring_mass=0.
							
						#Print relevant quantities
						fout_mav.write("%.4e %.4e %.4e %.4e\n"%(waiting_time,mean_mass,mean_chain_mass,mean_ring_mass))
						fout_nchain.write("%.4e %d %d\n"%(waiting_time,n_chains,n_rings))
						print(n_chains,n_rings)

			else:
				#Attempt cyclization

				k=random.choice(where(masses>0)[0])
				mk=masses[k]
				rmk=r_rate(mk)

				if(r_rate(mk)>r_max): 
					#Update max. rate estimate
					r_max=rmk

				else:
					if(random.rand()<r_rate(mk)/r_max):	
						#Cyclization is performed

						#Increment time
						waiting_time+=(1-Alpha)/(n_chains*rmk)

						#Update rings array
						rings.append(mk)
						masses[k]=0 #once the ring is formed, we remove the mass

						#Update relevant quantities
						n_rings+=1
						n_chains-=1		

						#Compute mean chain mass
						mean_mass=float((sum(masses[masses!=0])+sum(rings)))/(len(masses[masses!=0])+len(rings))
						mean_chain_mass=mean(masses[masses!=0])

						#Compute mean ring mass
						if(n_rings>0):
							mean_ring_mass=mean(rings)
						else:
							mean_ring_mass=0.	

						#Print relevant quantities
						fout_mav.write("%.4e %.4e %.4e %.4e\n"%(waiting_time,mean_mass,mean_chain_mass,mean_ring_mass))
						fout_nchain.write("%.4e %d %d\n"%(waiting_time,n_chains,n_rings))
						print(n_chains,n_rings)

