# smolu_dsmc
Direct Simulation Monte Carlo code (Pyhton) to simulate the Smoluchowski equation with fragmentation.

The code is explained in detail in Reference 1:

Ref. 1: "Fragmentation and entanglement limit vimentin intermediate filament assembly", Quang D. Tran, Valerio Sorichetti, Gerard Pehau-Arnaudet, Martin Lenz, Cécile Leduc (BiorXiv -- doi: https://doi.org/10.1101/2022.03.19.484978). Link: https://www.biorxiv.org/content/10.1101/2022.03.19.484978v3

There are two versions of the code, both of which implement a Direct Simulation Monte Carlo scheme to solve the Smoluchowski equation with fragmentation. The base version is smoluchowski_dsmc_frag.py. 

In the second version (smoluchowski_dsmc_frag_entanglement.py), a multiplicative term is added to the annealing rate in order to take into account the steric effect of entanglements (see Ref. 1 for a more detailed discussion).
