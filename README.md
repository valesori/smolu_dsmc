# smolu_dsmc
Direct Simulation Monte Carlo code (Pyhton) to simulate the Smoluchowski equation with fragmentation.

The codes are explained in detail in Ref. 1 and 2 (see below)

Three versions of the code are provided, both of which implement a Direct Simulation Monte Carlo scheme to solve the Smoluchowski equation with fragmentation.

smoluchowski_dsmc_frag.py mplements the classical Smoluchowski coagulation equation, by also taking into account the possibility of fragmentation (see Ref. 1).

smoluchowski_dsmc_frag_entanglement.py is a modified version of smoluchowski_dsmc_frag.py, in which a multiplicative term is added to the annealing rate in order to take into account the steric effect of entanglements (see Ref. 1).

smoluchowski_dsmc_rings.py simulates a modified version of the Smoluchowski equation, in the absence of fragmentation (irreversible coagulation), but taking into account the possibility of loop formation (see Ref. 2).

#######################################

References:

(1) "Fragmentation and entanglement limit vimentin intermediate filament assembly", Quang D. Tran, Valerio Sorichetti, Gerard Pehau-Arnaudet, Martin Lenz, Cécile Leduc (BiorXiv - DOI: https://doi.org/10.1101/2022.03.19.484978). Link: https://www.biorxiv.org/content/10.1101/2022.03.19.484978v3

(2) "Runaway Transition in Irreversible Polymer Condensation with Cyclisation", Maria Panoukidou, Simon Weir, Valerio Sorichetti, Yair Gutierrez Fosado, Martin Lenz, Davide Michieletto (arXiv - DOI: https://doi.org/10.48550/arXiv.2210.14010). Link: https://arxiv.org/abs/2210.14010. 
