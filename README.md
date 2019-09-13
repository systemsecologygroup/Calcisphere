# calcisphere
Calculate the energy requirements for calcium transport into vesicles vs over the plasma membrane in calcifying nannoplankton.

# General description
This code calculates the energy required to maintain calcium homeostasis in single celled organisms like calcifying nannoplankton. 
In order to maintain calcium homeostasis, a cell must export as many calcium ions as are entering the cytoplasm. Calcium enters 
the cell via channels that can be either opened or closed. The channels are relatively large pores that allow for a fast entry 
of calcium into the cell. The spike of increased intracellular calcium due to opening of the channels is used for intracellular 
signalling. However, in order to establish the required low intracellular calcium concentrations required for optimal cell 
functioning, cells use a combination of Na-Ca-exchangers, that have a relatively high maximal transport rate but a low affinity 
for intracellular calcium, and Ca-ATPases that due to a higher affinity for intracellular calcium are able to draw down 
intracellular calcium to micromolar concentrations. Since active ion transport requires energy in the form of ATP, we use this to 
estimate the energy required to establish a balance between passive influx and active export.

# Model description
In our model we assume that, in order to establish calcium homeostasis, the active ion transport equals passive calcium entry 
through the channels. The channels are simulated with a constant permeability that depends on the concentration gradient over 
the plasma membrane. Since intracellular calcium concentrations are usually in the micromolar range but extracellular calcium 
is around 10mM in the ocean, this gradient is driving the calcium influx.

Ca_out=10.0 # mM; extracellular calcium concentration in the ocean 10 mM=10e-3 mol L-1 

Ca_in=0.10e-3 # mM; intracellular calcium concentration


The ion flux for one calcium channel is 1pA, which corresponds to 3.0e6 divalent ions per second (Tsien 1983).

i=3.0e6 # Ca2+ s-1 channel-1

The density of calcium channels on the cell surface then determine the maximal flux of calcium ions entering the cytoplasm. It is generally 
believed that for coccolithophores to calcify, the calcium flux needs to be very high. However, the density of calcium channels on 
the plasma membrane is not known for calcifying nannoplankton. We therefore use the highest known density of calcium channels on 
snail axons (Tsien 1983) as a potential example to relate the observed calcium flux required for coccolithophore calcification to a possible 
calcium influx as observed in the animal kingdom.

N=1.0e12 # channels m-2

The maximum possible calcium flux with this density of calcium channels and the ion flux per channel is 2.3475e-14 mol s-1 given the surface 
area of the considered cell.

#r_cyt=10.0e-6 # diameter of calcispheres is around 20 micrometers
r_cyt=2.5e-6 # diameter of E. huxleyi is around 5 micrometers (Harvey et al 2015)
A_cyt=4.0*np.pi*r_cyt**2.0 # cell surface in m2

I_Ca=N*i # 1.8e20 Ca ions m-2 s-1
N_A=6.0221367e23 # ions per mol
J_Ca=I_Ca/N_A # 2.98e-4 mol Ca m-2 s-1 (as maximum flux)
F_Ca=J_Ca*A_cyt # 2.3475e-14 mol s-1

Since this flux is driven by the calcium gradient over the membrane, we can calculate the permeability of the membrane (PCa).

PCa=F_Ca/(Ca_bd-Ca_in) # because PCa is the rate that leads to this flux given the ion gradient

This calculation uses the calcium concentration directly at the cell surface (Ca_bd), i.e. the calcium concentration of the boundary layer around 
the cell and not the concentration in the open ocean. This concentration is calculated based on the diffusive flux of calcium around the 
cell as a consequence of calcium depletion due to calcification. For the calcification flux we use the precipitation of one coccolith with 22 fmol 
per one hour as estimated by (Holtz et al 2013). The original 
observation is from Paasche. The diffusion coefficient for calcium is set after Li and Gregory (1973).

Ca_bd=Ca_out-QCa/(4.0*np.pi*r_cyt*D_Ca)
QCa=6.11e-18 # 22 fmol h-1 -> 6.11e-18 mol s-1
D_Ca=7.93e-6/10000.0 # cm2 s-1 -> m2 s-1 (after Li and Gregory 1974)




