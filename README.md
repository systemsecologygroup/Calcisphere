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

The ion flux for one calcium channel is 1pA, which corresponds to 3.0e6 divalent ions per second (Tsien 1983).
i=3.0e6 # Ca2+ s-1 channel-1

