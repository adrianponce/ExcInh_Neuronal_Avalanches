## CONTENT ##

Codes to detect and analyze the neuronal avalanches, and to model them:

- _avalanche_analysis_EI.m_ : This function constructs the activity tensor, calculates the connected components at each time frame, detects the propagating clusters (neuronal avalanches), and calculate the durations and sizes of avalanches.

- _avalanche_timing_EI.m_ : This function computes the initiation and termination time of avalanches, together with the shape of the avalanches.

- _Run_EImodel.m_ : Runs the network model, detects neuronal avalanches, and performs their analysis.

- _Get_Connectivity_matrix.m_ : This function constructs the model’s connectivity

- _Gillespie_EImodel.m_ :  Simulates a stochastic Wilson-Cowan network model using the Gillespie algorithm. This function is based on Edward Wallace's code, publicly available here: https://github.com/ewallace/stochsimcode

- _SpikesToFluoresence.m_ : This function implements a phenomenological model that converts spike times  to synthetic fluorescence time series; as in Wei et al. (2020).

- _Get_NonSpatialAvalanches.m_ : Detection of neuronal avalanches using a spatially unconstraint definition.

- _neurons_type_and_coords.mat_ : example cell identities and coordinates.


To study the power-law distributions of neuronal avalanches, please download the MATLAB NCC (Neural Complexity and Criticality) Toolbox, from Marshall et al. (2016) https://doi.org/10.3389/fphys.2016.00250


**References**

Benayoun, M., Cowan, J.D., van Drongelen, W., and Wallace, E. (2010). Avalanches in a Stochastic Model of Spiking Neurons. PLoS Comput. Biol. 6, e1000846.
Marshall, N., Timme, N.M., Bennett, N., Ripp, M., Lautzenhiser, E., and Beggs, J.M. (2016). Analysis of power laws, shape collapses, and neural complexity: new techniques and MATLAB support via the NCC toolbox. Front. Physiol. 7: 250.
Wei, Z., Lin, B.-J., Chen, T.-W., Daie, K., Svoboda, K., and Druckmann, S. (2020). A comparison of neuronal population dynamics measured with calcium imaging and electrophysiology. PLoS Comput. Biol. 16(9): e1008198
