# Efficient-coding-in-a-spiking-predictive-coding-network
MATLAB code for a network filter model, consisting of 'filter-and-fire' neurons with both recurrent and feed-forward filters, that performs close to optimal tracking of its input. 

When using this code, please cite the following preprint: 

Zeldenrust, F., Gutkin, B., & Denève, S. (2019). Efficient and robust coding in heterogeneous recurrent networks. BioRxiv. https://doi.org/10.1101/804864


# About
This MATLAB toolbox creates and runs the filter network for efficient coding described in the paper above. A few examples scripts are given:
* In 'Example_run_model' it is shown how to create representing filters for each neuron, and then run the network. It first creates a set of basis-kernels (make_basisfunctions). Next, it creates a set of feed-forward filters as a random combination of these basisfunctions, and creates the corresponding recurrent filters for optimal coding (generate_filters). The network can now track any input signal (run_model). 
* In 'Example_run_network_paper' the different networks as used in the paper (parameter 'net' sets it to homogeneous, heterogeneous and type 1 & type 2) are pre-defined. A network is then created (make_kernels_network) and run.
* The following scripts will run the simulations for the figures of the paper:
  * Figure 3: run_model_if_sta_prc
  * Figure 4: run_model_redundancy
  * Figure 5: run_model_noise
  * Figure 6: run_model_correlations


# Use
see: Example_run_model and Example_run_network_paper.
