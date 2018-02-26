# Efficient-coding-in-a-spiking-predictive-coding-network
MATLAB code for a network filter model, consisting of 'filter-and-fire' neurons with both recurrent and feed-forward filters, that performs close to optimal tracking of its input. 

When using this code, please reference the following Cosyne abstract: 

Zeldenrust, F., Denève, S., & Gutkin, B. S. (2013). Matching encoding and decoding with spiking neurons. In Cosyne Abstracts 2013, Salt Lake City USA.

# About
This MATLAB toolbox creates a filter network. First, it creates a set of basis-kernels (make_basisfunctions). Next, it creates a set of feed-forward filters as a random combination of these basisfunctions, and creates the corresponding recurrent filters for optimal coding (generate_filters). The network can now track any input signal (run_model). 

# Use
see: Example_run_model
