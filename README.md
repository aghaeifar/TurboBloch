# bloch_simulator
Efficient and fast implementation of Bloch equation simulator for magnetic resonance imaging (MRI) sequences, supporting parallel transmission (pTx). 

MATLAB mex wrapper is provided and works with MATLAB R2018a and newer releases. Use command below to compile mex file:

mex bloch_sim.cpp bloch_sim_mex.cpp -R2018a
