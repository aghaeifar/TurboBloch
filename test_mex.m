clc
mex bloch_sim_mex.cpp bloch_sim.cpp -R2018a

%% test speed
gamma = 267522187.44;
flipangle = 20 * pi/180;
ntime = 400; % number of samples
dur = 2; % second
npos = 98952; % number of spatial positions

b1 = complex(ones(ntime,1) * flipangle/gamma/dur);
tp = ones(ntime,1) * dur/ntime;

gr = zeros(ntime, 3);
b0 = zeros(npos, 1);
pr = rand(npos, 3);

tic
bloch_sim_mex(b1, gr, tp, b0, pr);
toc
