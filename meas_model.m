function [X] = meas_model(theta,S,globals)

%global noisesd L NR

% rng('shuffle');

noisesd = globals.noisesd;
L = globals.L;
NR = globals.NR;

Nreal = normrnd(0,noisesd,NR,L);
Nimag = normrnd(0,noisesd,NR,L);

X = H(theta,globals)*S + Nreal + sqrt(-1)*Nimag;

