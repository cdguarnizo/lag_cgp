% Demo for Laguerre approximation of any order differential equation
clear
close all
clc

%%%%%%%%%%%%%%%%%
% Generate Data %
%%%%%%%%%%%%%%%%%
tini = 0;
tfin = 4;
%Generate first order ODE data
addpath(genpath('../toolbox'));

%MIMO parameters
K = 50; %Number of points over latent functions
N = 50; %NUmber of points over output values

%ODE parameters
params = [1, 1, 4];
sys = tf(1, params);

lengthscale = 0.8;
t = linspace(tini, tfin, K)';
epsilon = 1e-10;
rbfkern = kernCreate(t, 'rbf');

rng(1e6)
% Sample from the latent forces
rbfkern.inverseWidth = 2./(lengthscale^2);
Kuquq = kernCompute(rbfkern, t);
tempo = gsamp(zeros(K, 1), Kuquq + epsilon*eye(K), 1);
u = tempo';
% Generate output from input
f = lsim(sys, u, t, 'foh');
f = f + .05*randn(length(f),1);

save demo01_50 f u t sys