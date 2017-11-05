% Demo for Laguerre approximation of any order differential equation
clear all
clc

%%%%%%%%%%%%%%%%%
% Generate Data %
%%%%%%%%%%%%%%%%%
tini = 0;
tfin = 5;
%Generate first order ODE data
rng(4e6)

addpath(genpath('../../toolbox'));

%MIMO parameters
D = 2;
Q = 2;
Sdq = 1;

K = 50; %Number of points over latent functions
N = 50; %NUmber of points over output values

%ODE parameters
sys = cell(D, 1);
params = [1, 4, 9, 14, 8];
sys{1} = tf(1, params);
params = [1, 1];
sys{2} = tf(1, params);

% Sample from the latent forces
lengthscales = [0.8, 1.5];
t = linspace(tini, tfin, K)';
epsilon = 1e-10;
rbfkern = kernCreate(t, 'rbf');
uq = zeros(K, Q);
for q = 1:Q,
    rbfkern.inverseWidth = 2./(lengthscales(q)^2);
    Kuquq = kernCompute(rbfkern, t);
    tempo = gsamp(zeros(K, 1), Kuquq + epsilon*eye(K), 1);
    uq(:, q) = tempo';
end

f = zeros(N, D);
for d = 1:D,
    for q = 1:Q,
        f(:,d) = f(:,d) + lsim(sys{d}, uq(:,q)', t', 'foh');
    end
    f(:,d) = f(:,d) + .1*max(f(:,d))*randn(size(f(:,d)));
end

figure(1);
subplot(221)
plot(f(:,1));
subplot(222)
plot(f(:,2));
subplot(223)
plot(uq(:,1));
subplot(224)
plot(uq(:,2));

figure(2);
subplot(211)
impulseplot(sys{1},5)
subplot(212)
impulseplot(sys{2},5)