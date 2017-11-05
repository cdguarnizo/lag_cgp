% Demo for Laguerre approximation of first order differential equation
clear
close all
clc

addpath(genpath('../toolbox'));

%%%%%%%%%%%%%
% Load Data %
%%%%%%%%%%%%%
load('../datasets/demo01.mat')

y = [f(:); u(:)];
X.val = [t; t];
X.index = [ones(length(t),1); 2*ones(length(t),1)];

%%%%%%%%%%%%%%%%%%%%%%%%
% Train Laguerre Model %
%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../modelio', '../laguerre');

options.kernType = 'laguerre';
options.nlf = 1;
options.optimiser = 'scg';
options.tieOptions.tieIndices = true;
options.includeNoise = true;
options.isVarS = false;
options.kern.isArd = false;
options.kern.order = 10;
options.beta = 1e2;

Nr = 10;
hres = zeros(Nr,length(t));
gammat = zeros(1,Nr);
rng(1e6)
for k = 1:Nr,
    fprintf('Learning model %d\n',k)
    model = ioftcmgpCreate(X, y, 1, 1, options);
    model = ioftcmgpOptimise(model, 0, 100);
    gammat(k) = model.kern.gamma;
    hres(k,:) = sum(repmat(model.kern.c', 1,length(t)).*laguerreFuncEval(t(:)', model.kern.order, model.kern.gamma), 1);
end

hreso = impulse(sys,t);

save results_demo01 hres hreso f u t sys gammat

%%
%Plot results
load results_demo01

markerSize = 20;
markerWidth = 6;
markerType = 'k.';
lineWidth = 2;
fillColor = [0.8 0.8 0.8];
figure;

clf
% Mean and standard deviation
ym = mean(hres)';
yS = std(hres)';
fill([t; t(end:-1:1)], ...
    [ym; ym(end:-1:1)] ...
    + [yS; -yS(end:-1:1)], ...
    fillColor,'EdgeColor',fillColor)
hold on;
h = plot(t, ym, 'k--');
set(h, 'linewidth', lineWidth)

% Plot real function
p = plot(t, hreso, 'k-');
set(p, 'linewidth', lineWidth)

hold off;
grid on;
%legend('Deviation', 'Mean', 'True');
%title('Impulse response approximation')

matlab2tikz('Results01.tikz', 'showInfo', false, ...
        'parseStrings',false,'standalone', false, 'width', '\figurewidth', 'height', '\figureheight')

smse = zeros(size(hres,1),1);
rmse = zeros(size(hres,1),1);
for k = 1:10,
    %smse(k) = mysmse(hreso',hres(k,:));
    rmse(k) = sqrt(mean((hreso'-hres(k,:)).^2));
end