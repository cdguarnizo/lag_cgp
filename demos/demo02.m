% Demo for Laguerre approximation of first order differential equation
clear
clc

addpath(genpath('../toolbox'),'../modelio', '../laguerre');

%%%%%%%%%%%%%%%%%
%   Load Data   %
%%%%%%%%%%%%%%%%%
load('../datasets/mgp2o2i.mat')
N = length(t);
K = N;
y = [f(:); uq(:)];
X.val = [t; t; t; t];
X.index = [ones(N,1); 2*ones(N,1); 3*ones(N,1); 4*ones(N,1)];

ytrain = y;
ytrain([K+(30:40),2*K+(22:32)]) = [];
Xtrain.val = X.val;
Xtrain.val([K+(30:40),2*K+(22:32)]) = [];
Xtrain.index = X.index;
Xtrain.index([K+(30:40),2*K+(22:32)]) = [];

ytest = y([K+(30:40),2*K+(22:32)]);
Xtest.val = X.val([K+(30:40),2*K+(22:32)]);
Xtest.index = X.index([K+(30:40),2*K+(22:32)]);

%%%%%%%%%%%%%%%%%%%%%%%%
% Train Laguerre Model %
%%%%%%%%%%%%%%%%%%%%%%%%
options.kernType = 'laguerre';
options.nlf = 1;
options.optimiser = 'scg';
options.tieOptions.tieIndices = true;
options.includeNoise = true;
options.isVarS = false;
options.kern.isArd = false;
options.beta = 1e-2;
options.kern.order = 10;

rng(1e5)
model = ioftcmgpCreate(Xtrain, ytrain, 2, 2, options);
model = ioftcmgpOptimise(model, 1, 100);

[ym, ySd] = ioftcmgpPosterior(model, X);

%% Estimation of impulse response
hres = zeros(2,N);
td = t;
for d=1:2,
    hres(d,:) = sum(repmat(model.kern.c(1+(d-1)*(model.kern.order +1): d*(model.kern.order +1))', 1,length(td)).*laguerreFuncEval(td(:)', model.kern.order, model.kern.gamma(d)), 1);
end

save model2o2i X model ytrain Xtrain ytest Xtest ym ySd N t
%% Plot results
addpath(genpath('../toolbox'),'../modelio', '../laguerre');
load('model2o2i.mat')
close all
markerSize = 10;
markerWidth = 5;
markerType = 'k.';
lineWidth = 1.5;
fillColor = [0.8 0.8 0.8];
figure;
Names = cell(1,4);
Names{1} = 'Output 1';
Names{2} = 'Output 2';
Names{3} = 'Input 1';
Names{4} = 'Input 2';
smse = zeros(4,1);
nlpd = smse;
for d = 1:max(X.index),
    %subplot(2,2,d)
    figure(d);

    % GP prediction along all time inputs
    XPred = X.val(X.index == d);
    yPred = ym(X.index == d);
    yS = real(ySd(X.index == d));
    if ~all(yS == 0.),
    h = fill([XPred; XPred(end:-1:1)], ...
        [yPred; yPred(end:-1:1)] ...
        + 2*[yS; -yS], ...
        fillColor,'EdgeColor',fillColor);
    end
    hold on;
    h = plot(XPred, yPred, 'k-');
    
    % Plot training data
    XTr = Xtrain.val(Xtrain.index == d);
    yTr = ytrain(Xtrain.index == d);
    set(h, 'linewidth', lineWidth)
    p = plot(XTr, yTr, 'r.');
    set(p, 'markersize', 10, 'lineWidth', 1.5);
    indTe = find(Xtest.index == d);
    if ~isempty(indTe),
        XTe = Xtest.val(indTe);
        yTe = ytest(indTe);
        % Plot test data
        p2 = plot(XTe, yTe, 'b.');
        set(p2, 'markersize', 10, 'lineWidth', 1.5);
        
        indd = find(X.index == d);
        ini = find(XTe(1)==X.val(indd));
        fin = find(XTe(end)==X.val(indd));
        smse(d) = mysmse(yTe,yPred(ini:fin));
        nlpd(d) = mynlpd(yTe,yPred(ini:fin),yS(ini:fin).^2+1e-6);
    end
    hold off;
    %grid on;
    if d==5,
        legend('GP std','GP mean', 'Train data', 'Test data');
    end
    xlim([0,5])
    box on
    %title(Names{d})
end

hres = zeros(2,N);
h = hres;
irfsmse = zeros(2,1);
irfrmse = zeros(2,1);
td = t;
for d=1:2,
    hres(d,:) = sum(repmat(model.kern.c(1+(d-1)*(model.kern.order +1): d*(model.kern.order +1))', 1,length(td)).*laguerreFuncEval(td(:)', model.kern.order, model.kern.gamma(d)), 1);
    h(d,:) = impulse(sys{d},t);
    irfsmse(d) = mysmse(h(d,:),hres(d,:));
    irfrmse(d) = sqrt(mean((h(d,:)-hres(d,:)).^2));
    figure;
    hold on;
    plot(t,h(d,:),'-k','lineWidth', 1.5)
    plot(t,hres(d,:),'--k','lineWidth', 1.5)
    grid on;
end
