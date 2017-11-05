function model = ioftcmgpCreate(X, y, nout, nin, options)

% IOFTCMMGPCREATE
% Based on multigp code from GPmat. This function allows the creation of
% full GP models for a type of global kern.
% Inputs:
% X: Space/time data.
% y: Output data.
% Options: Definition of flag variables and type of inference.
% FTCMMGP

model.type = 'ioftcmgp';

model.isVarS = false;

%Each realization contain a cell of the output values
model.nout = nout;
model.nlf = nin;
model.kernType = options.kernType;
model.ndim = size(X.val,2);
model.optimiser = options.optimiser;
model.approx = 'ftc';

if isfield(options, 'scale') && ~isempty(options.scale)
    model.scale = options.scale;
else
    model.scale = ones(1, model.nout);
end
if isfield(options, 'bias') && ~isempty(options.bias)
    model.bias = options.bias;
else
    model.bias = zeros(1, model.nout);
end

model.tieIndices = options.tieOptions.tieIndices;
model.includeNoise = options.includeNoise;

% Set the dimension of each output in model, which will be useful for
% computing sizes and things like that 
%[param, names] = spmultimodelExtractParam(model)
model.outX = X;
model.y = y;
kern.tspace = zeros(1,2);

%Create index for multioutput data
model.index = zeros(2,model.nout);
start = 1;
for d = 1:model.nout,
    ind = X.index == d;
    maxX = max(X.val(ind));
    if kern.tspace(2) < maxX,
        kern.tspace(2) = maxX;
    end
    model.sizeX(d) = sum(ind);
    model.index(:,d) = [start; start + model.sizeX(d)-1];
    start = start + model.sizeX(d);
end

% Kern definition
% This is only for conv. cov. functions where the convoloution is approximated
if isfield(options, 'tsize'),
    kern.tsize = options.tsize;
else
    kern.tsize = 1000;
end
kern.type = options.kernType;
kern.inputDimension = model.ndim;
kern.options = options.kern;
kern.options.nlf = nin;
kern.options.nout = nout;
kern.options.tieOutputParams = true;
kern.options.approx = model.approx;
kern = kernParamInit(kern);
model.kern = kern;
model.kern.isVarS = options.isVarS;
model.kernType = kern.type;
model.nParams = model.kern.nParams;

% % Count number of parameters, adding lengthscale values for n realizations
% model.nParams = numParams - 1;
% 
% % Create Length-scales parameters for multi model
% if isfield(options, 'ls') && ~isempty(options.ls),
%     if size(options.ls, 2) == model.nrep,
%         model.ls = options.ls;
%     else
%         model.ls = options.ls*rand(1,model.nrep);
%     end
%     model.lsTransform =  optimiDefaultConstraint('positive');
%     model.nParams = model.nParams + model.nrep;
% end

% Set up a mean function if one is given.
if isfield(options, 'meanFunction') && ~isempty(options.meanFunction),
    % TODO: this is not implemented yet in Extract and Expand params
    % functions
    if isstruct(options.meanFunction),
        model.meanFunction = options.meanFunction;
    else
        if ~isempty(options.meanFunction),
            model.meanFunction = meanCreate(model.ndim, model.nout, X, y, options.meanFunctionOptions);
        end
    end
    model.nParams = model.nParams + model.meanFunction.nParams;
end

% Create noise model
if isfield(options, 'beta') && ~isempty(options.beta) && options.includeNoise,
    if size(options.beta,2) == model.nout,
        model.beta = options.beta;
    else
        model.beta = options.beta*ones(1,model.nout);
    end
    model.betaTransform =  optimiDefaultConstraint('positive');
    model.nParams = model.nParams + model.nout;
end

% Allocate temporal memory spaces
model.Kyyinv = [];
model.alpha = [];
model.logDet = 0;

model = ioftcmgpComputeKernel(model);

params = modelExtractParam(model);
model = modelExpandParam(model, params);