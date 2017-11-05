function kern = laguerreKernParamInit(kern)

% LAGUERREKERNPARAMINIT
%
% COPYRIGHT : Cristian Guarnizo, 2015
% MULTIGP

if isfield(kern, 'options') && isfield(kern.options, 'nout')
    kern.nout = kern.options.nout;    
else
    error('Number of outputs is required for this kernel')
end

if isfield(kern, 'options') && isfield(kern.options, 'nlf')
    kern.nlf = kern.options.nlf;    
else
    error('Number of latent forces is required for this kernel')
end

if isfield(kern, 'options') && isfield(kern.options, 'order')
    kern.order = kern.options.order;    
else
    error('Number of Laguerre approximation coefficients is required.')
end

if isfield(kern, 'options') && isfield(kern.options, 'approx')
    kern.approx = kern.options.approx;    
else
    error('Approximation method is required for this kernel')
end

if ~isfield(kern, 'inputDimension')
    warning('laguerreKernParamInit:noInputDimension', 'Input dimension has not been provided. Assuming is one.')
    kern.inputDimension = 1;
end

if isfield(kern, 'options') && isfield(kern.options, 'isArd') && kern.options.isArd
    kern.isArd = true;
    error('ARD for numerical integral kernels is not implemented yet.')
else
    kern.isArd = false;
    if isfield(kern, 'options') && isfield(kern.options, 'tieOutputParams')
        kern.tieOutputParams = kern.options.tieOutputParams;
        if kern.options.tieOutputParams
            kern.lq = 1. + rand(1, kern.nlf); % Length-scale
            kern.gamma = 0.5*ones(1, kern.nout); % Laguerre scale time parameter
            kern.c = rand(1, kern.nout*(kern.order + 1)); % coefficients
            kern.nParams = kern.nlf + kern.nout*(kern.order + 2);
        else
            kern.lq = 1. + rand(1, kern.nlf); %Length-scale
            kern.gamma = 0.5*ones(1, kern.nout); % scale time parameter
            kern.c = rand(1, kern.nout*(kern.order + 1)); % coefficients
            kern.nParams = kern.nlf + kern.nout*(kern.order + 2);
        end        
    else
        kern.tieOutputParams = true;
        kern.lq = 1. + rand(1, kern.nlf); %Length-scale
        kern.c = ones(1, kern.nout*(kern.order + 1)); % coefficients
        kern.gamma = ones(1, kern.nout); % scale time parameter
        kern.nParams = kern.nlf + kern.nout*(kern.order + 2);
    end
    kern.lfParamsTemplate = 1;
    kern.outParamsTemplate = 1;
end

if isfield(kern, 'options') && isfield(kern.options, 'isVarS') ...
        && kern.options.isVarS
    kern.transforms.index = 1:kern.nlf + kern.nout;
else
    %Only lengthscales and gamma values are positive
    kern.transforms.index = 1:kern.nlf + kern.nout;
end
kern.transforms.type = optimiDefaultConstraint('positive');
kern.isStationary = false;
kern.isNormalised = false;

params = laguerreKernExtractParam(kern);
kern = laguerreKernExpandParam(kern, params);