function g = ioftcmgpGradient(params, model)

% FTCMMGPGRADIENT 
% FTCMMGP

model = modelExpandParam(model, params);
model = ioftcmgpComputeKernel(model);
g = -ioftcmgpLowerBoundGradients(model);