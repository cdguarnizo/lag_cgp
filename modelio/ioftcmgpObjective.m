function f = ioftcmgpObjective(params, model)

% FTCMMGPOBJECTIVE Wrapper function for MODELOPTIMISE objective.
% FTCMMGP
model = ioftcmgpExpandParam(model, params);
model = ioftcmgpComputeKernel(model);
f = -ioftcmgpLowerBound(model);