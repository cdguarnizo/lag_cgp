function [f,g] = ftcmmgpObjectiveGradient(params, model)

model = ftcmmgpExpandParam(model, params);
model = ftcmmgpComputeKernel(model);
f = -ftcmmgpLowerBound(model);
g = -ftcmmgpLowerBoundGradients(model);