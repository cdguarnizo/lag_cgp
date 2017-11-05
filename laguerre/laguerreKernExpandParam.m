function kern = laguerreKernExpandParam(kern, params)

% LAGUERREKERNEXPANDPARAM
%
% COPYRIGHT : Cristian Guarnizo, 2015
% MULTIGP
% TODO: add sensitivities

nParamsLat = kern.nlf;
kern.lq = reshape(params(1:nParamsLat), 1, kern.nlf);
if kern.tieOutputParams,
    nParamsOut = kern.nout*(kern.order + 2);
    kern.gamma = reshape(params(nParamsLat+1:nParamsLat+kern.nout), 1, kern.nout);
    kern.c = reshape(params(nParamsLat+kern.nout+1:nParamsLat+nParamsOut), 1, kern.nout*(kern.order+1));
else
    error('Not implemented yet.')
end
% if ~(isfield(kern, 'isVarS') && kern.isVarS),
%     kern.sensitivity = reshape(params(nParamsLat+nParamsOut+1:end), kern.nout, kern.nlf);
% end