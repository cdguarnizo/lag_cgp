function model = ioftcmgpExpandParam(model, params)

% IBPMMGPEXPANDPARAM Expand the parameters into an FTCMMGP struct.
% IBPMMGP

paramPart = real(params);

if isfield(model, 'fix'),
    for i = 1:length(model.fix),
       paramPart(model.fix(i).index) = model.fix(i).value;
    end
end

% % Check if there is a Length-scale
% if isfield(model, 'ls') && ~isempty(model.ls),
%     startVal = 1;
%     endVal = model.nlf; %model.nrep;
%     fhandle = str2func([model.lsTransform 'Transform']);
%     model.ls = fhandle(paramPart(startVal:endVal), 'atox');
% end

startVal = 1; %endVal;
endVal = startVal + model.kern.nParams - 1;
kernParams = paramPart(startVal:endVal);
if length(kernParams) ~= model.kern.nParams,
    error('kern parameter vector is incorrect length');
end
model.kern = kernExpandParam(model.kern, kernParams);

% Check if there is a beta parameter.
if isfield(model, 'beta') && ~isempty(model.beta)
    startVal = endVal + 1;
    endVal = endVal + model.nout;
    fhandle = str2func([model.betaTransform 'Transform']);
    model.beta = fhandle(paramPart(startVal:endVal), 'atox');
end