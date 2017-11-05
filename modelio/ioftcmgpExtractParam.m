function [param, names] = ioftcmgpExtractParam(model)

% FTCMMGPEXTRACTPARAM Extract the parameters of an FTCMMGP struct.
% FTCMMGP

% Extract Length-scales values for multi model
% if isfield(model, 'ls') && ~isempty(model.ls),
%     fhandle = str2func([model.lsTransform 'Transform']);
%     lsParams = fhandle(model.ls, 'xtoa');    
%     if nargout>1,
%         lsParamNames = cell(model.nout,1);
%         for i = 1:length(lsParams)
%             lsParamNames{i} = ['LengthScale ' num2str(i)];
%         end
%     end
% else
%     lsParamNames = {};
%     lsParams =[];
% end

% Extract kernel parameters
if nargout > 1,
    [paramKern, namesKern] = kernExtractParam(model.kern);
    %paramKern = paramKern(2:end);
    %namesKern = namesKern{2:end};
else
    paramKern = kernExtractParam(model.kern);
    %paramKern = paramKern(2:end);
end

% Extract beta values (noise level)
if isfield(model, 'beta') && ~isempty(model.beta)
    fhandle = str2func([model.betaTransform 'Transform']);
    betaParams = fhandle(model.beta, 'xtoa');    
    if nargout>1
        betaParamNames = cell(model.nout,1);
        for i = 1:length(betaParams)
            betaParamNames{i} = ['Beta ' num2str(i)];
        end
    end
else
    betaParamNames = {};
    betaParams =[];
end

% Set parameters values in a vector
%param = [lsParams paramKern betaParams];
param = [paramKern betaParams];

% Fix the value of the parameters
if isfield(model, 'fix'),
    for i = 1:length(model.fix),
        param(model.fix(i).index) = model.fix(i).value;
    end
end

if nargout > 1,
    %names = {lsParamNames{:} namesKern{:}, betaParamNames{:}};
    names = {namesKern{:}, betaParamNames{:}};
end