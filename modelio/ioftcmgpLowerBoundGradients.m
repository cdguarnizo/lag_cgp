function gParam = ioftcmgpLowerBoundGradients(model)

% FTCMMLOWERBOUNDGRADIENTS
% FTCMMGP
% gBeta = [];
% if model.includeNoise,
%     gBeta = zeros(1,model.nout);
%     gBetat = zeros(1, model.nout);
%     gKernParam = zeros(1, model.nParams - model.nout);
% else
%     gKernParam = zeros(1, model.nParams);
% end

%gParam = zeros(1, model.nParams);
gParam = repmat(ioftcmgpLowerBound(model),1,model.nParams);
paramt = modelExtractParam(model);
for k = 1:model.nParams,
    p = paramt(k);
    paramt(k) = p - 1e-6;
    model = modelExpandParam(model, paramt);
    model = ioftcmgpComputeKernel(model);
    gParam(k) = (gParam(k) - ioftcmgpLowerBound(model))/1e-6;
    paramt(k) = p;
    model = modelExpandParam(model, paramt);
end

%Param
% Length-scales(nrep), KernParam(nParamsOut*D), Beta(D)

% For each n, input length scale is changed

% dLdKyy = model.alpha*model.alpha' - model.Kyyinv;
% if model.includeNoise,
%     indout = model.outX.index <= model.nout;
%     temp = diag(dLdKyy(indout, indout));
%     for d = 1:model.nout,
%         indout = model.outX.index == d;
%         gBetat(d) = (-1/model.beta(d)^2)*sum(temp(indout));
%     end
%     fhandle = str2func([model.betaTransform 'Transform']);
%     gBeta = gBeta + gBetat.*fhandle(model.beta, 'gradfact');
% end

% When we assume different lq for each realization
%model.kern.lq = model.ls(n);
%gKernParamt = kernGradient(model.kern, {model.outX{:,n}}', dLdKyy);
%gKernParam(n) = gKernParamt(1);
%gKernParam(model.nrep + 1 : end) = gKernParam(model.nrep + 1 : end) + gKernParamt(2:end);

%Here we need to add gradients from Kuu Kfu and Kff
%g = [glq, ggamma, gc]

%Kff gradients
%gKernParam = kernGradient(model.kern, convind(model.outX, model.nout), dLdKyy);

%gParam = .5*[gKernParam gBeta];

%     function g = gKuu()
%         indu = model.outX.index(model.outX.index > model.nout);
%         tempu = length(indu);
%         g = zeros(1, model.nlf);
%         for q = 1:model.nlf,
%             lq = model.kern.lq(q);
%             ind = indu(indu == q + model.nout);
%             
%             model.kern.lq(q) = lq + 1e-6;
%             temp = model.outX.val(model.outX.index == model.nout+q)/model.kern.lq(q);
%             Kuu = exp(-dist2(temp, temp));
%             
%             model.kern.lq(q) = lq - 1e-6;
%             temp = model.outX.val(model.outX.index == model.nout+q)/model.kern.lq(q);
%             Kuu = ((Kuu - exp(-dist2(temp, temp)))/2e-6).*dLdKyy(ind,ind);
%             g(q) = sum(Kuu(:));
%         end
%     end
% 
%     function g = gKfu()
%         
%     end

end