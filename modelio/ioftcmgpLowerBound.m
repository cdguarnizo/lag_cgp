function logp = ioftcmgpLowerBound(model)
% IOFTCMGPLOWERBOUND
% IOFTCMGP

log2pi = log(2*pi);
logp = - 0.5*(model.y'*model.alpha + length(model.outX.val)*log2pi + model.logDet);
