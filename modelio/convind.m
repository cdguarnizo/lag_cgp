function Xcell = convind(X, D)
% This function allows to convert values index scheme to cell
if nargin < 2,
    D = max(X.index);
end

Xcell = cell(D,1);
for d = 1:D,
    Xcell{d} = X.val(X.index == d);
end