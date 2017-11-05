function [tt, indt] = gridmaker(t, nsize)

dt = linspace(0, max(t), nsize)';
indt = logical([zeros(1, nsize),ones(1, length(t))]);
[tt, ind] = sort([dt; t]);
indt = indt(ind);
subt = find(indt);
bool = tt(indt) == tt(subt-1);
subt = subt(bool)-1;
tt(subt) = [];
indt(subt) = [];
indt = find(indt);