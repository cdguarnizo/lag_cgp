function l = lagurreFuncGradientEval(t, order, gamma)

gammat = gamma*t;
const = sqrt(2*gamma)*exp(-gammat);
%Calculate each L_i
L = zeros(order+1, length(t));
%Calcuate each Laguerre function
l = zeros(order+1, length(t));
power2gammat = zeros(order, length(t));
L(1, :) = 1/(2*gamma) - t; %L(0)
l(1, :) = const.*L(1, :);
for i = 1:order,
    power2gammat(i, :) = gammat.^i;
    L(i+1, :) = L(1, :)*(-1)^i;
    for k = 0:i-1,
        coef = ((1+2*(i-k))/(2*gamma) - t);
        L(i+1, :) = L(i+1, :) + ((-1)^k*factorial(i)*2^(i-k)/(factorial(k)*factorial(i-k)^2))*(coef.*power2gammat(i-k, :));
    end
    l(i+1, :) = L(i+1, :).*const;
end