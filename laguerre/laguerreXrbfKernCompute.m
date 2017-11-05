function Kfu = laguerreXrbfKernCompute(kern, t, t2)
% Calculate Kfu using Squared Exponential kernel as covariance function for
% Kuu, and calculate the impulse response filter values approximated by
% laguerre orthogonal functions.
% Inputs:
% t  Cell containing time input arrays of each output
% t2 Cell containing time input arrays of each latent functions

%we have g(1:nN) and Kuu(1:nN, 1:nN), also g(t) and Kuu(t, 1:nN)
%gdt, Kuudtdt, gt, Kuutdt
dimen = length(cell2mat(t));
flag_same_t =  false;
if nargin < 3,
    t2 = t;
    flag_same_t =  true;
    dimen2 = dimen;
else
    dimen2 = length(cell2mat(t2));
end

D = length(kern.gamma);
Q = length(kern.lq);
dt = linspace(kern.tspace(1), kern.tspace(2), kern.tsize);
deltax = dt(2)-dt(1);
%Calculation of approximated impulse response
gdt = cell(1,D);
gt = cell(1,D);
for d = 1:D,
    gdt{d} = sum(repmat(kern.c(1+(d-1)*(kern.order +1): d*(kern.order +1))', 1,length(dt)).*laguerreFuncEval(dt(:)', kern.order, kern.gamma(d)), 1); %TODO: change name of Volterra Eval
    gt{d} = sum(repmat(kern.c(1+(d-1)*(kern.order +1): d*(kern.order +1))', 1,length(t{d})).*laguerreFuncEval(t{d}(:)', kern.order, kern.gamma(d)), 1);
end

Kfu = cell(D,Q);
for q = 1:Q,
    % Create kuudtdt
    temp = dt(:)/kern.lq(q);
    t2temp = t2{q}/kern.lq(q);
    kuudtt2 = exp(-dist2(temp, t2{q}/(2.*kern.lq(q)) ));
    for d = 1:D,
        %Find portions
        %t time inputs for f, deltat grid difference
        portion = floor(t{d}(1:end)/deltax);
        % Create kuutdt
        Kuutt2 = exp(-dist2(t{d}/(2.*kern.lq(q)), t2temp));
        Kfutt2 = zeros(size(Kuutt2));
        for j = 2:length(t{d}),
            %Convolution and correction
            Kfutt2(j,:) = (deltax*gdt{d}(2:portion(j)))*kuudtt2(portion(j):-1:2,:);
            Kfutt2(j,:) = Kfutt2(j,:) + .5*((t{d}(j) - portion(j)*deltax)*gt{d}(j)*Kuutt2(1,:) + deltax*gt{d}(1)*Kuutt2(j,:));
        end
        %Adding Q matrices
        Kfu{d,q} = Kfutt2;
    end
end