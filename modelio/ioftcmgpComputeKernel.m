function model = ioftcmgpComputeKernel(model)

% FTCMMGPCOMPUTEKERNELS
% FTCMMGP

% y(t) = f(t) + e(t)
%       y |Kyy Kfu|
% KDD =   |       |
%       u |Kuf Kuu|

%TODO: change nlf for something else, there are no latent functions in this
%model

%%%%%%%%%%%%%%%%%%%
% Compute cov Kuu %
%%%%%%%%%%%%%%%%%%%
%Assign Block Diagonal Kuu
indu = model.outX.index(model.outX.index > model.nout) - model.nout;
tempu = length(indu);
Kuu = zeros(tempu, tempu);
for q = 1:model.nlf,
    ind = indu == q;
    temp = model.outX.val(model.outX.index == model.nout+q)/(model.kern.lq(q));
    Kuu(ind, ind) = exp(-dist2(temp, temp));
end
[row, col] = size(Kuu);
diagind = 1:row+1:row*col;
Kuu(diagind) = Kuu(diagind) + 1e-6;
% Kuu(diagind) = Kuu(diagind) + 1/model.beta(model.nout+1);

%%%%%%%%%%%%%%%%%%%
% Compute cov Kfu %
%%%%%%%%%%%%%%%%%%%
Kfu = cell(model.nout, model.nlf);
for q = 1:model.nlf,
    ind2 = model.outX.index == q + model.nout;
    tempu = model.outX.val(ind2)/(model.kern.lq(q));
    for d = 1:model.nout,
        ind = find(model.outX.index == d);
        [dt, inddt] = gridmaker(model.outX.val(ind), model.kern.tsize);
        gdt = sum(repmat(model.kern.c(1+(d-1)*(model.kern.order +1): d*(model.kern.order +1))', 1,length(dt)).*...
            laguerreFuncEval(dt(:)', model.kern.order, model.kern.gamma(d)), 1);
        gdtdel = diff(dt');
        gdtdel = .5*sum([0., gdtdel; gdtdel, 0.]);
        gdtdel = gdtdel.*gdt;

        kuudtt = exp(-dist2(dt/(model.kern.lq(q)), tempu));
        Kfu{d,q} = zeros(length(ind),length(tempu));
        for j = 2:length(inddt),
            %Convolution evaluation
            Kfu{d,q}(j,:) = gdtdel(1:inddt(j)-1)*kuudtt(inddt(j):-1:2,:) + ...
                (.5*(dt(inddt(j))-dt(inddt(j)-1))*gdt(inddt(j)))*kuudtt(1,:);
        end
    end
end
Kfu = cell2mat(Kfu);
%%%%%%%%%%%%%%%%%%%
% Compute cov Kyy %
%%%%%%%%%%%%%%%%%%%
fhandle = str2func([model.kernType 'KernCompute']);
Kyy = fhandle(model.kern, convind(model.outX, model.nout));
Kyy = .5*(Kyy + Kyy');
if model.includeNoise,
    % Generating noise vector
    noise = 1./model.beta;
    % Adding noise
    [row, col] = size(Kyy);
    diagind = 1:row+1:row*col;
    Kyy(diagind) = Kyy(diagind) + noise(model.outX.index(model.outX.index <= model.nout)');
end

%%%%%%%%%%%%%%%%%%%
% Compute cov Kgg %
%%%%%%%%%%%%%%%%%%%
Kgg = [Kyy Kfu; Kfu' Kuu];

% figure(1);
% ind = model.outX.index == 1;
% hres = sum(repmat(model.kern.c(1: model.kern.order +1)', 1,sum(ind)).*laguerreFuncEval(model.outX.val(ind)', model.kern.order, model.kern.gamma(1)), 1);
% plot(model.outX.val(ind), hres, 'r')

% Update alpha
L = jitChol(Kgg);
model.logDet = 2.*sum(log(diag(L)));
Linv = L\eye(size(L));
model.Kyyinv = Linv*Linv.'; %pdinv(KDD);
model.alpha =  model.Kyyinv*model.y;