function [ymean, ySd] = ioftcmgpPosterior(model, Xtest)

% FTCMMGPPOSTERIOR
% FTCMMGP


%             y   u  
%       ys |Kysy Kfsu|
% KDsD =   |         |
%       us |Kusf Kusu|

if nargout > 1,
    flag = true;
else
    flag = false;
end

%%%%%%%%%%%%%%%%%%%%%%
% Compute cov Ku_s u %
%%%%%%%%%%%%%%%%%%%%%%
indus = find(Xtest.index > model.nout);
if isempty(indus),
    Kusu = [];
else
    Xus = Xtest.val(indus);
    indus = Xtest.index(indus) - model.nout;
    uniq = unique(indus)';

    indu = model.outX.index > model.nout;
    Xu = model.outX.val(indu);
    indu = model.outX.index(indu) - model.nout;
    
    Kusu = zeros(length(indus), length(indu));
    for q = uniq,
        inds = indus == q;
        ind = indu == q;
        
        Kusu(inds, ind) = exp( -dist2( Xus(inds)/model.kern.lq(q),Xu(ind)/model.kern.lq(q) ) );
    end
end

%%%%%%%%%%%%%%%%%%%%
% Compute cov Kusf %
%%%%%%%%%%%%%%%%%%%%
if isempty(indus),
    Kusf = [];
else
    indf = model.outX.index(model.outX.index <= model.nout);
    
    Kusf = zeros(length(indus), length(indf));
    for q = uniq,
        indq = indus == q;
        tempu = Xus(indq)/model.kern.lq(q);
        for d = 1:model.nout,
            ind = find(indf == d);
            [dt, inddt] = gridmaker(model.outX.val(ind), model.kern.tsize);
            gdt = sum(repmat(model.kern.c(1+(d-1)*(model.kern.order +1): d*(model.kern.order +1))', 1,length(dt)).*...
                laguerreFuncEval(dt(:)', model.kern.order, model.kern.gamma(d)), 1);
            gdtdel = diff(dt');
            gdtdel = .5*sum([0., gdtdel; gdtdel, 0.]);
            gdtdel = gdtdel.*gdt;
            
            kuudtt = exp(-dist2(dt/(model.kern.lq(q)), tempu));
            Kfut = zeros(length(ind),length(tempu));
            for j = 2:length(inddt),
                %Convolution evaluation
                Kfut(j,:) = gdtdel(1:inddt(j)-1)*kuudtt(inddt(j):-1:2,:) + ...
                    (.5*(dt(inddt(j))-dt(inddt(j)-1))*gdt(inddt(j)))*kuudtt(1,:);
            end
            Kusf(indq, ind) = Kfut';
        end
        
    end
end

indfs = find(Xtest.index <= model.nout);
%%%%%%%%%%%%%%%%%%%%
% Compute cov Kfsu %
%%%%%%%%%%%%%%%%%%%%
if isempty(indfs),
    Kfsu = [];
else
    Xfs = Xtest.val(indfs);
    unif = unique(Xtest.index(indfs))';
    indfs = Xtest.index(indfs);
    Kfsu = zeros(length(indfs),length(indu));
    for q = 1:model.nlf,
        indq = indu == q;
        tempu = Xu(indq)/model.kern.lq(q);
        for d = unif,
            ind = find(indfs == d);
            [dt, inddt] = gridmaker(Xfs(ind), model.kern.tsize);
            gdt = sum(repmat(model.kern.c(1+(d-1)*(model.kern.order +1): d*(model.kern.order +1))', 1,length(dt)).*...
                laguerreFuncEval(dt(:)', model.kern.order, model.kern.gamma(d)), 1);
            gdtdel = diff(dt');
            gdtdel = .5*sum([0., gdtdel; gdtdel, 0.]);
            gdtdel = gdtdel.*gdt;
            
            kuudtt = exp(-dist2(dt/(model.kern.lq(q)), tempu));
            Kfut = zeros(length(ind),length(tempu));
            for j = 2:length(inddt),
                %Convolution evaluation
                Kfut(j,:) = gdtdel(1:inddt(j)-1)*kuudtt(inddt(j):-1:2,:) + ...
                    (.5*(dt(inddt(j))-dt(inddt(j)-1))*gdt(inddt(j)))*kuudtt(1,:);
            end
            Kfsu(ind,indq) = Kfut;
        end
    end
end

%%%%%%%%%%%%%%%%%%%
% Compute cov Kyy %
%%%%%%%%%%%%%%%%%%%
fhandle = str2func([model.kernType 'KernCompute']);
Kfsf = fhandle(model.kern, convind(Xtest, model.nout), convind(model.outX, model.nout));

Kgsg = [Kfsf,Kfsu; Kusf,Kusu];

ymean = Kgsg*model.alpha;
if flag,
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Compute cov Ku_s u_s %
    %%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(indus),
        Kusus = [];
    else
        Kusus = zeros(length(indus), length(indus));
        for q = uniq,
            inds = indus == q;
            a = Xus(inds)/model.kern.lq(q);
            Kusus(inds, inds) = exp( -dist2(a,a) );
        end
%         [row, col] = size(Kusus);
%         diagind = 1:row+1:row*col;
%         Kusus(diagind) = Kusus(diagind) + 1/model.beta(model.nout+1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    % Compute cov Kfsus %
    %%%%%%%%%%%%%%%%%%%%%
    if isempty(indfs) && isempty(indus),
        Kfsus = [];
    else
        Kfsus = zeros(length(indfs),length(indu));
        for q = uniq,
            indq = indus == q;
            tempu = Xus(indq)/model.kern.lq(q);
            for d = unif,
                ind = find(indfs == d);
                [dt, inddt] = gridmaker(Xfs(ind), model.kern.tsize);
                gdt = sum(repmat(model.kern.c(1+(d-1)*(model.kern.order +1): d*(model.kern.order +1))', 1,length(dt)).*...
                    laguerreFuncEval(dt(:)', model.kern.order, model.kern.gamma(d)), 1);
                gdtdel = diff(dt');
                gdtdel = .5*sum([0., gdtdel; gdtdel, 0.]);
                gdtdel = gdtdel.*gdt;
                
                kuudtt = exp(-dist2(dt/(model.kern.lq(q)), tempu));
                Kfut = zeros(length(ind),length(tempu));
                for j = 2:length(inddt),
                    %Convolution evaluation
                    Kfut(j,:) = gdtdel(1:inddt(j)-1)*kuudtt(inddt(j):-1:2,:) + ...
                        (.5*(dt(inddt(j))-dt(inddt(j)-1))*gdt(inddt(j)))*kuudtt(1,:);
                end
                Kfsus(ind,indq) = Kfut;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    % Compute cov Kfsfs %
    %%%%%%%%%%%%%%%%%%%%%
    fhandle = str2func([model.kernType 'KernCompute']);
    Kfsfs = fhandle(model.kern, convind(Xtest, model.nout));

    noise = 1./model.beta(:);
    
    [row, col] = size(Kfsfs);
    diagind = 1:row+1:row*col;
    Kfsfs(diagind) = Kfsfs(diagind) + noise(indfs)';
    
    Kgsgs = [Kfsfs,Kfsus; Kfsus',Kusus];
    
    ySd = real(sqrt(diag(Kgsgs) - diag(Kgsg*(model.Kyyinv*Kgsg'))));
%     for d = 1:model.nout,
%         ind = find(indfs == d);
%         ySd(ind) = ySd(ind) + sqrt(noise(d));
%     end
%     indus = find(Xtest.index > model.nout);
%     ySd(indus) = ySd(indus) + sqrt(noise(model.nout+1));
end