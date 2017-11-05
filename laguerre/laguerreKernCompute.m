function Kfft = laguerreKernCompute(kern, t, t2)
% Calculate Kff using Squared Exponential kernel as covariance function for
% Kuu, and calculate the impulse response filter values approximated by
% laguerre orthogonal functions.
% Inputs:
% t Cell containing time input arrays of each output

%we have g(1:nN) and Kuu(1:nN, 1:nN), also g(t) and Kuu(t, 1:nN)
%gdt, Kuudtdt, gt, Kuutdt

dimen = length(cell2mat(t));
flag_same_t = false;
if nargin < 3,
    flag_same_t = true;
    dimen2 = dimen;
else
    dimen2 = length(cell2mat(t2));
end

D = length(kern.gamma);
Q = length(kern.lq);
%Calculation of approximated impulse response
gdt = cell(1,D);
dt = cell(1,D);
inddt = cell(1,D);
gdtdel = cell(1,D);
for d = 1:D,
    [dt{d}, inddt{d}] = gridmaker(t{d}, kern.tsize);
    gdt{d} = sum(repmat(kern.c(1+(d-1)*(kern.order +1): d*(kern.order +1))', 1,length(dt{d})).*...
        laguerreFuncEval(dt{d}(:)', kern.order, kern.gamma(d)), 1);
    gdtdel{d} = diff(dt{d}');
    gdtdel{d} = .5*sum([0., gdtdel{d}; gdtdel{d}, 0.]);
    gdtdel{d} = gdtdel{d}.*gdt{d};
end

Kff = cell(D,D);
Kfft = zeros(dimen, dimen2);
if flag_same_t,
    for q = 1:Q,
        for d = 1:D,
            tempd = dt{d}(:)/kern.lq(q);
            for dp = d:D,
                % Create kuudtdt
                if dp ~= d,
                    tempdp = dt{dp}/kern.lq(q);
                    kuudtdt = exp(-dist2(tempd, tempdp));
                else
                    % Create kuudtdt
                    kuudtdt = exp(-dist2(tempd, tempd));
                end
                Kfutdt = zeros(length(t{d}),size(kuudtdt,2));
                for j = 2:length(t{d}),
                    %Convolution evaluation
                    Kfutdt(j,:) = gdtdel{d}(1:inddt{d}(j)-1)*kuudtdt(inddt{d}(j):-1:2,:) + ...
                        (.5*(dt{d}(inddt{d}(j))-dt{d}(inddt{d}(j)-1))*gdt{d}(inddt{d}(j)))*kuudtdt(1,:);
                end
                Kff{d, dp} = zeros(length(t{d}),length(t{dp}));
                
                for j = 2:length(t{dp}),
                    %Convolution evaluation
                    Kff{d, dp}(:, j) = Kfutdt(:,inddt{dp}(j):-1:2)*gdtdel{dp}(1:inddt{dp}(j)-1)' + ...
                        (.5*(dt{dp}(inddt{dp}(j))-dt{dp}(inddt{dp}(j)-1))*gdt{dp}(inddt{dp}(j)))*Kfutdt(:,1);
                end
                if d ~= dp,
                    %Filling symetric part of Kff
                    Kff{dp, d} = Kff{d, dp}.';
                end
            end
        end
        %Adding Q matrices
        Kfft = Kfft + cell2mat(Kff);
    end
else
    gdt2 = cell(1,D);
    dt2 = cell(1,D);
    inddt2 = cell(1,D);
    gdtdel2 = cell(1,D);
    for d = 1:D,
        [dt2{d}, inddt2{d}] = gridmaker(t2{d}, kern.tsize);
        gdt2{d} = sum(repmat(kern.c(1+(d-1)*(kern.order +1): d*(kern.order +1))', 1,length(dt2{d})).*laguerreFuncEval(dt2{d}(:)', kern.order, kern.gamma(d)), 1); %TODO: change name of Volterra Eval
        gdtdel2{d} = diff(dt2{d}');
        gdtdel2{d} = .5*sum([0., gdtdel2{d}; gdtdel2{d}, 0.]);
        gdtdel2{d} = gdtdel2{d}.*gdt2{d};
    end
    for q = 1:Q,
        for d = 1:D,
            tempd = dt{d}/kern.lq(q);
            for dp = 1:D,
                tempdp = dt2{dp}/kern.lq(q);
                kuudtdt = exp(-dist2(tempd, tempdp));
                
                Kfutdt = zeros(length(t{d}),size(kuudtdt,2));
                for j = 2:length(t{d}),
                    %Convolution evaluation
                    Kfutdt(j,:) = gdtdel{d}(1:inddt{d}(j)-1)*kuudtdt(inddt{d}(j):-1:2,:) + ...
                        (.5*(dt{d}(inddt{d}(j))-dt{d}(inddt{d}(j)-1))*gdt{d}(inddt{d}(j)))*kuudtdt(1,:);
                end
                Kff{d, dp} = zeros(length(t{d}),length(t2{dp}));
                for j = 2:length(t2{dp}),
                    %Convolution evaluation
                    Kff{d, dp}(:, j) = Kfutdt(:,inddt2{dp}(j):-1:2)*gdtdel2{dp}(1:inddt2{dp}(j)-1)' + ...
                        (.5*(dt2{dp}(inddt2{dp}(j))-dt2{dp}(inddt2{dp}(j)-1))*gdt2{dp}(inddt2{dp}(j)))*Kfutdt(:,1);
                end
            end
        end
        %Adding Q matrices
        Kfft = Kfft + cell2mat(Kff);
    end
end