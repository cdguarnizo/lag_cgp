function g = laguerreKernGradient(kern, t, dL_dKff)
% Calculate the gradient of Kff w.r.t kernel hyperarameters:
% legnthscale, Lagguerre coefficients and Laguerre parameter gamma.
% Inputs:
% kern Structure containing information about kernel.
% t Cell containing time input arrays of each output.
% dL_dKff Gradient of likelihood wrt to Kff

%we have g(1:nN) and Kuu(1:nN, 1:nN), also g(t) and Kuu(t, 1:nN)
%gdt, Kuudtdt, gt, Kuutdt

D = length(kern.gamma);
Q = length(kern.lq);
dimen = length(cell2mat(t));
glq = zeros(1,Q);
gc = zeros(1,D*(kern.order+1));
ggamma = zeros(1, D);

% dt = linspace(kern.tspace(1), kern.tspace(2), kern.tsize);
% deltax = dt(2)-dt(1);
% %Calculation of approximated impulse response
% gdt = cell(1,D);
% glidt = cell(1,D);
% gt = cell(1,D);
% glit = cell(1,D);
% for d = 1:D,
%     glidt{d} = laguerreFuncEval(dt(:)', kern.order, kern.gamma(d));
%     gdt{d} = sum(repmat(kern.c(1+(d-1)*(kern.order +1): d*(kern.order +1))', 1,length(dt)).*glidt{d}, 1);
%     %gdt{d} = zeros(1, length(dt));
%     %for k=1:kern.order+1,
%     %    gdt{d} = gdt{d} + kern.c(k+(d-1)*(kern.order+1))*glidt{d}(k,:);
%     %end
%     %gdt{d} = sum(bsxfun(@times, kern.c(1+(d-1)*(kern.order+1):d*(kern.order+1))', glidt{d}), 1);
%     glit{d} = laguerreFuncEval(t{d}(:)', kern.order, kern.gamma(d));
%     gt{d} = sum(repmat(kern.c(1+(d-1)*(kern.order +1): d*(kern.order +1))', 1,length(t{d})).*glit{d}, 1);
%     %gt{d} = zeros(1, length(t{d}));
%     %for k=1:kern.order+1,
%     %    gt{d} = gt{d} + kern.c(k+(d-1)*(kern.order+1))*glit{d}(k,:);
%     %end
%     %gt{d} = sum(bsxfun(@times, kern.c(1+(d-1)*(kern.order+1):d*(kern.order+1))', glit{d}), 1);
% end
% 
% Kff = cell(D,D);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Calculation of gradient wrt lq %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for q = 1:Q,
%     % Create kuudtdt
%     temp = dt(:)/kern.lq(q);
%     dist = dist2(temp, temp);
%     kuudtdt = (2/kern.lq(q))*(dist.*exp(-dist));
%     for d = 1:D,
%         %Find portions
%         %t time inputs for f, deltat grid difference
%         portion = floor(t{d}(1:end)/deltax)+1;
%         % Create kuutdt
%         dist = dist2(t{d}/kern.lq(q), temp);
%         Kuutdt = (2/kern.lq(q))*(dist.*exp(-dist));
%         Kfutdt = zeros(size(Kuutdt));
%         for j = 2:length(t{d}),
%             %Convolution and correction
%             if t{d}(j) > dt(portion(j)),
%                 %Convolution and correction
%                 Kfutdt(j,:) = (deltax*gdt{d}(2:portion(j)))*kuudtdt(portion(j):-1:2,:);
%                 %Kfutdt(j,:) = Kfutdt(j,:) + .5*((t{d}(j) - portion(j)*deltax)*gt{d}(j)*Kuutdt(1,:) + deltax*gt{d}(1)*Kuutdt(j,:));
%                 Kfutdt(j,:) = Kfutdt(j,:) + .5*((t{d}(j) - dt(portion(j)))*gt{d}(j)*Kuutdt(1,:) + deltax*gt{d}(1)*Kuutdt(j,:));
%             else %Then is equal
%                 %Convolution and correction
%                 Kfutdt(j,:) = (deltax*gdt{d}(2:portion(j)-1))*kuudtdt(portion(j)-1:-1:2,:);
%                 %Kfutdt(j,:) = Kfutdt(j,:) + .5*((t{d}(j) - portion(j)*deltax)*gt{d}(j)*Kuutdt(1,:) + deltax*gt{d}(1)*Kuutdt(j,:));
%                 Kfutdt(j,:) = Kfutdt(j,:) + .5*(deltax*(gt{d}(j)*Kuutdt(1,:) + gt{d}(1)*Kuutdt(j,:)));
%             end
%         end
%         for dp = d:D,
%             %t time inputs for f, deltat grid difference
%             portion = floor(t{dp}(1:end)/deltax)+1;
%             %Take portion up to t_j for convolution
%             Kff{d, dp} = zeros(length(t{d}),length(t{dp}));
%             for j = 2:length(t{dp}),
%                 if t{dp}(j) > dt(portion(j)),
%                     %Here the last value Kfu(t_d,t_dp) is approximated
%                     Kfutt = (Kfutdt(:, portion(j)+1) - Kfutdt(:, portion(j)))...
%                         /(dt(portion(j)+1) - dt(portion(j)))*t{dp}(j) +  Kfutdt(:, portion(j));
%                     %Convolution and correction
%                     Kff{d, dp}(:, j) = Kfutdt(:, 2:portion(j))*(deltax*gdt{dp}(portion(j):-1:2)')...
%                         + .5*(deltax*gt{dp}(1)*Kfutt + (t{dp}(j) - dt(portion(j)))*gt{dp}(j)*Kfutdt(:, 1));
%                 else
%                     %Convolution and correction
%                     Kff{d, dp}(:, j) = Kfutdt(:, 2:portion(j)-1)*(deltax*gdt{dp}(portion(j)-1:-1:2)')...
%                         + .5*(deltax*(gt{dp}(1)*Kfutdt(:, portion(j)) + gt{dp}(j)*Kfutdt(:, 1)));
%                 end
%             end
%             if d ~= dp,
%                 %Filling symetric part of Kff
%                 Kff{dp, d} = Kff{d, dp}.';
%             end
%         end
%     end
%     glq(q) = sum(sum(dL_dKff.*cell2mat(Kff)));
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Calculation of gradient wrt c_i %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% index = zeros(2, D);
% start = 1;
% for d = 1:D,
%     fin = start + length(t{d});
%     index(:, d) = [start; fin-1];
%     start = fin;
% end
% for i = 1:kern.order+1,
%     Kfft = zeros(dimen,dimen);
%     for q = 1:Q,
%         % Create kuudtdt
%         temp = dt(:)/kern.lq(q);
%         kuudtdt = exp(-dist2(temp, temp));
%         for d = 1:D,
%             %Find portions
%             %t time inputs for f, deltat grid difference
%             portion = floor(t{d}(1:end)/deltax)+1;
%             % Create kuutdt
%             Kuutdt = exp(-dist2(t{d}/kern.lq(q), temp));
%             Kfutdt = zeros(size(Kuutdt));
%             for j = 2:length(t{d}),
%                 if t{d}(j) > dt(portion(j)),
%                     %Convolution and correction
%                     Kfutdt(j,:) = (deltax*gdt{d}(2:portion(j)))*kuudtdt(portion(j):-1:2,:);
%                     %Kfutdt(j,:) = Kfutdt(j,:) + .5*((t{d}(j) - portion(j)*deltax)*gt{d}(j)*Kuutdt(1,:) + deltax*gt{d}(1)*Kuutdt(j,:));
%                     Kfutdt(j,:) = Kfutdt(j,:) + .5*((t{d}(j) - dt(portion(j)))*gt{d}(j)*Kuutdt(1,:) + deltax*gt{d}(1)*Kuutdt(j,:));
%                 else %Then is equal
%                     %Convolution and correction
%                     Kfutdt(j,:) = (deltax*gdt{d}(2:portion(j)-1))*kuudtdt(portion(j)-1:-1:2,:);
%                     %Kfutdt(j,:) = Kfutdt(j,:) + .5*((t{d}(j) - portion(j)*deltax)*gt{d}(j)*Kuutdt(1,:) + deltax*gt{d}(1)*Kuutdt(j,:));
%                     Kfutdt(j,:) = Kfutdt(j,:) + .5*(deltax*(gt{d}(j)*Kuutdt(1,:) + gt{d}(1)*Kuutdt(j,:)));
%                 end
%             end
%             for dp = 1:D,
%                 %t time inputs for f, deltat grid difference
%                 portion = floor(t{dp}(1:end)/deltax)+1;
%                 %Take portion up to t_j for convolution
%                 Kff{d, dp} = zeros(length(t{d}),length(t{dp}));
%                 for j = 2:length(t{dp}),
%                     if t{dp}(j) > dt(portion(j)),
%                         %Here the last value Kfu(t_d,t_dp) is approximated
%                         Kfutt = (Kfutdt(:, portion(j)+1) - Kfutdt(:, portion(j)))...
%                             /(dt(portion(j)+1) - dt(portion(j)))*t{dp}(j) +  Kfutdt(:, portion(j));
%                         %Convolution and correction
%                         Kff{d, dp}(:, j) = Kfutdt(:, 2:portion(j))*(deltax*glidt{dp}(i, portion(j):-1:2)')...
%                             + .5*(deltax*glit{dp}(i,1)*Kfutt + (t{dp}(j) - dt(portion(j)))*glit{dp}(i,j)*Kfutdt(:, 1));
%                     else
%                         %Convolution and correction
%                         Kff{d, dp}(:, j) = Kfutdt(:, 2:portion(j)-1)*(deltax*glidt{dp}(i, portion(j)-1:-1:2)')...
%                             + .5*(deltax*(glit{dp}(i,1)*Kfutdt(:, portion(j)) + glit{dp}(i,j)*Kfutdt(:, 1)));
%                     end
%                 end
%             end
%         end
%         Kfft = Kfft + cell2mat(Kff);
%     end
%     for d = 1:D,
%         gc(i + (d-1)*(kern.order+1)) = 2*sum(sum(dL_dKff(:,index(1,d):index(2,d)).*Kfft(:,index(1,d):index(2,d))));
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Calculation of gradient wrt gamma_d %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ggdt = cell(1, D);
% ggt = cell(1, D);
% for d = 1:D,
%     ggdt{d} = sum(repmat(kern.c(1+(d-1)*(kern.order+1):d*(kern.order+1))'...
%         ,1,length(dt)).*laguerreFuncGradientEval(dt(:)', kern.order, kern.gamma(d)), 1);
%     ggt{d} = sum(repmat(kern.c(1+(d-1)*(kern.order+1):d*(kern.order+1))'...
%         ,1,length(t{d})).*laguerreFuncGradientEval(t{d}(:)', kern.order, kern.gamma(d)), 1);
% end
% Kfft = zeros(dimen,dimen);
% 
% for q = 1:Q,
%     % Create kuudtdt
%     temp = dt(:)/kern.lq(q);
%     kuudtdt = exp(-dist2(temp, temp));
%     for d = 1:D,
%         %Find portions
%         %t time inputs for f, deltat grid difference
%         portion = floor(t{d}(1:end)/deltax)+1;
%         % Create kuutdt
%         Kuutdt = exp(-dist2(t{d}/kern.lq(q), temp));
%         Kfutdt = zeros(size(Kuutdt));
%         for j = 2:length(t{d}),
%             if t{d}(j) > dt(portion(j)),
%                 %Convolution and correction
%                 Kfutdt(j,:) = (deltax*gdt{d}(2:portion(j)))*kuudtdt(portion(j):-1:2,:);
%                 %Kfutdt(j,:) = Kfutdt(j,:) + .5*((t{d}(j) - portion(j)*deltax)*gt{d}(j)*Kuutdt(1,:) + deltax*gt{d}(1)*Kuutdt(j,:));
%                 Kfutdt(j,:) = Kfutdt(j,:) + .5*((t{d}(j) - dt(portion(j)))*gt{d}(j)*Kuutdt(1,:) + deltax*gt{d}(1)*Kuutdt(j,:));
%             else %Then is equal
%                 %Convolution and correction
%                 Kfutdt(j,:) = (deltax*gdt{d}(2:portion(j)-1))*kuudtdt(portion(j)-1:-1:2,:);
%                 %Kfutdt(j,:) = Kfutdt(j,:) + .5*((t{d}(j) - portion(j)*deltax)*gt{d}(j)*Kuutdt(1,:) + deltax*gt{d}(1)*Kuutdt(j,:));
%                 Kfutdt(j,:) = Kfutdt(j,:) + .5*(deltax*(gt{d}(j)*Kuutdt(1,:) + gt{d}(1)*Kuutdt(j,:)));
%             end
%         end
%         for dp = 1:D,
%             %t time inputs for f, deltat grid difference
%             portion = floor(t{dp}(1:end)/deltax)+1;
%             %Take portion up to t_j for convolution
%             Kff{d, dp} = zeros(length(t{d}),length(t{dp}));
%             for j = 2:length(t{dp}),
% %                 %Here the last value Kfu(t_d,t_dp) is approximated
% %                 Kfutt = (Kfutdt(:, portion(j)+1) - Kfutdt(:, portion(j)))...
% %                     /(dt(portion(j)+1) - dt(portion(j)))*t{dp}(j) +  Kfutdt(:, portion(j));
% %                 %Convolution and correction
% %                 Kff{d, dp}(:, j) = Kfutdt(:, 2:portion(j))*(deltax*ggdt{dp}(portion(j):-1:2)')...
% %                     + .5*(deltax*ggt{dp}(1)*Kfutt + (t{dp}(j) - portion(j)*deltax)*ggt{dp}(j)*Kfutdt(:, 1));
%                 if t{dp}(j) > dt(portion(j)),
%                     %Here the last value Kfu(t_d,t_dp) is approximated
%                     Kfutt = (Kfutdt(:, portion(j)+1) - Kfutdt(:, portion(j)))...
%                         /(dt(portion(j)+1) - dt(portion(j)))*t{dp}(j) +  Kfutdt(:, portion(j));
%                     %Convolution and correction
%                     Kff{d, dp}(:, j) = Kfutdt(:, 2:portion(j))*(deltax*ggdt{dp}(portion(j):-1:2)')...
%                         + .5*(deltax*ggt{dp}(1)*Kfutt + (t{dp}(j) - dt(portion(j)))*ggt{dp}(j)*Kfutdt(:, 1));
%                 else
%                     %Convolution and correction
%                     Kff{d, dp}(:, j) = Kfutdt(:, 2:portion(j)-1)*(deltax*ggdt{dp}(portion(j)-1:-1:2)')...
%                         + .5*(deltax*(ggt{dp}(1)*Kfutdt(:, portion(j)) + ggt{dp}(j)*Kfutdt(:, 1)));
%                 end
%             end
%         end
%     end
%     Kfft = Kfft + cell2mat(Kff);
% end
% % Sample from the latent forces
% for d = 1:D,
%     ggamma(d) = 2*sum(sum(dL_dKff(:,index(1,d):index(2,d)).*Kfft(:,index(1,d):index(2,d))));
% end

%Lengthscales
for q=1:Q,
    lq = kern.lq(q);
    kern.lq(q) = lq + 1e-6;
    K=laguerreKernCompute(kern, t);
    kern.lq(q) = lq - 1e-6;
    K = (K - laguerreKernCompute(kern, t))/2e-6;
    glq(q) = sum(sum(dL_dKff.*K));
    kern.lq(q) = lq;
end

index = zeros(2, D);
start = 1;
for d = 1:D,
    fin = start + length(t{d});
    index(:, d) = [start; fin-1];
    start = fin;
end

%Gamma
for d = 1:D,
    gam = kern.gamma(d);
    kern.gamma(d) = gam + 1e-6;
    K=laguerreKernCompute(kern, t);
    kern.gamma(d) = gam - 1e-6;
    K = (K - laguerreKernCompute(kern, t))/2e-6;
    ggamma(d) = sum(sum(dL_dKff.*K));
    kern.gamma(d) = gam;
end

%Coefficients
for i = 1:kern.order+1,
    for d = 1:D,
        c = kern.c(i+(d-1)*(kern.order+1));
        kern.c(i+(d-1)*(kern.order+1)) = c + 1e-6;
        K = laguerreKernCompute(kern, t);
        kern.c(i+(d-1)*(kern.order+1)) = c - 1e-6;
        K = (K - laguerreKernCompute(kern, t))/2e-6;
        gc(i + (d-1)*(kern.order+1)) = sum(sum(dL_dKff.*K));
        kern.c(i+(d-1)*(kern.order+1)) = c;
    end
end


g = [glq, ggamma, gc];