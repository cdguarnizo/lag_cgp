function [params, names] = laguerreKernExtractParam(kern)

% LAGUERREKERNEXTRACTPARAM
%
% COPYRIGHT : Cristian Guarnizo, 2015
% MULTIGP
% TODO: add sensitivities variables

if isfield(kern, 'isVarS') && kern.isVarS,
    params = [kern.lq(:)' kern.gamma(:)' kern.c(:)'];
    if nargout > 1,
        lsLatNames = cell(1, kern.nlf); % lenght scales
        for i=1:kern.nlf,
            force = num2str(i);
            lsLatNames{i} = ['length scale latent: force ' force '.'];
        end
        OutNames = cell(1, kern.nout*(kern.order+1));
        for i=1:kern.nout,
            output = num2str(i);
            OutNames{i} = ['gamma output: output ' output '.'];
            for j = 1:kern.order+1,
                coeff = num2str(j);
                OutNames{kern.nout + (j-1)*(kern.order+1)} = ['coeff ', coeff ,' output: output ' output '.'];
            end
        end
        
    end
    names = [lsLatNames(:)' OutNames(:)'];
else
    params = [kern.lq(:)' kern.gamma(:)' kern.c(:)'];
    if nargout > 1,
        %sensitivityNames = cell(kern.nout, kern.nlf);
        lsLatNames = cell(1, kern.nlf);
        for i=1:kern.nlf,
            force = num2str(i);
            lsLatNames{i} = ['length scale latent: force ' force '.'];
        end

        OutNames = cell(1, kern.nout*(kern.order+1));
        for i=1:kern.nout,
            output = num2str(i);
            OutNames{i} = ['gamma output: output ' output '.'];
            for j = 1:kern.order+1,
                coeff = num2str(j);
                OutNames{kern.nout + (j-1)*(kern.order+1)} = ['coeff ', coeff ,' output: output ' output '.'];
            end
        end
        names = [lsLatNames(:)' OutNames(:)'];
    end
end