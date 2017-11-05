function plotpostmgp(X, y, XTest, yPred, ySd)

[D, nrep] = size(X);

if ~iscell(y),
    yTemp = cell(D,nrep);
    %obtain row sizes from X
    for n = 1:nrep,
        fin = 0;
        for d = 1:D,
            delta = size(X{d, n});
            ini = fin + 1;
            fin = fin + delta;
            yTemp{d,n} = y{n}(ini:fin);
        end
    end
    y = yTemp;
end

if ~iscell(yPred{1,1}),
    yTemp = cell(D,nrep);
    %obtain row sizes from X
    for n = 1:nrep,
        fin = 0;
        for d = 1:D,
            delta = size(XTest{d, n},1);
            ini = fin + 1;
            fin = fin + delta;
            yTemp{d,n} = yPred{n}(ini:fin);
        end
    end
    yPred = yTemp;
end

if ~iscell(ySd{1,1}),
    yTemp = cell(D,nrep);
    %obtain row sizes from X
    for n = 1:nrep,
        fin = 0;
        for d = 1:D,
            delta = size(XTest{d, n},1);
            ini = fin + 1;
            fin = fin + delta;
            yTemp{d,n} = ySd{n}(ini:fin);
            
        end
    end
    ySd = yTemp;
end

markerSize = 20;
markerWidth = 6;
markerType = 'k.';
lineWidth = 2;
fillColor = [0.8 0.8 0.8];
c = 0;
for n=1:nrep,
    for d = 1:D,
        c = c + 1;

        hold off;
        figure(c);
        clf

        fill([XTest{d,n}; XTest{d,n}(end:-1:1)], ...
            [yPred{d,n}; yPred{d,n}(end:-1:1)] ...
            + 2*[ySd{d,n}; -ySd{d,n}], ...
            fillColor,'EdgeColor',fillColor)
        hold on;
        h = plot(XTest{d,n}, yPred{d,n}, 'k-');
        set(h, 'linewidth', lineWidth)
        p = plot(X{d,n}, y{d,n}, markerType);
        set(p, 'markersize', markerSize, 'lineWidth', markerWidth);
    end
end