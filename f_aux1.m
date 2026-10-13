function [basis] = f_aux1(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N  = p+1;
dh = cell(1, N);
fh = cell(1, N);
switch p
    case 0
        xp = [1/2]; %#ok<NBRAK2> 
    case 1
        xp = [0, 1];
    case 2
        xp = [0, 1, 1/2];
    case 3
        xp = [0, 1, 1/3, 2/3];
    case 4
        xp = [0, 1, 1/4, 1/2, 3/4];
    otherwise
        return;
end
%--------------------------------------------------------------------------
syms x;
for i = 1:N
    fi = sym(1);
    for j = 1:N
        if j ~= i
            fi = fi.*(x-xp(1, j))./(xp(1, i)-xp(1, j));
        end
    end
    di = simplify(diff(fi, 1, x));
    if di == 0
        dh{1, i} = @(x) zeros(size(x)); 
        fh{1, i} = @(x)  ones(size(x));
    else
        if diff(di, 1, x) == 0
            dh{1, i} = @(x) ones(size(x)).*double(di);
        else
            dh{1, i} = matlabFunction(di, 'Vars', x);
        end
        fh{1, i} = matlabFunction(fi, 'Vars', x);
    end
end
%--------------------------------------------------------------------------
Qs  = quadGaussLegendre(2*p+2, 'Domain', [0, 1]);
QsP = Qs.Points;
W_  = Qs.Weights;
R   = size(Qs.Weights, 1);
%--------------------------------------------------------------------------
basis.dh  = dh;
basis.fh  = fh;
basis.xp  = xp;
basis.QsP = QsP;
basis.W_  = W_;
basis.N   = N;
basis.R   = R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end