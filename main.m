%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
%--------------------------------------------------------------------------
dirs = [...
    "linspecer", ... % https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap
    "QuadTools", ... % https://www.mathworks.com/matlabcentral/fileexchange/72378-quadtools
    "spinit/src", ... % https://www.mathworks.com/matlabcentral/fileexchange/30293-spinit?s_tid=srchtitle
    "debug"];
for i = 1:numel(dirs)
    folder = dirs(1, i);
    if ~isfolder(folder)
        mkdir(folder);
    end
    addpath(folder);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA:
%--------------------------------------------------------------------------
K     = [10, 20, 40, 80];
p     = [1, 2, 3, 4];
tol   = 10*eps;
Xl    =-0.25;
Xr    = 1.25;
dX    = Xr-Xl;
h     = dX./K;
%--------------------------------------------------------------------------
nu    = 1;
c     = 0.5;
arg   = @(x) pi.*x;
Uh    = @(x) sin(arg(x));
dUh   = @(x) pi.*cos(arg(x));
dAdUh = @(x)-pi.^2.*(sin(arg(x)).^2+nu.*sin(arg(x))-cos(arg(x)).^2);
fh    = @(x) pi.^2.*((nu+sin(pi.*x)).*sin(pi.*x)-(1+c).*cos(pi.*x).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP:
%--------------------------------------------------------------------------
flag  = false;
nh    = size(K , 2);
np    = size(p , 2);
ec    = cell(nh, 2);
l2    = cell(1 , 2);
for i = 1:2
    for j = 1:nh
        ec{j, i} = zeros(K(1, j), np);
    end
    l2{1, i} = zeros(nh, np);
    if i == 2
        flag = true;
    end
    for j = 1:nh
        for k = 1:np
            %--------------------------------------------------------------
            Xvt   = Xl:dX/K(1, j):Xr;
            basis = f_aux1(p(1, k));
            mesh  = f_aux2(basis, Xvt);
            intg  = f_aux3(basis, mesh, fh);
            N     = basis.N;
            U_dof = zeros(K(1, j), N);
            delta = realmax;
            %--------------------------------------------------------------
            while delta > tol
                C_dof  = nu+U_dof;
                g1     = Res(basis, mesh, intg, flag, nu, c, U_dof, Uh, dUh, dAdUh);
                dU     = g1.dU;
                dUf    = g1.dUf;
                g2     = Jac(basis, mesh, intg, flag, nu, c, U_dof, Uh, dUh, dU, dUf); % this "dU" already contains the lifting contribution!
                R1dof  = g1.R_dof;
                R2dof  = g2.R_dof;
                J_mat  = g2.MAT;
                B_dof  = R1dof;%-R2dof;
                B_doft = B_dof.';
                B1doft = B_doft(:);
                X1doft = J_mat\B1doft;
                X1dof  = reshape(X1doft, [], K(1, j)).'; % \delta u
                U_dof  = U_dof-X1dof;
                delta  = norm(X1dof, Inf);
            end
            [ec{j, i}(:, k), ...
                l2{1, i}(j, k)] = enorm(basis, mesh, intg, U_dof, Uh);
            %--------------------------------------------------------------
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT:
%--------------------------------------------------------------------------
fprintf("  N   |   Y          N             Y\n-------------   -------------------------\n");
for k = 1:np
    for j = 2:nh
        fprintf("%.3f | %.3f   %.5e | %.5e\n", ...
            log(l2{1, 1}(j, k)./l2{1, 1}(j-1, k))./log(h(1, j)./h(1, j-1)), ...
            log(l2{1, 2}(j, k)./l2{1, 2}(j-1, k))./log(h(1, j)./h(1, j-1)), l2{1, 1}(j, k), l2{1, 2}(j, k));
    end
    fprintf("-------------   -------------------------\n");
end
%Plot1(mesh, ec, l2, 0, "Toy (1).pdf");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%