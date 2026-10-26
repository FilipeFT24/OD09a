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
%     = 4*2.^(0:4);
K     = 4*2.^(1:4);
p     = 1:4;
tol   = 100*eps;
Xl    = 0;
Xr    = 1;
dX    = Xr-Xl;
h     = dX./K;
%--------------------------------------------------------------------------
nu    = 1;
c     = 0.5;
arg   = @(x) pi.*x;
Uh    = @(x) sin(arg(x));
dUh   = @(x) pi.*cos(arg(x));
dAdUh = @(x)-pi.^2.*(sin(arg(x)).^2+nu.*sin(arg(x))-cos(arg(x)).^2);
Fh    = @(x) pi.^2.*((nu+sin(pi.*x)).*sin(pi.*x)-(1+c).*cos(pi.*x).^2);
%--------------------------------------------------------------------------
Ph    = @(x) sin(arg(x));
dPh   = @(x) pi.*cos(arg(x));
Gh    = @(x) pi.^2.*(2.*c.*(cos(arg(x)).^2-sin(arg(x)).^2)+(nu+sin(arg(x))).*sin(arg(x)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP:
%--------------------------------------------------------------------------
flag  = false;
nh    = size(K, 2);
np    = size(p, 2);
Fd    = cell(1, 2);  % primal  error distribution
Fda   = cell(1, 2);  % adjoint error distribution
l2    = cell(1, 2);  % primal  error in the        l2-norm
h1    = cell(1, 2);  % primal  error in the broken h1-norm
l2a   = cell(1, 2);  % adjoint error in the        l2-norm
h1a   = cell(1, 2);  % adjoint error in the broken h1-norm
FUNC  = cell(1, np);
VISU  = cell(1, np);
for i = 1:2
    for j = 1:nh
        for k = 1:np
            Fd{1, i}{j, k} = zeros(K(1, j), p(1, k)+1);
        end
    end
    l2{1, i} = zeros(nh, np);
    if i == 2
        flag = true;
    end
    for j = 1:nh
        for k = 1:np
            %--------------------------------------------------------------
            Xvt   = Xl:dX/K(1, j):Xr;
            aux   = func_aux(p(1, k), Xvt, nu, Fh, Gh, Uh, dUh);
            N     = aux.N;
            det   = aux.det;
            mass  = aux.mass;
            U_dof = zeros(K(1, j), N);
            delta = realmax;
            if j == 1
                FUNC{1, k} = aux.fh;
                VISU{1, k} = aux.visu;
            end
            %--------------------------------------------------------------
            while delta > tol
                g1     = Res(aux, flag, nu, c, U_dof, dUh, dAdUh);
                g2     = Jac(aux, flag, nu, c, U_dof, g1.dU, g1.dUf); % this "dU" already contains the lifting contribution!
                R_dof  = g1.R_dof;
                J_mat  = g2.MAT;
                B_dof  = R_dof;
                B_doft = B_dof.';
                B1doft = B_doft(:);
                X1doft = J_mat\B1doft;
                X1dof  = reshape(X1doft, [], K(1, j)).'; % \delta u
                U_dof  = U_dof-X1dof;
                delta  = norm(X1dof, Inf);
            end
            [Fd{1, i}{j, k}, l2{1, i}(j, k), h1{1, i}(j, k)] = ...
                enorm(aux, U_dof, Uh, dUh);
            %--------------------------------------------------------------
            g1     = Res(aux, flag, nu, c, U_dof, dUh, dAdUh);
            g2     = Jac(aux, flag, nu, c, U_dof, g1.dU, g1.dUf);
            J_mat  = g2.MAT;
            G_dof  = aux.Gkdof;
            G_doft = G_dof.';
            G1doft = G_doft(:);
            P_doft = J_mat.'\G1doft;
            P_dof  = reshape(P_doft, [], K(1, j)).';
            [Fda{1, i}{j, k}, l2a{1, i}(j, k), h1a{1, i}(j, k)] = ...
                enorm(aux, P_dof, Ph, dPh);
            %--------------------------------------------------------------
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT:
%--------------------------------------------------------------------------
fprintf("______________________________________________________________________________________________\n");
fprintf("______________________________________________________________________________________________\n");
fprintf("                                        PRIMAL PROBLEM\n");
fprintf("----------------------------------------------------------------------------------------------\n");
fprintf("                   l2-norm                                          h1-norm (broken)          \n");
fprintf("---------------------------------------------    ---------------------------------------------\n");
fprintf("   N   |    Y           N      |      Y             N   |    Y           N      |      Y      \n---------------   ---------------------------    ---------------   ---------------------------\n");
for k = 1:np
    for j = 1:nh
        if j == 1
            fprintf("                  %.6e | %.6e                      %.6e | %.6e\n", l2{1, 1}(j, k), l2{1, 2}(j, k), h1{1, 1}(j, k), h1{1, 2}(j, k));
        else
            fprintf("%.4f | %.4f   %.6e | %.6e    %.4f | %.4f   %.6e | %.6e\n", ...
                log(l2{1, 1}(j, k)./l2{1, 1}(j-1, k))./log(h(1, j)./h(1, j-1)), ...
                log(l2{1, 2}(j, k)./l2{1, 2}(j-1, k))./log(h(1, j)./h(1, j-1)), l2{1, 1}(j, k), l2{1, 2}(j, k), ...
                log(h1{1, 1}(j, k)./h1{1, 1}(j-1, k))./log(h(1, j)./h(1, j-1)), ...
                log(h1{1, 2}(j, k)./h1{1, 2}(j-1, k))./log(h(1, j)./h(1, j-1)), h1{1, 1}(j, k), h1{1, 2}(j, k));
        end
    end
    fprintf("---------------   ---------------------------    ---------------   ---------------------------\n");
end
fprintf("\n");
fprintf("______________________________________________________________________________________________\n");
fprintf("______________________________________________________________________________________________\n");
fprintf("                                       ADJOINT PROBLEM\n");
fprintf("----------------------------------------------------------------------------------------------\n");
fprintf("                   l2-norm                                          h1-norm (broken)          \n");
fprintf("---------------------------------------------    ---------------------------------------------\n");
fprintf("   N   |    Y           N      |      Y             N   |    Y           N      |      Y      \n---------------   ---------------------------    ---------------   ---------------------------\n");
for k = 1:np
    for j = 1:nh
        if j == 1
            fprintf("                  %.6e | %.6e                      %.6e | %.6e\n", l2a{1, 1}(j, k), l2a{1, 2}(j, k), h1a{1, 1}(j, k), h1a{1, 2}(j, k));
        else
            if h1a{1, 1}(j-1, k) >= h1a{1, 1}(j, k)
                fprintf("%.4f | %.4f   %.6e | %.6e    %.4f | %.4f   %.6e | %.6e\n", ...
                    log(l2a{1, 1}(j, k)./l2a{1, 1}(j-1, k))./log(h(1, j)./h(1, j-1)), ...
                    log(l2a{1, 2}(j, k)./l2a{1, 2}(j-1, k))./log(h(1, j)./h(1, j-1)), l2a{1, 1}(j, k), l2a{1, 2}(j, k), ...
                    log(h1a{1, 1}(j, k)./h1a{1, 1}(j-1, k))./log(h(1, j)./h(1, j-1)), ...
                    log(h1a{1, 2}(j, k)./h1a{1, 2}(j-1, k))./log(h(1, j)./h(1, j-1)), h1a{1, 1}(j, k), h1a{1, 2}(j, k));
            else
                fprintf("%.4f | %.4f   %.6e | %.6e   %.4f | %.4f   %.6e | %.6e\n", ...
                    log(l2a{1, 1}(j, k)./l2a{1, 1}(j-1, k))./log(h(1, j)./h(1, j-1)), ...
                    log(l2a{1, 2}(j, k)./l2a{1, 2}(j-1, k))./log(h(1, j)./h(1, j-1)), l2a{1, 1}(j, k), l2a{1, 2}(j, k), ...
                    log(h1a{1, 1}(j, k)./h1a{1, 1}(j-1, k))./log(h(1, j)./h(1, j-1)), ...
                    log(h1a{1, 2}(j, k)./h1a{1, 2}(j-1, k))./log(h(1, j)./h(1, j-1)), h1a{1, 1}(j, k), h1a{1, 2}(j, k));
            end
        end
    end
    fprintf("---------------   ---------------------------    ---------------   ---------------------------\n");
end
%Plot1(Xl, Xr, Fd, Fda, h, l2, l2a, h1, h1a, FUNC, VISU, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%