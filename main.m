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
flag  = 1; % 0: w/o lifting, 1: w/ lifting
p     = 1;
%{
h     = 0.01;
Xvt   = 0:h:1;
% Xvt(1, 2) = 0.01;
% Xvt(1, 3) = 0.072;
% Xvt(1, 4) = 0.123;
% Xvt(1, 7) = 0.265;

%}
h     = 0.005;
Xvt   =-0.25:h:1.25;


nu    = 1;
c     = 0.5;
arg   = @(x) pi.*x;
Uh    = @(x) sin(arg(x));
dUh   = @(x) pi.*cos(arg(x));
dAdUh = @(x)-pi.^2.*(sin(arg(x)).^2+nu.*sin(arg(x))-cos(arg(x)).^2);
fh    = @(x) pi.^2.*((nu+sin(pi.*x)).*sin(pi.*x)-(1+c).*cos(pi.*x).^2);
%--------------------------------------------------------------------------
basis = f_aux1(p);
mesh  = f_aux2(basis, Xvt);
intg  = f_aux3(basis, mesh, fh);
visu  = intg.visu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP:
%--------------------------------------------------------------------------
delta = realmax;
tol   = 1e-12;
%--------------------------------------------------------------------------
N     = basis.N;
K     = mesh.K;
Ke    = K+1;
Xl    = Xvt(1,   1);
Xr    = Xvt(1, end);
%--------------------------------------------------------------------------
U_dof = zeros(K, N);
ulb1  =  Uh(Xl);
urb2  =  Uh(Xr);
dlb1  = dUh(Xl);
drb2  = dUh(Xr);
%{
rng("default");
U_dof = rand(size(U_dof));
ulb1  = 1.50;
urb2  = 2.25;
dlb1  = 2.50;
drb2  = 4.00;
%}
%--------------------------------------------------------------------------
while delta > tol
    %----------------------------------------------------------------------
    C_dof  = nu+U_dof;
    g1     = Res(basis, mesh, intg, flag, nu, c, U_dof, ulb1, urb2, dlb1, drb2, dUh, dAdUh);
    dU     = g1.dU;
    dUf    = g1.dUf;
    g2     = Jac(basis, mesh, intg, flag, nu, c, U_dof, ulb1, urb2, dlb1, drb2, dU, dUf); % this "dU" already contains the lifting contribution!
    %----------------------------------------------------------------------
    R1dof  = g1.R_dof;
    R2dof  = g2.R_dof;
    J_mat  = g2.MAT;
    B_dof  = R1dof;%-R2dof;
    B_doft = B_dof.';
    B1doft = B_doft(:);
    X1doft = J_mat\B1doft;
    X1dof  = reshape(X1doft, [], K).'; % \delta u
    U_dof  = U_dof-X1dof;
    delta  = norm(X1dof, Inf);
    %----------------------------------------------------------------------

%     lhs    = load("LHS.mat"); LHS = lhs.J;
%     rhs    = load("RHS.mat"); RHS = rhs.R;
%     new    = load("NEW.mat"); SOL = new.U;
%     diff_l = J_mat -LHS; max(abs(diff_l), [], 'all')
%     diff_r = B1doft-RHS; max(abs(diff_r), [], 'all')
%     U_doft = U_dof.';
%     U1doft = U_doft(:);
%     diff_u = U1doft-SOL; max(abs(diff_u), [], 'all')

end
%--------------------------------------------------------------------------
Plot1(basis, mesh, intg, U_dof, Uh, 0, "Toy (1).pdf");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%