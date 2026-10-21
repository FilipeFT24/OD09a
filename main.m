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
p     = 2;
h     = 0.01;
Xvt   = 0:h:1;

Xvt(1, 2) = 0.01;
Xvt(1, 3) = 0.072;
Xvt(1, 4) = 0.123;
Xvt(1, 7) = 0.265;

nu    = 1;
c     = 0.5;
arg   = @(x) pi.*x;
uh    = @(x) sin(arg(x));
d1uh  = @(x) pi.*cos(arg(x));
cd2uh = @(x)-pi.^2.*(sin(arg(x)).^2+nu.*sin(arg(x))-cos(arg(x)).^2);
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
K     = mesh.K;
Ke    = K+1;
Xl    = Xvt(1,   1);
Xr    = Xvt(1, end);

ulb1  =  1.5;%uh(0);
urb2  = 2.25;%uh(1);
ulb1  =   uh(Xl);
urb2  =   uh(Xr);
dlb1  = d1uh(Xl);
drb2  = d1uh(Xr);
U_dof = zeros(K, p+1);
%--------------------------------------------------------------------------

rng("default");
U_dof = rand(size(U_dof));
%U_dof = inittype(0, uh, mesh.xdc, intg.xqc, intg.Wbf_m);
%U_dof = zeros(K, p+1);

while delta > tol
    %----------------------------------------------------------------------
    C_dof    = nu+U_dof;
    g1       = Res(basis, mesh, intg, flag, nu, c, U_dof, ulb1, urb2, dlb1, drb2, uh, d1uh, cd2uh);
    d1Ul     = g1.d1Ul;
    d1Ux     = g1.d1Ux;
    d2U_dofx = g1.d2U_dofx;
    g2       = Jac(basis, mesh, intg, flag, nu, c, U_dof, ulb1, urb2, dlb1, drb2, d1Ul, d1Ux, d2U_dofx); % this "d1UL" already contains the lifting contribution!
    %----------------------------------------------------------------------
    R1dof    = g1.R_dof;
    R2dof    = g2.R_dof;
    J_mat    = g2.MAT;
    B_dof    = R1dof+R2dof;
    B_doft   = B_dof.';
    B1doft   = B_doft(:);
    X1doft   = J_mat\B1doft;
    X1dof    = reshape(X1doft, [], K).'; % \delta u
    U_dof    = U_dof-X1dof;

    figure;
    hold on;
    plot(U_dof*visu, '-b');

    delta  = norm(X1dof, inf);
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%