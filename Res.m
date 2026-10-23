function [g] = Res(basis, mesh, intg, flag, nu, c, U_dof, ulb1, urb2, dlb1, drb2, dUh, dAdUh)
%--------------------------------------------------------------------------
det      = mesh.det;
bf_      = intg.bf_;
Wbf      = intg.Wbf;
%--------------------------------------------------------------------------
A_dof    = nu+U_dof;
alb1     = nu+ulb1;
arb2     = nu+urb2;
X1       = func_dU  (basis, mesh, intg,        U_dof, flag,             ulb1, urb2,               dUh);
dAdU_dof = func_dAdU(basis, mesh, intg, A_dof, U_dof, flag, alb1, arb2, ulb1, urb2, dlb1, drb2, dAdUh);
%--------------------------------------------------------------------------
dU_dof   = X1.dU_doft;
F_dof    = intg.F_dof;
F        = F_dof*bf_;
dU       = dU_dof*bf_;
dAdU     = dAdU_dof*bf_;
dU2      = dU.^2;
R        =-dAdU-c.*dU2-F;
R_dof    = R*Wbf.*det;
%--------------------------------------------------------------------------
dU_doff  = X1.dU_doff;
dUf      = dU_doff*bf_;
g.dU     = dU;
g.dUf    = dUf;   % aux. variable
g.R_dof  = R_dof; % residual


% res1 = load("res.mat");
% res  = res1.R;
% diffres = res-R_dof;
% 
% max(abs(diffres), [], 'all')

%--------------------------------------------------------------------------
end