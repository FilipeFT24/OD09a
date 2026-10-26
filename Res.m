function [g] = Res(g, flag, nu, c, U_dof, dUh, dAdUh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_dof    = nu+U_dof;
X1       = func_dU  (g,        U_dof, flag, dUh);
dAdU_dof = func_dAdU(g, A_dof, U_dof,       dAdUh);
%--------------------------------------------------------------------------
det      = g.det;
bf_      = g.bf_;
Wbf      = g.Wbf;
F_dof    = g.F_dof;
dU_dof   = X1.dU_doft;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end