function [d1U_dof] = f_loop(mesh, intg, flag, c, nu, U_dof, Ulb1, Urb2, dUlb1, dUrb2, uh, duh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
det      = mesh.det;
bf_      = intg.bf_;
Wbf      = intg.Wbf;
F_dof    = intg.F_dof;
%--------------------------------------------------------------------------
A_dof    = nu+U_dof;
Alb1     = nu+Ulb1;
Arb2     = nu+Urb2;
d1U_dof  = G1(mesh, intg, flag,        U_dof,             Ulb1, Urb2);
Ad2U_dof = G2(mesh, intg, flag, A_dof, U_dof, Alb1, Arb2, Ulb1, Urb2, dUlb1, dUrb2);
%--------------------------------------------------------------------------
%{
F        =   F_dof*bf_;
d1U      = d1U_dof*bf_;
Ad2U     = Ad2U_dof*bf_;
d1U2     = d1U.^2;
R        = Ad2U+c.*d1U2+F;
%}
Ad2U     = Ad2U_dof*bf_;
R        = Ad2U;
R_dof    =-R*Wbf.*det;

val = load("this.mat");
val1 = val.res;

diff1 = val1+Ad2U_dof;%-R_dof;
max(abs(diff1), [], 'all')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end