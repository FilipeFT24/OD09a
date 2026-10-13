function [d1U_dof_h] = f_loop(basis, mesh, intg, flag, c, nu, U_dof, ulb1, urb2, dlb1, drb2, uh, d1uh, ad2uh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
det       = mesh.det;
bf_       = intg.bf_;
Wbf       = intg.Wbf;
F_dof_h   = intg.F_dof;
%--------------------------------------------------------------------------
C_dof     = nu+U_dof;
Clb1      = nu+ulb1;
Crb2      = nu+urb2;
d1U_dof_h = G1 (       mesh, intg, flag,        U_dof,             ulb1, urb2);
d2U_dofCh = G2C(basis, mesh, intg, flag, C_dof, U_dof, Clb1, Crb2, ulb1, urb2, dlb1, drb2);
%--------------------------------------------------------------------------
%{
d1U_dof_a = inittype(1,  d1uh, mesh.xdc, intg.xqc, intg.Wbf_m);
d2U_dofCa = inittype(1, ad2uh, mesh.xdc, intg.xqc, intg.Wbf_m);
figure;
subplot(1, 2, 1);
hold on;
plot(d1U_dof_h*intg.visu, '-r');
plot(d1U_dof_a*intg.visu, 'og');
subplot(1, 2, 2);
hold on;
plot(d2U_dofCh*intg.visu, '-r');
plot(d2U_dofCa*intg.visu, 'og');
%}
%
F        =   F_dof_h*bf_;
d1U      = d1U_dof_h*bf_;
d2UC     = d2U_dofCh*bf_;
d1U2     = d1U.^2;
R        =-d2UC+c.*d1U2+F;
R_dof    =-R*Wbf.*det;

% val = load("this.mat");
% val1 = val.res;
% 
% diff1 = val1+d2U_dof;%-R_dof;
% max(abs(diff1), [], 'all')

figure;
hold on;
plot(R_dof*intg.visu);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end