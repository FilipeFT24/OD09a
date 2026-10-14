function [d1U_dof_h] = f_loop(basis, mesh, intg, flag, c, nu, U_dof, ulb1, urb2, dlb1, drb2, uh, d1uh, cd2uh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_dof = nu+U_dof;
clb1  = nu+ulb1;
crb2  = nu+urb2;
R_dof = Res(basis, mesh, intg, flag, c, C_dof, U_dof, clb1, crb2, ulb1, urb2, dlb1, drb2, uh, d1uh, cd2uh);
%--------------------------------------------------------------------------
J_dof = Jac(basis, mesh, intg, flag, c, C_dof, U_dof);


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