function [e_c, l2, h1] = enorm(basis, mesh, intg, U_dof, Uh, dUh)
%--------------------------------------------------------------------------
W_     = basis.W_;
det    = mesh.det;
Xq     = intg.Xq;
bd_    = intg.bd_;
bf_    = intg.bf_;
Wbf_m  = intg.Wbf_m;
%--------------------------------------------------------------------------
Ddq    = U_dof*bd_./det;
Udq    = U_dof*bf_;
Daq    = dUh(Xq);
Uaq    =  Uh(Xq);
e_q1   = Udq-Uaq;
e_q2   = Ddq-Daq;
e2c1   = e_q1.^2*W_.*det(:, 1);
e2c2   = e_q2.^2*W_.*det(:, 1);
l2     = sqrt(sum(e2c1, 1));
h1     = sqrt(sum(e2c1+e2c2, 1));
%--------------------------------------------------------------------------
e_qabs = abs(e_q1);
e_c    = e_qabs*Wbf_m;
%--------------------------------------------------------------------------
end