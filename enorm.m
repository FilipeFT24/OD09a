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
e2q1   = e_q1.^2;
e2q2   = e_q2.^2;
e2q1p2 = e2q1+e2q2;
e2c1   = e2q1  *W_.*det;
e2c1p2 = e2q1p2*W_.*det;
l2     = sqrt(sum(e2c1  , 1));
h1     = sqrt(sum(e2c1p2, 1));
%--------------------------------------------------------------------------
e_qabs = abs(e_q1);
e_c    = e_qabs*Wbf_m;
%--------------------------------------------------------------------------
end