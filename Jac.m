function [g] = Jac(basis, mesh, intg, flag, nu, c, U_dof, Uh, dUh, dU, dUf) %#ok<INUSD>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
N            = basis.N;
K            = mesh.K;
det          = mesh.det;
TNN          = mesh.TNN;
IJd_         = mesh.IJd_;
Mat          = mesh.Jac;
R            = intg.R;
bd_          = intg.bd_;
DFr          = intg.DFr;
FDr          = intg.FDr;
NN           = N*N;
val          = zeros(TNN, 1);
rhs          = zeros(K, N);
%--------------------------------------------------------------------------
A_dof        = nu+U_dof;
dUc          = U_dof*bd_./det;
dU_t         = dU.';
dUct         = dUc.';
d1Ult1       = dU_t(:);
d1Uct1       = dUct(:);
Jacd_RK_NN   = DFr.*d1Uct1-2.*c.*FDr.*d1Ult1;
Jacd_NN_K    = sum(reshape(Jacd_RK_NN.', NN, R, []), 2);
Jacd_NN_K_1  = Jacd_NN_K(:);
val(IJd_, 1) = Jacd_NN_K_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
Xl   = mesh.Xl;
Xr   = mesh.Xr;
ulb1 =  Uh(Xl);
urb2 =  Uh(Xr);
dlb1 = dUh(Xl);
drb2 = dUh(Xr);
alb1 = nu+ulb1;
arb2 = nu+urb2;
%--------------------------------------------------------------------------
if flag
    %----------------------------------------------------------------------
    Km1      = K-1;
    IJdl     = mesh.IJdl; % ll
    IJol     = mesh.IJol; % lr
    IJor     = mesh.IJor; % rl
    IJdr     = mesh.IJdr; % rr
    IJb1     = mesh.IJb1; % rrb1
    IJb2     = mesh.IJb2; % llb2
    Wbf      = intg.Wbf;
    bfLll_it = intg.bfLll_it;
    bfLlr_it = intg.bfLlr_it;
    bfLrl_it = intg.bfLrl_it;
    bfLrr_it = intg.bfLrr_it;
    bfLllb2t = intg.bfLllb2t;
    bfLlrb2t = intg.bfLlrb2t;
    bfLrlb1t = intg.bfLrlb1t;
    bfLrrb1t = intg.bfLrrb1t;
    %----------------------------------------------------------------------
    J1Lll_it = zeros(N, N, Km1);
    J1Llr_it = zeros(N, N, Km1);
    J1Lrl_it = zeros(N, N, Km1);
    J1Lrr_it = zeros(N, N, Km1);
    for i = 1:Km1
        J1Lll_it(:, :, i) = bfLll_it*(-2.*c.*dU_t(:, i  ).*Wbf); % det already factored in! (transpose)
        J1Llr_it(:, :, i) = bfLlr_it*(-2.*c.*dU_t(:, i  ).*Wbf); % =
        J1Lrl_it(:, :, i) = bfLrl_it*(-2.*c.*dU_t(:, i+1).*Wbf); % =
        J1Lrr_it(:, :, i) = bfLrr_it*(-2.*c.*dU_t(:, i+1).*Wbf); % =
    end
    J1Lll_it     = reshape(J1Lll_it, [NN, Km1]);
    J1Llr_it     = reshape(J1Llr_it, [NN, Km1]);
    J1Lrl_it     = reshape(J1Lrl_it, [NN, Km1]);
    J1Lrr_it     = reshape(J1Lrr_it, [NN, Km1]);
    J1Lrrb1t     = bfLrrb1t*(-2.*c.*dU_t(:, 1).*Wbf);
    J1Lllb2t     = bfLllb2t*(-2.*c.*dU_t(:, K).*Wbf);
    VAL1dil_NN_1 = J1Lll_it(:);
    VAL1oil_NN_1 = J1Llr_it(:);
    VAL1oir_NN_1 = J1Lrl_it(:);
    VAL1dir_NN_1 = J1Lrr_it(:);
    VAL1db1_NN_1 = J1Lrrb1t(:);
    VAL1db2_NN_1 = J1Lllb2t(:);
    val(IJdl, 1) = val(IJdl, 1)-VAL1dil_NN_1;
    val(IJol, 1) =-VAL1oil_NN_1;
    val(IJor, 1) =-VAL1oir_NN_1;
    val(IJdr, 1) = val(IJdr, 1)-VAL1dir_NN_1;
    val(IJb1, 1) = val(IJb1, 1)-VAL1db1_NN_1;
    val(IJb2, 1) = val(IJb2, 1)-VAL1db2_NN_1;
    rhs(1, :)    =-(bfLrlb1t.*ulb1)*(-2.*c.*dU_t(:, 1).*Wbf);
    rhs(K, :)    =-(bfLlrb2t.*urb2)*(-2.*c.*dU_t(:, K).*Wbf);
    %----------------------------------------------------------------------
    %{
    aux          = zeros(TNN);
    aux(IJdl, 1) = aux(IJdl, 1)-VAL1dil_NN_1;
    aux(IJol, 1) =-VAL1oil_NN_1;
    aux(IJor, 1) =-VAL1oir_NN_1;
    aux(IJdr, 1) = aux(IJdr, 1)-VAL1dir_NN_1;
    aux(IJb1, 1) = aux(IJb1, 1)-VAL1db1_NN_1;
    aux(IJb2, 1) = aux(IJb2, 1)-VAL1db2_NN_1;
    AUX          = Mat(aux);
    B_doft       = U_dof.';
    B1doft       = B_doft(:);
    X1doft       = AUX*B1doft;
    X1dofI       = reshape(X1doft, [], K).'+rhs;
    X1dofE       =-dUf.*(-2.*c.*dU)*Wbf.*det;
    diffl1       = X1dofI-X1dofE; max(abs(diffl1), [], 'all')
    %}
    %----------------------------------------------------------------------
end
%--------------------------------------------------------------------------
X_dAdUjac = func_dAdUjac(basis, mesh, intg, A_dof, U_dof, alb1, arb2, ulb1, urb2, dlb1, drb2);
val       = val-X_dAdUjac.lhs;
rhs       = rhs-X_dAdUjac.rhs;
g.MAT     = Mat(val);
g.R_dof   = rhs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end