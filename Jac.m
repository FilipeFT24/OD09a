function [g] = Jac(basis, mesh, intg, flag, nu, c, U_dof, ulb1, urb2, dlb1, drb2, dU, dUf) %#ok<INUSD>
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
bf_          = intg.bf_;
DDr          = intg.DDr;
DFr          = intg.DFr;
FDr          = intg.FDr;
NN           = N*N;
val          = zeros(TNN, 1);
rhs          = zeros(K, N);
%--------------------------------------------------------------------------
A_dof        = nu+U_dof;
dUc          = U_dof*bd_./det;
A            = A_dof*bf_;
dU_t         = dU.';
dUct         = dUc.';
At           = A.';
d1Ult1       = dU_t(:);
d1Uct1       = dUct(:);
At1          =   At(:);
Jacd_RK_NN   = DFr.*d1Uct1+DDr.*At1-2.*c.*FDr.*d1Ult1;
Jacd_NN_K    = sum(reshape(Jacd_RK_NN.', NN, R, []), 2);
Jacd_NN_K_1  = Jacd_NN_K(:);
val(IJd_, 1) = Jacd_NN_K_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
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
    val(IJdl)    = val(IJdl)-VAL1dil_NN_1;
    val(IJol)    =-VAL1oil_NN_1;
    val(IJor)    =-VAL1oir_NN_1;
    val(IJdr)    = val(IJdr)-VAL1dir_NN_1;
    val(IJb1)    = val(IJb1)-VAL1db1_NN_1;
    val(IJb2)    = val(IJb2)-VAL1db2_NN_1;
    rhs(1, :)    =-(bfLrlb1t.*ulb1)*(-2.*c.*dU_t(:, 1).*Wbf);
    rhs(K, :)    =-(bfLlrb2t.*urb2)*(-2.*c.*dU_t(:, K).*Wbf);
    %----------------------------------------------------------------------
    %{
    aux          = zeros(TNN);
    aux(IJdl)    = aux(IJdl)-VAL1dil_NN_1;
    aux(IJol)    =-VAL1oil_NN_1;
    aux(IJor)    =-VAL1oir_NN_1;
    aux(IJdr)    = aux(IJdr)-VAL1dir_NN_1;
    aux(IJb1)    = aux(IJb1)-VAL1db1_NN_1;
    aux(IJb2)    = aux(IJb2)-VAL1db2_NN_1;
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
alb1      = nu+ulb1;
arb2      = nu+urb2;
X_dAdUjac = func_dAdUjac(basis, mesh, intg, A_dof, U_dof, flag, alb1, arb2, ulb1, urb2, dlb1, drb2);
val       = val+X_dAdUjac.lhs;
rhs       = rhs+X_dAdUjac.rhs;
g.MAT     = Mat(val);
g.B_dof   = rhs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dia = load("mat_d.mat"); matd = dia.Jd;
low = load("mat_l.mat"); matl = low.Jlo;
upp = load("mat_u.mat"); matu = upp.Jup;

soll = zeros(N, N, K);
solu = zeros(N, N, K);
sold = zeros(N, N, K);
MAT =  Mat(val);
for i = 1:K
    sold(:, :, i) = MAT((i-1)*N+(1:N), (i-1)*N+(1:N));
    if i > 1
        soll(:, :, i) = MAT((i-1)*N+(1:N), (i-2)*N+(1:N));
    end
    if i < K
        solu(:, :, i) = MAT((i-1)*N+(1:N), (i)*N+(1:N));
    end
end
A = full(MAT);

diffd = sold-matd; max(abs(diffd), [], 'all')
diffu = solu-matu; max(abs(diffu), [], 'all')
diffl = soll-matl; max(abs(diffl), [], 'all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end