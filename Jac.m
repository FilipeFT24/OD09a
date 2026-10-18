function [g] = Jac(basis, mesh, intg, flag, nu, c, U_dof, ulb1, urb2, dlb1, drb2, d1Ul, d1Ux, d2U_dofx) %#ok<INUSD>
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
C_dof        = nu+U_dof;
d1Uc         = U_dof*bd_./det;
C            = C_dof*bf_;
d1Ult        = d1Ul.';
d1Uct        = d1Uc.';
Ct           = C.';
d1Ult1       = d1Ult(:);
d1Uct1       = d1Uct(:);
Ct1          =    Ct(:);
Jacd_RK_NN   = DFr.*d1Uct1+DDr.*Ct1-2.*c.*FDr.*d1Ult1;
Jacd_NN_K    = sum(reshape(Jacd_RK_NN.', NN, R, []), 2);
Jacd_NN_K_1  = Jacd_NN_K(:);
val(IJd_, 1) = Jacd_NN_K_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
K       = mesh.K;
Km1     = K-1;
IJdl    = mesh.IJdl; % ll
IJol    = mesh.IJol; % lr
IJor    = mesh.IJor; % rl
IJdr    = mesh.IJdr; % rr
IJb1    = mesh.IJb1; % rrb1
IJb2    = mesh.IJb2; % llb2
bDl     = intg.bDl;
bDr     = intg.bDr;
bfl     = intg.bfl;
bfr     = intg.bfr;
bDl_    = intg.bDl_;
bDr_    = intg.bDr_;
bfl_    = intg.bfl_;
bfr_    = intg.bfr_;
dl_i    = zeros(Km1, 1);
dr_i    = zeros(Km1, 1);
for i = 1:Km1
    dl_i(i, 1) = bDl(i  , :)*U_dof(i  , :).';
    dr_i(i, 1) = bDr(i+1, :)*U_dof(i+1, :).';
end
drb1    = bDr(1, :)*U_dof(1, :).';
dlb2    = bDl(K, :)*U_dof(K, :).';
ul_i    = U_dof(1:Km1, :)*bfl_;
ur_i    = U_dof(2:K  , :)*bfr_;
urb1    = U_dof(    1, :)*bfr_;
ulb2    = U_dof(    K, :)*bfl_;
nupul_i = nu+ul_i;
nupur_i = nu+ur_i;
nupurb1 = nu+urb1;
nupulb2 = nu+ulb2;
nupulb1 = nu+ulb1;
nupurb2 = nu+urb2;
%--------------------------------------------------------------------------
if flag
    %----------------------------------------------------------------------
    Wbf      = intg.Wbf;
    bfLll_it = intg.bfLll_it;
    bfLlr_it = intg.bfLlr_it;
    bfLrl_it = intg.bfLrl_it;
    bfLrr_it = intg.bfLrr_it;
    bfLllb2t = intg.bfLllb2t;
    bfLlrb2t = intg.bfLlrb2t;
    bfLrlb1t = intg.bfLrlb1t;
    bfLrrb1t = intg.bfLrrb1t;
    J1Lll_it = zeros(N, N, Km1);
    J1Llr_it = zeros(N, N, Km1);
    J1Lrl_it = zeros(N, N, Km1);
    J1Lrr_it = zeros(N, N, Km1);
    for i = 1:Km1
        J1Lll_it(:, :, i) = bfLll_it*(-2.*c.*d1Ult(:, i  ).*Wbf); % det already factored in! (transpose)
        J1Llr_it(:, :, i) = bfLlr_it*(-2.*c.*d1Ult(:, i  ).*Wbf); % =
        J1Lrl_it(:, :, i) = bfLrl_it*(-2.*c.*d1Ult(:, i+1).*Wbf); % =
        J1Lrr_it(:, :, i) = bfLrr_it*(-2.*c.*d1Ult(:, i+1).*Wbf); % =
    end
    J1Lll_it     = reshape(J1Lll_it, [NN, Km1]);
    J1Llr_it     = reshape(J1Llr_it, [NN, Km1]);
    J1Lrl_it     = reshape(J1Lrl_it, [NN, Km1]);
    J1Lrr_it     = reshape(J1Lrr_it, [NN, Km1]);
    J1Lrrb1t     = bfLrrb1t*(-2.*c.*d1Ult(:, 1).*Wbf);
    J1Lllb2t     = bfLllb2t*(-2.*c.*d1Ult(:, K).*Wbf);
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
    rhs(1, :)    = (bfLrlb1t.*ulb1)*(-2.*c.*d1Ult(:, 1).*Wbf.*det(1, 1));
    rhs(K, :)    = (bfLlrb2t.*urb2)*(-2.*c.*d1Ult(:, K).*Wbf.*det(K, 1));
    %----------------------------------------------------------------------
    %{
    aux          = zeros(TNN);
    aux(IJdl)    = aux(IJdl)-VAL1dil_NN_1;
    aux(IJol)    =-VALoil_NN_1;
    aux(IJor)    =-VALoir_NN_1;
    aux(IJdr)    = aux(IJdr)-VAL1oil_NN_1;
    aux(IJb1)    = aux(IJb1)-VAL1oir_NN_1;
    aux(IJb2)    = aux(IJb2)-VAL1oir_NN_1;
    AUX          = Mat(aux);
    B_doft       = U_dof.';
    B1doft       = B_doft(:);
    X1doft       = AUX*B1doft;
    X1dofI       = reshape(X1doft, [], K).'-rhs;
    X1dofE       =-2.*c.*d1Ul.*d1Ux*Wbf.*det;
    diffl1       = X1dofI-X1dofE; max(abs(diffl1), [], 'all')
    %}
    %----------------------------------------------------------------------
end
%--------------------------------------------------------------------------

val = val.*0;
rhs = rhs.*0;

J2aLll_it = zeros(N, N, Km1);
J2aLlr_it = zeros(N, N, Km1);
J2aLrl_it = zeros(N, N, Km1);
J2aLrr_it = zeros(N, N, Km1);
J2bLll_it = zeros(N, N, Km1);
J2bLlr_it = zeros(N, N, Km1);
J2bLrl_it = zeros(N, N, Km1);
J2bLrr_it = zeros(N, N, Km1);
J2cLll_it = zeros(N, N, Km1);
J2cLlr_it = zeros(N, N, Km1);
J2cLrl_it = zeros(N, N, Km1);
J2cLrr_it = zeros(N, N, Km1);
J2bLrrb1t = zeros(N, N);
J2bLllb2t = zeros(N, N);
if flag
    for i = 1:Km1
        J2bLll_it(:, :, i) = bfl_*(1./2.*bDl(i  , :).*(nupul_i(i, 1)+(ul_i(i, 1)-ur_i(i, 1)))); % (transpose)
        J2bLlr_it(:, :, i) =-bfr_*(1./2.*bDl(i  , :).*(nupul_i(i, 1)));                         % =
        J2bLrl_it(:, :, i) = bfl_*(1./2.*bDr(i+1, :).*(nupur_i(i, 1)));                         % =
        J2bLrr_it(:, :, i) =-bfr_*(1./2.*bDr(i+1, :).*(nupur_i(i, 1)+(ur_i(i, 1)-ul_i(i, 1)))); % =
    end
    J2bLrrb1t =-bfr_*(bDr(1, :).*(nupurb1+(urb1-ulb1))); % =
    J2bLllb2t = bfl_*(bDl(K, :).*(nupulb2+(ulb2-urb2))); % =
end
for i = 1:Km1
    J2aLll_it(:, :, i) = (dl_i(i, 1).*bfl_+nupul_i(i, 1).*bDl_(:, i  ))./2*bfl; % (transpose)
    J2aLlr_it(:, :, i) = (dr_i(i, 1).*bfr_+nupur_i(i, 1).*bDr_(:, i+1))./2*bfl; % =
    J2aLrl_it(:, :, i) =-(dl_i(i, 1).*bfl_+nupul_i(i, 1).*bDl_(:, i  ))./2*bfr; % =
    J2aLrr_it(:, :, i) =-(dr_i(i, 1).*bfr_+nupur_i(i, 1).*bDr_(:, i+1))./2*bfr; % =
end
J2aLrrb1t    =-(drb1.*bfr_+nupurb1.*bDr_(:, 1))./2*bfr; % =
J2aLllb2t    = (dlb2.*bfl_+nupulb2.*bDl_(:, K))./2*bfl; % =
J2Lll_it     = J2aLll_it+J2bLll_it+J2cLll_it;
J2Llr_it     = J2aLlr_it+J2bLlr_it+J2cLlr_it;
J2Lrl_it     = J2aLrl_it+J2bLrl_it+J2cLrl_it;
J2Lrr_it     = J2aLrr_it+J2bLrr_it+J2cLrr_it;
J2Lrrb1t     = J2aLrrb1t+J2bLrrb1t;
J2Lllb2t     = J2aLllb2t+J2bLllb2t;
J2Lll_it     = reshape(J2Lll_it, [NN, Km1]);
J2Llr_it     = reshape(J2Llr_it, [NN, Km1]);
J2Lrl_it     = reshape(J2Lrl_it, [NN, Km1]);
J2Lrr_it     = reshape(J2Lrr_it, [NN, Km1]);
VAL2dil_NN_1 = J2Lll_it(:);
VAL2oil_NN_1 = J2Llr_it(:);
VAL2oir_NN_1 = J2Lrl_it(:);
VAL2dir_NN_1 = J2Lrr_it(:);
VAL2db1_NN_1 = J2Lrrb1t(:);
VAL2db2_NN_1 = J2Lllb2t(:);
val(IJdl)    = val(IJdl)-VAL2dil_NN_1;
val(IJol)    = val(IJol)-VAL2oil_NN_1;
val(IJor)    = val(IJor)-VAL2oir_NN_1;
val(IJdr)    = val(IJdr)-VAL2dir_NN_1;
val(IJb1)    = val(IJb1)-VAL2db1_NN_1;
val(IJb2)    = val(IJb2)-VAL2db2_NN_1;
rhs(1, :)    = rhs(1, :)-(bfr.*(dlb1.*ulb1+nupulb1.*dlb1)./2-bDr(1, :).*nupurb1.*ulb1); % <-o
rhs(K, :)    = rhs(K, :)+(bfl.*(drb2.*urb2+nupurb2.*drb2)./2-bDl(K, :).*nupulb2.*urb2); % o->
%
aux          = zeros(TNN);
aux(IJdl)    = aux(IJdl)+VAL2dil_NN_1;
aux(IJol)    = VAL2oil_NN_1;
aux(IJor)    = VAL2oir_NN_1;
aux(IJdr)    = aux(IJdr)+VAL2dir_NN_1;
aux(IJb1)    = aux(IJb1)+VAL2db1_NN_1;
aux(IJb2)    = aux(IJb2)+VAL2db2_NN_1;
AUX          = Mat(aux);
B_doft       = U_dof.';
B1doft       = B_doft(:);
X1doft       = AUX*B1doft;
X1dofI       = reshape(X1doft, [], K).'+rhs;
X1dofE       = d2U_dofx;
mass         = intg.mass;
X1dofE       = X1dofE*mass.*det;
diff_2       = X1dofI-X1dofE; max(abs(diff_2), [], 'all')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE JACOBIAN AND ASSIGN TO...
%--------------------------------------------------------------------------
g.MAT   = Mat(val);
g.R_dof = rhs;


% dia = load("mat_d.mat"); matd = dia.Jd;
% low = load("mat_l.mat"); matl = low.Jlo;
% upp = load("mat_u.mat"); matu = upp.Jup;
% 
% soll = zeros(N, N, K);
% solu = zeros(N, N, K);
% sold = zeros(N, N, K);
% MAT = g.MAT;
% for i = 1:K
%     sold(:, :, i) = MAT((i-1)*N+(1:N), (i-1)*N+(1:N));
%     if i > 1
%         soll(:, :, i) = MAT((i-1)*N+(1:N), (i-2)*N+(1:N));
%     end
%     if i < K
%         solu(:, :, i) = MAT((i-1)*N+(1:N), (i)*N+(1:N));
%     end
% end
% A = full(MAT);
% 
% diffd = sold-matd; max(abs(diffd), [], 'all')
% diffu = solu-matu; max(abs(diffu), [], 'all')
% diffl = soll-matl; max(abs(diffl), [], 'all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end