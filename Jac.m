function [MAT] = Jac(basis, mesh, intg, flag, c, C_dof, U_dof, ulb1, urb2, d1Ul, d1Ux) %#ok<INUSD> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
N            = basis.N;
NN           = N*N;
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
val          = zeros(TNN, 1);
%--------------------------------------------------------------------------
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
K    = mesh.K;
IJdl = mesh.IJdl; % ll
IJol = mesh.IJol; % lr
IJor = mesh.IJor; % rl
IJdr = mesh.IJdr; % rr
IJb1 = mesh.IJb1; % rrb1
IJb2 = mesh.IJb2; % llb2
Km1  = K-1;
%--------------------------------------------------------------------------
if flag
    %----------------------------------------------------------------------
    Wbf         = intg.Wbf;
    bfLll_it    = intg.bfLll_it;
    bfLlr_it    = intg.bfLlr_it;
    bfLrl_it    = intg.bfLrl_it;
    bfLrr_it    = intg.bfLrr_it;
    bfLrrb1t    = intg.bfLrrb1t; % rr
    bfLrlb1t    = intg.bfLrlb1t; % rl
    bfLllb2t    = intg.bfLllb2t; % ll
    bfLlrb2t    = intg.bfLlrb2t; % lr
    JLll_it     = zeros(N, N, Km1);
    JLlr_it     = zeros(N, N, Km1);
    JLrl_it     = zeros(N, N, Km1);
    JLrr_it     = zeros(N, N, Km1);
    for i = 1:Km1
        JLll_it(:, :, i) =-bfLll_it(:, :, i)*(-2.*c.*d1Ult(:, i  ).*Wbf.*det(i  , 1));
        JLlr_it(:, :, i) =-bfLlr_it(:, :, i)*(-2.*c.*d1Ult(:, i  ).*Wbf.*det(i  , 1));
        JLrl_it(:, :, i) =-bfLrl_it(:, :, i)*(-2.*c.*d1Ult(:, i+1).*Wbf.*det(i+1, 1));
        JLrr_it(:, :, i) =-bfLrr_it(:, :, i)*(-2.*c.*d1Ult(:, i+1).*Wbf.*det(i+1, 1));
    end
    JLrrb1t     = bfLrrb1t*(-2.*c.*d1Ult(:, 1).*Wbf.*det(1, 1));
    JLllb2t     =-bfLllb2t*(-2.*c.*d1Ult(:, K).*Wbf.*det(K, 1));
    JLll_it     = reshape(JLll_it, [NN, Km1]);
    JLlr_it     = reshape(JLlr_it, [NN, Km1]);
    JLrl_it     = reshape(JLrl_it, [NN, Km1]);
    JLrr_it     = reshape(JLrr_it, [NN, Km1]);
    VALdil_NN_1 = JLll_it(:);
    VALoil_NN_1 = JLlr_it(:);
    VALoir_NN_1 = JLrl_it(:);
    VALdir_NN_1 = JLrr_it(:);
    VALdb1_NN_1 = JLrrb1t(:);
    VALdb2_NN_1 = JLllb2t(:);
    val(IJdl)   = val(IJdl)+VALdil_NN_1;
    val(IJol)   = VALoil_NN_1;
    val(IJor)   = VALoir_NN_1;
    val(IJdr)   = val(IJdr)+VALdir_NN_1;
    val(IJb1)   = val(IJb1)+VALdb1_NN_1;
    val(IJb2)   = val(IJb2)+VALdb2_NN_1;
    %----------------------------------------------------------------------
    %
    aux         = zeros(TNN);
    aux(IJdl)   = aux(IJdl)+VALdil_NN_1;
    aux(IJol)   = VALoil_NN_1;
    aux(IJor)   = VALoir_NN_1;
    aux(IJdr)   = aux(IJdr)+VALdir_NN_1;
    aux(IJb1)   = aux(IJb1)+VALdb1_NN_1;
    aux(IJb2)   = aux(IJb2)+VALdb2_NN_1;
    AUX         = Mat(aux);
    RHSL        = zeros(K, N);
    RHSL(1, :)  = ulb1.*bfLrlb1t*(-2.*c.*d1Ult(:, 1).*Wbf.*det(1, 1));
    RHSL(K, :)  = urb2.*bfLlrb2t*(-2.*c.*d1Ult(:, K).*Wbf.*det(K, 1));
    B_doft      = U_dof.';
    B1doft      = B_doft(:);
    X1doft      = AUX*B1doft;
    X1dofI      = reshape(X1doft, [], K).'+RHSL;
    X1dofE      =-2.*c.*d1Ul.*d1Ux*Wbf.*det;
    diffl1      = X1dofI-X1dofE; max(abs(diffl1), [], 'all')
    %}
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE JACOBIAN:
%--------------------------------------------------------------------------
MAT = Mat(val);


% dia = load("mat_d.mat"); matd = dia.Jd;
% low = load("mat_l.mat"); matl = low.Jlo;
% upp = load("mat_u.mat"); matu = upp.Jup;
% 
% soll = zeros(N, N, K);
% solu = zeros(N, N, K);
% sold = zeros(N, N, K);
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