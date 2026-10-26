function [g] = func_aux(p, Xvt, nu, Fh, Gh, Uh, dUh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASIS:
%--------------------------------------------------------------------------
N  = p+1;
dh = cell(1, N);
fh = cell(1, N);
switch p
    case 0
        xp = [1/2]; %#ok<NBRAK2>
    case 1
        xp = [0, 1];
    case 2
        xp = [0, 1, 1/2];
    case 3
        xp = [0, 1, 1/3, 2/3];
    case 4
        xp = [0, 1, 1/4, 1/2, 3/4];
    otherwise
        return;
end
%--------------------------------------------------------------------------
syms x;
for i = 1:N
    fi = sym(1);
    for j = 1:N
        if j ~= i
            fi = fi.*(x-xp(1, j))./(xp(1, i)-xp(1, j));
        end
    end
    di = simplify(diff(fi, 1, x));
    if di == 0
        dh{1, i} = @(x) zeros(size(x));
        fh{1, i} = @(x)  ones(size(x));
    else
        if diff(di, 1, x) == 0
            dh{1, i} = @(x) ones(size(x)).*double(di);
        else
            dh{1, i} = matlabFunction(di, 'Vars', x);
        end
        fh{1, i} = matlabFunction(fi, 'Vars', x);
    end
end
%--------------------------------------------------------------------------
Qs   = quadGaussLegendre(2*p+2, 'Domain', [0, 1]);
QsP  = Qs.Points;
W_   = Qs.Weights;
R    = size(Qs.Weights, 1);
g.dh = dh;
g.fh = fh;
g.W_ = W_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH:
%--------------------------------------------------------------------------
NN   = N*N;
Xv   = Xvt.';
Ke   = size(Xv, 1);
Kem1 = Ke-1;
K    = Kem1;
%{
Xc   = (Xv(2:Ke, 1)+Xv(1:Kem1, 1))./2;
%}
det  = (Xv(2:Ke, 1)-Xv(1:Kem1, 1));
V0K  = [1:Kem1; 2:Ke].';
V0V  = zeros(K, 2);
xdc  = zeros(K, N);
for k = 1:2
    V0V(:, k) = Xv(V0K(:, k), :);
end
for i = 1:N
    xdc(:, i) = det(:, 1).*xp(1, i)+V0V(:, 1);
end
%--------------------------------------------------------------------------
Km1          = K-1;
Km2          = K-2;
T            = 3*K-2;
Tm1          = T-1;
KN           = K*N;
TNN          = T*NN;
Km1NN        = Km1*NN;
im           = zeros(T  , 1);
jm           = zeros(T  , 1);
XT           = zeros(TNN, 2);
im(  1:2, 1) = [1  , 1];
jm(  1:2, 1) = [1  , 2];
im(Tm1:T, 1) = [K  , K];
jm(Tm1:T, 1) = [Km1, K];
for i = 1:Km2
    im(2+3*(i-1)+1:2+3*i, 1) = i+1;
    jm(2+3*(i-1)+1:2+3*i, 1) = i:i+2;
end
for i = 1:T
    XT((i-1)*NN+(1:NN), :) = ...
        [repelem((im(i, 1)-1)*N+(1:N), N); repmat((jm(i, 1)-1)*N+(1:N), 1, N)].';
end
XT_i   = XT(:, 1);
XT_j   = XT(:, 2);
IJd_   = reshape(((3+3*((1:K  ).'-2))*NN+(1:NN)).', [], 1); % diag.
IJor   = reshape(((2+3*((2:K  ).'-2))*NN+(1:NN)).', [], 1); % diag. lower (off-diag. blocks of r cell)
IJol   = reshape(((4+3*((1:Km1).'-2))*NN+(1:NN)).', [], 1); % diag. upper (off-diag. blocks of l cell)
IJdl   = IJd_( 0*NN+(1:Km1NN), 1);
IJdr   = IJd_( 1*NN+(1:Km1NN), 1);
IJb1   = IJd_( 0*NN+(1:NN)   , 1);
IJb2   = IJd_(Km1NN+(1:NN)   , 1);
%--------------------------------------------------------------------------
g.K    = K;
g.det  = det;
g.KN   = KN;
g.TNN  = TNN;
g.IJd_ = IJd_; % cell
g.IJdl = IJdl; % fint ll
g.IJdr = IJdr; % face rr
g.IJol = IJol; % face lr
g.IJor = IJor; % face rl
g.IJb1 = IJb1; % bnd1
g.IJb2 = IJb2; % bnd2
%{
g.Jaux = spinit(XT_i, XT_j, [KN, KN]);
g.Jac  = spinit(XT_i, XT_j, [KN, KN]);
%}
g.XT_i = XT_i;
g.XT_j = XT_j;
g.Jaux = sparse(XT_i, XT_j, zeros(TNN, 1), KN, KN);
g.Jac  = sparse(XT_i, XT_j, zeros(TNN, 1), KN, KN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
Xq = zeros(K, R);
for i = 1:K
    for r = 1:R
        Xq(i, r) = det(i, 1).*QsP(r, 1)+V0V(i, 1);
    end
end
bf = zeros(R, N);
bd = zeros(R, N);
for i = 1:N
    bd(:, i) = dh{1, i}(QsP);
    bf(:, i) = fh{1, i}(QsP);
end
bd_     = bd.';
bf_     = bf.';
Wbd     = W_.*bd;
Wbf     = W_.*bf;
Wbd_    = Wbd.';
Wbf_    = Wbf.';
mass    = bf_*Wbf;
Wbd_m   = Wbd/mass;
Wbf_m   = Wbf/mass;
F_dof   = inittype(1, Fh, xdc, Xq, Wbf_m);
G_dof   = inittype(1, Gh, xdc, Xq, Wbf_m);
G       = G_dof*bf_;
Gkdof   = G*Wbf.*det;
visu    = sum(Wbf_, 2)./sum(W_, 1);
%--------------------------------------------------------------------------
g.K     = K;
g.N     = N;
g.R     = R;
g.Xq    = Xq;
g.bd_   = bd_;
g.bf_   = bf_;
g.Wbf   = Wbf;
g.Wbd_m = Wbd_m;
g.Wbf_m = Wbf_m;
g.F_dof = F_dof;
g.Gkdof = Gkdof;
g.visu  = visu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
bdl = zeros(1, N);
bdr = zeros(1, N);
bfl = zeros(1, N);
bfr = zeros(1, N);
for i = 1:N
    bdl(1, i) = dh{1, i}(1);
    bdr(1, i) = dh{1, i}(0);
    bfr(1, i) = fh{1, i}(0);
    bfl(1, i) = fh{1, i}(1);
end
bfl_  = bfl.';
bfr_  = bfr.';
bdl_m = bdl/mass;
bdr_m = bdr/mass;
bfl_m = bfl/mass;
bfr_m = bfr/mass;
bDl_M = zeros(K, N);
bDr_M = zeros(K, N);
bfl_M = zeros(K, N);
bfr_M = zeros(K, N);
for i = 1:K
    bDl_M(i, :) = bdl_m./det(i, 1).^2;
    bDr_M(i, :) = bdr_m./det(i, 1).^2;
    bfl_M(i, :) = bfl_m./det(i, 1);
    bfr_M(i, :) = bfr_m./det(i, 1);
end
ffll    = bfl_*bfl;
fflr    = bfl_*bfr;
ffrl    = bfr_*bfl;
ffrr    = bfr_*bfr;
ffll_m  = mass\ffll;
fflr_m  = mass\fflr;
ffrl_m  = mass\ffrl;
ffrr_m  = mass\ffrr;
Lll_m_i = ffll_m./2; % +      o->
Llr_m_i =-fflr_m./2; % -      o->
Lrl_m_i = ffrl_m./2; % + = -- <-o
Lrr_m_i =-ffrr_m./2; % -      <-o
Lll_mb2 = ffll_m;    % +      o->
Llr_mb2 =-mass\bfl_; % -      o->
Lrl_mb1 = mass\bfr_; % + = -- <-o
Lrr_mb1 =-ffrr_m;    % -      <-o
Lll_m_i =-Lll_m_i;   % (w/ minus sign: R([u]) = -{v}[u])
Llr_m_i =-Llr_m_i;   % =
Lrl_m_i =-Lrl_m_i;   % =
Lrr_m_i =-Lrr_m_i;   % =
Lrl_mb1 =-Lrl_mb1;   % =
Lrr_mb1 =-Lrr_mb1;   % =
Llr_mb2 =-Llr_mb2;   % =
Lll_mb2 =-Lll_mb2;   % =
%--------------------------------------------------------------------------
bdl_          = bdl.';
bdr_          = bdr.';
bfl_          = bfl.';
bfr_          = bfr.';
bDl           = bdl./det;
bDr           = bdr./det;
bDl_          = bDl.';
bDr_          = bDr.';
bflLll_i_m    = bfl*Lll_m_i;
bflLlr_i_m    = bfl*Llr_m_i;
bfrLrl_i_m    = bfr*Lrl_m_i;
bfrLrr_i_m    = bfr*Lrr_m_i;
bflLll_i_M    = zeros(Km1, N);
bflLlr_i_M    = zeros(Km1, N);
bfrLrl_i_M    = zeros(Km1, N);
bfrLrr_i_M    = zeros(Km1, N);
for i = 1:Km1
    bflLll_i_M(i, :) = bflLll_i_m./det(i  , 1); % alLllul (1)
    bflLlr_i_M(i, :) = bflLlr_i_m./det(i  , 1); % alLlrur (3)
    bfrLrl_i_M(i, :) = bfrLrl_i_m./det(i+1, 1); % arLrlul (2)
    bfrLrr_i_M(i, :) = bfrLrr_i_m./det(i+1, 1); % arLrrur (4)
end
bflLll_i_Mt   = bflLll_i_M.';
bfrLrl_i_Mt   = bfrLrl_i_M.';
bflLlr_i_Mt   = bflLlr_i_M.';
bfrLrr_i_Mt   = bfrLrr_i_M.';
bfrLrlb1_m    = bfr*Lrl_mb1;
bfrLrrb1_m    = bfr*Lrr_mb1;
bflLllb2_m    = bfl*Lll_mb2;
bflLlrb2_m    = bfl*Llr_mb2;
bfrLrlb1_M    = bfrLrlb1_m./det(1, 1);
bfrLrrb1_M    = bfrLrrb1_m./det(1, 1);
bflLllb2_M    = bflLllb2_m./det(K, 1);
bflLlrb2_M    = bflLlrb2_m./det(K, 1);
bflLllb2_Mt   = bflLllb2_M.';
bflLlrb2_Mt   = bflLlrb2_M.';
bfrLrlb1_Mt   = bfrLrlb1_M.';
bfrLrrb1_Mt   = bfrLrrb1_M.';
eta           = 2;
g.bDl         = bDl;
g.bDr         = bDr;
g.bfl         = bfl;
g.bfr         = bfr;
g.bdl_        = bdl_;
g.bdr_        = bdr_;
g.bfl_        = bfl_;
g.bfr_        = bfr_;
g.bDl_        = bDl_;
g.bDr_        = bDr_;
g.bDl_M       = bDl_M;
g.bDr_M       = bDr_M;
g.bfl_M       = bfl_M;
g.bfr_M       = bfr_M;
g.eta         = eta;
g.bflLll_i_Mt = bflLll_i_Mt;
g.bflLlr_i_Mt = bflLlr_i_Mt;
g.bfrLrl_i_Mt = bfrLrl_i_Mt;
g.bfrLrr_i_Mt = bfrLrr_i_Mt;
g.bflLllb2_Mt = bflLllb2_Mt;
g.bflLlrb2_Mt = bflLlrb2_Mt;
g.bfrLrlb1_Mt = bfrLrlb1_Mt;
g.bfrLrrb1_Mt = bfrLrrb1_Mt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jac.m:
%--------------------------------------------------------------------------
RK  = R*K;
dfr = zeros(N, N, R);
fdr = zeros(N, N, R);
DFr = zeros(RK, NN);
FDr = zeros(RK, NN);
for r = 1:R
    dfr(:, :, r) = Wbd_(:, r)*bf(r, :); % wd_*f
    fdr(:, :, r) = Wbf_(:, r)*bd(r, :); % wf_*d
end
for i = 1:K
    for r = 1:R
        DFr((i-1)*R+r, :) = reshape(dfr(:, :, r).', NN, []);
        FDr((i-1)*R+r, :) = reshape(fdr(:, :, r).', NN, []);
    end
end
g.DFr      = DFr;
g.FDr      = FDr;
%--------------------------------------------------------------------------
g.ffll_m   = ffll_m; % DEBUG
g.fflr_m   = fflr_m; % =
g.ffrl_m   = ffrl_m; % =
g.ffrr_m   = ffrr_m; % =
bfffll_m   = bf*ffll_m;
bffflr_m   = bf*fflr_m;
bfffrl_m   = bf*ffrl_m;
bfffrr_m   = bf*ffrr_m;
bfLll_it   = 1./2.*bfffll_m.'; % +      o->
bfLlr_it   =-1./2.*bffflr_m.'; % -      o->
bfLrl_it   = 1./2.*bfffrl_m.'; % + = -- <-o
bfLrr_it   =-1./2.*bfffrr_m.'; % -      <-o
bfLllb2t   =       bfffll_m.'; % +      o->
bfLlrb2t   =      -bfl_m*bf_;  % -      o->
bfLrlb1t   =       bfr_m*bf_;  % + = -- <-o
bfLrrb1t   =      -bfffrr_m.'; % -      <-o
g.bfLll_it = bfLll_it;
g.bfLlr_it = bfLlr_it;
g.bfLrl_it = bfLrl_it;
g.bfLrr_it = bfLrr_it;
g.bfLrrb1t = bfLrrb1t; % rr
g.bfLrlb1t = bfLrlb1t; % rl
g.bfLllb2t = bfLllb2t; % ll
g.bfLlrb2t = bfLlrb2t; % lr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% func_dAdUjac.m:
%--------------------------------------------------------------------------
DDr = zeros(RK, NN);
ddr = zeros(N, N, R);
for r = 1:R
    ddr(:, :, r) = Wbd_(:, r)*bd(r, :); % wd_*d
end
for i = 1:K
    for r = 1:R
        DDr((i-1)*R+r, :) = reshape(ddr(:, :, r).', NN, [])./det(i, 1);
    end
end
g.DDr      = DDr;
g.mass     = mass; % DEBUG
%--------------------------------------------------------------------------
FDll_i     = zeros(N, N, Km1);
FDlr_i     = zeros(N, N, Km1);
FDrl_i     = zeros(N, N, Km1);
FDrr_i     = zeros(N, N, Km1);
DFll_i     = zeros(N, N, Km1);
DFlr_i     = zeros(N, N, Km1);
DFrl_i     = zeros(N, N, Km1);
DFrr_i     = zeros(N, N, Km1);
PPll_li    = zeros(N, N, Km1);
PPll_ri    = zeros(N, N, Km1);
PPlr_li    = zeros(N, N, Km1);
PPlr_ri    = zeros(N, N, Km1);
PPrl_li    = zeros(N, N, Km1);
PPrl_ri    = zeros(N, N, Km1);
PPrr_li    = zeros(N, N, Km1);
PPrr_ri    = zeros(N, N, Km1);
for i = 1:Km1
    FDll_i (:, :, i) = bfl_*bDl(i  , :)./2;      % +      o-> (CON.)
    FDlr_i (:, :, i) = bfl_*bDr(i+1, :)./2;      % +      o->
    FDrl_i (:, :, i) =-bfr_*bDl(i  , :)./2;      % -      <-o
    FDrr_i (:, :, i) =-bfr_*bDr(i+1, :)./2;      % -      <-o
    DFll_i (:, :, i) = bDl_(:, i  )./2*bfl;      % +      o-> (SYM.)
    DFlr_i (:, :, i) =-bDl_(:, i  )./2*bfr;      % -      <-o
    DFrl_i (:, :, i) = bDr_(:, i+1)./2*bfl;      % + = -- <-o
    DFrr_i (:, :, i) =-bDr_(:, i+1)./2*bfr;      % -      <-o
    PPll_li(:, :, i) = bfl_*bflLll_i_M(i, :)./2; % +      o-> (PEN.)
    PPll_ri(:, :, i) = bfl_*bfrLrl_i_M(i, :)./2; % +      o->
    PPlr_li(:, :, i) = bfl_*bflLlr_i_M(i, :)./2; % +      o->
    PPlr_ri(:, :, i) = bfl_*bfrLrr_i_M(i, :)./2; % +      o->
    PPrl_li(:, :, i) =-bfr_*bflLll_i_M(i, :)./2; % -      <-o
    PPrl_ri(:, :, i) =-bfr_*bfrLrl_i_M(i, :)./2; % -      <-o
    PPrr_li(:, :, i) =-bfr_*bflLlr_i_M(i, :)./2; % -      <-o
    PPrr_ri(:, :, i) =-bfr_*bfrLrr_i_M(i, :)./2; % -      <-o
end
FDb1       =-bfr_*bDr(1, :)./2; % - <-o (CON.)
FDb2       = bfl_*bDl(K, :)./2; % + o->
DFb1       =-bDr_(:, 1)*bfr;    % - <-o (SYM.)
DFb2       = bDl_(:, K)*bfl;    % + o->
PPb1       =-bfr_*bfrLrrb1_M;   % - <-o (PEN.)
PPb2       = bfl_*bflLllb2_M;   % + o->
Jmll_li    = PPll_li.*eta+FDll_i+DFll_i;
Jmll_ri    = PPll_ri.*eta;
Jmlr_li    = PPlr_li.*eta+DFlr_i;
Jmlr_ri    = PPlr_ri.*eta+FDlr_i;
Jmrl_li    = PPrl_li.*eta+FDrl_i;
Jmrl_ri    = PPrl_ri.*eta+DFrl_i;
Jmrr_li    = PPrr_li.*eta;
Jmrr_ri    = PPrr_ri.*eta+FDrr_i+DFrr_i;
Jmb1       = FDb1+DFb1+PPb1.*eta;
Jmb2       = FDb2+DFb2+PPb2.*eta;
Jmb1t      = Jmb1.';
Jmb2t      = Jmb2.';
Jmb1t_1    = Jmb1t(:);
Jmb2t_1    = Jmb2t(:);
g.Jmll_li  = reshape(pagetranspose(Jmll_li), [NN, Km1]);
g.Jmll_ri  = reshape(pagetranspose(Jmll_ri), [NN, Km1]);
g.Jmlr_li  = reshape(pagetranspose(Jmlr_li), [NN, Km1]);
g.Jmlr_ri  = reshape(pagetranspose(Jmlr_ri), [NN, Km1]);
g.Jmrl_li  = reshape(pagetranspose(Jmrl_li), [NN, Km1]);
g.Jmrl_ri  = reshape(pagetranspose(Jmrl_ri), [NN, Km1]);
g.Jmrr_li  = reshape(pagetranspose(Jmrr_li), [NN, Km1]);
g.Jmrr_ri  = reshape(pagetranspose(Jmrr_ri), [NN, Km1]);
g.Jmb1     = Jmb1t_1;
g.Jmb2     = Jmb2t_1;
%--------------------------------------------------------------------------
ffll_to2   = ffll.'./2;
fflr_to2   = fflr.'./2;
ffrl_to2   = ffrl.'./2;
ffrr_to2   = ffrr.'./2;
ffll_to2_1 = ffll_to2(:);
fflr_to2_1 = fflr_to2(:);
ffrl_to2_1 = ffrl_to2(:);
ffrr_to2_1 = ffrr_to2(:);
Xdull_i    = PPll_li.*eta+DFll_i;
Xdulr_i    = PPlr_ri.*eta;
Xdurl_i    = PPrl_li.*eta;
Xdurr_i    = PPrr_ri.*eta+DFrr_i;
Xdub1      = DFb1+PPb1.*eta;
Xdub2      = DFb2+PPb2.*eta;
Xdub1t     = Xdub1.';
Xdub2t     = Xdub2.';
Xdub1t_1   = Xdub1t(:);
Xdub2t_1   = Xdub2t(:);
g.FFll     = ffll_to2_1; % + o->
g.FFlr     = fflr_to2_1; % + o->
g.FFrl     =-ffrl_to2_1; % - <-o
g.FFrr     =-ffrr_to2_1; % - <-o
g.Xdull_i  = reshape(pagetranspose(Xdull_i), [NN, Km1]);
g.Xdulr_i  = reshape(pagetranspose(Xdulr_i), [NN, Km1]);
g.Xdurl_i  = reshape(pagetranspose(Xdurl_i), [NN, Km1]);
g.Xdurr_i  = reshape(pagetranspose(Xdurr_i), [NN, Km1]);
g.FFb1     =-ffrr_to2_1; % - <-o
g.FFb2     = ffll_to2_1; % + o->
g.Xdub1    = Xdub1t_1;
g.Xdub2    = Xdub2t_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMON VARS:
%--------------------------------------------------------------------------
Xl     = Xv( 1, 1);
Xr     = Xv(Ke, 1);
ulb1   =  Uh(Xl);
urb2   =  Uh(Xr);
dlb1   = dUh(Xl);
drb2   = dUh(Xr);
alb1   = nu+ulb1;
arb2   = nu+urb2;
g.dlb1 = dlb1;
g.drb2 = drb2;
g.ulb1 = ulb1;
g.urb2 = urb2;
g.alb1 = alb1;
g.arb2 = arb2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end