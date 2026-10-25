function [intg] = f_aux3(basis, mesh, Fh, Gh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
K   = mesh.K;
det = mesh.det;
V0T = mesh.V0T;
xdc = mesh.xdc;
dh  = basis.dh;
fh  = basis.fh;
QsP = basis.QsP;
W_  = basis.W_;
N   = basis.N;
R   = basis.R;
Km1 = K-1;
NN  = N*N;
RK  = R*K;
Xq  = zeros(K, R);
for i = 1:K
    for r = 1:R
        Xq(i, r) = det(i, 1).*QsP(r, 1)+V0T(i, 1);
    end
end
%--------------------------------------------------------------------------
bf  = zeros(R, N);
bd  = zeros(R, N);
for i = 1:N
    bd(:, i) = dh{1, i}(QsP);
    bf(:, i) = fh{1, i}(QsP);
end
bd_        = bd.';
bf_        = bf.';
Wbd        = W_.*bd;
Wbf        = W_.*bf;
Wbd_       = Wbd.';
Wbf_       = Wbf.';
mass       = bf_*Wbf;
Wbd_m      = Wbd/mass;
Wbf_m      = Wbf/mass;
intg.R     = R;
intg.Xq    = Xq;
intg.bd_   = bd_;
intg.bf_   = bf_;
intg.Wbf   = Wbf;
intg.Wbd_m = Wbd_m;
intg.Wbf_m = Wbf_m;
%--------------------------------------------------------------------------
F_dof      = inittype(1, Fh, xdc, Xq, Wbf_m);
G_dof      = inittype(1, Gh, xdc, Xq, Wbf_m);
G          = G_dof*bf_;
Gkdof      = G*Wbf.*det;
visu       = sum(Wbf_, 2)./sum(W_, 1);
intg.F_dof = F_dof;
intg.Gkdof = Gkdof;
intg.visu  = visu;
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
bfl_   = bfl.';
bfr_   = bfr.';
bdl_m  = bdl/mass;
bdr_m  = bdr/mass;
bfl_m  = bfl/mass;
bfr_m  = bfr/mass;
bDl_M  = zeros(K, N);
bDr_M  = zeros(K, N);
bfl_M  = zeros(K, N);
bfr_M  = zeros(K, N);
for i  = 1:K
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
Lll_m_i = 1./2.*ffll_m; % +      o->
Llr_m_i =-1./2.*fflr_m; % -      o->
Lrl_m_i = 1./2.*ffrl_m; % + = -- <-o
Lrr_m_i =-1./2.*ffrr_m; % -      <-o
Lll_mb2 =       ffll_m; % +      o->
Llr_mb2 =   -mass\bfl_; % -      o->
Lrl_mb1 =    mass\bfr_; % + = -- <-o
Lrr_mb1 =      -ffrr_m; % -      <-o

% Lll_m_i =-Lll_m_i; % (w/ minus sign: R([u]) = -{v}[u])
% Llr_m_i =-1./2.*fflr_m; % -      o->
% Lrl_m_i = 1./2.*ffrl_m; % + = -- <-o
% Lrr_m_i =-1./2.*ffrr_m; % -      <-o
% Lll_mb2 =       ffll_m; % +      o->
% Llr_mb2 =   -mass\bfl_; % -      o->
% Lrl_mb1 =    mass\bfr_; % + = -- <-o
% Lrr_mb1 =      -ffrr_m; % -      <-o

%--------------------------------------------------------------------------
bdl_             = bdl.';
bdr_             = bdr.';
bfl_             = bfl.';
bfr_             = bfr.';
bDl              = bdl./det;
bDr              = bdr./det;
bDl_             = bDl.';
bDr_             = bDr.';
bflLll_i_m       = bfl*Lll_m_i;
bflLlr_i_m       = bfl*Llr_m_i;
bfrLrl_i_m       = bfr*Lrl_m_i;
bfrLrr_i_m       = bfr*Lrr_m_i;
bflLll_i_M       = zeros(Km1, N);
bflLlr_i_M       = zeros(Km1, N);
bfrLrl_i_M       = zeros(Km1, N);
bfrLrr_i_M       = zeros(Km1, N);
for i = 1:Km1
    bflLll_i_M(i, :) = bflLll_i_m./det(i  , 1); % alLllul (1)
    bflLlr_i_M(i, :) = bflLlr_i_m./det(i  , 1); % alLlrur (3)
    bfrLrl_i_M(i, :) = bfrLrl_i_m./det(i+1, 1); % arLrlul (2)
    bfrLrr_i_M(i, :) = bfrLrr_i_m./det(i+1, 1); % arLrrur (4)
end
bflLll_i_Mt      = bflLll_i_M.';
bfrLrl_i_Mt      = bfrLrl_i_M.';
bflLlr_i_Mt      = bflLlr_i_M.';
bfrLrr_i_Mt      = bfrLrr_i_M.';
bfrLrlb1_m       = bfr*Lrl_mb1;
bfrLrrb1_m       = bfr*Lrr_mb1;
bflLllb2_m       = bfl*Lll_mb2;
bflLlrb2_m       = bfl*Llr_mb2;
bfrLrlb1_M       = bfrLrlb1_m./det(1, 1);
bfrLrrb1_M       = bfrLrrb1_m./det(1, 1);
bflLllb2_M       = bflLllb2_m./det(K, 1);
bflLlrb2_M       = bflLlrb2_m./det(K, 1);
bflLllb2_Mt      = bflLllb2_M.';
bflLlrb2_Mt      = bflLlrb2_M.';
bfrLrlb1_Mt      = bfrLrlb1_M.';
bfrLrrb1_Mt      = bfrLrrb1_M.';
eta              = 2;
intg.bDl         = bDl;
intg.bDr         = bDr;
intg.bfl         = bfl;
intg.bfr         = bfr;
intg.bdl_        = bdl_;
intg.bdr_        = bdr_;
intg.bfl_        = bfl_;
intg.bfr_        = bfr_;
intg.bDl_        = bDl_;
intg.bDr_        = bDr_;
intg.bDl_M       = bDl_M;
intg.bDr_M       = bDr_M;
intg.bfl_M       = bfl_M;
intg.bfr_M       = bfr_M;
intg.eta         = eta;
intg.bflLll_i_Mt = bflLll_i_Mt;
intg.bflLlr_i_Mt = bflLlr_i_Mt;
intg.bfrLrl_i_Mt = bfrLrl_i_Mt;
intg.bfrLrr_i_Mt = bfrLrr_i_Mt;
intg.bflLllb2_Mt = bflLllb2_Mt;
intg.bflLlrb2_Mt = bflLlrb2_Mt;
intg.bfrLrlb1_Mt = bfrLrlb1_Mt;
intg.bfrLrrb1_Mt = bfrLrrb1_Mt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jac.m:
%--------------------------------------------------------------------------
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
intg.DFr      = DFr;
intg.FDr      = FDr;
%--------------------------------------------------------------------------
%{
intg.ffll_m   = ffll_m; % DEBUG
intg.fflr_m   = fflr_m;
intg.ffrl_m   = ffrl_m;
intg.ffrr_m   = ffrr_m;
%}
bfffll_m      = bf*ffll_m;
bffflr_m      = bf*fflr_m;
bfffrl_m      = bf*ffrl_m;
bfffrr_m      = bf*ffrr_m;
bfLll_it      = 1./2.*bfffll_m.'; % +      o->
bfLlr_it      =-1./2.*bffflr_m.'; % -      o->
bfLrl_it      = 1./2.*bfffrl_m.'; % + = -- <-o
bfLrr_it      =-1./2.*bfffrr_m.'; % -      <-o
bfLllb2t      =       bfffll_m.'; % +      o->
bfLlrb2t      =      -bfl_m*bf_;  % -      o->
bfLrlb1t      =       bfr_m*bf_;  % + = -- <-o
bfLrrb1t      =      -bfffrr_m.'; % -      <-o
intg.bfLll_it = bfLll_it;
intg.bfLlr_it = bfLlr_it;
intg.bfLrl_it = bfLrl_it;
intg.bfLrr_it = bfLrr_it;
intg.bfLrrb1t = bfLrrb1t; % rr
intg.bfLrlb1t = bfLrlb1t; % rl
intg.bfLllb2t = bfLllb2t; % ll
intg.bfLlrb2t = bfLlrb2t; % lr
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
intg.DDr     = DDr;
%{
intg.mass    = mass; % DEBUG
%}
%--------------------------------------------------------------------------
FDll_i       = zeros(N, N, Km1);
FDlr_i       = zeros(N, N, Km1);
FDrl_i       = zeros(N, N, Km1);
FDrr_i       = zeros(N, N, Km1);
DFll_i       = zeros(N, N, Km1);
DFlr_i       = zeros(N, N, Km1);
DFrl_i       = zeros(N, N, Km1);
DFrr_i       = zeros(N, N, Km1);
PPll_li      = zeros(N, N, Km1);
PPll_ri      = zeros(N, N, Km1);
PPlr_li      = zeros(N, N, Km1);
PPlr_ri      = zeros(N, N, Km1);
PPrl_li      = zeros(N, N, Km1);
PPrl_ri      = zeros(N, N, Km1);
PPrr_li      = zeros(N, N, Km1);
PPrr_ri      = zeros(N, N, Km1);
for i = 1:Km1
    FDll_i (:, :, i) = bfl_*bDl(i  , :)./2;      % +      o-> (CON.)
    FDlr_i (:, :, i) = bfl_*bDr(i+1, :)./2;      % +      o->
    FDrl_i (:, :, i) =-bfr_*bDl(i  , :)./2;      % -      <-o
    FDrr_i (:, :, i) =-bfr_*bDr(i+1, :)./2;      % -      <-o
    DFll_i (:, :, i) = bDl_(:, i  )./2*bfl;      % +      o-> (SYM.)
    DFlr_i (:, :, i) =-bDl_(:, i  )./2*bfr;      % -      <-o
    DFrl_i (:, :, i) = bDr_(:, i+1)./2*bfl;      % + = -- <-o
    DFrr_i (:, :, i) =-bDr_(:, i+1)./2*bfr;      % -      <-o
    PPll_li(:, :, i) =-bfl_*bflLll_i_M(i, :)./2; % +      o-> (PEN.)
    PPll_ri(:, :, i) =-bfl_*bfrLrl_i_M(i, :)./2; % +      o-> (w/ minus sign: R([u]) = -{v}[u])
    PPlr_li(:, :, i) =-bfl_*bflLlr_i_M(i, :)./2; % +      o-> (=)
    PPlr_ri(:, :, i) =-bfl_*bfrLrr_i_M(i, :)./2; % +      o-> (=)
    PPrl_li(:, :, i) = bfr_*bflLll_i_M(i, :)./2; % -      <-o (=)
    PPrl_ri(:, :, i) = bfr_*bfrLrl_i_M(i, :)./2; % -      <-o (=)
    PPrr_li(:, :, i) = bfr_*bflLlr_i_M(i, :)./2; % -      <-o (=)
    PPrr_ri(:, :, i) = bfr_*bfrLrr_i_M(i, :)./2; % -      <-o (=)
end
FDb1         =-bfr_*bDr(1, :)./2; % - <-o      (CON.)
FDb2         = bfl_*bDl(K, :)./2; % + o->
DFb1         =-bDr_(:, 1)*bfr;    % - <-o      (SYM.)
DFb2         = bDl_(:, K)*bfl;    % + o->
PPb1         = bfr_*bfrLrrb1_M;   % + = -- <-o (PEN.) (w/ minus sign: R([u]) = -{v}[u])
PPb2         =-bfl_*bflLllb2_M;   % -      o->        (=)
Jmll_li      = PPll_li.*eta+FDll_i+DFll_i;
Jmll_ri      = PPll_ri.*eta;
Jmlr_li      = PPlr_li.*eta+DFlr_i;
Jmlr_ri      = PPlr_ri.*eta+FDlr_i;
Jmrl_li      = PPrl_li.*eta+FDrl_i;
Jmrl_ri      = PPrl_ri.*eta+DFrl_i;
Jmrr_li      = PPrr_li.*eta;
Jmrr_ri      = PPrr_ri.*eta+FDrr_i+DFrr_i;
Jmb1         = FDb1+DFb1+PPb1.*eta;
Jmb2         = FDb2+DFb2+PPb2.*eta;
Jmb1t        = Jmb1.';
Jmb2t        = Jmb2.';
Jmb1t_1      = Jmb1t(:);
Jmb2t_1      = Jmb2t(:);
intg.Jmll_li = reshape(pagetranspose(Jmll_li), [NN, Km1]);
intg.Jmll_ri = reshape(pagetranspose(Jmll_ri), [NN, Km1]);
intg.Jmlr_li = reshape(pagetranspose(Jmlr_li), [NN, Km1]);
intg.Jmlr_ri = reshape(pagetranspose(Jmlr_ri), [NN, Km1]);
intg.Jmrl_li = reshape(pagetranspose(Jmrl_li), [NN, Km1]);
intg.Jmrl_ri = reshape(pagetranspose(Jmrl_ri), [NN, Km1]);
intg.Jmrr_li = reshape(pagetranspose(Jmrr_li), [NN, Km1]);
intg.Jmrr_ri = reshape(pagetranspose(Jmrr_ri), [NN, Km1]);
intg.Jmb1    = Jmb1t_1;
intg.Jmb2    = Jmb2t_1;
%--------------------------------------------------------------------------
ffll_to2     = ffll.'./2;
fflr_to2     = fflr.'./2;
ffrl_to2     = ffrl.'./2;
ffrr_to2     = ffrr.'./2;
ffll_to2_1   = ffll_to2(:);
fflr_to2_1   = fflr_to2(:);
ffrl_to2_1   = ffrl_to2(:);
ffrr_to2_1   = ffrr_to2(:);
Xdull_i      = PPll_li.*eta+DFll_i;
Xdulr_i      = PPlr_ri.*eta;
Xdurl_i      = PPrl_li.*eta;
Xdurr_i      = PPrr_ri.*eta+DFrr_i;
Xdub1        = DFb1+PPb1.*eta;
Xdub2        = DFb2+PPb2.*eta;
Xdub1t       = Xdub1.';
Xdub2t       = Xdub2.';
Xdub1t_1     = Xdub1t(:);
Xdub2t_1     = Xdub2t(:);
intg.FFll    = ffll_to2_1; % + o->
intg.FFlr    = fflr_to2_1; % + o->
intg.FFrl    =-ffrl_to2_1; % - <-o
intg.FFrr    =-ffrr_to2_1; % - <-o
intg.Xdull_i = reshape(pagetranspose(Xdull_i), [NN, Km1]);
intg.Xdulr_i = reshape(pagetranspose(Xdulr_i), [NN, Km1]);
intg.Xdurl_i = reshape(pagetranspose(Xdurl_i), [NN, Km1]);
intg.Xdurr_i = reshape(pagetranspose(Xdurr_i), [NN, Km1]);
intg.FFb1    =-ffrr_to2_1; % - <-o
intg.FFb2    = ffll_to2_1; % + o->
intg.Xdub1   = Xdub1t_1;
intg.Xdub2   = Xdub2t_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end