function [intg] = f_aux3(basis, mesh, Fh)
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
Xq  = zeros(K, R);
for i = 1:K
    for r = 1:R
        Xq(i, r) = det(i, 1).*QsP(r, 1)+V0T(i, 1);
    end
end
%--------------------------------------------------------------------------
bf    = zeros(R, N);
bd    = zeros(R, N);
for i = 1:N
    bd(:, i) = dh{1, i}(QsP);
    bf(:, i) = fh{1, i}(QsP);
end
bd_   = bd.';
bf_   = bf.';
Wbd   = W_.*bd;
Wbf   = W_.*bf;
Wbd_  = Wbd.';
Wbf_  = Wbf.';
mass  = bf_*Wbf;
Wbd_m = Wbd/mass;
Wbf_m = Wbf/mass;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL JACOBIAN:
%--------------------------------------------------------------------------
ddr = zeros(N, N, R);
dfr = zeros(N, N, R);
fdr = zeros(N, N, R);
NN  = N*N;
RK  = R*K;
DDr = zeros(RK, NN);
DFr = zeros(RK, NN);
FDr = zeros(RK, NN);
for r = 1:R
    ddr(:, :, r) = Wbd_(:, r)*bd(r, :); % wd_*d
    dfr(:, :, r) = Wbd_(:, r)*bf(r, :); % wd_*f
    fdr(:, :, r) = Wbf_(:, r)*bd(r, :); % wf_*d
end
for i = 1:K
    for r = 1:R
        DDr((i-1)*R+r, :) = reshape(ddr(:, :, r).', NN, [])./det(i, 1);
        DFr((i-1)*R+r, :) = reshape(dfr(:, :, r).', NN, []);
        FDr((i-1)*R+r, :) = reshape(fdr(:, :, r).', NN, []);
    end
end
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
ffll   = bfl_*bfl;
fflr   = bfl_*bfr;
ffrl   = bfr_*bfl;
ffrr   = bfr_*bfr;
ffll_m = mass\ffll;
fflr_m = mass\fflr;
ffrl_m = mass\ffrl;
ffrr_m = mass\ffrr;
%{
ffll   = zeros(N, N);
fflr   = zeros(N, N);
ffrl   = zeros(N, N);
ffrr   = zeros(N, N);
for j = 1:N
    ffll(:, j) = bfl.*bfl(:, j);
    fflr(:, j) = bfl.*bfr(:, j);
    ffrl(:, j) = bfr.*bfl(:, j);
    ffrr(:, j) = bfr.*bfr(:, j);
end
%}
%--------------------------------------------------------------------------
Lll_m_i    = 1./2.*ffll_m; % +      o->
Llr_m_i    =-1./2.*fflr_m; % -      o->
Lrl_m_i    = 1./2.*ffrl_m; % + = -- <-o
Lrr_m_i    =-1./2.*ffrr_m; % -      <-o
Lll_mb2    =       ffll_m; % +      o->
Llr_mb2    =   -mass\bfl_; % -      o->
Lrl_mb1    =    mass\bfr_; % + = -- <-o
Lrr_mb1    =      -ffrr_m; % -      <-o
bflLll_i_m = bfl*Lll_m_i;
bflLlr_i_m = bfl*Llr_m_i;
bfrLrl_i_m = bfr*Lrl_m_i;
bfrLrr_i_m = bfr*Lrr_m_i;
bflLllb2_m = bfl*Lll_mb2;
bflLlrb2_m = bfl*Llr_mb2;
bfrLrlb1_m = bfr*Lrl_mb1;
bfrLrrb1_m = bfr*Lrr_mb1;
bflLll_i_M = zeros(Km1, N);
bflLlr_i_M = zeros(Km1, N);
bfrLrl_i_M = zeros(Km1, N);
bfrLrr_i_M = zeros(Km1, N);
for i = 1:Km1
    bflLll_i_M(i, :) = bflLll_i_m./det(i  , 1); % alLllul (1)
    bflLlr_i_M(i, :) = bflLlr_i_m./det(i  , 1); % alLlrur (3)
    bfrLrl_i_M(i, :) = bfrLrl_i_m./det(i+1, 1); % arLrlul (2)
    bfrLrr_i_M(i, :) = bfrLrr_i_m./det(i+1, 1); % arLrrur (4)
end
bflLllb2_M = bflLllb2_m./det(K, 1);
bflLlrb2_M = bflLlrb2_m./det(K, 1);
bfrLrlb1_M = bfrLrlb1_m./det(1, 1);
bfrLrrb1_M = bfrLrrb1_m./det(1, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL/FACE JACOBIAN (det. already factored in!):
%--------------------------------------------------------------------------
bfffll_m  = bf*ffll_m;
bffflr_m  = bf*fflr_m;
bfffrl_m  = bf*ffrl_m;
bfffrr_m  = bf*ffrr_m;
bfLll_it  = 1./2.*bfffll_m.'; % +      o->
bfLlr_it  =-1./2.*bffflr_m.'; % -      o->
bfLrl_it  = 1./2.*bfffrl_m.'; % + = -- <-o
bfLrr_it  =-1./2.*bfffrr_m.'; % -      <-o
bfLllb2t  =       bfffll_m.'; % +      o->
bfLlrb2t  =      -bfl_m*bf_;  % -      o->
bfLrlb1t  =       bfr_m*bf_;  % + = -- <-o
bfLrrb1t  =      -bfffrr_m.'; % -      <-o
%{
penll_l_i = zeros(N, N, Km1);
penll_r_i = zeros(N, N, Km1);
penlr_l_i = zeros(N, N, Km1);
penlr_r_i = zeros(N, N, Km1);
penrl_l_i = zeros(N, N, Km1);
penrl_r_i = zeros(N, N, Km1);
penrr_l_i = zeros(N, N, Km1);
penrr_r_i = zeros(N, N, Km1);
for i = 1:Km1
    penll_l_i(:, :, i) = bfl_*bflLll_i_M(i, :);
    penll_r_i(:, :, i) = bfl_*bfrLrl_i_M(i, :);
    penlr_l_i(:, :, i) = bfl_*bflLlr_i_M(i, :);
    penlr_r_i(:, :, i) = bfl_*bfrLrr_i_M(i, :);
    penrl_l_i(:, :, i) = bfr_*bflLll_i_M(i, :);
    penrl_r_i(:, :, i) = bfr_*bfrLrl_i_M(i, :);
    penrr_l_i(:, :, i) = bfr_*bflLlr_i_M(i, :);
    penrr_r_i(:, :, i) = bfr_*bfrLrr_i_M(i, :);
end
%}
bfl_i_M   = zeros(Km1, N);
bfr_i_M   = zeros(Km1, N);
for i = 1:Km1
    bfl_i_M(i, :) = bfl_m./det(i, 1);
    bfr_i_M(i, :) = bfr_m./det(i, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGN TO (CELL)...
%--------------------------------------------------------------------------
Wbf_             = Wbf.';
dc_m             = bd_*Wbf_m;
intg.R           = R;
intg.Xq          = Xq;
intg.bd          = bd;
intg.bf          = bf;
intg.bd_         = bd_;
intg.bf_         = bf_;
intg.Wbf         = Wbf;
intg.Wbf_        = Wbf_;
intg.mass        = mass;
intg.Wbd_m       = Wbd_m;
intg.Wbf_m       = Wbf_m;
intg.dc_m        = dc_m;
%--------------------------------------------------------------------------
intg.DDr         = DDr;
intg.DFr         = DFr;
intg.FDr         = FDr;
intg.visu        = sum(Wbf_, 2)./sum(W_, 1);
intg.F_dof       = inittype(1, Fh, xdc, Xq, Wbf_m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGN TO (FACE)...
%--------------------------------------------------------------------------
bdl_             = bdl.';
bdr_             = bdr.';
bfl_             = bfl.';
bfr_             = bfr.';
bDl              = bdl./det;
bDr              = bdr./det;
bDl_             = bDl.';
bDr_             = bDr.';
bflLll_i_Mt      = bflLll_i_M.';
bfrLrl_i_Mt      = bfrLrl_i_M.';
bflLlr_i_Mt      = bflLlr_i_M.';
bfrLrr_i_Mt      = bfrLrr_i_M.';
bflLllb2_Mt      = bflLllb2_M.';
bflLlrb2_Mt      = bflLlrb2_M.';
bfrLrlb1_Mt      = bfrLrlb1_M.';
bfrLrrb1_Mt      = bfrLrrb1_M.';
bfl_i_Mt         = bfl_i_M.';
bfr_i_Mt         = bfr_i_M.';
eta              = 2;
%--------------------------------------------------------------------------
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
intg.bfl_i_Mt    = bfl_i_Mt;
intg.bfr_i_Mt    = bfr_i_Mt;
%--------------------------------------------------------------------------
intg.ffll_m      = ffll_m;
intg.fflr_m      = fflr_m;
intg.ffrl_m      = ffrl_m;
intg.ffrr_m      = ffrr_m;
intg.bfLll_it    = bfLll_it;
intg.bfLlr_it    = bfLlr_it;
intg.bfLrl_it    = bfLrl_it;
intg.bfLrr_it    = bfLrr_it;
intg.bfLrrb1t    = bfLrrb1t; % rr
intg.bfLrlb1t    = bfLrlb1t; % rl
intg.bfLllb2t    = bfLllb2t; % ll
intg.bfLlrb2t    = bfLlrb2t; % lr
%{
intg.penll_l_i   = penll_l_i;
intg.penll_r_i   = penll_r_i;
intg.penlr_l_i   = penlr_l_i;
intg.penlr_r_i   = penlr_r_i;
intg.penrl_l_i   = penrl_l_i;
intg.penrl_r_i   = penrl_r_i;
intg.penrr_l_i   = penrr_l_i;
intg.penrr_r_i   = penrr_r_i;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end