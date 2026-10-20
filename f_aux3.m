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
xqc = zeros(K, R);
for i = 1:K
    for r = 1:R
        xqc(i, r) = det(i, 1).*QsP(r, 1)+V0T(i, 1);
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
ffll = zeros(N, N);
fflr = zeros(N, N);
ffrl = zeros(N, N);
ffrr = zeros(N, N);
for j = 1:N
    ffll(:, j) = bfl.*bfl(:, j);
    fflr(:, j) = bfl.*bfr(:, j);
    ffrl(:, j) = bfr.*bfl(:, j);
    ffrr(:, j) = bfr.*bfr(:, j);
end
%}
%--------------------------------------------------------------------------
Lll_m     = 1./2.*ffll_m; % +      o->
Llr_m     =-1./2.*fflr_m; % -      o->
Lrl_m     = 1./2.*ffrl_m; % + = -- <-o
Lrr_m     =-1./2.*ffrr_m; % -      <-o
bflLll_m  = bfl*Lll_m;
bflLlr_m  = bfl*Llr_m;
bfrLrl_m  = bfr*Lrl_m;
bfrLrr_m  = bfr*Lrr_m;
bflLll_Mi = zeros(Km1, N);
bflLlr_Mi = zeros(Km1, N);
bfrLrl_Mi = zeros(Km1, N);
bfrLrr_Mi = zeros(Km1, N);
%{
penl_ll   = zeros(Km1, N); % vl*(alLllul)/2
penl_lr   = zeros(Km1, N); % vl*(arLrlul)/2
penl_rl   = zeros(Km1, N); % vl*(alLlrur)/2
penl_rr   = zeros(Km1, N); % vl*(arLrrur)/2
penr_ll   = zeros(Km1, N); % vr*(alLllul)/2
penr_lr   = zeros(Km1, N); % vr*(arLrlul)/2
penr_rl   = zeros(Km1, N); % vr*(alLlrur)/2
penr_rr   = zeros(Km1, N); % vr*(arLrrur)/2
%}
for i = 1:Km1
    bflLll_Mi(i, :) = bflLll_m./det(i  , 1); % alLllul (1)
    bflLlr_Mi(i, :) = bflLlr_m./det(i  , 1); % alLlrur (3)
    bfrLrl_Mi(i, :) = bfrLrl_m./det(i+1, 1); % arLrlul (2)
    bfrLrr_Mi(i, :) = bfrLrr_m./det(i+1, 1); % arLrrur (4)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL/FACE JACOBIAN:
%--------------------------------------------------------------------------
bfffll_m = bf*ffll_m;
bffflr_m = bf*fflr_m;
bfffrl_m = bf*ffrl_m;
bfffrr_m = bf*ffrr_m;
bfrb1_M  = bfr_M(1, :);
bflb2_M  = bfl_M(K, :);
bfLll_it = 1./2.*bfffll_m.';  % +      o->
bfLlr_it =-1./2.*bffflr_m.';  % -      o->
bfLrl_it = 1./2.*bfffrl_m.';  % + = -- <-o
bfLrr_it =-1./2.*bfffrr_m.';  % -      <-o
bfLllb2t =       bfffll_m.';  % +      o->
bfLlrb2t =      -bflb2_M*bf_; % -      o->
bfLrlb1t =       bfrb1_M*bf_; % + = -- <-o % isto devia ser _m nao _M %%%%%%%%%%%%%%%%%%%%%%%%% CHECK
bfLrrb1t =      -bfffrr_m.';  % -      <-o
%{
bfLll_it = zeros(N, R, Km1);
bfLlr_it = zeros(N, R, Km1);
bfLrl_it = zeros(N, R, Km1);
bfLrr_it = zeros(N, R, Km1);
for i = 1:Km1
    bfLll_it(:, :, i) = 1./2.*bfffll_m.'./det(i  , 1); % +      o->
    bfLlr_it(:, :, i) =-1./2.*bffflr_m.'./det(i  , 1); % -      o->
    bfLrl_it(:, :, i) = 1./2.*bfffrl_m.'./det(i+1, 1); % + = -- <-o
    bfLrr_it(:, :, i) =-1./2.*bfffrr_m.'./det(i+1, 1); % -      <-o
end
bfLllb2t = bfffll_m.'./det(K, 1); % +      o->
bfLlrb2t =-bflb2_M*bf_;           % -      o->
bfLrlb1t = bfrb1_M*bf_;           % + = -- <-o
bfLrrb1t =-bfffrr_m.'./det(1, 1); % -      <-o
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGN TO (CELL)...
%--------------------------------------------------------------------------
Wbf_          = Wbf.';
dc_m          = bd_*Wbf_m;
intg.R        = R;
intg.xqc      = xqc;
intg.bd       = bd;
intg.bf       = bf;
intg.bd_      = bd_;
intg.bf_      = bf_;
intg.Wbf      = Wbf;
intg.Wbf_     = Wbf_;
intg.mass     = mass;
intg.Wbd_m    = Wbd_m;
intg.Wbf_m    = Wbf_m;
intg.dc_m     = dc_m;
%--------------------------------------------------------------------------
intg.DDr      = DDr;
intg.DFr      = DFr;
intg.FDr      = FDr;
intg.visu     = sum(Wbf_, 2)./sum(W_, 1);
intg.F_dof    = inittype(1, Fh, xdc, xqc, Wbf_m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGN TO (FACE)...
%--------------------------------------------------------------------------
bDl            = bdl./det;
bDr            = bdr./det;
bDl_           = bDl.';
bDr_           = bDr.';
bdl_           = bdl.';
bdr_           = bdr.';
bfl_           = bfl.';
bfr_           = bfr.';
eta            = 2;
intg.bDl       = bDl;
intg.bDr       = bDr;
intg.bfl       = bfl;
intg.bfr       = bfr;
intg.bDl_      = bDl_;
intg.bDr_      = bDr_;
intg.bdl_      = bdl_;
intg.bdr_      = bdr_;
intg.bfl_      = bfl_;
intg.bfr_      = bfr_;
intg.bDl_M     = bDl_M;
intg.bDr_M     = bDr_M;
intg.bfl_M     = bfl_M;
intg.bfr_M     = bfr_M;
intg.eta       = eta;
intg.bflLll_Mi = bflLll_Mi;
intg.bfrLrl_Mi = bfrLrl_Mi;
intg.bflLlr_Mi = bflLlr_Mi;
intg.bfrLrr_Mi = bfrLrr_Mi;
%--------------------------------------------------------------------------
%{
Lll_it         = zeros(N, N, Km1);
Llr_it         = zeros(N, N, Km1);
Lrl_it         = zeros(N, N, Km1);
Lrr_it         = zeros(N, N, Km1);
for i = 1:Km1
    Lll_it(:, :, i) = 1./2.*ffll_mt./det(i  , 1); % +      o->
    Llr_it(:, :, i) =-1./2.*fflr_mt./det(i  , 1); % -      o->
    Lrl_it(:, :, i) = 1./2.*ffrl_mt./det(i+1, 1); % + = -- <-o
    Lrr_it(:, :, i) =-1./2.*ffrr_mt./det(i+1, 1); % -      <-o
end
Lrrb1t         = ffrr_mt./det(1, 1); % + = -- <-o
Lllb2t         = ffll_mt./det(K, 1); % +      o->
TNN            = mesh.TNN;
IJdl           = mesh.IJdl; % ll
IJol           = mesh.IJol; % lr
IJor           = mesh.IJor; % rl
IJdr           = mesh.IJdr; % rr
IJb1           = mesh.IJb1; % rrb1
IJb2           = mesh.IJb2; % llb2
JAUX           = mesh.Jaux; % aux. matrix (for debugging)
VALL           = zeros(TNN, 1);
Lll_it_res     = reshape(Lll_it, [NN, Km1]);
Llr_it_res     = reshape(Llr_it, [NN, Km1]);
Lrl_it_res     = reshape(Lrl_it, [NN, Km1]);
Lrr_it_res     = reshape(Lrr_it, [NN, Km1]);
VALLdil_NN_1   = Lll_it_res(:);
VALLoil_NN_1   = Llr_it_res(:);
VALLoir_NN_1   = Lrl_it_res(:);
VALLdir_NN_1   = Lrr_it_res(:);
VALLdb1_NN_1   = Lrrb1t    (:);
VALLdb2_NN_1   = Lllb2t    (:);
VALL(IJdl)     = VALL(IJdl)-VALLdil_NN_1;
VALL(IJol)     =-VALLoil_NN_1;
VALL(IJor)     =-VALLoir_NN_1;
VALL(IJdr)     = VALL(IJdr)-VALLdir_NN_1;
VALL(IJb1)     = VALL(IJb1)-VALLdb1_NN_1;
VALL(IJb2)     = VALL(IJb2)-VALLdb2_NN_1;
MATL           = JAUX(VALL);
intg.ffll_mt   = ffll_mt; % DEBUG
intg.fflr_mt   = fflr_mt;
intg.ffrl_mt   = ffrl_mt;
intg.ffrr_mt   = ffrr_mt;
intg.MATL      = MATL;
%}
intg.bfLll_it  = bfLll_it;
intg.bfLlr_it  = bfLlr_it;
intg.bfLrl_it  = bfLrl_it;
intg.bfLrr_it  = bfLrr_it;
intg.bfLrrb1t  = bfLrrb1t; % rr
intg.bfLrlb1t  = bfLrlb1t; % rl
intg.bfLllb2t  = bfLllb2t; % ll
intg.bfLlrb2t  = bfLlrb2t; % lr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end