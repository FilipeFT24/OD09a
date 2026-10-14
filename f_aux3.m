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
mass  = bf_*Wbf;
Wbd_m = Wbd/mass;
Wbf_m = Wbf/mass;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL JAC.:
%--------------------------------------------------------------------------
ddr_m = zeros(N, N, R);
dfr_m = zeros(N, N, R);
fdr_m = zeros(N, N, R);
NN    = N*N;
RK    = R*K;
DDr_M = zeros(RK, NN);
DFr_M = zeros(RK, NN);
FDr_M = zeros(RK, NN);
for r = 1:R
    ddr_m(:, :, r) = Wbd_m(r, :).'*bd(r, :); % d_m*d
    dfr_m(:, :, r) = Wbd_m(r, :).'*bf(r, :); % d_m*f
    fdr_m(:, :, r) = Wbf_m(r, :).'*bd(r, :); % f_m*d
end
for i = 1:K
    for r = 1:R
        DDr_M((i-1)*R+r, :) = reshape(ddr_m(:, :, r).', NN, [])./det(i, 1).^3;
        DFr_M((i-1)*R+r, :) = reshape(dfr_m(:, :, r).', NN, [])./det(i, 1).^2;
        FDr_M((i-1)*R+r, :) = reshape(fdr_m(:, :, r).', NN, [])./det(i, 1).^2;
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
%--------------------------------------------------------------------------
%{
ffll = zeros(N, N);
fflr = zeros(N, N);
ffrl = zeros(N, N);
ffrr = zeros(N, N);
for j = 1:N
    ffll(:, j) = bfl.*bfl(:, j); % E (<-)
    fflr(:, j) = bfl.*bfr(:, j); % E (<-)
    ffrl(:, j) = bfr.*bfl(:, j); % W (->)
    ffrr(:, j) = bfr.*bfr(:, j); % W (->)
end
%}
ffll     = bfl_*bfl;
fflr     = bfl_*bfr;
ffrl     = bfr_*bfl;
ffrr     = bfr_*bfr;
ffll3_m  = mass\ffll;
fflr3_m  = mass\fflr;
ffrl3_m  =-mass\ffrl;
ffrr3_m  =-mass\ffrr;
bfl_ll   = bfl*ffll3_m./2;
bfl_lr   =-bfl*fflr3_m./2;
bfr_rl   =-bfr*ffrl3_m./2;
bfr_rr   = bfr*ffrr3_m./2;
Km1      = K-1;
bflLll_i = zeros(Km1, N);
bflLlr_i = zeros(Km1, N);
bfrLrl_i = zeros(Km1, N);
bfrLrr_i = zeros(Km1, N);
for i = 1:Km1
    bflLll_i(i, :) = bfl_ll./det(i  , 1);
    bflLlr_i(i, :) = bfl_lr./det(i  , 1);
    bfrLrl_i(i, :) = bfr_rl./det(i+1, 1);
    bfrLrr_i(i, :) = bfr_rr./det(i+1, 1);
end
bflLll_b = bfl_ll./det(K, 1);
bflLlr_b = bfl_lr./det(K, 1);
bfrLrl_b = bfr_rl./det(1, 1);
bfrLrr_b = bfr_rr./det(1, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE JAC.:
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGN TO (CELL)...
%--------------------------------------------------------------------------
Wbf_          = Wbf.';
dc_m          = bd_*Wbf_m;
intg.R        = R;
intg.xqc      = xqc;
intg.bd       = bd;
intg.bd_      = bd_;
intg.bf_      = bf_;
intg.Wbf      = Wbf;
intg.mass     = mass;
intg.Wbd_m    = Wbd_m;
intg.Wbf_m    = Wbf_m;
intg.dc_m     = dc_m;
%--------------------------------------------------------------------------
intg.DDr_M    = DDr_M;
intg.DFr_M    = DFr_M;
intg.FDr_M    = FDr_M;
intg.visu     = sum(Wbf_, 2)./sum(W_, 1);
intg.F_dof    = inittype(1, Fh, xdc, xqc, Wbf_m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGN TO (FACE)...
%--------------------------------------------------------------------------
bDl           = bdl./det;
bDr           = bdr./det;
bfl_          = bfl.';
bfr_          = bfr.';
eta           = 2;
intg.bDl      = bDl;
intg.bDr      = bDr;
intg.bfl      = bfl;
intg.bfr      = bfr;
intg.bfl_     = bfl_;
intg.bfr_     = bfr_;
intg.bDl_M    = bDl_M;
intg.bDr_M    = bDr_M;
intg.bfl_M    = bfl_M;
intg.bfr_M    = bfr_M;
intg.eta      = eta;
intg.bflLll_i = bflLll_i;
intg.bflLlr_i = bflLlr_i;
intg.bfrLrl_i = bfrLrl_i;
intg.bfrLrr_i = bfrLrr_i;
intg.bflLll_b = bflLll_b;
intg.bflLlr_b = bflLlr_b;
intg.bfrLrl_b = bfrLrl_b;
intg.bfrLrr_b = bfrLrr_b;
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end