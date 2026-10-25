function [X] = func_dAdUjac(basis, mesh, intg, A_dof, U_dof, alb1, arb2, ulb1, urb2, dlb1, drb2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPLICIT:
%--------------------------------------------------------------------------
N            = basis.N;
K            = mesh.K;
det          = mesh.det;
TNN          = mesh.TNN;
IJd_         = mesh.IJd_;
IJdl         = mesh.IJdl; % ll
IJol         = mesh.IJol; % lr
IJor         = mesh.IJor; % rl
IJdr         = mesh.IJdr; % rr
IJb1         = mesh.IJb1; % rrb1
IJb2         = mesh.IJb2; % llb2
Km1          = K-1;
NN           = N*N;
val          = zeros(TNN, 1);
%--------------------------------------------------------------------------
bf_          = intg.bf_;
DDr          = intg.DDr;
R            = intg.R;
bDl          = intg.bDl;
bDr          = intg.bDr;
bfl          = intg.bfl;
bfr          = intg.bfr;
bDl_         = intg.bDl_;
bDr_         = intg.bDr_;
bdl_         = intg.bdl_;
bdr_         = intg.bdr_;
bfl_         = intg.bfl_;
bfr_         = intg.bfr_;
bflLlrb2_Mt  = intg.bflLlrb2_Mt;
bfrLrlb1_Mt  = intg.bfrLrlb1_Mt;
Jmll_li      = intg.Jmll_li;
Jmll_ri      = intg.Jmll_ri;
Jmlr_li      = intg.Jmlr_li;
Jmlr_ri      = intg.Jmlr_ri;
Jmrl_li      = intg.Jmrl_li;
Jmrl_ri      = intg.Jmrl_ri;
Jmrr_li      = intg.Jmrr_li;
Jmrr_ri      = intg.Jmrr_ri;
Jmb1         = intg.Jmb1;
Jmb2         = intg.Jmb2;
FFll         = intg.FFll;
FFlr         = intg.FFlr;
FFrl         = intg.FFrl;
FFrr         = intg.FFrr;
Xdull_i      = intg.Xdull_i;
Xdulr_i      = intg.Xdulr_i;
Xdurl_i      = intg.Xdurl_i;
Xdurr_i      = intg.Xdurr_i;
FFb1         = intg.FFb1;
FFb2         = intg.FFb2;
Xdub1        = intg.Xdub1;
Xdub2        = intg.Xdub2;
eta          = intg.eta;
%--------------------------------------------------------------------------
A            = A_dof*bf_;
At           = A.';
At1          = At(:);
VAL2_RK_NN   =-DDr.*At1;
VAL2_NN_K    = sum(reshape(VAL2_RK_NN.', NN, R, []), 2);
VAL2_NN_K_1  = VAL2_NN_K(:);
val(IJd_, 1) = VAL2_NN_K_1;
%--------------------------------------------------------------------------
al_i         = A_dof(1:Km1, :)*bfl_;
ar_i         = A_dof(2:K  , :)*bfr_;
arb1         = A_dof(1    , :)*bfr_;
alb2         = A_dof(  K  , :)*bfl_;
dl_i         = U_dof(1:Km1, :)*bdl_./det(1:Km1, 1); % bD = bd/det
dr_i         = U_dof(2:K  , :)*bdr_./det(2:K  , 1); % =
ul_i         = U_dof(1:Km1, :)*bfl_;
ur_i         = U_dof(2:K  , :)*bfr_;
al_it        = al_i.';
ar_it        = ar_i.';
jac1ll_it    = Jmll_li.*al_it+Jmll_ri.*ar_it; % [v]{adu}+{adv}[u]+pen
jac1lr_it    = Jmlr_li.*al_it+Jmlr_ri.*ar_it; % =
jac1rl_it    = Jmrl_li.*al_it+Jmrl_ri.*ar_it; % =
jac1rr_it    = Jmrr_li.*al_it+Jmrr_ri.*ar_it; % =
jac1b1_t     = Jmb1.*arb1;
jac1b2_t     = Jmb2.*alb2;
%--------------------------------------------------------------------------
dl_it        = dl_i.';
dr_it        = dr_i.';
du_i         = ul_i-ur_i;
dulit        = du_i.';
durit        =-dulit;
jac2ll_it    = FFll.*dl_it+Xdull_i.*dulit; % [v]{dau}+dv[u]+pen
jac2lr_it    = FFlr.*dr_it+Xdulr_i.*durit; % =
jac2rl_it    = FFrl.*dl_it+Xdurl_i.*dulit; % =
jac2rr_it    = FFrr.*dr_it+Xdurr_i.*durit; % =
jac_ll_it    = jac1ll_it+jac2ll_it;
jac_lr_it    = jac1lr_it+jac2lr_it;
jac_rl_it    = jac1rl_it+jac2rl_it;
jac_rr_it    = jac1rr_it+jac2rr_it;
%--------------------------------------------------------------------------
drb1         = U_dof(1, :)*bDr_(:, 1);
dlb2         = U_dof(K, :)*bDl_(:, K);
urb1         = U_dof(1, :)*bfr_;
ulb2         = U_dof(K, :)*bfl_;
dub1         = urb1-ulb1;
dub2         = ulb2-urb2;
jac2b1_t     = FFb1.*drb1+Xdub1.*dub1;
jac2b2_t     = FFb2.*dlb2+Xdub2.*dub2;
jac_b1_t     = jac1b1_t+jac2b1_t;
jac_b2_t     = jac1b2_t+jac2b2_t;
%--------------------------------------------------------------------------
VAL2dil_NN_1 = jac_ll_it(:);
VAL2oil_NN_1 = jac_lr_it(:);
VAL2oir_NN_1 = jac_rl_it(:);
VAL2dir_NN_1 = jac_rr_it(:);
VAL2db1_NN_1 = jac_b1_t (:);
VAL2db2_NN_1 = jac_b2_t (:);
val(IJdl, 1) = val(IJdl, 1)+VAL2dil_NN_1;
val(IJol, 1) = VAL2oil_NN_1;
val(IJor, 1) = VAL2oir_NN_1;
val(IJdr, 1) = val(IJdr, 1)+VAL2dir_NN_1;
val(IJb1, 1) = val(IJb1, 1)+VAL2db1_NN_1;
val(IJb2, 1) = val(IJb2, 1)+VAL2db2_NN_1;
%--------------------------------------------------------------------------
rhs          = zeros(K, N);
plb1         = ulb1.*bfrLrlb1_Mt;
prb2         = urb2.*bflLlrb2_Mt;
rhs1b1       =-(bfr.*(alb1.*dlb1)./2-bDr(1, :).*arb1.*ulb1-bfr.*(arb1.*plb1).*eta); % <-o
rhs1b2       = (bfl.*(arb2.*drb2)./2-bDl(K, :).*alb2.*urb2-bfl.*(alb2.*prb2).*eta); % o->
rhs2b1       =-(bfr.*(dlb1.*ulb1)./2); % <-o
rhs2b2       = (bfl.*(drb2.*urb2)./2); % o->
rhs(1, :)    = rhs1b1+rhs2b1;
rhs(K, :)    = rhs1b2+rhs2b2;
X.lhs        = val;
X.rhs        = rhs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPLICIT:
%--------------------------------------------------------------------------
%{
bd_          = intg.bd_;
Wbd_m        = intg.Wbd_m;
bDl_M        = intg.bDl_M;
bDr_M        = intg.bDr_M;
bfl_M        = intg.bfl_M;
bfr_M        = intg.bfr_M;
bDrb1_M      = bDr_M(1, :);
bDlb2_M      = bDl_M(K, :);
bfrb1_M      = bfr_M(1, :);
bflb2_M      = bfl_M(K, :);
bflLll_i_Mt  = intg.bflLll_i_Mt;
bflLlr_i_Mt  = intg.bflLlr_i_Mt;
bfrLrl_i_Mt  = intg.bfrLrl_i_Mt;
bfrLrr_i_Mt  = intg.bfrLrr_i_Mt;
bflLllb2_Mt  = intg.bflLllb2_Mt;
bfrLrrb1_Mt  = intg.bfrLrrb1_Mt;
%--------------------------------------------------------------------------
plli         = zeros(Km1, 1);
plri         = zeros(Km1, 1);
prli         = zeros(Km1, 1);
prri         = zeros(Km1, 1);
for i = 1:Km1
    plli(i, :) = U_dof(i  , :)*bflLll_i_Mt(:, i);
    plri(i, :) = U_dof(i+1, :)*bflLlr_i_Mt(:, i);
    prli(i, :) = U_dof(i  , :)*bfrLrl_i_Mt(:, i);
    prri(i, :) = U_dof(i+1, :)*bfrLrr_i_Mt(:, i);
end
pl_i         = plli+plri;
pr_i         = prli+prri;
du_i         = ul_i-ur_i;
duli         = du_i;
duri         =-duli;
aduo2_i      = (al_i.*dl_i+ar_i.*dr_i)./2;
duuo2_i      = (dr_i.*ur_i+dl_i.*ul_i)./2;
ap_o2_i      = (al_i.*pl_i+ar_i.*pr_i+du_i.*(plli-prri))./2; % (alLll+arLrl)*ul+(alLlr+arLrr)*ur+(ul-ur)*(Lll-Lrr) = al*(Lll+Llr)+ar*(Lrl+Lrr)+(ul-ur)*(Lll-Lrr)
dAdU_doflia  = bfl_M(1:Km1, :).*(aduo2_i+duuo2_i);
dAdU_dofria  = bfr_M(2:K  , :).*(aduo2_i+duuo2_i);
dAdU_doflib  = bDl_M(1:Km1, :)./2.*(al_i+ul_i).*duli;
dAdU_dofrib  = bDr_M(2:K  , :)./2.*(ar_i+ur_i).*duri;
dAdU_doflic  = bfl_M(1:Km1, :).*ap_o2_i;
dAdU_dofric  = bfr_M(2:K  , :).*ap_o2_i;
%--------------------------------------------------------------------------
prb1         = U_dof(1, :)*bfrLrrb1_Mt;
plb2         = U_dof(K, :)*bflLllb2_Mt;
p_b1         = prb1+plb1;
p_b2         = plb2+prb2;
dub1         = urb1-ulb1;
dub2         = ulb2-urb2;
ap_b1        = (arb1.*p_b1+dub1.*prb1);
ap_b2        = (alb2.*p_b2+dub2.*plb2);
aduo2b1      = (arb1.*drb1+alb1.*dlb1)./2;
aduo2b2      = (alb2.*dlb2+arb2.*drb2)./2;
duuo2b1      = (drb1.*urb1+dlb1.*ulb1)./2;
duuo2b2      = (dlb2.*ulb2+drb2.*urb2)./2;
dAdU_dofb1a  = bfrb1_M.*(aduo2b1+duuo2b1);
dAdU_dofb2a  = bflb2_M.*(aduo2b2+duuo2b2);
dAdU_dofb1b  = bDrb1_M.*(arb1+urb1).*(urb1-ulb1);
dAdU_dofb2b  = bDlb2_M.*(alb2+ulb2).*(ulb2-urb2);
dAdU_dofb1c  = bfrb1_M.*(ap_b1);
dAdU_dofb2c  = bflb2_M.*(ap_b2);
%--------------------------------------------------------------------------
A                      = A_dof*bf_; % nu+U
dU                     = U_dof*bd_;
AdU                    = A.*dU;
dAdU_doffjac           =-AdU*Wbd_m./det.^2;
dAdU_doffjac(1:Km1, :) = dAdU_doffjac(1:Km1, :)+(dAdU_doflia+dAdU_doflib-dAdU_doflic.*eta); % o-> (con+sym+pen) pen term is [v]{ar([u])}, where r([u]) = -{v}[u]) => minus sign
dAdU_doffjac(2:K  , :) = dAdU_doffjac(2:K  , :)-(dAdU_dofria+dAdU_dofrib-dAdU_dofric.*eta); % <-o (con+sym+pen)
dAdU_doffjac(1    , :) = dAdU_doffjac(1    , :)-(dAdU_dofb1a+dAdU_dofb1b-dAdU_dofb1c.*eta); % =
dAdU_doffjac(  K  , :) = dAdU_doffjac(  K  , :)+(dAdU_dofb2a+dAdU_dofb2b-dAdU_dofb2c.*eta); % =
%--------------------------------------------------------------------------
aux          = zeros(TNN);
aux(IJd_, 1) = VAL2_NN_K_1;
aux(IJdl, 1) = aux(IJdl, 1)+VAL2dil_NN_1;
aux(IJol, 1) = VAL2oil_NN_1;
aux(IJor, 1) = VAL2oir_NN_1;
aux(IJdr, 1) = aux(IJdr, 1)+VAL2dir_NN_1;
aux(IJb1, 1) = aux(IJb1, 1)+VAL2db1_NN_1;
aux(IJb2, 1) = aux(IJb2, 1)+VAL2db2_NN_1;
Mat          = mesh.Jaux;
AUX          = Mat(aux);
B_doft       = U_dof.';
B1doft       = B_doft(:);
X1doft       = AUX*B1doft;
X1dofI       = reshape(X1doft, [], K).'+rhs;
X1dofE       = dAdU_doffjac;
mass         = intg.mass;
X1dofE       = X1dofE*mass.*det;
diff_2       = X1dofI-X1dofE; max(abs(diff_2), [], 'all')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end