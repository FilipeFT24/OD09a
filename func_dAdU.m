function [dAdU_doft] = func_dAdU(basis, mesh, intg, A_dof, U_dof, alb1, arb2, ulb1, urb2, dlb1, drb2, funch) %#ok<INUSD>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
det         = mesh.det;
bd_         = intg.bd_;
bf_         = intg.bf_;
Wbd_m       = intg.Wbd_m;
%--------------------------------------------------------------------------
A           = A_dof*bf_; % nu+U
dU          = U_dof*bd_;
AdU         = A.*dU;
dAdU_doft   =-AdU*Wbd_m./det.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
N           = basis.N;
K           = mesh.K;
Km1         = K-1;
bdl_        = intg.bdl_;
bdr_        = intg.bdr_;
bfl_        = intg.bfl_;
bfr_        = intg.bfr_;
bDl_        = intg.bDl_;
bDr_        = intg.bDr_;
bDl_M       = intg.bDl_M;
bDr_M       = intg.bDr_M;
bfl_M       = intg.bfl_M;
bfr_M       = intg.bfr_M;
bflLll_i_Mt = intg.bflLll_i_Mt;
bflLlr_i_Mt = intg.bflLlr_i_Mt;
bfrLrl_i_Mt = intg.bfrLrl_i_Mt;
bfrLrr_i_Mt = intg.bfrLrr_i_Mt;
bflLllb2_Mt = intg.bflLllb2_Mt;
bflLlrb2_Mt = intg.bflLlrb2_Mt;
bfrLrlb1_Mt = intg.bfrLrlb1_Mt;
bfrLrrb1_Mt = intg.bfrLrrb1_Mt;
eta         = intg.eta;
dAdU_doff   = zeros(K, N);
%--------------------------------------------------------------------------
bDrb1_M     = bDr_M(1, :);
bDlb2_M     = bDl_M(K, :);
bfrb1_M     = bfr_M(1, :);
bflb2_M     = bfl_M(K, :);
al_i        = A_dof(1:Km1, :)*bfl_;
ar_i        = A_dof(2:K  , :)*bfr_;
ul_i        = U_dof(1:Km1, :)*bfl_;
ur_i        = U_dof(2:K  , :)*bfr_;
dl_i        = U_dof(1:Km1, :)*bdl_./det(1:Km1, 1); % bD = bd/det
dr_i        = U_dof(2:K  , :)*bdr_./det(2:K  , 1); % =
arb1        = A_dof(1    , :)*bfr_;
alb2        = A_dof(  K  , :)*bfl_;
urb1        = U_dof(1    , :)*bfr_;
ulb2        = U_dof(  K  , :)*bfl_;
drb1        = U_dof(1    , :)*bDr_(:, 1);
dlb2        = U_dof(  K  , :)*bDl_(:, K);
plli        = zeros(Km1, 1);
plri        = zeros(Km1, 1);
prli        = zeros(Km1, 1);
prri        = zeros(Km1, 1);
for i = 1:Km1
    plli(i, :) = U_dof(i  , :)*bflLll_i_Mt(:, i);
    plri(i, :) = U_dof(i+1, :)*bflLlr_i_Mt(:, i);
    prli(i, :) = U_dof(i  , :)*bfrLrl_i_Mt(:, i);
    prri(i, :) = U_dof(i+1, :)*bfrLrr_i_Mt(:, i);
end
prb1        = U_dof(1    , :)*bfrLrrb1_Mt;
plb1        = ulb1.*bfrLrlb1_Mt;
prb2        = urb2.*bflLlrb2_Mt;
plb2        = U_dof(  K  , :)*bflLllb2_Mt;
%--------------------------------------------------------------------------
aduo2_i     = (al_i.*dl_i+ar_i.*dr_i)./2;
aduo2b1     = (arb1.*drb1+alb1.*dlb1)./2;
aduo2b2     = (alb2.*dlb2+arb2.*drb2)./2;
%--------------------------------------------------------------------------
pl_i        = plli+plri;
pr_i        = prli+prri;
ap_o2_i     = (al_i.*pl_i+ar_i.*pr_i)./2; % ((alLll+arLrl)*ul+(alLlr+arLrr)*ur)/2 = (al*(Lll+Llr)+ar*(Lrl+Lrr))/2
p_b1        = prb1+plb1;
p_b2        = plb2+prb2;
ap_b1       = arb1.*p_b1;
ap_b2       = alb2.*p_b2;
%--------------------------------------------------------------------------
dAdU_doflia = bfl_M(1:Km1, :).*aduo2_i;
dAdU_dofria = bfr_M(2:K  , :).*aduo2_i;
dAdU_doflib = 1./2.*bDl_M(1:Km1, :).*al_i.*(ul_i-ur_i);
dAdU_dofrib = 1./2.*bDr_M(2:K  , :).*ar_i.*(ur_i-ul_i);
dAdU_doflic = bfl_M(1:Km1, :).*ap_o2_i;
dAdU_dofric = bfr_M(2:K  , :).*ap_o2_i;
%--------------------------------------------------------------------------
dAdU_dofb1a = bfrb1_M.*aduo2b1;
dAdU_dofb2a = bflb2_M.*aduo2b2;
dAdU_dofb1b = bDrb1_M.*arb1.*(urb1-ulb1);
dAdU_dofb2b = bDlb2_M.*alb2.*(ulb2-urb2);
dAdU_dofb1c = bfrb1_M.*ap_b1;
dAdU_dofb2c = bflb2_M.*ap_b2;
%--------------------------------------------------------------------------
dAdU_doff(1:Km1, :) = dAdU_doff(1:Km1, :)+(dAdU_doflia+dAdU_doflib-dAdU_doflic.*eta); % o-> (con+sym+pen) (pen term is [v]{ar([u])}, where r([u]) = -{v}[u]) => minus sign
dAdU_doff(2:K  , :) = dAdU_doff(2:K  , :)-(dAdU_dofria+dAdU_dofrib-dAdU_dofric.*eta); % <-o (con+sym+pen)
dAdU_doff(1    , :) = dAdU_doff(1    , :)-(dAdU_dofb1a+dAdU_dofb1b-dAdU_dofb1c.*eta); % =
dAdU_doff(  K  , :) = dAdU_doff(  K  , :)+(dAdU_dofb2a+dAdU_dofb2b-dAdU_dofb2c.*eta); % =
dAdU_doft           = dAdU_doft+dAdU_doff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK:
%--------------------------------------------------------------------------
%{
visu      = intg.visu;
xdc       = mesh.xdc;
xqc       = intg.xqc;
Wbf_m     = intg.Wbf_m;
dAdU_dofa = inittype(1, funch, xdc, xqc, Wbf_m); % funch = function handle
figure;
hold on;
plot(dAdU_doft*visu, '-r');
plot(dAdU_dofa*visu, 'og');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end