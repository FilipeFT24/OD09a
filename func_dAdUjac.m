function [X] = func_dAdUjac(basis, mesh, intg, A_dof, U_dof, flag, alb1, arb2, ulb1, urb2, dlb1, drb2)
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
bflLll_i_Mt  = intg.bflLll_i_Mt;
bflLlr_i_Mt  = intg.bflLlr_i_Mt;
bfrLrl_i_Mt  = intg.bfrLrl_i_Mt;
bfrLrr_i_Mt  = intg.bfrLrr_i_Mt;
bflLllb2_Mt  = intg.bflLllb2_Mt;
bflLlrb2_Mt  = intg.bflLlrb2_Mt;
bfrLrlb1_Mt  = intg.bfrLrlb1_Mt;
bfrLrrb1_Mt  = intg.bfrLrrb1_Mt;
eta          = intg.eta;
Km1          = K-1;
NN           = N*N;
val          = zeros(TNN, 1);
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
ul_i         = U_dof(1:Km1, :)*bfl_;
ur_i         = U_dof(2:K  , :)*bfr_;
dl_i         = U_dof(1:Km1, :)*bdl_./det(1:Km1, 1); % bD = bd/det
dr_i         = U_dof(2:K  , :)*bdr_./det(2:K  , 1); % =
arb1         = A_dof(1    , :)*bfr_;
alb2         = A_dof(  K  , :)*bfl_;
urb1         = U_dof(1    , :)*bfr_;
ulb2         = U_dof(  K  , :)*bfl_;
drb1         = U_dof(1    , :)*bDr_(:, 1);
dlb2         = U_dof(  K  , :)*bDl_(:, K);
plb1         = ulb1.*bfrLrlb1_Mt;
prb2         = urb2.*bflLlrb2_Mt;
%--------------------------------------------------------------------------
J2aLll_it    = zeros(N, N, Km1);
J2aLlr_it    = zeros(N, N, Km1);
J2aLrl_it    = zeros(N, N, Km1);
J2aLrr_it    = zeros(N, N, Km1);
J2bLll_it    = zeros(N, N, Km1);
J2bLlr_it    = zeros(N, N, Km1);
J2bLrl_it    = zeros(N, N, Km1);
J2bLrr_it    = zeros(N, N, Km1);
J2cLll_it    = zeros(N, N, Km1);
J2cLlr_it    = zeros(N, N, Km1);
J2cLrl_it    = zeros(N, N, Km1);
J2cLrr_it    = zeros(N, N, Km1);
J2bLrrb1t    = zeros(N, N);
J2bLllb2t    = zeros(N, N);
if flag
    for i = 1:Km1
        J2bLll_it(:, :, i) = bfl_*(1./2.*bDl(i  , :).*(al_i(i, 1)+(ul_i(i, 1)-ur_i(i, 1)))); % (transpose): {adv}[u]+dv[u]
        J2bLlr_it(:, :, i) =-bfr_*(1./2.*bDl(i  , :).*(al_i(i, 1)));                         % =
        J2bLrl_it(:, :, i) = bfl_*(1./2.*bDr(i+1, :).*(ar_i(i, 1)));                         % =
        J2bLrr_it(:, :, i) =-bfr_*(1./2.*bDr(i+1, :).*(ar_i(i, 1)+(ur_i(i, 1)-ul_i(i, 1)))); % =
    end
    J2bLrrb1t =-bfr_*(bDr(1, :).*(arb1+(urb1-ulb1))); % =
    J2bLllb2t = bfl_*(bDl(K, :).*(alb2+(ulb2-urb2))); % =
end
for i = 1:Km1
    J2aLll_it(:, :, i) = (dl_i(i, 1).*bfl_+al_i(i, 1).*bDl_(:, i  ))./2*bfl; % (transpose): [v]{dau+adu}
    J2aLlr_it(:, :, i) = (dr_i(i, 1).*bfr_+ar_i(i, 1).*bDr_(:, i+1))./2*bfl; % =
    J2aLrl_it(:, :, i) =-(dl_i(i, 1).*bfl_+al_i(i, 1).*bDl_(:, i  ))./2*bfr; % =
    J2aLrr_it(:, :, i) =-(dr_i(i, 1).*bfr_+ar_i(i, 1).*bDr_(:, i+1))./2*bfr; % =
    J2cLll_it(:, :, i) = (bflLll_i_Mt(:, i).*(ul_i(i, 1)-ur_i(i, 1)+al_i(i, 1))+ar_i(i, 1).*bfrLrl_i_Mt(:, i))./2*bfl; % (transpose)
    J2cLlr_it(:, :, i) = (bfrLrr_i_Mt(:, i).*(ur_i(i, 1)-ul_i(i, 1)+ar_i(i, 1))+al_i(i, 1).*bflLlr_i_Mt(:, i))./2*bfl; % =
    J2cLrl_it(:, :, i) =-(bflLll_i_Mt(:, i).*(ul_i(i, 1)-ur_i(i, 1)+al_i(i, 1))+ar_i(i, 1).*bfrLrl_i_Mt(:, i))./2*bfr; % =
    J2cLrr_it(:, :, i) =-(bfrLrr_i_Mt(:, i).*(ur_i(i, 1)-ul_i(i, 1)+ar_i(i, 1))+al_i(i, 1).*bflLlr_i_Mt(:, i))./2*bfr; % =
end
J2aLrrb1t    =-(drb1.*bfr_+arb1.*bDr_(:, 1))./2*bfr; % (transpose)
J2aLllb2t    = (dlb2.*bfl_+alb2.*bDl_(:, K))./2*bfl; % =
J2cLrrb1t    =-(bfrLrrb1_Mt.*(urb1-ulb1+arb1))*bfr;  % =
J2cLllb2t    = (bflLllb2_Mt.*(ulb2-urb2+alb2))*bfl;  % =
%--------------------------------------------------------------------------
J2Lll_it     = J2aLll_it+J2bLll_it-J2cLll_it.*eta; % pen term is [v]{ar([u])}, where r([u]) = -{v}[u]) => minus sign
J2Llr_it     = J2aLlr_it+J2bLlr_it-J2cLlr_it.*eta; % =
J2Lrl_it     = J2aLrl_it+J2bLrl_it-J2cLrl_it.*eta; % =
J2Lrr_it     = J2aLrr_it+J2bLrr_it-J2cLrr_it.*eta; % =
J2Lrrb1t     = J2aLrrb1t+J2bLrrb1t-J2cLrrb1t.*eta;
J2Lllb2t     = J2aLllb2t+J2bLllb2t-J2cLllb2t.*eta;
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
val(IJdl, 1) = val(IJdl, 1)+VAL2dil_NN_1;
val(IJol, 1) = VAL2oil_NN_1;
val(IJor, 1) = VAL2oir_NN_1;
val(IJdr, 1) = val(IJdr, 1)+VAL2dir_NN_1;
val(IJb1, 1) = val(IJb1, 1)+VAL2db1_NN_1;
val(IJb2, 1) = val(IJb2, 1)+VAL2db2_NN_1;
rhs          = zeros(K, N);
rhs(1, :)    =-(bfr.*(dlb1.*ulb1+alb1.*dlb1)./2-bDr(1, :).*arb1.*ulb1-bfr.*(arb1.*plb1).*eta); % <-o
rhs(K, :)    = (bfl.*(drb2.*urb2+arb2.*drb2)./2-bDl(K, :).*alb2.*urb2-bfl.*(alb2.*prb2).*eta); % o->
X.lhs        = val;
X.rhs        = rhs;
%}
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
D_dof        = U_dof(1:Km1, :)-U_dof(2:K, :);
plli         = zeros(Km1, 1);
plri         = zeros(Km1, 1);
prli         = zeros(Km1, 1);
prri         = zeros(Km1, 1);
xlli         = zeros(Km1, 1);
xrri         = zeros(Km1, 1);
for i = 1:Km1
    plli(i, :) = U_dof(i  , :)*bflLll_i_Mt(:, i);
    plri(i, :) = U_dof(i+1, :)*bflLlr_i_Mt(:, i);
    prli(i, :) = U_dof(i  , :)*bfrLrl_i_Mt(:, i);
    prri(i, :) = U_dof(i+1, :)*bfrLrr_i_Mt(:, i);
    xlli(i, :) = D_dof(i  , :)*bflLll_i_Mt(:, i); % ul-ur
    xrri(i, :) =-D_dof(i  , :)*bfrLrr_i_Mt(:, i); % ur-ul
end
prb1         = U_dof(1, :)*bfrLrrb1_Mt;
plb2         = U_dof(K, :)*bflLllb2_Mt;
%--------------------------------------------------------------------------
aduo2_i      = (al_i.*dl_i+ar_i.*dr_i)./2;
aduo2b1      = (arb1.*drb1+alb1.*dlb1)./2;
aduo2b2      = (alb2.*dlb2+arb2.*drb2)./2;
duuo2_i      = (dr_i.*ur_i+dl_i.*ul_i)./2;
duuo2b1      = (drb1.*urb1+dlb1.*ulb1)./2;
duuo2b2      = (dlb2.*ulb2+drb2.*urb2)./2;
pl_i         = plli+plri;
pr_i         = prli+prri;
p_b1         = prb1+plb1;
p_b2         = plb2+prb2;
du_i         = ul_i-ur_i;
dub1         = urb1-ulb1;
dub2         = ulb2-urb2;
ap_o2_i      = (al_i.*pl_i+ar_i.*pr_i+du_i.*(plli-prri))./2; % (alLll+arLrl)*ul+(alLlr+arLrr)*ur+(ul-ur)*(Lll-Lrr) = al*(Lll+Llr)+ar*(Lrl+Lrr)+(ul-ur)*(Lll-Lrr)
ap_b1        = (arb1.*p_b1+dub1.*prb1);
ap_b2        = (alb2.*p_b2+dub2.*plb2);
%--------------------------------------------------------------------------
dAdU_doflib  = zeros(Km1, N);
dAdU_dofrib  = zeros(Km1, N);
if flag
    dAdU_doflib = 1./2.*bDl_M(1:Km1, :).*(al_i+ul_i).*(ul_i-ur_i);
    dAdU_dofrib = 1./2.*bDr_M(2:K  , :).*(ar_i+ur_i).*(ur_i-ul_i);
end
dAdU_doflia  = bfl_M(1:Km1, :).*(aduo2_i+duuo2_i);
dAdU_dofria  = bfr_M(2:K  , :).*(aduo2_i+duuo2_i);
dAdU_doflic  = bfl_M(1:Km1, :).*(ap_o2_i);
dAdU_dofric  = bfr_M(2:K  , :).*(ap_o2_i);
%--------------------------------------------------------------------------
dAdU_dofb1b  = zeros(1, N);
dAdU_dofb2b  = zeros(1, N);
if flag
    dAdU_dofb1b = bDrb1_M.*(arb1+urb1).*(urb1-ulb1);
    dAdU_dofb2b = bDlb2_M.*(alb2+ulb2).*(ulb2-urb2);
end
dAdU_dofb1a  = bfrb1_M.*(aduo2b1+duuo2b1);
dAdU_dofb2a  = bflb2_M.*(aduo2b2+duuo2b2);
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