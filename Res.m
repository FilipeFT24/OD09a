function [g] = Res(basis, mesh, intg, flag, nu, c, U_dof, ulb1, urb2, dlb1, drb2, uh, d1uh, cd2uh) %#ok<INUSD> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
N        = basis.N;
K        = mesh.K;
det      = mesh.det;
bd_      = intg.bd_;
bf_      = intg.bf_;
Wbf      = intg.Wbf;
Wbd_m    = intg.Wbd_m;
dc_m     = intg.dc_m;
d1U_dofx = zeros(K, N); % aux. variable (using for debugging later on)
%--------------------------------------------------------------------------
d1U      = U_dof*bd_;
U        = U_dof*bf_;
nupU     = nu+U;
nupUd1U  = nupU.*d1U;
d1U_dofl = U_dof*dc_m./det;
d2U_dofC =-nupUd1U*Wbd_m./det.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
Km1         = K-1;
bDl_        = intg.bDl_;
bDr_        = intg.bDr_;
bdl_        = intg.bdl_;
bdr_        = intg.bdr_;
bfl_        = intg.bfl_;
bfr_        = intg.bfr_;
bDl_M       = intg.bDl_M;
bDr_M       = intg.bDr_M;
bfl_M       = intg.bfl_M;
bfr_M       = intg.bfr_M;
bf_L_l_i_Mt = intg.bf_L_l_i_Mt;
bf_L_r_i_Mt = intg.bf_L_r_i_Mt;
bflLllb2_Mt = intg.bflLllb2_Mt;
bflLlrb2_Mt = intg.bflLlrb2_Mt;
bfrLrlb1_Mt = intg.bfrLrlb1_Mt;
bfrLrrb1_Mt = intg.bfrLrrb1_Mt;
eta         = intg.eta;
bDrb1_M     = bDr_M(1, :);
bDlb2_M     = bDl_M(K, :);
bfrb1_M     = bfr_M(1, :);
bflb2_M     = bfl_M(K, :);
bfl_i_M     = bfl_M(1:Km1, :);
bfr_i_M     = bfr_M(2:K  , :);
%--------------------------------------------------------------------------
ul_i        = U_dof(1:Km1, :)*bfl_;
ur_i        = U_dof(2:K  , :)*bfr_;
dl_i        = U_dof(1:Km1, :)*bdl_./det(1:Km1, 1); % bD = bd/det
dr_i        = U_dof(2:K  , :)*bdr_./det(2:K  , 1); % =
urb1        = U_dof(1    , :)*bfr_;
ulb2        = U_dof(  K  , :)*bfl_;
drb1        = U_dof(1    , :)*bDr_(:, 1);
dlb2        = U_dof(  K  , :)*bDl_(:, K);
p_li        = zeros(Km1, 1);
p_ri        = zeros(Km1, 1);
for i = 1:Km1
    p_li(i, :) = U_dof(i  , :)*bf_L_l_i_Mt(:, i);
    p_ri(i, :) = U_dof(i+1, :)*bf_L_r_i_Mt(:, i);
end
prb1        = U_dof(1    , :)*bfrLrrb1_Mt;
plb1        = ulb1.*bfrLrlb1_Mt;
prb2        = urb2.*bflLlrb2_Mt;
plb2        = U_dof(  K  , :)*bflLllb2_Mt;
nupul_i     = nu+ul_i;
nupur_i     = nu+ur_i;
nupurb1     = nu+urb1;
nupulb2     = nu+ulb2;
nupulb1     = nu+ulb1;
nupurb2     = nu+urb2;
nupuduo2_i  = (nupul_i.*dl_i+nupur_i.*dr_i)./2;
nupuduo2b1  = (nupurb1.*drb1+nupulb1.*dlb1)./2;
nupuduo2b2  = (nupulb2.*dlb2+nupurb2.*drb2)./2;
nupup_o2_i  = (nupul_i.*p_li+nupur_i.*p_ri)./2; % ((alLll+arLrl)*ul+(alLlr+arLrr)*ur)/2
nupup_o2b1  = (nupurb1.*prb1+nupulb1.*plb1);
nupup_o2b2  = (nupulb2.*plb2+nupurb2.*prb2);
%--------------------------------------------------------------------------
if flag
    %----------------------------------------------------------------------
    ujil               = ul_i-ur_i;
    ujir               =-ujil;
    ujb1               = urb1-ulb1;
    ujb2               = ulb2-urb2;
    rhsiL              = 1./2.*bfl_i_M.*ujil; % o->
    rhsiR              =-1./2.*bfr_i_M.*ujir; % <-o
    rhsb1              =-bfrb1_M.*ujb1; % <-o
    rhsb2              = bflb2_M.*ujb2; % o->
    d1U_dofx(1:Km1, :) = d1U_dofx(1:Km1, :)-rhsiL;
    d1U_dofx(2:K  , :) = d1U_dofx(2:K  , :)-rhsiR;
    d1U_dofx(    1, :) = d1U_dofx(    1, :)-rhsb1;
    d1U_dofx(    K, :) = d1U_dofx(    K, :)-rhsb2;
    d1U_dofl           = ...
        d1U_dofl+d1U_dofx;
    %----------------------------------------------------------------------
    %{
    ffll_mt            = intg.ffll_mt;
    fflr_mt            = intg.fflr_mt;
    ffrl_mt            = intg.ffrl_mt;
    ffrr_mt            = intg.ffrr_mt;
    U_dofl_i           = U_dof(1:Km1, :);
    U_dofr_i           = U_dof(2:K  , :);
    detl_i             = det  (1:Km1, 1);
    detr_i             = det  (2:K  , 1);
    rhsiLI             = 1./2.*(U_dofl_i*ffll_mt-U_dofr_i*fflr_mt)./detl_i;
    rhsiRI             =-1./2.*(U_dofr_i*ffrr_mt-U_dofl_i*ffrl_mt)./detr_i;
    rhsb1I             =-U_dof(1, :)*ffrr_mt./det(1, 1)+bfrb1_M.*ulb1;
    rhsb2I             = U_dof(K, :)*ffll_mt./det(K, 1)-bflb2_M.*urb2;
    diffl1             = rhsiL-rhsiLI; max(abs(diffl1), [], 'all')
    diffr1             = rhsiR-rhsiRI; max(abs(diffr1), [], 'all')
    diffl2             = rhsb1-rhsb1I; max(abs(diffl2), [], 'all')
    diffr2             = rhsb2-rhsb2I; max(abs(diffr2), [], 'all')
    %
    lift1              = zeros(K, N);
    lift1(1:Km1, :)    = lift1(1:Km1, :)-rhsiLI;
    lift1(2:K  , :)    = lift1(2:K  , :)-rhsiRI;
    lift1(    1, :)    = lift1(    1, :)-rhsb1I;
    lift1(    K, :)    = lift1(    K, :)-rhsb2I;
    MATL               = intg.MATL;
    RHSL               = zeros(K, N);
    RHSL(1, :)         = bfrb1_M.*ulb1;
    RHSL(K, :)         = bflb2_M.*urb2;
    U_doft             = U_dof.';
    U                  = U_doft(:);
    lift1I_1           = MATL*U;
    lift1I             = reshape(lift1I_1, [], K).'+RHSL;
    diff_3             = lift1-lift1I; max(abs(diff_3), [], 'all')
    %}
    %----------------------------------------------------------------------
end
%--------------------------------------------------------------------------
%
d2U_dof1ib = zeros(Km1, N);
d2U_dof2ib = zeros(Km1, N);
if flag
    d2U_dof1ib = 1./2.*bDr_M(2:K  , :).*nupur_i.*(ur_i-ul_i);
    d2U_dof2ib = 1./2.*bDl_M(1:Km1, :).*nupul_i.*(ul_i-ur_i);
end
d2U_dof1ia = bfr_M(2:K  , :).*nupuduo2_i;
d2U_dof2ia = bfl_M(1:Km1, :).*nupuduo2_i;
d2U_dof1ic = bfr_M(2:K  , :).*nupup_o2_i.*eta;
d2U_dof2ic = bfl_M(1:Km1, :).*nupup_o2_i.*eta;
%--------------------------------------------------------------------------
d2U_dofb1b = zeros(1, N);
d2U_dofb2b = zeros(1, N);
if flag
    d2U_dofb1b = bDrb1_M.*nupurb1.*(urb1-ulb1);
    d2U_dofb2b = bDlb2_M.*nupulb2.*(ulb2-urb2);
end
d2U_dofb1a = bfrb1_M.*nupuduo2b1;
d2U_dofb2a = bflb2_M.*nupuduo2b2;
d2U_dofb1c = bfrb1_M.*nupup_o2b1.*eta;
d2U_dofb2c = bflb2_M.*nupup_o2b2.*eta;
%--------------------------------------------------------------------------
d2U_dofC(2:K  , :) = d2U_dofC(2:K  , :)-(d2U_dof1ia+d2U_dof1ib+d2U_dof1ic); % <-o (con+sym+pen)
d2U_dofC(1:Km1, :) = d2U_dofC(1:Km1, :)+(d2U_dof2ia+d2U_dof2ib+d2U_dof2ic); % o-> (con+sym+pen)
d2U_dofC(1    , :) = d2U_dofC(1    , :)-(d2U_dofb1a+d2U_dofb1b+d2U_dofb1c); % =
d2U_dofC(  K  , :) = d2U_dofC(  K  , :)+(d2U_dofb2a+d2U_dofb2b+d2U_dofb2c); % =
%--------------------------------------------------------------------------
%{
d1Uadof_ = inittype(1,  d1uh, mesh.xdc, intg.xqc, intg.Wbf_m);
d2UadofC = inittype(1, cd2uh, mesh.xdc, intg.xqc, intg.Wbf_m);
figure;
subplot(1, 2, 1);
hold on;
plot(d1U_dofl*intg.visu, '-r');
plot(d1Uadof_*intg.visu, 'og');
subplot(1, 2, 2);
hold on;
plot(d2U_dofC*intg.visu, '-r');
plot(d2UadofC*intg.visu, 'og');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE RESIDUAL AND ASSIGN TO...
%--------------------------------------------------------------------------
F_dof_  = intg.F_dof;
F       =   F_dof_*bf_;
d1Ul    = d1U_dofl*bf_;
d2UC    = d2U_dofC*bf_;
d1Ul2   = d1Ul.^2;
R       = d2UC-c.*d1Ul2-F;
R_dof   = R*Wbf.*det;
%--------------------------------------------------------------------------
g.d1Ul  = d1Ul;  % aux. variable
g.R_dof = R_dof; % residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG JAC.:
%--------------------------------------------------------------------------
%
d2U_dofx   = zeros(K  , N);
d2U_dof1ib = zeros(Km1, N);
d2U_dof2ib = zeros(Km1, N);
if flag
    d2U_dof1ib = 1./2.*bDr_M(2:K  , :).*(ur_i-ul_i).*(nupur_i+ur_i);
    d2U_dof2ib = 1./2.*bDl_M(1:Km1, :).*(ul_i-ur_i).*(nupul_i+ul_i);
end
d2U_dof1ia = bfr_M(2:K  , :).*((dr_i.*ur_i+dl_i.*ul_i)./2+nupuduo2_i);
d2U_dof2ia = bfl_M(1:Km1, :).*((dl_i.*ul_i+dr_i.*ur_i)./2+nupuduo2_i);
d2U_dof1ic = bfr_M(2:K  , :).*nupup_o2_i.*eta;
d2U_dof2ic = bfl_M(1:Km1, :).*nupup_o2_i.*eta;
%--------------------------------------------------------------------------
d2U_dofb1b = zeros(1, N);
d2U_dofb2b = zeros(1, N);
if flag
    d2U_dofb1b = bDrb1_M.*(urb1-ulb1).*(nupurb1+urb1);
    d2U_dofb2b = bDlb2_M.*(ulb2-urb2).*(nupulb2+ulb2);
end
d2U_dofb1a = bfrb1_M.*((drb1.*urb1+dlb1.*ulb1)./2+nupuduo2b1);
d2U_dofb2a = bflb2_M.*((dlb2.*ulb2+drb2.*urb2)./2+nupuduo2b2);
d2U_dofb1c = bfrb1_M.*nupup_o2b1.*eta;
d2U_dofb2c = bflb2_M.*nupup_o2b2.*eta;
%--------------------------------------------------------------------------
d2U_dofx(2:K  , :) = d2U_dofx(2:K  , :)-(d2U_dof1ia+d2U_dof1ib+d2U_dof1ic); % <-o (con+sym+pen)
d2U_dofx(1:Km1, :) = d2U_dofx(1:Km1, :)+(d2U_dof2ia+d2U_dof2ib+d2U_dof2ic); % o-> (con+sym+pen)
d2U_dofx(1    , :) = d2U_dofx(1    , :)-(d2U_dofb1a+d2U_dofb1b+d2U_dofb1c); % =
d2U_dofx(  K  , :) = d2U_dofx(  K  , :)+(d2U_dofb2a+d2U_dofb2b+d2U_dofb2c); % =
%--------------------------------------------------------------------------
g.d1Ux             = d1U_dofx*bf_; % aux. variable (using for debugging later on)
g.d2U_dofx         = d2U_dofx;     % =
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end