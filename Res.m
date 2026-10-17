function [g] = Res(basis, mesh, intg, flag, c, C_dof, U_dof, clb1, crb2, ulb1, urb2, dlb1, drb2, uh, d1uh, cd2uh) %#ok<INUSD> 
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
C        = C_dof*bf_;
d1U      = U_dof*bd_;
Cd1U     = C.*d1U;
d1U_dofl = U_dof*dc_m./det;
d2U_dofC =-Cd1U*Wbd_m./det.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
Km1      = K-1;
bDl      = intg.bDl;
bDr      = intg.bDr;
bfl_     = intg.bfl_;
bfr_     = intg.bfr_;
bDl_M    = intg.bDl_M;
bDr_M    = intg.bDr_M;
bfl_M    = intg.bfl_M;
bfr_M    = intg.bfr_M;
bflLll_i = intg.bflLll_i;
bflLlr_i = intg.bflLlr_i;
bfrLrl_i = intg.bfrLrl_i;
bfrLrr_i = intg.bfrLrr_i;
bflLll_b = intg.bflLll_b;
bflLlr_b = intg.bflLlr_b;
bfrLrl_b = intg.bfrLrl_b;
bfrLrr_b = intg.bfrLrr_b;
eta      = intg.eta;
bDrb1_M  = bDr_M(1, :);
bDlb2_M  = bDl_M(K, :);
bfrb1_M  = bfr_M(1, :);
bflb2_M  = bfl_M(K, :);
bfl_i_M  = bfl_M(1:Km1, :);
bfr_i_M  = bfr_M(2:K  , :);
%--------------------------------------------------------------------------
c_li     = C_dof(1:Km1, :)*bfl_;
c_ri     = C_dof(2:K  , :)*bfr_;
crb1     = C_dof(    1, :)*bfr_;
clb2     = C_dof(    K, :)*bfl_;
ul_i     = U_dof(1:Km1, :)*bfl_;
ur_i     = U_dof(2:K  , :)*bfr_;
urb1     = U_dof(    1, :)*bfr_;
ulb2     = U_dof(    K, :)*bfl_;
dl_i     = zeros(Km1, 1);
dr_i     = zeros(Km1, 1);
plli     = zeros(Km1, 1);
plri     = zeros(Km1, 1);
prli     = zeros(Km1, 1);
prri     = zeros(Km1, 1);
for i = 1:Km1
    dl_i(i, 1) = bDl     (i  , :)*U_dof(i  , :).';
    dr_i(i, 1) = bDr     (i+1, :)*U_dof(i+1, :).';
    plli(i, 1) = bflLll_i(i  , :)*U_dof(i  , :).';
    plri(i, 1) = bflLlr_i(i  , :)*U_dof(i+1, :).';
    prli(i, 1) = bfrLrl_i(i  , :)*U_dof(i  , :).';
    prri(i, 1) = bfrLrr_i(i  , :)*U_dof(i+1, :).';
end
drb1     = bDr(1, :)*U_dof(1, :).';
dlb2     = bDl(K, :)*U_dof(K, :).';
al_i     = (plli+plri)./2;
ar_i     = (prli+prri)./2;
arb1     = (bfrLrr_b*U_dof(1, :).'+0)./2; % fix w/ bc
alb1     = (bflLlr_b*U_dof(1, :).'+0)./2;
alb2     = (bflLll_b*U_dof(K, :).'+0)./2;
arb2     = (bfrLrl_b*U_dof(K, :).'+0)./2;
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
    RHSL(1, :)         =-bfrb1_M.*ulb1;
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
for i = 1:K
    if i > 1
        if ~flag
            d2U_dof1ib = zeros(1, N);
        else
            d2U_dof1ib = 1./2.*bDr_M(i, :).*c_ri(i-1, 1).*(ur_i(i-1, 1)-ul_i(i-1, 1));
        end
        d2U_dof1ia     = bfr_M(i, :).*(c_ri(i-1, 1).*dr_i(i-1, 1)+c_li(i-1, 1).*dl_i(i-1, 1))./2;
        d2U_dof1ic     = bfr_M(i, :).*(c_ri(i-1, 1).*ar_i(i-1, 1)+c_li(i-1, 1).*al_i(i-1, 1));
        d2U_dofC(i, :) = ...
            d2U_dofC(i, :)-(d2U_dof1ia+d2U_dof1ib+d2U_dof1ic); % <-o (con+sym+pen)
    end
    if i < K
        if ~flag
            d2U_dof2ib = zeros(1, N);
        else
            d2U_dof2ib = 1./2.*bDl_M(i, :).*c_li(i, 1).*(ul_i(i, 1)-ur_i(i, 1));
        end
        d2U_dof2ia     = bfl_M(i, :).*(c_li(i, 1).*dl_i(i, 1)+c_ri(i, 1).*dr_i(i, 1))./2;
        d2U_dof2ic     = bfl_M(i, :).*(c_li(i, 1).*al_i(i, 1)+c_ri(i, 1).*ar_i(i, 1));
        d2U_dofC(i, :) = ...
            d2U_dofC(i, :)+(d2U_dof2ia+d2U_dof2ib+d2U_dof2ic); % o-> (con+sym+pen)
    end
end
%--------------------------------------------------------------------------
if ~flag
    d2U_dofb1b = zeros(1, N);
    d2U_dofb2b = zeros(1, N);
else
    d2U_dofb1b = bDrb1_M.*crb1.*(urb1-ulb1);
    d2U_dofb2b = bDlb2_M.*clb2.*(ulb2-urb2);
end
d2U_dofb1a     = bfrb1_M.*(crb1.*drb1+clb1.*dlb1)./2;
d2U_dofb2a     = bflb2_M.*(clb2.*dlb2+crb2.*drb2)./2;
d2U_dofb1c     = bfrb1_M.*(crb1.*arb1+clb1.*alb1).*eta;
d2U_dofb2c     = bflb2_M.*(clb2.*alb2+crb2.*arb2).*eta;
d2U_dofC(1, :) = ...
    d2U_dofC(1, :)-(d2U_dofb1a+d2U_dofb1b+d2U_dofb1c); % <-o (con+sym+pen)
d2U_dofC(K, :) = ...
    d2U_dofC(K, :)+(d2U_dofb2a+d2U_dofb2b+d2U_dofb2c); % o-> (con+sym+pen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE RESIDUAL AND ASSIGN TO...
%--------------------------------------------------------------------------
F_dof_  = intg.F_dof;
F       =   F_dof_*bf_;
d1Ul    = d1U_dofl*bf_;
d1Ux    = d1U_dofx*bf_;
d2UC    = d2U_dofC*bf_;
d1Ul2   = d1Ul.^2;
R       = d2UC-c.*d1Ul2-F;
R_dof   = R*Wbf.*det;
%--------------------------------------------------------------------------
g.d1Ul  = d1Ul;  % aux. variable
g.d1Ux  = d1Ux;  % aux. variable (using for debugging later on)
g.R_dof = R_dof; % residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK:
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
end