function [d2U_dof] = G2C(basis, mesh, intg, flag, C_dof, U_dof, clb1, crb2, ulb1, urb2, dlb1, drb2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
det      = mesh.det;
bd_      = intg.bd_;
bf_      = intg.bf_;
Wbd_m    = intg.Wbd_m;
%--------------------------------------------------------------------------
C        = C_dof*bf_;
d1U      = U_dof*bd_;
Cd1U     = C.*d1U;
d2U_dof  =-Cd1U*Wbd_m./det.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
K        = mesh.K;
N        = basis.N;
Km1      = K-1;
bDl      = intg.bDl;
bDr      = intg.bDr;
bfl_     = intg.bfl_;
bfr_     = intg.bfr_;
bDl_M    = intg.bDl_M;
bDr_M    = intg.bDr_M;
bfl_M    = intg.bfl_M;
bfr_M    = intg.bfr_M;
bDrb1_M  = bDr_M(1, :);
bDlb2_M  = bDl_M(K, :);
bfrb1_M  = bfr_M(1, :);
bflb2_M  = bfl_M(K, :);
bflLll_i = intg.bflLll_i;
bflLlr_i = intg.bflLlr_i;
bfrLrl_i = intg.bfrLrl_i;
bfrLrr_i = intg.bfrLrr_i;
bflLll_b = intg.bflLll_b;
bflLlr_b = intg.bflLlr_b;
bfrLrl_b = intg.bfrLrl_b;
bfrLrr_b = intg.bfrLrr_b;
eta      = intg.eta;
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
for i = 1:K
    if i > 1
        if ~flag
            d2U_dof1ia = zeros(1, N);
        else
            d2U_dof1ia = 1./2.*c_ri(i-1, 1).*bDr_M(i, :).*(ur_i(i-1, 1)-ul_i(i-1, 1));
        end
        d2U_dof1ib    = bfr_M(i, :).*(c_ri(i-1, 1).*dr_i(i-1, 1)+c_li(i-1, 1).*dl_i(i-1, 1))./2; % con
        d2U_dof1ic    = bfr_M(i, :).*(c_ri(i-1, 1).*ar_i(i-1, 1)+c_li(i-1, 1).*al_i(i-1, 1));    % pen
        d2U_dof(i, :) = ...
            d2U_dof(i, :)-(d2U_dof1ia+d2U_dof1ib+d2U_dof1ic); % <-o
    end
    if i < K
        if ~flag
            d2U_dof2ia = zeros(1, N);
        else
            d2U_dof2ia = 1./2.*c_li(i, 1).*bDl_M(i, :).*(ul_i(i, 1)-ur_i(i, 1));
        end
        d2U_dof2ib    = bfl_M(i, :).*(c_li(i, 1).*dl_i(i, 1)+c_ri(i, 1).*dr_i(i, 1))./2; % con
        d2U_dof2ic    = bfl_M(i, :).*(c_li(i, 1).*al_i(i, 1)+c_ri(i, 1).*ar_i(i, 1));    % pen
        d2U_dof(i, :) = ...
            d2U_dof(i, :)+(d2U_dof2ia+d2U_dof2ib+d2U_dof2ic); % o->
    end
end
%--------------------------------------------------------------------------
if ~flag
    d2U_dofb1a = zeros(1, N);
    d2U_dofb2a = zeros(1, N);
else
    d2U_dofb1a = crb1.*bDrb1_M.*(urb1-ulb1); % lift: {v} = v, not v/2!
    d2U_dofb2a = clb2.*bDlb2_M.*(ulb2-urb2); % =
end
d2U_dofb1b    = bfrb1_M.*(crb1.*drb1+clb1.*dlb1)./2;   % con
d2U_dofb2b    = bflb2_M.*(clb2.*dlb2+crb2.*drb2)./2;   % =
d2U_dofb1c    = bfrb1_M.*(crb1.*arb1+clb1.*alb1).*eta; % pen
d2U_dofb2c    = bflb2_M.*(clb2.*alb2+crb2.*arb2).*eta; % =
d2U_dof(1, :) = ...
    d2U_dof(1, :)-(d2U_dofb1a+d2U_dofb1b+d2U_dofb1c);
d2U_dof(K, :) = ...
    d2U_dof(K, :)+(d2U_dofb2a+d2U_dofb2b+d2U_dofb2c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end