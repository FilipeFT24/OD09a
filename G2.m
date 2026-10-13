function [d2U_dof] = G2(mesh, intg, flag, A_dof, U_dof, Alb1, Arb2, flb1, frb2, dlb1, drb2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
det      = mesh.det;
bd_      = intg.bd_;
bf_      = intg.bf_;
Wbd_m    = intg.Wbd_m;
%--------------------------------------------------------------------------
A        = A_dof*bf_;
d1U      = U_dof*bd_;
Ad1U     = A.*d1U;
d2U_dof  =-Ad1U*Wbd_m./det.^2;
%
d2U_dof  = d2U_dof.*0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
K        = mesh.K;
Km1      = K-1;
Km2      = K-2;
bDl      = intg.bDl;
bDr      = intg.bDr;

bDrb1    = bDr(1, :);
bDlb2    = bDl(K, :);


bfl      = intg.bfl;
bfr      = intg.bfr;
bfl_     = intg.bfl_;
bfr_     = intg.bfr_;
eta      = intg.eta;
bflLll_i = intg.bflLll_i;
bflLlr_i = intg.bflLlr_i;
bfrLrl_i = intg.bfrLrl_i;
bfrLrr_i = intg.bfrLrr_i;
bflLll_b = intg.bflLll_b;
bflLlr_b = intg.bflLlr_b;
bfrLrl_b = intg.bfrLrl_b;
bfrLrr_b = intg.bfrLrr_b;
%--------------------------------------------------------------------------
A_li     = A_dof(1:Km1, :)*bfl_;
A_ri     = A_dof(2:K  , :)*bfr_;
Arb1     = A_dof(    1, :)*bfr_;
Alb2     = A_dof(    K, :)*bfl_;
fl_i     = U_dof(1:Km1, :)*bfl_;
fr_i     = U_dof(2:K  , :)*bfr_;
frb1     = U_dof(    1, :)*bfr_;
flb2     = U_dof(    K, :)*bfl_;
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
pl_i     = (plli+plri)./2;
pr_i     = (prli+prri)./2;
drb1     = bDr(1, :)*U_dof(1, :).';
dlb2     = bDl(K, :)*U_dof(K, :).';
arb1     = (bfrLrr_b*U_dof(1, :).'+0)./2; % fix w/ bc
alb1     = (bflLlr_b*U_dof(1, :).'+0)./2;
alb2     = (bflLll_b*U_dof(K, :).'+0)./2;
arb2     = (bfrLrl_b*U_dof(K, :).'+0)./2;
%--------------------------------------------------------------------------
for i = 1:Km2
    %{
    d2U_dof(i+1, :) = d2U_dof(i+1, :)+...
        -(A_ri(i  , 1).*dr_i(i  , 1)+A_li(i  , 1).*dl_i(i  , 1)).*bfr./2 ...
        +(A_li(i+1, 1).*dl_i(i+1, 1)+A_ri(i+1, 1).*dr_i(i+1, 1)).*bfl./2;
    d2U_dof(i+1, :) = d2U_dof(i+1, :)-eta.*(...
        -(A_ri(i  , 1).*pr_i(i  , :)+A_li(i  , 1).*pl_i(i  , :)).*bfr ...
        +(A_li(i+1, 1).*pl_i(i+1, :)+A_ri(i+1, 1).*pr_i(i+1, :)).*bfl);
    if flag
        d2U_dof(i+1, :) = d2U_dof(i+1, :)+...
            -(A_ri(i  , 1).*fr_i(i  , 1)-A_ri(i  , 1).*fl_i(i  , 1)).*bDr(i+1, :)./2 ...
            +(A_li(i+1, 1).*fl_i(i+1, 1)-A_li(i+1, 1).*fr_i(i+1, 1)).*bDl(i+1, :)./2;
    end
    %}


    
    d2U_dof(i+1, :) = d2U_dof(i+1, :)+...
        -1./2.*A_ri(i  , 1).*bDr(i+1, :).*(fr_i(i  , 1)-fl_i(i  , 1)) ...
        +1./2.*A_li(i+1, 1).*bDl(i+1, :).*(fl_i(i+1, 1)-fr_i(i+1, 1));
end


d2U_dof(1, :) = d2U_dof(1, :)-Arb1.*bDrb1.*(frb1-flb1)+1./2.*A_li(  1, 1).*bDl(1, :).*(fl_i(  1, 1)-fr_i(  1, 1));
d2U_dof(K, :) = d2U_dof(K, :)+Alb2.*bDlb2.*(flb2-frb2)-1./2.*A_ri(Km1, 1).*bDr(K, :).*(fr_i(Km1, 1)-fl_i(Km1, 1));


%{
d2U_dof(1, :) = d2U_dof(1, :)-(Arb1.*drb1+Alb1.*dlb1).*bfr;
d2U_dof(K, :) = d2U_dof(K, :)+(Alb2.*dlb2+Arb2.*drb2).*bfl;
d2U_dof(1, :) = d2U_dof(1, :)-(Arb1.*arb1+Alb1.*alb1).*bfr.*eta;
d2U_dof(K, :) = d2U_dof(K, :)+(Alb2.*alb2+Arb2.*arb2).*bfl.*eta;
if flag
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end