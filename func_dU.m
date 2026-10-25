function [X] = func_dU(basis, mesh, intg, U_dof, flag, ulb1, urb2, dUh) %#ok<INUSD> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
det     = mesh.det;
bd_     = intg.bd_;
Wbf_m   = intg.Wbf_m;
%--------------------------------------------------------------------------
dU      = U_dof*bd_;
dU_doft = dU*Wbf_m./det;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
N       = basis.N;
K       = mesh.K;
dU_doff = zeros(K, N);
%--------------------------------------------------------------------------
if flag
    %----------------------------------------------------------------------
    Km1               = K-1;
    dU_doff           = zeros(K, N);
    bfl_              = intg.bfl_;
    bfr_              = intg.bfr_;
    bfl_M             = intg.bfl_M;
    bfr_M             = intg.bfr_M;
    bfl_i_M           = bfl_M(1:Km1, :);
    bfr_i_M           = bfr_M(2:K  , :);
    bfrb1_M           = bfr_M(1, :);
    bflb2_M           = bfl_M(K, :);
    %----------------------------------------------------------------------
    ul_i              = U_dof(1:Km1, :)*bfl_;
    ur_i              = U_dof(2:K  , :)*bfr_;
    urb1              = U_dof(1    , :)*bfr_;
    ulb2              = U_dof(  K  , :)*bfl_;
    ujil              = ul_i-ur_i;
    ujir              =-ujil;
    ujb1              = urb1-ulb1;
    ujb2              = ulb2-urb2;
    dU_dofli          = 1./2.*bfl_i_M.*ujil; % o->
    dU_dofri          =-1./2.*bfr_i_M.*ujir; % <-o
    dU_dofb1          =      -bfrb1_M.*ujb1; % <-o
    dU_dofb2          =       bflb2_M.*ujb2; % o->
    dU_doff(1:Km1, :) = dU_doff(1:Km1, :)+dU_dofli;
    dU_doff(2:K  , :) = dU_doff(2:K  , :)+dU_dofri;
    dU_doff(    1, :) = dU_doff(    1, :)+dU_dofb1;
    dU_doff(    K, :) = dU_doff(    K, :)+dU_dofb2;
    dU_doft           = ...
        dU_doft-dU_doff;
    %----------------------------------------------------------------------
end
%--------------------------------------------------------------------------
X.dU_doff = dU_doff; % lift (for debugging)
X.dU_doft = dU_doft;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK #1:
%--------------------------------------------------------------------------
%{
ffll_m        = intg.ffll_m;
fflr_m        = intg.fflr_m;
ffrl_m        = intg.ffrl_m;
ffrr_m        = intg.ffrr_m;
ffll_mt       = ffll_m.';
fflr_mt       = fflr_m.';
ffrl_mt       = ffrl_m.';
ffrr_mt       = ffrr_m.';
U_dofli       = U_dof(1:Km1, :);
U_dofri       = U_dof(2:K  , :);
detli         = det  (1:Km1, 1);
detri         = det  (2:K  , 1);
dU_dofliI     = 1./2.*(U_dofli*ffll_mt-U_dofri*fflr_mt)./detli;
dU_dofriI     =-1./2.*(U_dofri*ffrr_mt-U_dofli*ffrl_mt)./detri;
dU_dofb1I     =-U_dof(1, :)*ffrr_mt./det(1, 1)+bfrb1_M.*ulb1; % i-e w/ - <-o
dU_dofb2I     = U_dof(K, :)*ffll_mt./det(K, 1)-bflb2_M.*urb2; % =   w/ + o->
diffl1        = dU_dofli-dU_dofliI; max(abs(diffl1), [], 'all')
diffr1        = dU_dofri-dU_dofriI; max(abs(diffr1), [], 'all')
diffl2        = dU_dofb1-dU_dofb1I; max(abs(diffl2), [], 'all')
diffr2        = dU_dofb2-dU_dofb2I; max(abs(diffr2), [], 'all')
%--------------------------------------------------------------------------
NN            = N*N;
TNN           = mesh.TNN;
IJdl          = mesh.IJdl; % ll
IJol          = mesh.IJol; % lr
IJor          = mesh.IJor; % rl
IJdr          = mesh.IJdr; % rr
IJb1          = mesh.IJb1; % rrb1
IJb2          = mesh.IJb2; % llb2
JAUX          = mesh.Jaux; % aux. matrix (for debugging)
Lll_it        = zeros(N, N, Km1);
Llr_it        = zeros(N, N, Km1);
Lrl_it        = zeros(N, N, Km1);
Lrr_it        = zeros(N, N, Km1);
RHSL          = zeros(K, N);
for i = 1:Km1
    Lll_it(:, :, i) = 1./2.*ffll_mt./det(i  , 1); % +      o->
    Llr_it(:, :, i) =-1./2.*fflr_mt./det(i  , 1); % -      o->
    Lrl_it(:, :, i) = 1./2.*ffrl_mt./det(i+1, 1); % + = -- <-o
    Lrr_it(:, :, i) =-1./2.*ffrr_mt./det(i+1, 1); % -      <-o
end
Lrrb1t        =-ffrr_mt./det(1, 1); % - <-o
Lllb2t        = ffll_mt./det(K, 1); % + o->
VALL          = zeros(TNN, 1);
Lll_it_res    = reshape(Lll_it, [NN, Km1]);
Llr_it_res    = reshape(Llr_it, [NN, Km1]);
Lrl_it_res    = reshape(Lrl_it, [NN, Km1]);
Lrr_it_res    = reshape(Lrr_it, [NN, Km1]);
VALLdil_NN_1  = Lll_it_res(:);
VALLoil_NN_1  = Llr_it_res(:);
VALLoir_NN_1  = Lrl_it_res(:);
VALLdir_NN_1  = Lrr_it_res(:);
VALLdb1_NN_1  = Lrrb1t    (:);
VALLdb2_NN_1  = Lllb2t    (:);
VALL(IJdl, 1) = VALL(IJdl, 1)+VALLdil_NN_1;
VALL(IJol, 1) = VALLoil_NN_1;
VALL(IJor, 1) = VALLoir_NN_1;
VALL(IJdr, 1) = VALL(IJdr, 1)+VALLdir_NN_1;
VALL(IJb1, 1) = VALL(IJb1, 1)+VALLdb1_NN_1;
VALL(IJb2, 1) = VALL(IJb2, 1)+VALLdb2_NN_1;
MATL          = JAUX(VALL);
RHSL(1, :)    = bfrb1_M.*ulb1; % + = -- <-o
RHSL(K, :)    =-bflb2_M.*urb2; % -      o->
U_doft        = U_dof.';
U             = U_doft(:);
lift1I_1      = MATL*U;
lift1I        = reshape(lift1I_1, [], K).'+RHSL;
diff_3        = dU_doff-lift1I; max(abs(diff_3), [], 'all')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK #2:
%--------------------------------------------------------------------------
%{
visu    = intg.visu;
xdc     = mesh.xdc;
xqc     = intg.xqc;
Wbf_m   = intg.Wbf_m;
dU_dofa = inittype(1, funch, xdc, xqc, Wbf_m); % funch = function handle
figure;
hold on;
plot(dU_doft*visu, '-r');
plot(dU_dofa*visu, 'og');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end