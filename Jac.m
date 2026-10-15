function [MAT] = Jac(basis, mesh, intg, flag, c, C_dof, U_dof, d1Ul)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
N            = basis.N;
R            = intg.R;
NN           = N*N;
det          = mesh.det;
TNN          = mesh.TNN;
IJd_         = mesh.IJd_;
Mat          = mesh.Jac;
bd_          = intg.bd_;
bf_          = intg.bf_;
DDr          = intg.DDr;
DFr          = intg.DFr;
FDr          = intg.FDr;
val          = zeros(TNN, 1);
%--------------------------------------------------------------------------
d1Uc         = U_dof*bd_./det;
C            = C_dof*bf_;
d1Ult        = d1Ul.';
d1Uct        = d1Uc.';
Ct           = C.';
d1Ult1       = d1Ult(:);
d1Uct1       = d1Uct(:);
Ct1          =    Ct(:);
Jacd_RK_NN   =-DFr.*d1Uct1-DDr.*Ct1+2.*c.*FDr.*d1Ult1;
Jacd_NN_K    = sum(reshape(Jacd_RK_NN.', NN, R, []), 2);
Jacd_NN_K_1  = Jacd_NN_K(:);
val(IJd_, 1) = Jacd_NN_K_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
K     = mesh.K;
Km1   = K-1;
IJdl  = mesh.IJdl; % ll
IJol  = mesh.IJol; % lr
IJor  = mesh.IJor; % rl
IJdr  = mesh.IJdr; % rr
bfl_  = intg.bfl_;
bfr_  = intg.bfr_;
ul_i  = U_dof(1:Km1, :)*bfl_;
ur_i  = U_dof(2:K  , :)*bfr_;
ul_it = ul_i.';
ur_it = ur_i.';
%--------------------------------------------------------------------------
if flag
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    val = val.*0;



    Jacll_i   =-1./2.*FFll.*ll_it; % (-   = -) o-> (int)
    Jaclr_i   = 1./2.*FFlr.*lr_it; % (--  = +) o-> (ext) 
    Jacrl_i   =-1./2.*FFrl.*rl_it; % (--- = -) <-o (ext)
    Jacrr_i   = 1./2.*FFrr.*rr_it; % (--  = +) <-o (int)
    Jacll_i_1 = Jacll_i(:);
    Jaclr_i_1 = Jaclr_i(:);
    Jacrl_i_1 = Jacrl_i(:);
    Jacrr_i_1 = Jacrr_i(:);
    val(IJdl) = val(IJdl)+Jacll_i_1;
    val(IJol) = Jaclr_i_1;
    val(IJor) = Jacrl_i_1;
    val(IJdr) = val(IJdr)+Jacrr_i_1;
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE JACOBIAN:
%--------------------------------------------------------------------------
MAT = Mat(val);

JACf = full(MAT);
K    = mesh.K;
soll = zeros(N, N, Km1);
solu = zeros(N, N, Km1);
mass = intg.mass;
minv = inv(mass);
for i = 1:K
    if i > 1
        soll(:, :, i) =-JACf((i-1)*N+(1:N), (i-2)*N+(1:N));
    end
    if i < K
        solu(:, :, i) =-JACf((i-1)*N+(1:N), (i)*N+(1:N));
    end
end
val   = load("this.mat");
val1  = val.jac;

diff1 = val1-soll;
diffM = max(abs(diff1), [], 'all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end