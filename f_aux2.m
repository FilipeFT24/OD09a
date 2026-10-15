function [mesh] = f_aux2(basis, Xvt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N    = basis.N;
NN   = N*N;
xp   = basis.xp;
Xv   = Xvt.';
Ke   = size(Xv, 1);
Kem1 = Ke-1;
K    = Kem1;
Xc   = (Xv(2:Ke, 1)+Xv(1:Kem1, 1))./2;
det  = (Xv(2:Ke, 1)-Xv(1:Kem1, 1));
V0K  = [1:Kem1; 2:Ke].';
V0V  = zeros(K, 2);
xdc  = zeros(K, N);
for k = 1:2
    V0V(:, k) = Xv(V0K(:, k), :);
end
for i = 1:N
    xdc(:, i) = det(:, 1).*xp(1, i)+V0V(:, 1);
end
%--------------------------------------------------------------------------
Km1          = K-1;
Km2          = K-2;
T            = 3*K-2;
Tm1          = T-1;
KN           = K*N;
TNN          = T*NN;
Km1NN        = Km1*NN;
im           = zeros(T  , 1);
jm           = zeros(T  , 1);
XT           = zeros(TNN, 2);
im(  1:2, 1) = [1  , 1];
jm(  1:2, 1) = [1  , 2];
im(Tm1:T, 1) = [K  , K];
jm(Tm1:T, 1) = [Km1, K];
for i = 1:Km2
    im(2+3*(i-1)+1:2+3*i, 1) = i+1;
    jm(2+3*(i-1)+1:2+3*i, 1) = i:i+2;
end
for i = 1:T
    XT((i-1)*NN+(1:NN), :) = ...
        [repelem((im(i, 1)-1)*N+(1:N), N); repmat((jm(i, 1)-1)*N+(1:N), 1, N)].';
end
XT_i = XT(:, 1);
XT_j = XT(:, 2);
IJd_ = reshape(((3+3*((1:K  ).'-2))*NN+(1:NN)).', [], 1); % diag.
IJor = reshape(((2+3*((2:K  ).'-2))*NN+(1:NN)).', [], 1); % diag. lower (off-diag. blocks of r cell)
IJol = reshape(((4+3*((1:Km1).'-2))*NN+(1:NN)).', [], 1); % diag. upper (off-diag. blocks of l cell)
IJdl = IJd_( 0*NN+(1:Km1NN), 1);
IJdr = IJd_( 1*NN+(1:Km1NN), 1);
IJb1 = IJd_( 0*NN+(1:NN)   , 1);
IJb2 = IJd_(Km1NN+(1:NN)   , 1);
%--------------------------------------------------------------------------
mesh.K    = K;
mesh.Xc   = Xc;
mesh.Xv   = Xv;
mesh.det  = det;
mesh.V0T  = V0V;
mesh.xdc  = xdc;
mesh.TNN  = TNN;
mesh.IJd_ = IJd_; % cell
mesh.IJdl = IJdl; % fint ll
mesh.IJdr = IJdr; % face rr
mesh.IJol = IJol; % face lr
mesh.IJor = IJor; % face rl
mesh.IJb1 = IJb1; % bnd1
mesh.IJb2 = IJb2; % bnd2
mesh.Jaux = spinit(XT_i, XT_j, [KN, KN]);
mesh.Jac  = spinit(XT_i, XT_j, [KN, KN]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end