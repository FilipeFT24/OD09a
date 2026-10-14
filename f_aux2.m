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
    XT((i-1)*NN+(1:NN), :) = [repelem((im(i, 1)-1)*N+(1:N), N); repmat((jm(i, 1)-1)*N+(1:N), 1, N)].';
end
XTi = XT(:, 1);
XTj = XT(:, 2);
IJd = reshape(((3+3*((1:K  ).'-2))*NN+(1:NN)).', [], 1); % diag.
IJl = reshape(((2+3*((2:K  ).'-2))*NN+(1:NN)).', [], 1); % diag. lower (off-diag. blocks)
IJu = reshape(((4+3*((1:Km1).'-2))*NN+(1:NN)).', [], 1); % diag. upper (off-diag. blocks)
%--------------------------------------------------------------------------
mesh.K   = K;
mesh.Xc  = Xc;
mesh.Xv  = Xv;
mesh.det = det;
mesh.V0T = V0V;
mesh.xdc = xdc;
mesh.TNN = TNN;
mesh.IJd = IJd;
mesh.IJl = IJl;
mesh.IJu = IJu;
mesh.Jac = spinit(XTi, XTj, [KN, KN]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end