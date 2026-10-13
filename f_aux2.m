function [mesh] = f_aux2(basis, Xvt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xp  = basis.xp;
N   = basis.N;
Xv  = Xvt.';
Ke  = size(Xv, 1);
K   = Ke-1;
Xc  = (Xv(2:Ke, 1)+Xv(1:Ke-1, 1))./2;
det = (Xv(2:Ke, 1)-Xv(1:Ke-1, 1));
V0K = [1:Ke-1; 2:Ke].';
V0V = zeros(K, 2);
xdc = zeros(K, N);
for k = 1:2
    V0V(:, k) = Xv(V0K(:, k), :);
end
for i = 1:N
    xdc(:, i) = det(:, 1).*xp(1, i)+V0V(:, 1);
end
%--------------------------------------------------------------------------
mesh.K   = K;
mesh.Xv  = Xv;
mesh.Xc  = Xc;
mesh.det = det;
mesh.V0T = V0V;
mesh.xdc = xdc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end