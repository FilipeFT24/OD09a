function [] = Plot1(Xl, Xr, Fd, l2, bfunc, vfunc, export, fid) %#ok<INUSL> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = fieldnames(get(groot, 'factory'));
fc = find(contains(fn, 'Interpreter'));
for i = 1:size(fc, 1)
    set(groot, strrep(fn{fc(i, 1)}, 'factory', 'default'), 'latex');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color       = linspecer(9, 'Qualitative');
flag        = [1, 0];
nh          = size(l2{1, 1}, 1);
np          = size(l2{1, 1}, 2);
kh          = nh;
%--------------------------------------------------------------------------
xlim        = [Xl, Xr];
xtick       = 0.00:0.25:1.00;
xticklabel  = ["$\hphantom{-}0.00$", "$0.25$", "$0.50$", "$0.75$", "$1.00$"];
xtitle      = "$x$";
ylim1       = 10.^[-15, -3];
ytick1      = 10.^(-15:3:-3);
yticklabel1 = ["$\hspace{0.65ex}10^{-15}$", "$10^{-12}$", "$10^{-09}$", "$10^{-06}$", "$10^{-03}$"];
ytitle1     = "$|e_{c}|$";
layer       = "bottom";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    %----------------------------------------------------------------------
    %{
    dX = Xr-Xl;
    K1 = size(Fd{1, 1}{kh, 1}, 1);
    V2 = K1*2+1;
    K2 = V2-1;
    X1 = Xl:dX/K1:Xr;
    X2 = Xl:dX/K2:Xr;
    Y2 = cell(1, 2);
    for i = 1:2
        Y2{1, i} = zeros(V2, np);
    end 
    for j = 1:V2
        x2j = X2(1, j);
        c1j = min(sum(x2j >= X1, 2), K1);
        x2g = (x2j-X1(1, c1j))./(X1(1, c1j+1)-X1(1, c1j));
        for i = 1:2
            for k = 1:np
                N = size(bfunc{1, k}, 2);
                for l = 1:N
                    Y2{1, i}(j, k) = ...
                        Y2{1, i}(j, k)+bfunc{1, k}{1, l}(x2g)*Fd{1, i}{kh, k}(c1j, l);
                end
            end
        end
    end
    %}
    dX   = Xr-Xl;
    K2   = size(Fd{1, 1}{kh, 1}, 1);
    dx   = dX/K2;
    dxo2 = dx/2;
    X2   = Xl+dxo2:dx:Xr-dxo2;
    Y2   = cell(1, 2);
    for i = 1:2
        Y2{1, i} = zeros(K2, np);
        for k = 1:np
            Y2{1, i}(:, k) = Fd{1, i}{kh, k}*vfunc{1, k};
        end
    end
    %}
    %----------------------------------------------------------------------
    f1 = figure('Color', 'w', 'Windowstate', 'maximized', 'Renderer', 'painters');
    hold on;
    for i = 1:np
        plot(X2, Y2{1, 1}(:, i), '-', 'Color', color(i, :), 'LineWidth', 4);
        plot(X2, Y2{1, 2}(:, i), ':', 'Color', color(i, :), 'LineWidth', 4);
    end
    set(gca, 'Yscale', 'log');
    Plot1NoExp(...
        xlim, ...
        ylim1, ...
        xtick, ...
        ytick1, ...
        xticklabel, ...
        yticklabel1, ...
        layer, ...
        xtitle, ...
        ytitle1, ...
        [], ...
        [], [], [])
    if export
        set(f1, 'Units', 'inches', 'PaperSize', 19*[1, 0.935], 'PaperPosition', 20*[0, 0, 1, 0.935]); print(f1, fid(1, 1), '-dpdf', '-r600');
    end
    %----------------------------------------------------------------------
end
if flag(1, 2)
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end