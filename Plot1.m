function [] = Plot1(Xl, Xr, Fd, Fda, h, l2, l2a, h1, h1a, FUNC, VISU, export) %#ok<INUSL> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = fieldnames(get(groot, 'factory'));
fc = find(contains(fn, 'Interpreter'));
for i = 1:size(fc, 1)
    set(groot, strrep(fn{fc(i, 1)}, 'factory', 'default'), 'latex');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color        = linspecer(9, 'Qualitative');
flag         = [0, 1];
nh           = size(l2{1, 1}, 1);
np           = size(l2{1, 1}, 2);
kh           = nh;
fid1         = ["Toy (1.1).pdf", "Toy (1.2).pdf"];
fid2         = ["Toy (2.1)-N.pdf", "Toy (2.2)-N.pdf", "Toy (2.3)-N.pdf", "Toy (2.4)-N.pdf", "Toy (2.1)-Y.pdf", "Toy (2.2)-Y.pdf", "Toy (2.3)-Y.pdf", "Toy (2.4)-Y.pdf"];
%--------------------------------------------------------------------------
xlim1        = [Xl, Xr];
xlim2        = [0, 4];
xtick1       = 0.00:0.25:1.00;
xtick2       = 0:1:4;
xticklabel1  = ["$\hphantom{-}0.00$", "$0.25$", "$0.50$", "$0.75$", "$1.00$"];
xticklabel2  = ["$0$", "$1$", "$2$", "$3$", "$4$"];
xtitle1      = "$x$";
xtitle2      = "$log_{2}(h_{0}/h)$";
ylim11       = 10.^[-14, -2];
ylim12       = 10.^[-13, -1];
ylim2        = 10.^[-6, 0];
ytick11      = 10.^(-14:3:-2);
ytick12      = 10.^(-13:3:-1);
ytick2       = 10.^(-6:2:0);
yticklabel11 = ["$\hspace{0.65ex}10^{-14}$", "$10^{-11}$", "$10^{-08}$", "$10^{-05}$", "$10^{-02}$"];
yticklabel12 = ["$\hspace{0.65ex}10^{-13}$", "$10^{-10}$", "$10^{-07}$", "$10^{-04}$", "$10^{-01}$"];
yticklabel2  = ["$\hspace{0.65ex}10^{-06}$", "$10^{-04}$", "$10^{-02}$", "$10^{+00}$"];
ytitle1      = "$|e_{c}|$";
ytitle2      = "$|e_{c}|/|e_{c_{0}}|$";
legendlabel  = ["$k=1$", "$k=2$", "$k=3$", "$k=4$"];
legendloc    = "northeast";
layer        = "bottom";
h0           = 1e-03;
f1           = {@(x) h0.*2.^(-2.*x), @(x) h0.*2.^(-3.*x), @(x) h0.*2.^(-4.*x), @(x) h0.*2.^(-5.*x)}; % l2-norm
%{
f2           = {@(x) h0.*2.^(-1.*x), @(x) h0.*2.^(-2.*x), @(x) h0.*2.^(-3.*x), @(x) h0.*2.^(-4.*x)}; % h1-form
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    %----------------------------------------------------------------------
    Y2  = cell(1, 2);
    Y2a = cell(1, 2);
    %{
    dX  = Xr-Xl;
    K1  = size(Fd{1, 1}{kh, 1}, 1);
    V2  = K1*5+1;
    K2  = V2-1;
    X1  = Xl:dX/K1:Xr;
    X2  = Xl:dX/K2:Xr;
    for i = 1:2
        Y2 {1, i} = zeros(V2, np);
        Y2a{1, i} = zeros(V2, np);
    end 
    for j = 1:V2
        x2j = X2(1, j);
        c1j = min(sum(x2j >= X1, 2), K1);
        x2g = (x2j-X1(1, c1j))./(X1(1, c1j+1)-X1(1, c1j));
        for i = 1:2
            for k = 1:np
                N = size(FUNC{1, k}, 2);
                for l = 1:N
                    Y2 {1, i}(j, k) = Y2 {1, i}(j, k)+FUNC{1, k}{1, l}(x2g)*Fd{1, i}{kh, k}(c1j, l);
                    Y2a{1, i}(j, k) = Y2a{1, i}(j, k)+FUNC{1, k}{1, l}(x2g)*Fd{1, i}{kh, k}(c1j, l);
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
    for i = 1:2
        Y2 {1, i} = zeros(K2, np);
        Y2a{1, i} = zeros(K2, np);
        for k = 1:np
            Y2 {1, i}(:, k) = Fd {1, i}{kh, k}*VISU{1, k};
            Y2a{1, i}(:, k) = Fda{1, i}{kh, k}*VISU{1, k};
        end
    end
    %}
    %----------------------------------------------------------------------
    if ~export
        f11 = figure('Color', 'w', 'Windowstate', 'maximized', 'Renderer', 'painters');
    else
        f11 = figure('Color', 'w', 'Windowstate', 'maximized', 'Renderer', 'painters', 'MenuBar', 'none', 'Toolbar', 'none');
    end
    hold on;
    for i = 1:np
        plot(X2, Y2{1, 1}(:, i), '-', 'Color', color(i, :), 'LineWidth', 4);
        plot(X2, Y2{1, 2}(:, i), ':', 'Color', color(i, :), 'LineWidth', 4);
    end
    set(gca, 'Yscale', 'log');
    Plot1NoExp(...
        xlim1, ...
        ylim11, ...
        xtick1, ...
        ytick11, ...
        xticklabel1, ...
        yticklabel11, ...
        layer, ...
        xtitle1, ...
        ytitle1, ...
        [], ...
        [], [], [])
    if export
        set(f11, 'Units', 'inches', 'PaperSize', 19*[1, 0.935], 'PaperPosition', 20*[0, 0, 1, 0.935]); print(f11, fid1(1, 1), '-dpdf', '-r600');
    end
    %----------------------------------------------------------------------
    if ~export
        f12 = figure('Color', 'w', 'Windowstate', 'maximized', 'Renderer', 'painters');
    else
        f12 = figure('Color', 'w', 'Windowstate', 'maximized', 'Renderer', 'painters', 'MenuBar', 'none', 'Toolbar', 'none');
    end
    hold on;
    hold on;
    for i = 1:np
        plot(X2, Y2a{1, 1}(:, i), '-', 'Color', color(i, :), 'LineWidth', 4);
        plot(X2, Y2a{1, 2}(:, i), ':', 'Color', color(i, :), 'LineWidth', 4);
    end
    set(gca, 'Yscale', 'log');
    Plot1NoExp(...
        xlim1, ...
        ylim12, ...
        xtick1, ...
        ytick12, ...
        xticklabel1, ...
        yticklabel12, ...
        layer, ...
        xtitle1, ...
        ytitle1, ...
        [], ...
        [], [], [])
    if export
        set(f12, 'Units', 'inches', 'PaperSize', 19*[1, 0.935], 'PaperPosition', 20*[0, 0, 1, 0.935]); print(f12, fid1(1, 2), '-dpdf', '-r600');
    end
    %----------------------------------------------------------------------
end
if flag(1, 2)
    %----------------------------------------------------------------------
    x1 = 1;
    x2 = 2;
    dx = 0.01;
    for i = 1:2
        if ~export
            fi = figure('Color', 'w', 'Windowstate', 'maximized', 'Renderer', 'painters');
        else
            fi = figure('Color', 'w', 'Windowstate', 'maximized', 'Renderer', 'painters', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        Mat = [log2(h(1, 1)./h).', l2{1, i}./l2{1, i}(1, :)];
        for j = 1:np
            fplot(@(x) f1{1, j}(x), "-", 'Color', 'k', 'LineWidth', 0.5);
            y1 = f1{1, j}(x1);
            y2 = f1{1, j}(x2);
            plot(...
                [x1, x2, x2, x1], ...
                [y1, y2, y1, y1], 'Color', 'k', 'LineWidth', 0.5);
            text(x2+dx, sqrt(y1*y2), num2str(1+j), 'FontSize', 20, 'VerticalAlignment', 'middle'); % hard-coded for L2-norm
        end
        Ph = cell(1, np);
        for k = 1:np
            Ph{1, k} = plot(Mat(:, 1), Mat(:, 1+k), '-o', 'Color', color(k, :), 'LineWidth', 4, 'MarkerFacecolor', 'w', 'Markersize', 10);
        end
        set(gca, 'Yscale', 'log', 'GridLineStyle', '-', 'GridAlpha', 0.100, 'MinorGridAlpha', 0.075); grid(gca, 'on');
        Plot1NoExp(...
            xlim2, ...
            ylim2, ...
            xtick2, ...
            ytick2, ...
            xticklabel2, ...
            yticklabel2, ...
            layer, ...
            xtitle2, ...
            ytitle2, ...
            Ph, ...
            legendlabel, ...
            legendloc, ...
            1);
        if export
            set(fi, 'Units', 'inches', 'PaperSize', 19*[1, 0.935], 'PaperPosition', 20*[0, 0, 1, 0.935]); print(fi, fid2(1, i), '-dpdf', '-r600');
        end
    end
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end