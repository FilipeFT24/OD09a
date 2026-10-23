function [] = Plot1(basis, mesh, intg, U_dof, Uh, export, fid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = fieldnames(get(groot, 'factory'));
fc = find(contains(fn, 'Interpreter'));
for i = 1:size(fc, 1)
    set(groot, strrep(fn{fc(i, 1)}, 'factory', 'default'), 'latex');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color       = linspecer(9, 'Qualitative');
flag        = [0, 1];
%--------------------------------------------------------------------------
K           = mesh.K;
Xc          = mesh.Xc;
xlim        = [Xc(1, 1), Xc(K, 1)];
ylim1       = [-0.05, 1.05];
ylim2       = 10.^[-14, -6];
xtick       = 0.00:0.25:1.00;
ytick1      = 0.00:0.25:1.00;
ytick2      = 10.^(-14:2:-6);
xticklabel  = ["$\hphantom{-}0.00$", "$0.25$", "$0.50$", "$0.75$", "$1.00$"];
yticklabel1 = ["$0.00$", "$0.25$", "$0.50$", "$0.75$", "$1.00$"];
yticklabel2 = ["$\hspace{0.65ex}10^{-14}$", "$10^{-12}$", "$10^{-10}$", "$10^{-08}$", "$10^{-06}$"];
xtitle      = "$x$";
ytitle1     = "$u$";
ytitle2     = "$e_{c}$";
layer       = "bottom";
%--------------------------------------------------------------------------
W_          = basis.W_;
det         = mesh.det;
Xq          = intg.Xq;
bf_         = intg.bf_;
visu        = intg.visu;
%--------------------------------------------------------------------------
Udc         = U_dof*visu;
Udq         = U_dof*bf_;
Uaq         = Uh(Xq);
e_q         = Udq-Uaq;
e_c         = e_q.^2*W_.*det(:, 1);
en2         = sqrt(sum(e_c, 1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    f1 = figure('Color', 'w', 'Windowstate', 'maximized', 'Renderer', 'painters');
    hold on;
    Ph = plot(Xc, Udc, '-', 'Color', color(2, :), 'LineWidth', 4);
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
        Ph, ...
        [], [], [])
    if export
        set(f1, 'Units', 'inches', 'PaperSize', 19*[1, 0.935], 'PaperPosition', 20*[0, 0, 1, 0.935]); print(f1, fid(1, 1), '-dpdf', '-r600');
    end
end
if flag(1, 2)
    f2 = figure('Color', 'w', 'Windowstate', 'maximized', 'Renderer', 'painters');
    hold on;
    Ph = plot(Xc, e_c, '-', 'Color', color(2, :), 'LineWidth', 4);
    Plot1NoExp(...
        xlim, ...
        ylim2, ...
        xtick, ...
        ytick2, ...
        xticklabel, ...
        yticklabel2, ...
        layer, ...
        xtitle, ...
        ytitle1, ...
        Ph, ...
        [], [], [])
    set(gca, 'Yscale', 'log');
    if export
        set(f2, 'Units', 'inches', 'PaperSize', 19*[1, 0.935], 'PaperPosition', 20*[0, 0, 1, 0.935]); print(f2, fid(1, 2), '-dpdf', '-r600');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end