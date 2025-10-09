%% Allocation grid (players × sub-frequencies) — clean under x-axis
% Title outside the axes; f_k labels and column constraints separated with extra bottom space.

% ===== EDIT THESE =====
Nplayers = 3;
Nsubs    = 3;
Pmax     = [10 10 10];    % per-player budgets (rows)
Pcap     = [10 10 10];    % per-sub-frequency caps (cols)
units    = 'W';
mainTitle = '3 Players – 3 Sub-Frequency Channel';
fsBase    = 14;           % base font size
% ======================

% Normalize vectors if scalars provided
if isscalar(Pmax), Pmax = repmat(Pmax,1,Nplayers); end
if isscalar(Pcap), Pcap = repmat(Pcap,1,Nsubs);    end

% Grid geometry
cellW = 1; cellH = 1; W = Nsubs*cellW; H = Nplayers*cellH;

% Figure + layout (title outside axes so it can't collide)
fig = figure('Color','w','Position',[100 100 880 700]);
tlo = tiledlayout(fig,1,1,'TileSpacing','compact','Padding','compact');
sgtitle(tlo, ['\bf ' mainTitle], 'Interpreter','latex','FontSize',fsBase+2);

ax = nexttile(tlo,1); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');

% Margins (extra bottom room to separate label rows)
leftPad = 0.80; rightPad = 2.60; topPad = 0.40; bottomPad = 2.20;
set(ax,'XLim',[-leftPad, W + rightPad], ...
       'YLim',[-bottomPad, H + topPad]);

% Draw grid
for r = 0:Nplayers, plot([0 W],[r r],'k-','LineWidth',1.2); end
for c = 0:Nsubs,    plot([c c],[0 H],'k-','LineWidth',1.2); end

la = {'Interpreter','latex'};  % LaTeX everywhere

% Cell labels P_{ik}
for i = 1:Nplayers
    for k = 1:Nsubs
        xc = (k-0.5)*cellW; yc = H - (i-0.5)*cellH;
        text(xc,yc,sprintf('$P_{%d%d}$',i,k), 'FontSize',fsBase, ...
            'HorizontalAlignment','center','VerticalAlignment','middle', la{:});
    end
end

% Row labels (Player i)
for i = 1:Nplayers
    yc = H - (i-0.5)*cellH;
    text(-0.35,yc, sprintf('\\textbf{Player %d}',i), 'FontSize',fsBase, ...
        'HorizontalAlignment','right','VerticalAlignment','middle', la{:});
end

% Row constraints on the right
for i = 1:Nplayers
    yc  = H - (i-0.5)*cellH;
    txt = sprintf('$\\sum_k P_{%d k} \\le P_{%d}^{\\max} = %g\\,\\mathrm{%s}$', i, i, Pmax(i), units);
    text(W+0.20, yc, txt, 'FontSize',fsBase-1, ...
        'HorizontalAlignment','left','VerticalAlignment','middle', la{:});
end

% --- UNDER X-AXIS LABELS (separated rows) ---
% 1) f_k labels (close to grid)
yF = -0.25;
for k = 1:Nsubs
    text((k-0.5)*cellW, yF, sprintf('$f_{%d}$',k), 'FontSize',fsBase, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', la{:});
end

% 2) Column constraints (further down, TWO LINES to reduce width)
yCols = -1.15;
for k = 1:Nsubs
    xc   = (k-0.5)*cellW;
    line1 = sprintf('$\\sum_i P_{i %d} \\le \\bar{P}_{%d}$', k, k);
    line2 = sprintf('$=\\,%g\\,\\mathrm{%s}$', Pcap(k), units);
    text(xc, yCols, {line1; line2}, 'FontSize',fsBase-1, ...
        'HorizontalAlignment','center','VerticalAlignment','top', la{:});
end

% 3) Legend/note even lower
yNote = -1.85;
noteTxt = ['$\sum_k P_{ik}\!\le\!P_i^{\max}$ (row budgets),  ', ...
           '$\sum_i P_{ik}\!\le\!\bar P_k$ (sub-channel caps)'];
text(W/2, yNote, noteTxt, 'FontSize',fsBase-1, ...
     'HorizontalAlignment','center','VerticalAlignment','top', la{:});

% Export (vector + PNG)
exportgraphics(fig,'allocation_grid_stacked_labels.pdf','ContentType','vector');
exportgraphics(fig,'allocation_grid_stacked_labels.png','Resolution',300);

