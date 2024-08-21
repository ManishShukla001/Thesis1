% Read the data
sigmab = readmatrix('D:\Dec\LCS Code\All Results\Forward\NFSigma1_810.txt');

% Create figure
f = figure('Visible', 'on', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Data preprocessing
sigmab_filtered = sigmab;
sigmab_filtered(sigmab_filtered < 0.00005) = 0.00005;
sigmab_filtered(sigmab_filtered > 1.5) = 1.5;
sigmab_filtered(sigmab_filtered < 0.00008 | sigmab_filtered > 0.00014) = 0.00005;

% Determine colorbar range automatically
min_val = min(sigmab_filtered(:));
max_val = max(sigmab_filtered(:));

% Create the plot
p = pcolor(XXb, YYb, sigmab_filtered);
axis equal tight
colorbar
shading interp
lighting phong
camlight('left')

% Customize colorbar
hcb = colorbar;
caxis([min_val max_val])
set(get(hcb, 'Title'), 'String', 'FTLE', 'FontSize', 15, 'FontWeight', 'bold', 'Color', 'k')
colormap('jet')

% Enhance plot appearance
set(gca, 'FontSize', 12)
xlabel('X', 'FontSize', 14)
ylabel('Y', 'FontSize', 14)
title('LCS Results', 'FontSize', 16)

% Optional: Add grid lines
grid on

% add scatter plot
% hold on
% s = scatter(filx, fily, 'filled');
% s.LineWidth = 0.1;
% s.MarkerEdgeColor = 'w';
% s.MarkerFaceColor = [1 1 1];

% Save the figure as a high-resolution image
print(f, 'LCS_Results', '-dpng', '-r300')
