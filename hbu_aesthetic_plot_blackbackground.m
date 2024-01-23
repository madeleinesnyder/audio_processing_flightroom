x = linspace(0, 2*pi, 1000); % More points for a smoother gradient
y = sin(x);

figure();
set(gcf, 'Color', 'k');
hold on;

colors = [linspace(0, 1, length(x))', linspace(0, 0.5, length(x))', linspace(1, 0, length(x))'];

for i = 1:length(x)-1
    plot(x(i:i+1), y(i:i+1), 'Color', colors(i,:), 'LineWidth', 2);
end

xlabel("x axis",'FontWeight','bold');
ylabel("y axis",'FontWeight','bold');

set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'GridColor', 'w');

set(gca, 'FontName', 'Arial'); % Changes the font of the axis labels and tick labels
set(get(gca, 'XLabel'), 'FontName', 'Arial'); % Changes the font of the X-axis label
set(get(gca, 'YLabel'), 'FontName', 'Arial'); % Changes the font of the Y-axis label
set(get(gca, 'Title'), 'FontName', 'Arial'); % Changes the font of the title

set(gca, 'FontSize', 12); % Changes the font size of the axis labels and tick labels
set(get(gca, 'XLabel'), 'FontSize', 14); % Changes the font size of the X-axis label
set(get(gca, 'YLabel'), 'FontSize', 14); % Changes the font size of the Y-axis label
set(get(gca, 'Title'), 'FontSize', 16); % Changes the font size of the title

set(get(gca, 'Title'), 'FontWeight', 'bold');

