% Plot (2D space)
function NematicPlot(u1,u2,u3, X, Y, i)
N = (length(u1(:)))^(1/2);
% u1 = reshape(U(1:end/3), [N, N]);
% u2 = reshape(U(end/3+1:2*end/3), [N, N]);
% u3 = reshape(U(2*end/3+1:end), [N, N]);

figure;
Z = 0*X; % mesh grid
I = round((N-1)/4+1 :1 :(N-1)/4*3);
q=quiver3(X(I,I), Y(I,I), Z(I,I), u1(I,I), u2(I,I), u3(I,I), 'ShowArrowHead','on'); % Plot Line fields
q.Color = [0.52941,0.80784,0.98039]; % line color
q.LineWidth = 2;  % linewidth
q.AutoScaleFactor = 1; % ScaleFactor
axis equal
xticks(-0.5:0.2:0.5);  % 设置 X 轴的刻度
yticks(-0.5:0.2:0.5);  % 设置 Y 轴的刻度
% view(0,90);
view(-40,40);
set(gca, 'FontSize', 8);  % 设置坐标轴刻度字体大小
% set (gca,'XTickLabel', []); set (gca,'YTickLabel', []);
set(gca, 'Box', 'off');

%folderPath = '/Users/youyou/Documents/MATLAB/IMEX-RK-unit-constraint/LM2_CN';
%fileName = sprintf('u_%d.fig', i); % 动态生成文件名
%fullFilePath = fullfile(folderPath, fileName);
%saveas(gcf, fileName); % 保存当前图像