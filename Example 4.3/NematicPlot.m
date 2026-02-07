% % Plot (3D space)
% function NematicPlot(U, g, X, i)
% % g_left = [g{1,1}, g{1,2}, g{1,3}];
% % g_right = [g{2,1}, g{2,2}, g{2,3}];
% % g_bottom = [g{3,1},...;
% % g_top = [g{4,1},...;
% % g_front = [g{5,1},...;
% % g_back = [g{6,1},...;
% N = length(X);
% u1 = zeros(N+2, N+2, N+2);
% u2 = zeros(N+2, N+2, N+2);
% u3 = zeros(N+2, N+2, N+2);
% u1(2:end-1, 2:end-1, 2:end-1) = reshape(U(1:end/3), [N, N, N]);
% u2(2:end-1, 2:end-1, 2:end-1) = reshape(U(end/3+1:2*end/3), [N, N, N]);
% u3(2:end-1, 2:end-1, 2:end-1) = reshape(U(2*end/3+1:end), [N, N, N]);
% 
% h = 1/(N+1);
% x = 0: h: 1;  y = x;  z = x;
% [X, Y, Z] = meshgrid(x, y, z);
% 
% u1(1, 2:end-1, 2:end-1) = g{6,1};  u2(1, 2:end-1, 2:end-1) = g{6,2};  u3(1, 2:end-1, 2:end-1) = g{6,3};
% u1(end, 2:end-1, 2:end-1) = g{5,1};  u2(end, 2:end-1, 2:end-1) = g{5,2};  u3(end, 2:end-1, 2:end-1) = g{5,3};
% u1(2:end-1, 2:end-1, 1) = g{3,1};  u2(2:end-1, 2:end-1, 1) = g{3,2};  u3(2:end-1, 2:end-1, 1) = g{3,3};
% u1(2:end-1, 2:end-1, end) = g{4,1};  u2(2:end-1, 2:end-1, end) = g{4,2};  u3(2:end-1, 2:end-1, end) = g{4,3};
% u1(2:end-1, end, 2:end-1) = g{2,1};  u2(2:end-1, end, 2:end-1) = g{2,2};  u3(2:end-1, end, 2:end-1) = g{2,3};
% u1(2:end-1, 1, 2:end-1) = g{1,1};  u2(2:end-1, 1, 2:end-1) = g{1,2};  u3(2:end-1, 1, 2:end-1) = g{1,3};
% 
% 
% % 定义 x 和 y 的索引
% I1 = 1:2:N+2;  % x方向上每隔2个点取一个
% I2 = 1:2:N+2;  % y方向上每隔2个点取一个
% 
% % 选择 z 方向上的五个层
% I3 = [1, (N+1)/4+1, 2*(N+1)/4+1, 3*(N+1)/4+1, N+2];
% % 定义每一层的颜色
% colors = [
%     0, 0.5, 0;  % 深绿色
%     1, 0, 0;    % 红色      
%     0, 0, 1;    % 蓝色
%     1, 0.5, 0;  % 橙色
%     0.7, 0.5, 1 % 浅紫色
% ];
% 
% 
% % 绘制每一层并设置颜色
% hold on;
% for k = 1:length(I3)
%     q = quiver3(X(I1, I2, I3(k)), Y(I1, I2, I3(k)), Z(I1, I2, I3(k)), ...
%                 u1(I1, I2, I3(k)), u2(I1, I2, I3(k)), u3(I1, I2, I3(k)), ...
%                 'ShowArrowHead', 'on');
%     q.Color = colors(k, :);  % 设置当前层的颜色
% end
% hold off;% q.Color = [0.52941,0.80784,0.98039]; % line color
% q.LineWidth = 2;  % linewidth
% q.AutoScaleFactor = 1; % ScaleFactor
% view(20,13);
% axis equal
% set(gca, 'LineStyleOrder', {'--'}); % 设置虚线样式
% % 设置隐藏的线条为虚线
% hidden off;
% 
% % 添加网格线
% grid on;
% % set(gca, 'BoxStyle', 'full');
% % xlabel('y','Fontsize',16);
% % ylabel('x','Fontsize',16);
% % zlabel('z','Fontsize',16);
% 
% 
% axis([-0.1 1.1 -0.1 1.1 -0.1 1.1]);  % 设置 x, y, z 坐标轴范围为 [0, 1]
% xticks(0:0.5:1);  % 设置 x 轴刻度
% yticks(0:0.5:1);  % 设置 y 轴刻度
% zticks(0:0.5:1);  % 设置 z 轴刻度
% set(gca, 'Box', 'on');
% 
% folderPath = 'D:\Work relevant\Matlab\IMEX-RK-unit-constraint\Dirichlet_3D\ex1_figure';
% fileName = sprintf('u_%d.fig', i); % 动态生成文件名
% fullFilePath = fullfile(folderPath, fileName);
% saveas(gcf, fullFilePath); % 保存当前图像



% Plot (3D space)
function NematicPlot(U, g, X, i)
% g_left = [g{1,1}, g{1,2}, g{1,3}];
% g_right = [g{2,1}, g{2,2}, g{2,3}];
% g_bottom = [g{3,1},...;
% g_top = [g{4,1},...;
% g_front = [g{5,1},...;
% g_back = [g{6,1},...;
N = length(X);
u1 = zeros(N+2, N+2, N+2);
u2 = zeros(N+2, N+2, N+2);
u3 = zeros(N+2, N+2, N+2);
u1(2:end-1, 2:end-1, 2:end-1) = reshape(U(1:end/3), [N, N, N]);
u2(2:end-1, 2:end-1, 2:end-1) = reshape(U(end/3+1:2*end/3), [N, N, N]);
u3(2:end-1, 2:end-1, 2:end-1) = reshape(U(2*end/3+1:end), [N, N, N]);

h = 1/(N+1);
x = 0: h: 1;  y = x;  z = x;
[X, Y, Z] = meshgrid(x, y, z);

u1(1, 2:end-1, 2:end-1) = g{6,1};  u2(1, 2:end-1, 2:end-1) = g{6,2};  u3(1, 2:end-1, 2:end-1) = g{6,3};
u1(end, 2:end-1, 2:end-1) = g{5,1};  u2(end, 2:end-1, 2:end-1) = g{5,2};  u3(end, 2:end-1, 2:end-1) = g{5,3};
u1(2:end-1, 2:end-1, 1) = g{3,1};  u2(2:end-1, 2:end-1, 1) = g{3,2};  u3(2:end-1, 2:end-1, 1) = g{3,3};
u1(2:end-1, 2:end-1, end) = g{4,1};  u2(2:end-1, 2:end-1, end) = g{4,2};  u3(2:end-1, 2:end-1, end) = g{4,3};
u1(2:end-1, end, 2:end-1) = g{2,1};  u2(2:end-1, end, 2:end-1) = g{2,2};  u3(2:end-1, end, 2:end-1) = g{2,3};
u1(2:end-1, 1, 2:end-1) = g{1,1};  u2(2:end-1, 1, 2:end-1) = g{1,2};  u3(2:end-1, 1, 2:end-1) = g{1,3};


% 定义 x 和 y 的索引
% I1 = 1:2:N+2;  % x方向上每隔2个点取一个
% I2 = 1:2:N+2;  % y方向上每隔2个点取一个
I1 = 1:1:N+2;  % x方向上每隔2个点取一个
I2 = 1:1:N+2;  % y方向上每隔2个点取一个

% 选择 z 方向上的五个层
% I3 = [1, (N+1)/4+1, 2*(N+1)/4+1, 3*(N+1)/4+1, N+2];
I3 =  2*(N+1)/4+1;
% 定义每一层的颜色
colors = [
%     0, 0.5, 0;  % 深绿色
%     1, 0, 0;    % 红色      
    0, 0, 1;    % 蓝色
%     1, 0.5, 0;  % 橙色
%     0.7, 0.5, 1 % 浅紫色
];


% 绘制每一层并设置颜色
hold on;
for k = 1:length(I3)
    q = quiver3(X(I1, I2, I3(k)), Y(I1, I2, I3(k)), Z(I1, I2, I3(k)), ...
                u1(I1, I2, I3(k)), u2(I1, I2, I3(k)), u3(I1, I2, I3(k)), ...
                'ShowArrowHead', 'on');
    q.Color = colors(k, :);  % 设置当前层的颜色
end
hold off;% q.Color = [0.52941,0.80784,0.98039]; % line color
q.LineWidth = 2;  % linewidth
q.AutoScaleFactor = 1; % ScaleFactor
% view(20,13);
view(0,90);
axis equal
set(gca, 'LineStyleOrder', {'--'}); % 设置虚线样式
% 设置隐藏的线条为虚线
hidden off;

% 添加网格线
grid on;
% set(gca, 'BoxStyle', 'full');
% xlabel('y','Fontsize',16);
% ylabel('x','Fontsize',16);
% zlabel('z','Fontsize',16);

axis([-0.1 1.1 -0.1 1.1 0.4 0.6]);  % 设置 x, y, z 坐标轴范围为 [0, 1]
% axis([-0.1 1.1 -0.1 1.1 -0.1 1.1]);  % 设置 x, y, z 坐标轴范围为 [0, 1]
xticks(0:0.5:1);  % 设置 x 轴刻度
yticks(0:0.5:1);  % 设置 y 轴刻度
zticks(0:0.5:1);  % 设置 z 轴刻度
set(gca, 'Box', 'on');

folderPath = 'D:\Work relevant\Matlab\IMEX-RK-unit-constraint\Dirichlet_3D\ex1_figure';
fileName = sprintf('u_%d.fig', i); % 动态生成文件名
fullFilePath = fullfile(folderPath, fileName);
saveas(gcf, fullFilePath); % 保存当前图像

