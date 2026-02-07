% % Plot (3D space)
% function NematicPlot(U, X, i)
% N = length(X);
% h = 1/(N-1);
% x = 0: h: 1;  y = x;  z = x;
% [X, Y, Z] = meshgrid(x, y, z);
% 
% u1 = zeros(N, N, N);
% u2 = zeros(N, N, N);
% u3 = zeros(N, N, N);
% 
% u1(:,:,2:end-1) = reshape(U(1:end/3), [N, N, N-2]);
% u2(:,:,2:end-1) = reshape(U(end/3+1:2*end/3), [N, N, N-2]);
% u3(:,:,2:end-1) = reshape(U(2*end/3+1:end), [N, N, N-2]);
% 
% u1(:,:,1) = ones(N, N);
% u2(:,:,end) = ones(N, N);
% 
% % 定义 x 和 y 的索引
% I1 = 1:2:N;  % x方向上每隔2个点取一个
% I2 = 1:2:N;  % y方向上每隔2个点取一个
% 
% % 选择 z 方向上的五个层
% I3 = [1, (N-1)/4+1, 2*(N-1)/4+1, 3*(N-1)/4+1, N];
% 
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
% axis([0 1 0 1 0 1]);  % 设置 x, y, z 坐标轴范围为 [0, 1]
% xticks(0:0.5:1);  % 设置 x 轴刻度
% yticks(0:0.5:1);  % 设置 y 轴刻度
% zticks(0:0.5:1);  % 设置 z 轴刻度
% set(gca, 'Box', 'on');
% 
% folderPath = 'D:\Work relevant\Matlab\IMEX-RK-unit-constraint\Dirichlet_Neumann_3D\ex1_figure';
% fileName = sprintf('u_%d.fig', i); % 动态生成文件名
% fullFilePath = fullfile(folderPath, fileName);
% saveas(gcf, fullFilePath); % 保存当前图像
% 


% Plot (3D space)
function NematicPlot(U, X, i)
N = length(X);
h = 1/(N-1);
x = 0: h: 1;  y = x;  z = x;
[X, Y, Z] = meshgrid(x, y, z);

u1 = zeros(N, N, N);
u2 = zeros(N, N, N);
u3 = zeros(N, N, N);

u1(:,:,2:end-1) = reshape(U(1:end/3), [N, N, N-2]);
u2(:,:,2:end-1) = reshape(U(end/3+1:2*end/3), [N, N, N-2]);
u3(:,:,2:end-1) = reshape(U(2*end/3+1:end), [N, N, N-2]);

u1(:,:,1) = ones(N, N);
u2(:,:,end) = ones(N, N);

% % 定义 x 和 y 的索引
% I1 = (N-1)/2 + 1;  
% I2 = (N-1)/2 + 1; 
% 
% % 选择 z 方向
% I3 = 1:1:N;

% % 定义每一层的颜色
% colors = [
%     0, 0.5, 0;  % 深绿色
%     1, 0, 0;    % 红色      
%     0, 0, 1;    % 蓝色
%     1, 0.5, 0;  % 橙色
%     0.7, 0.5, 1 % 浅紫色
% ];


% % 绘制每一层并设置颜色
% hold on;
% % for k = 1:length(I3)
%     q = quiver3(X(I1, I2, I3), Y(I1, I2, I3), Z(I1, I2, I3), ...
%                 u1(I1, I2, I3), u2(I1, I2, I3), u3(I1, I2, I3), ...
%                 'ShowArrowHead', 'on');
% %     q.Color = colors(k, :);  % 设置当前层的颜色
% % end
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
% axis([0 1 0 1 0 1]);  % 设置 x, y, z 坐标轴范围为 [0, 1]
% xticks(0:0.5:1);  % 设置 x 轴刻度
% yticks(0:0.5:1);  % 设置 y 轴刻度
% zticks(0:0.5:1);  % 设置 z 轴刻度
% set(gca, 'Box', 'on');

figure;
c = 0.05; % vector ScaleFactor
% z = 0;

% 绘制圆柱体
for k = 13
    for j = 13
        for l = 1:N
        a = [x(k), y(j), z(l)]; 
        b = [x(k)+c*u1(k,j,l), y(j)+c*u2(k,j,l), z(l)+c*u3(k,j,l)];
        myplotcylinder(a, b);
        hold on
        end
    end
end

% for k = 12
%     for j = 12
%         for l = 1:N
%         a = [x(k), y(j), z(l)]; 
%         b = [x(k)+c*u1(k,j,l), y(j)+c*u2(k,j,l), z(l)+c*u3(k,j,l)];
%         myplotcylinder(a, b);
%         hold on
%         end
%     end
% end
% 
% for k = 14
%     for j = 14
%         for l = 1:N
%         a = [x(k), y(j), z(l)]; 
%         b = [x(k)+c*u1(k,j,l), y(j)+c*u2(k,j,l), z(l)+c*u3(k,j,l)];
%         myplotcylinder(a, b);
%         hold on
%         end
%     end
% end

% 光照设置
light('Position', [0.5, 0.5, 2], 'Style', 'infinite'); % 光源放置在顶部偏侧
lighting gouraud; % 光照模型为Gouraud，增强光影效果
material shiny;   % 设置为 shiny 模式，增强光亮

% 坐标轴设置
axis equal
% xticks(0:0.5:1);  % 设置 X 轴的刻度
% yticks(0:0.5:1);  % 设置 Y 轴的刻度
% zticks(0:0.5:1);   % 设置 Z 轴的刻度

% 设置观察视角
view(140, 30); % 改变视角，az = 135°, el = 45°
%set(gca, 'ZTickLabel', []); % 移除刻度标签
set(gca, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []); % 移除刻度标签
zlabel('z','Fontsize',16);
%set(gca, 'Box', 'on');
xlim([0.4,0.6]); ylim([0.4,0.6]); zlim([-0.01,1.01]);
folderPath = 'ex1_figure';
fileName = sprintf('u_%d.fig', i); % 动态生成文件名
fullFilePath = fullfile(folderPath, fileName);
saveas(gcf, fullFilePath); % 保存当前图像
end



