% LL_HeatFlow_Scheme.m(仅对保模1加了拉格朗日乘子)
% 基于论文 arXiv:2206.02882 对方程 (2.9) 的数值格式实现
% 求解方程: m_t = Laplace(m) + lambda * m
% 约束: |m| = 1
% 边界: Zero Neumann on [-0.5, 0.5]^2

clear; clc; close all;

    for i=1:5

    %% ================= 1. Parameters =================
    Nx = 64; Ny = 64;h=1/Nx;
    Lx = 1; Ly = 1;
    hx = Lx/Nx; hy = Ly/Ny;
    
    dt = 128e-5/2^i; 
    T  = 0.01024; 
    Nt = round(T/dt);

% 网格生成 (节点中心, 范围 [-0.5, 0.5])
x = linspace(-0.5, 0.5, Nx+1);
y = linspace(-0.5, 0.5, Ny+1);
[X, Y] = meshgrid(x, y);

%% ================= 2. 算子矩阵构建 =================
% 构建支持 Neumann 边界的拉普拉斯矩阵 L
[L, Dx, Dy] = get_neumann_laplacian_matrix(Nx, Ny, hx, hy);

N_dof = numel(X);
I = speye(N_dof);

% Step 1 的左端项矩阵 A = (I/dt - 0.5*Delta)
% 方程: (m_tilde - m)/dt = 0.5*L*(m_tilde + m) + lambda*m
% 整理: (I/dt - 0.5*L) * m_tilde = (I/dt + 0.5*L) * m + lambda*m
A = (1/dt)*I - 0.5*L;

%% ================= 3. 初始条件 =================
    m1 = 0.3 * sin(pi*X) .* sin(pi*Y);
    m2 = 0.3 * sin(3*pi*X) .* sin(pi*Y);
    m3 = 1.0 + 0.2 * cos(2*pi*X) .* cos(2*pi*Y);

[m1, m2, m3] = normalize_m(m1, m2, m3);

% --- [新增] 定义梯形积分权重 W ---
    % 内部点权重=1，边界点=0.5，角点=0.25
    wx = ones(1, Nx+1); wx(1) = 0.5; wx(end) = 0.5;
    wy = ones(Ny+1, 1); wy(1) = 0.5; wy(end) = 0.5;
    W = wy * wx; % 2D 权重矩阵     
    
    % 初始 Lambda = |grad m|^2
    Lm1 = reshape(L * m1(:), size(m1));
    Lm2 = reshape(L * m2(:), size(m2));
    Lm3 = reshape(L * m3(:), size(m3));
    
    lambda = -sum(sum( W.*(m1.*Lm1 + m2.*Lm2 + m3.*Lm3) )) * hx * hy;

Energy = zeros(Nt+1, 1);
Energy(1) = discrete_energy(m1,m2,m3,hx,hy,L,W);

fprintf('Simulating Heat Flow (Eq 2.9)... dt=%.1e\n', dt);

%% ================= 4. 时间演化 =================
tic;
for n = 1:Nt
    %% --- Step 1: Predictor ---
    % 求解 tilde_m
    % RHS = m/dt + 0.5*L*m + lambda*m
    
    % 计算 Laplacian 项
    Lm1 = reshape(L * m1(:), size(m1));
    Lm2 = reshape(L * m2(:), size(m2));
    Lm3 = reshape(L * m3(:), size(m3));
    
    % 显式非线性源项 lambda^n * m^n
    src1 = lambda .* m1;
    src2 = lambda .* m2;
    src3 = lambda .* m3;
    
    rhs1 = m1(:)/dt + 0.5*Lm1(:) + src1(:);
    rhs2 = m2(:)/dt + 0.5*Lm2(:) + src2(:);
    rhs3 = m3(:)/dt + 0.5*Lm3(:) + src3(:);
    
    % 求解线性方程组
    m1_tilde = reshape(A \ rhs1, size(m1));
    m2_tilde = reshape(A \ rhs2, size(m2));
    m3_tilde = reshape(A \ rhs3, size(m3));
    
    %% --- Step 2: Corrector ---
    % 公式: (m^{n+1} - tilde_m)/dt = 0.5 * (lambda^{n+1} * m^{n+1} - lambda^n * m^n)
    % 移项: m^{n+1} * (1 - 0.5*dt*lambda^{n+1}) = tilde_m - 0.5*dt*lambda^n*m^n
    % 令 R = tilde_m - 0.5*dt*lambda^n*m^n
    % 则 m^{n+1} 方向与 R 相同。由于 |m^{n+1}|=1，故 m^{n+1} = R / |R|
    
    term_old1 = 0.5 * dt * lambda .* m1;
    term_old2 = 0.5 * dt * lambda .* m2;
    term_old3 = 0.5 * dt * lambda .* m3;
    
    R1 = m1_tilde - term_old1;
    R2 = m2_tilde - term_old2;
    R3 = m3_tilde - term_old3;
    
    norm_R = sqrt(R1.^2 + R2.^2 + R3.^2);
    
    % 更新 m^{n+1}
    m1 = R1 ./ norm_R;
    m2 = R2 ./ norm_R;
    m3 = R3 ./ norm_R;
    
    % 更新 lambda^{n+1} (根据模长关系导出)
    % |R| = 1 - 0.5*dt*lambda^{n+1}  => lambda^{n+1} = 2*(1 - |R|) / dt
    lambda = 2 * (1 - norm_R) / dt;
    
    %% --- 记录能量 ---
    if mod(n, 1) == 0
        % [dm1x, dm1y] = calc_grad_neumann(m1, hx, hy);
        % [dm2x, dm2y] = calc_grad_neumann(m2, hx, hy);
        % [dm3x, dm3y] = calc_grad_neumann(m3, hx, hy);
        % grad_sq_new = dm1x.^2 + dm1y.^2 + dm2x.^2 + dm2y.^2 + dm3x.^2 + dm3y.^2;
        Energy(n+1) = discrete_energy(m1,m2,m3,hx,hy,L,W);
        
        % 简单的进度条
        if mod(n, 1) == 0
            fprintf('Step %d/%d, Energy = %.5f\n', n, Nt, Energy(n+1));
        end
    end
end
    toc;

    load("data_1/BDF4_T_0p01024.mat");
    error_L2 =  (m1_ex-m1).^2 + (m2_ex-m2).^2 + (m3_ex-m3).^2;
    Error_L2(i) =  h^2*sqrt( sum(error_L2(:)) );
    Error_inf(i) = ( norm(m1(:)-m1_ex(:),inf)+norm(m2(:)-m2_ex(:),inf)+norm(m3(:)-m3_ex(:),inf) )/3;

    end

    order_L2_1 = log(Error_L2(1:4)./Error_L2(2:5))/log(2)
    order_L2_2 = log((Error_L2(1:3)-Error_L2(2:4))./(Error_L2(2:4)-Error_L2(3:5)))/log(2)

    order_inf_1 = log(Error_inf(1:4)./Error_inf(2:5))/log(2)
    order_inf_2 = log((Error_inf(1:3)-Error_inf(2:4))./(Error_inf(2:4)-Error_inf(3:5)))/log(2)

%% ================= 5. 绘图 =================
figure;
subplot(1,2,1);
plot(0:dt:T, Energy(1:1:end), 'LineWidth', 1.5);
xlabel('Time'); ylabel('Energy'); title('Dirichlet Energy Decay');
grid on;

subplot(1,2,2);
step = 2;
quiver(X(1:step:end, 1:step:end), Y(1:step:end, 1:step:end), ...
       m1(1:step:end, 1:step:end), m2(1:step:end, 1:step:end));
axis equal; xlim([-0.5 0.5]); ylim([-0.5 0.5]);
title(['Solution at T=' num2str(T)]);
drawnow;

%% ================= 辅助函数 =================

function [m1, m2, m3] = normalize_m(m1, m2, m3)
    nn = sqrt(m1.^2 + m2.^2 + m3.^2);
    m1 = m1 ./ nn; m2 = m2 ./ nn; m3 = m3 ./ nn;
end

function [gx, gy] = calc_grad_neumann(f, hx, hy)
    % 二阶精度中心差分梯度，适配 Zero Neumann
    [Ny, Nx] = size(f);
    gx = zeros(Ny, Nx); gy = zeros(Ny, Nx);
    
    % x方向: 内部用中心差分
    gx(:, 2:end-1) = (f(:, 3:end) - f(:, 1:end-2)) / (2*hx);
    % x方向: 边界 Neumann => 导数为0
    gx(:, 1) = 0; gx(:, end) = 0;
    
    % y方向: 内部用中心差分
    gy(2:end-1, :) = (f(3:end, :) - f(1:end-2, :)) / (2*hy);
    % y方向: 边界 Neumann => 导数为0
    gy(1, :) = 0; gy(end, :) = 0;
end

function [L, Dx, Dy] = get_neumann_laplacian_matrix(Nx, Ny, hx, hy)
    % 构建基于 kron 的稀疏拉普拉斯矩阵，零纽曼边界条件
    % 使用 Ghost Point 方法: u_{-1} = u_{1} => u_{xx} = 2(u_{1}-u_{0})/h^2
    
    % 1D Matrix construction for x
    e = ones(Nx+1, 1);
    Lx_1d = spdiags([e -2*e e], -1:1, Nx+1, Nx+1);
    % Neumann fix: 2*u_{inner} - 2*u_{boundary}
    Lx_1d(1,2) = 2; Lx_1d(Nx+1, Nx) = 2; 
    Lx_1d = Lx_1d / hx^2;

    % 1D Matrix construction for y
    e = ones(Ny+1, 1);
    Ly_1d = spdiags([e -2*e e], -1:1, Ny+1, Ny+1);
    % Neumann fix
    Ly_1d(1,2) = 2; Ly_1d(Ny+1, Ny) = 2; 
    Ly_1d = Ly_1d / hy^2;

    Ix = speye(Nx+1);
    Iy = speye(Ny+1);

    % 2D Laplacian = Lx kron Iy + Ix kron Ly
    L = kron(Iy, Lx_1d) + kron(Ly_1d, Ix);
    
    % 仅返回空值，因为梯度我们在外部函数 calc_grad_neumann 中处理
    Dx = []; Dy = [];
end