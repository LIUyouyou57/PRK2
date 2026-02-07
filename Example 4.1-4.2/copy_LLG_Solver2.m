% function HeatFlow_Solver() 
    % (对保模1和保能量耗散都加了拉格朗日乘子)
    % 基于 arXiv:2206.02882 求解方程 (2.9)
    % 区域: [-1/2, 1/2]^2
    % BC: Zero Neumann
    
    clear; clc; close all;

    for i=1:5

    %% ================= 1. Parameters =================
    Nx = 64; Ny = 64; h=1/Nx;
    Lx = 1; Ly = 1;
    hx = Lx/Nx; hy = Ly/Ny;
    
    dt = 128e-5/2^i; 
    T  = 0.01024; 
    Nt = round(T/dt);
    
    %% ================= 2. Grid & Operators =================
    % 使用节点网格 (Node-based) 以匹配 Neumann 边界处理
    x = linspace(-Lx/2, Lx/2, Nx+1);
    y = linspace(-Ly/2, Ly/2, Ny+1);
    [X, Y] = meshgrid(x, y);
    
    % 生成二阶精度 Neumann 拉普拉斯矩阵
    [L, Dx, Dy] = get_neumann_laplacian_matrix(Nx, Ny, hx, hy);
    
    N_dof = (Nx+1)*(Ny+1);
    I_mat = speye(N_dof);
    
    % Step 1 隐式矩阵: (I/dt - 0.5*Delta)
    A_step1 = (1/dt)*I_mat - 0.5 * L;

    %% ================= 3. Initial Condition =================
    % % 之前的涡旋初值
    % LL = sqrt(X.^2 + Y.^2);
    % AA = (1 - 2*LL).^4;
    % BB = AA.^2 + LL.^2;
    % mask = (LL <= 0.5);
    % 
    % % 避免除零
    % BB(BB==0) = 1; 
    % 
    % m1 = mask.*(2*X.*AA./BB);
    % m2 = mask.*(2*Y.*AA./BB);
    % m3 = mask.*((AA.^2 - LL.^2)./BB) - (~mask).*(ones(size(LL)));

    m1 = 0.3 * sin(pi*X) .* sin(pi*Y);
    m2 = 0.3 * sin(3*pi*X) .* sin(pi*Y);
    m3 = 1.0 + 0.2 * cos(2*pi*X) .* cos(2*pi*Y);
    
    [m1, m2, m3] = normalize_m(m1, m2, m3);
    
    % 初始 Lambda = |grad m|^2
    Lm1 = reshape(L * m1(:), size(m1));
    Lm2 = reshape(L * m2(:), size(m2));
    Lm3 = reshape(L * m3(:), size(m3));
    
    lambda = -sum(sum( m1.*Lm1 + m2.*Lm2 + m3.*Lm3 )) * hx * hy;
    % 初始能量
    Energy = zeros(Nt+1, 1);
    Energy(1) = calc_energy(m1, m2, m3, hx, hy, L);
    
    fprintf('Start Simulation... Nx=%d, dt=%.1e\n', Nx, dt);

    %% ================= 4. Time Loop =================
    tic;
    for n = 1:Nt
        %% --- Step 1: Predictor ---
        % Formula: (m_tilde - m)/dt = 0.5*Delta(m_tilde + m) + lambda*m
        % LHS: (1/dt - 0.5*L) * m_tilde
        % RHS: (1/dt + 0.5*L) * m + lambda*m
        
        Lap_m1 = reshape(L * m1(:), size(m1));
        Lap_m2 = reshape(L * m2(:), size(m2));
        Lap_m3 = reshape(L * m3(:), size(m3));
        
        rhs1 = m1(:)/dt + 0.5*Lap_m1(:) + lambda(:).*m1(:);
        rhs2 = m2(:)/dt + 0.5*Lap_m2(:) + lambda(:).*m2(:);
        rhs3 = m3(:)/dt + 0.5*Lap_m3(:) + lambda(:).*m3(:);
        
        m1_tilde = reshape(A_step1 \ rhs1, size(m1));
        m2_tilde = reshape(A_step1 \ rhs2, size(m2));
        m3_tilde = reshape(A_step1 \ rhs3, size(m3));

        %% --- Step 2: Corrector ---
        % R = m_tilde - 0.5*dt*lambda*m
        % m_hat = R / |R|
        % lambda_new = 2/dt * (1 - |R|)
        
        factor = 0.5 * dt * lambda;
        R1 = m1_tilde - factor .* m1;
        R2 = m2_tilde - factor .* m2;
        R3 = m3_tilde - factor .* m3;
        
        norm_R = sqrt(R1.^2 + R2.^2 + R3.^2);
        
        m1_hat = R1 ./ norm_R;
        m2_hat = R2 ./ norm_R;
        m3_hat = R3 ./ norm_R;
        
        lambda = (2/dt) * (1 - norm_R);

        %% --- Step 3: Energy Dissipation (Scalar xi) ---
        % Calculate Dissipation term
        % vector V = 0.5 * (m_hat + m_n)
        v1 = 0.5 * (m1_hat + m1);
        v2 = 0.5 * (m2_hat + m2);
        v3 = 0.5 * (m3_hat + m3);
        
        Lv1 = reshape(L*v1(:), size(v1));
        Lv2 = reshape(L*v2(:), size(v2));
        Lv3 = reshape(L*v3(:), size(v3));
        
        % Cross Product: V x Delta V
        c1 = v2.*Lv3 - v3.*Lv2;
        c2 = v3.*Lv1 - v1.*Lv3;
        c3 = v1.*Lv2 - v2.*Lv1;
        
        Diss = sum(sum(c1.^2 + c2.^2 + c3.^2)) * hx * hy;
        
        % Target Energy
        E_old = Energy(n);
        E_target = E_old - dt * Diss;
        
        % Solve for xi using Secant Method
        % Func(xi) = Energy(normalize(m_hat + xi)) - E_target
        func_xi = @(xi) energy_wrapper(xi, m1_hat, m2_hat, m3_hat, hx, hy, L) - E_target;
        
        xi0 = 0;
        xi1 =  -dt^2;% Guess
        f0 = func_xi(xi0);
        f1 = func_xi(xi1);
        
        iter = 0;
        while abs(f1)>1e-12
            iter = iter + 1;
            xi_new = xi1 - f1*(xi1-xi0)/(f1-f0);
            xi0 = xi1; f0 = f1;
            xi1 = xi_new; f1 = func_xi(xi1);
            xi_sol = xi1;
        end

        
        %% --- Update ---
        t1 = m1_hat + xi_sol;
        t2 = m2_hat + xi_sol;
        t3 = m3_hat + xi_sol;
        nt = sqrt(t1.^2 + t2.^2 + t3.^2);
        
        m1 = t1 ./ nt;
        m2 = t2 ./ nt;
        m3 = t3 ./ nt;
        
        Energy(n+1) = calc_energy(m1, m2, m3, hx, hy, L);
        
        if mod(n, 1) == 0
            fprintf('Step %d/%d, E=%.5e, |xi|=%.2e, iter number=%d\n', n, Nt, Energy(n+1), abs(xi_sol), iter);
        end
    end
    toc;

    load("data_1/LLG_Solver2_N_128_T_0p01024.mat");
    error_L2 =  (m1_ex-m1).^2 + (m2_ex-m2).^2 + (m3_ex-m3).^2;
    Error_L2(i) =  h^2*sqrt( sum(error_L2(:)) );
    Error_inf(i) = ( norm(m1(:)-m1_ex(:),inf)+norm(m2(:)-m2_ex(:),inf)+norm(m3(:)-m3_ex(:),inf) )/3;

    end

    order_L2_1 = log(Error_L2(1:4)./Error_L2(2:5))/log(2)
    order_L2_2 = log((Error_L2(1:3)-Error_L2(2:4))./(Error_L2(2:4)-Error_L2(3:5)))/log(2)

    order_inf_1 = log(Error_inf(1:4)./Error_inf(2:5))/log(2)
    order_inf_2 = log((Error_inf(1:3)-Error_inf(2:4))./(Error_inf(2:4)-Error_inf(3:5)))/log(2)
%     %% ================= 5. Plots =================
%     figure;
%     plot((0:Nt)*dt, Energy, 'LineWidth', 1.5);
%     xlabel('t'); ylabel('Energy'); title('Energy Decay (Heat Flow)');
%     grid on;
% 
%     figure;
%     subplot(1,2,1);
%     quiver(X(1:2:end,1:2:end), Y(1:2:end,1:2:end), m1(1:2:end,1:2:end), m2(1:2:end,1:2:end));
%     title('Final Spin Field (m1, m2)'); axis equal; xlim([-0.5 0.5]); ylim([-0.5 0.5]);
% 
%     subplot(1,2,2);
%     imagesc(x, y, m3); colorbar;
%     title('Final m3 component'); axis equal; set(gca,'YDir','normal');
% % end

%% ================= Helper Functions =================

function [m1, m2, m3] = normalize_m(m1, m2, m3)
    n = sqrt(m1.^2 + m2.^2 + m3.^2);
    m1 = m1./n; m2 = m2./n; m3 = m3./n;
end

function [gx, gy] = calc_grad_neumann(f, hx, hy)
    % 二阶精度 Neumann 梯度计算 (边界法向导数为0)
    [Ny, Nx] = size(f);
    gx = zeros(Ny, Nx); gy = zeros(Ny, Nx);
    
    % X方向中心差分
    gx(:, 2:end-1) = (f(:, 3:end) - f(:, 1:end-2)) / (2*hx);
    gx(:, 1) = 0; gx(:, end) = 0; % Neumann BC
    
    % Y方向中心差分
    gy(2:end-1, :) = (f(3:end, :) - f(1:end-2, :)) / (2*hy);
    gy(1, :) = 0; gy(end, :) = 0; % Neumann BC
end

% function E = calc_energy(m1, m2, m3, hx, hy)
%     [dx1, dy1] = calc_grad_neumann(m1, hx, hy);
%     [dx2, dy2] = calc_grad_neumann(m2, hx, hy);
%     [dx3, dy3] = calc_grad_neumann(m3, hx, hy);
% 
%     dens = dx1.^2 + dy1.^2 + dx2.^2 + dy2.^2 + dx3.^2 + dy3.^2;
%     E = 0.5 * sum(dens(:)) * hx * hy;
% end

% 修改 calc_energy 函数，使其接受 L 矩阵
function E = calc_energy(m1, m2, m3, hx, hy, L) 
    % 注意：需要在调用处增加 L 参数
    Lm1 = reshape(L * m1(:), size(m1));
    Lm2 = reshape(L * m2(:), size(m2));
    Lm3 = reshape(L * m3(:), size(m3));
    
    total = sum(sum( m1.*Lm1 + m2.*Lm2 + m3.*Lm3 ));
    E = -0.5 * total * hx * hy;
end

% 相应地修改 energy_wrapper
function E = energy_wrapper(xi, m1h, m2h, m3h, hx, hy, L)
    t1 = m1h + xi; t2 = m2h + xi; t3 = m3h + xi;
    nt = sqrt(t1.^2 + t2.^2 + t3.^2);
    % 传入 L
    E = calc_energy(t1./nt, t2./nt, t3./nt, hx, hy, L);
end

% function E = energy_wrapper(xi, m1h, m2h, m3h, hx, hy)
%     t1 = m1h + xi; t2 = m2h + xi; t3 = m3h + xi;
%     nt = sqrt(t1.^2 + t2.^2 + t3.^2);
%     E = calc_energy(t1./nt, t2./nt, t3./nt, hx, hy);
% end

function [L, Dx, Dy] = get_neumann_laplacian_matrix(Nx, Ny, hx, hy)
    % 构建 Neumann 边界条件的离散拉普拉斯矩阵
    % 采用幽灵点法：(u_{i+1} - 2u_i + u_{i-1}) / h^2
    % 边界处 u_{-1} = u_{1} => (2u_{1} - 2u_{0}) / h^2
    
    Nx1 = Nx + 1; Ny1 = Ny + 1;
    
    % X 算子
    ex = ones(Nx1, 1);
    Lx_1d = spdiags([ex -2*ex ex], -1:1, Nx1, Nx1);
    Lx_1d(1, 2) = 2; Lx_1d(Nx1, Nx1-1) = 2;
    Lx_1d = Lx_1d / hx^2;
    
    % Y 算子
    ey = ones(Ny1, 1);
    Ly_1d = spdiags([ey -2*ey ey], -1:1, Ny1, Ny1);
    Ly_1d(1, 2) = 2; Ly_1d(Ny1, Ny1-1) = 2;
    Ly_1d = Ly_1d / hy^2;
    
    Ix = speye(Nx1); Iy = speye(Ny1);
    
    L = kron(Iy, Lx_1d) + kron(Ly_1d, Ix);
    
    Dx = []; Dy = []; % 用 calc_grad_neumann 代替矩阵形式，更灵活
end