% function LL_Scheme_Type1_Modified()
    % 基于论文 arXiv:2206.02882 修改版 Type-I 格式
    % 求解 Landau-Lifshitz 方程 (5.4)
    % 边界: Zero Neumann
    % 区域: [-1/2, 1/2]^2
    
    clear; clc; close all;

    for i=1:5

    %% ================= 1. Parameters =================
    Nx = 64; Ny = 64; h=1/Nx;
    Lx = 1; Ly = 1;
    hx = Lx/Nx; hy = Ly/Ny;
    
    dt = 128e-5/2^i; 
    T  = 0.01024; 
    Nt = round(T/dt);
    
    gamma = 1.0;
    beta  = -1; 
    S     = 0;  % Stabilization constant S > 0

    %% ================= Grid =================
    x = linspace(-0.5, 0.5, Nx+1);
    y = linspace(-0.5, 0.5, Ny+1);
    [X,Y] = meshgrid(x,y);
    
    %% ================= Operators =================
    % 生成 Neumann 拉普拉斯矩阵
    [L, ~, ~] = NeumannLaplacian2D_Matrix(Nx, Ny, hx, hy);
    N = numel(X);
    I = speye(N);
    
    % Step 1 Linear Operator Matrix (Left Hand Side)
    % Eq: (1/dt - (S+gamma)/2 * Delta) * m_tilde_new
    A_step1 = (1/dt)*I - 0.5*(S+gamma)*L;
    
    % Linear Operator part for Right Hand Side
    % Eq: ... + (S+gamma)/2 * Delta * m_tilde_old
    B_step1 = 0.5*(S+gamma)*L;

    %% ================= Initial Condition =================
    % LL_r = sqrt(X.^2 + Y.^2);
    % AA = (1 - 2*LL_r).^4;
    % BB = AA.^2 + LL_r.^2;
    % mask = (LL_r <= 0.5);
    % 
    % m1 = mask.*(2*X.*AA./(BB + eps)); 
    % m2 = mask.*(2*Y.*AA./(BB + eps));
    % m3 = mask.*((AA.^2 - LL_r.^2)./(BB + eps)) - (~mask).*(ones(size(LL_r)));
 
    m1 = 0.3 * sin(pi*X) .* sin(pi*Y);
    m2 = 0.3 * sin(3*pi*X) .* sin(pi*Y);
    m3 = 1.0 + 0.2 * cos(2*pi*X) .* cos(2*pi*Y);

    [m1, m2, m3] = normalize_m(m1, m2, m3);
    
    % 初始化变量
    % m: 物理时间步的值 (n)
    % mt: predictor的值 (\tilde{m})
    
    % Initial m^0
    m1_n = m1; m2_n = m2; m3_n = m3;
    
    % Initial \tilde{m}^0 (assumption: predictor matches initial state)
    mt1_n = m1; mt2_n = m2; mt3_n = m3;
    
    % History for extrapolation (\tilde{m}^{n-1})
    mt1_nm1 = mt1_n; mt2_nm1 = mt2_n; mt3_nm1 = mt3_n;

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

    Energy = zeros(Nt+1,1);
    Energy(1) = calculate_energy(m1_n, m2_n, m3_n, hx, hy, L, W);
    
    fprintf('Start Time Stepping... S=%.1f, dt=%.1e\n', S, dt);

    %% ================= Time Loop =================
    tic;
    for n = 1:Nt
        
        %% === Step 1: Gauss-Seidel Predictor ===
        % 计算外推项 \tilde{m}^\dagger
        % \tilde{m}^\dagger = 1.5 * \tilde{m}^n - 0.5 * \tilde{m}^{n-1}
        if n == 1
            mt1_dag = mt1_n; mt2_dag = mt2_n; mt3_dag = mt3_n;
        else
            mt1_dag = 1.5*mt1_n - 0.5*mt1_nm1;
            mt2_dag = 1.5*mt2_n - 0.5*mt2_nm1;
            mt3_dag = 1.5*mt3_n - 0.5*mt3_nm1;
        end
        
        Lap_mt1_dag = reshape(L*mt1_dag(:), size(m1));
        Lap_mt2_dag = reshape(L*mt2_dag(:), size(m2));
        Lap_mt3_dag = reshape(L*mt3_dag(:), size(m3));
        
        % 公共项: (S+gamma)/2 * Delta * \tilde{m}^n
        Diff_mt1_n = reshape(B_step1*mt1_n(:), size(m1));
        Diff_mt2_n = reshape(B_step1*mt2_n(:), size(m2));
        Diff_mt3_n = reshape(B_step1*mt3_n(:), size(m3));

        % --- Solve \tilde{m}_1^{n+1} ---
        % RHS = m1^n/dt + Diff_mt1_n + gamma*lambda^n*m1^n 
        %       - beta*(mt2_dag*Lap_mt3_dag - mt3_dag*Lap_mt2_dag) - S*Lap_mt1_dag
        
        Cross1 = mt2_dag.*Lap_mt3_dag - mt3_dag.*Lap_mt2_dag;
        RHS1 = m1_n(:)/dt + Diff_mt1_n(:) ...
             + gamma * lambda(:) .* m1_n(:) ...
             - beta * Cross1(:) ...
             - S * Lap_mt1_dag(:);
         
        mt1_np1 = reshape(A_step1 \ RHS1, size(m1));
        
        % Prepare half step for Gauss-Seidel: \tilde{m}_1^{n+1/2}
        mt1_half = 0.5 * (mt1_np1 + mt1_n);
        Lap_mt1_half = reshape(L*mt1_half(:), size(m1));
        
        % --- Solve \tilde{m}_2^{n+1} ---
        % RHS includes: beta*(mt1_half*Lap_mt3_dag - mt3_dag*Lap_mt1_half)
        Cross2 = mt1_half.*Lap_mt3_dag - mt3_dag.*Lap_mt1_half;
        RHS2 = m2_n(:)/dt + Diff_mt2_n(:) ...
             + gamma * lambda(:) .* m2_n(:) ...
             + beta * Cross2(:) ...
             - S * Lap_mt2_dag(:);
        
        mt2_np1 = reshape(A_step1 \ RHS2, size(m2));
        
        % Prepare half step: \tilde{m}_2^{n+1/2}
        mt2_half = 0.5 * (mt2_np1 + mt2_n);
        Lap_mt2_half = reshape(L*mt2_half(:), size(m2));
        
        % --- Solve \tilde{m}_3^{n+1} ---
        % RHS includes: -beta*(mt1_half*Lap_mt2_half - mt2_half*Lap_mt1_half)
        Cross3 = mt1_half.*Lap_mt2_half - mt2_half.*Lap_mt1_half;
        RHS3 = m3_n(:)/dt + Diff_mt3_n(:) ...
             + gamma * lambda(:) .* m3_n(:) ...
             - beta * Cross3(:) ...
             - S * Lap_mt3_dag(:);
         
        mt3_np1 = reshape(A_step1 \ RHS3, size(m3));

        %% === Step 2: Corrector ===
        % Solve (\hat{m}, \lambda^{n+1})
        % Formula: v = \tilde{m}^{n+1} - dt*gamma/2 * lambda^n * m^n
        %          \lambda^{n+1} = (1 - |v|) / (dt*gamma/2)
        %          \hat{m} = v / |v|
        % 注意：这里 \tilde{m}^{n+1} 是 Step 1 算出来的 mt_np1
        
        factor = 0.5 * dt * gamma;
        v1 = mt1_np1 - factor * lambda .* m1_n;
        v2 = mt2_np1 - factor * lambda .* m2_n;
        v3 = mt3_np1 - factor * lambda .* m3_n;
        
        nv = sqrt(v1.^2 + v2.^2 + v3.^2);
        
        % Update lambda for next step n+1
        lambda_next = (1 - nv) ./ factor;
        
        mhat1 = v1 ./ nv;
        mhat2 = v2 ./ nv;
        mhat3 = v3 ./ nv;
        
        %% === Step 3: Preserving Energy Dissipation ===
        % Solve scalar xi
        % m^{n+1} = (mhat + xi) / |mhat + xi|
        % E^{n+1} - E^n = -dt * || (mhat + m^n)/2 x Delta((mhat + m^n)/2) ||^2
        
        m_bar1 = 0.5 * (mhat1 + m1_n);
        m_bar2 = 0.5 * (mhat2 + m2_n);
        m_bar3 = 0.5 * (mhat3 + m3_n);
        
        L_bar1 = reshape(L*m_bar1(:), size(m1));
        L_bar2 = reshape(L*m_bar2(:), size(m2));
        L_bar3 = reshape(L*m_bar3(:), size(m3));
        
        c1 = m_bar2.*L_bar3 - m_bar3.*L_bar2;
        c2 = m_bar3.*L_bar1 - m_bar1.*L_bar3;
        c3 = m_bar1.*L_bar2 - m_bar2.*L_bar1;
        
        Diss = sum(sum(W.*(c1.^2 + c2.^2 + c3.^2))) * hx * hy;
        
        E_target = Energy(n) - dt * gamma * Diss;
        
        % Secant method for xi
        % Func(xi) = Energy(m(xi)) - E_target
        Func = @(xi_val) energy_xi_wrapper(xi_val, mhat1, mhat2, mhat3, hx, hy,L,W) - E_target;
        
        xi0 = 0; xi1 = -dt^2;
        F0 = Func(xi0); F1 = Func(xi1);
        
        iter = 0;
        while abs(F1)>1e-12
            iter = iter + 1;
            xi2 = xi1 - F1 * (xi1 - xi0) / (F1 - F0);
            xi0 = xi1; F0 = F1;
            xi1 = xi2; F1 = Func(xi1);
            xi_sol = xi1;
        end
        
        % Update final m^{n+1}
        temp1 = mhat1 + xi_sol;
        temp2 = mhat2 + xi_sol;
        temp3 = mhat3 + xi_sol;
        nn = sqrt(temp1.^2 + temp2.^2 + temp3.^2);
        
        m1_next = temp1 ./ nn;
        m2_next = temp2 ./ nn;
        m3_next = temp3 ./ nn;
        
        %% === Update History Variables ===
        % Update Predictor History
        mt1_nm1 = mt1_n; mt2_nm1 = mt2_n; mt3_nm1 = mt3_n; % \tilde{m}^{n-1} <- \tilde{m}^n
        mt1_n   = mt1_np1; mt2_n   = mt2_np1; mt3_n   = mt3_np1; % \tilde{m}^n   <- \tilde{m}^{n+1}
        
        % Update Physical Variable
        m1_n = m1_next; m2_n = m2_next; m3_n = m3_next; % m^n <- m^{n+1}
        
        % Update Lambda (calculated in Step 2)
        lambda = lambda_next;
        
        Energy(n+1) = calculate_energy(m1_n, m2_n, m3_n, hx, hy, L, W);
        
        if mod(n, 1) == 0
            fprintf('Step %d/%d, E=%.5e, |xi|=%.2e, iter number=%d\n', n, Nt, Energy(n+1), abs(xi_sol), iter);
        end
    end
    toc;

    load("data_1/LL_BDF4_T_0p01024.mat");
    error_L2 =  (m1_ex-m1_n).^2 + (m2_ex-m2_n).^2 + (m3_ex-m3_n).^2;
    Error_L2(i) =  h^2*sqrt( sum(error_L2(:)) );
    Error_inf(i) = ( norm(m1_n(:)-m1_ex(:),inf)+norm(m2_n(:)-m2_ex(:),inf)+norm(m3_n(:)-m3_ex(:),inf) )/3;

    end

    order_L2_1 = log(Error_L2(1:4)./Error_L2(2:5))/log(2)
    order_L2_2 = log((Error_L2(1:3)-Error_L2(2:4))./(Error_L2(2:4)-Error_L2(3:5)))/log(2)

    order_inf_1 = log(Error_inf(1:4)./Error_inf(2:5))/log(2)
    order_inf_2 = log((Error_inf(1:3)-Error_inf(2:4))./(Error_inf(2:4)-Error_inf(3:5)))/log(2)


    %% ================= 5. Plots =================
    figure;
    plot((0:Nt)*dt, Energy, 'LineWidth', 1.5);
    xlabel('t'); ylabel('Energy'); title('Energy Decay (Heat Flow)');
    grid on;
    
    figure;
    subplot(1,2,1);
    quiver(X(1:2:end,1:2:end), Y(1:2:end,1:2:end), m1_n(1:2:end,1:2:end), m2_n(1:2:end,1:2:end));
    title('Final Spin Field (m1, m2)'); axis equal; xlim([-0.5 0.5]); ylim([-0.5 0.5]);
    
    subplot(1,2,2);
    imagesc(x, y, m3); colorbar;
    title('Final m3 component'); axis equal; set(gca,'YDir','normal');
% end

%% ================= Helper Functions =================

function [m1, m2, m3] = normalize_m(m1, m2, m3)
    mod_m = sqrt(m1.^2 + m2.^2 + m3.^2);
    m1 = m1 ./ mod_m;
    m2 = m2 ./ mod_m;
    m3 = m3 ./ mod_m;
end

function [L, Dx, Dy] = NeumannLaplacian2D_Matrix(Nx, Ny, hx, hy)
    % 1D Matrix construction for x
    e = ones(Nx+1, 1);
    Lx_1d = spdiags([e -2*e e], -1:1, Nx+1, Nx+1);
    Lx_1d(1,2) = 2; Lx_1d(Nx+1, Nx) = 2; % Neumann Fix: u_{-1}=u_1
    Lx_1d = Lx_1d / hx^2;

    % 1D Matrix construction for y
    e = ones(Ny+1, 1);
    Ly_1d = spdiags([e -2*e e], -1:1, Ny+1, Ny+1);
    Ly_1d(1,2) = 2; Ly_1d(Ny+1, Ny) = 2; % Neumann Fix
    Ly_1d = Ly_1d / hy^2;

    Ix = speye(Nx+1);
    Iy = speye(Ny+1);

    % 2D Laplacian = Lx kron Iy + Ix kron Ly
    L = kron(Iy, Lx_1d) + kron(Ly_1d, Ix);
    
    Dx = []; Dy = []; % Not needed for matrix mult
end

function [gx, gy] = numerical_grad_neumann(m, hx, hy)
    % Gradient calculation respecting Zero Neumann BC (Gradient = 0 at boundary normal)
    [Ny, Nx] = size(m);
    gx = zeros(Ny, Nx);
    gy = zeros(Ny, Nx);
    
    % X-direction
    gx(:, 2:end-1) = (m(:, 3:end) - m(:, 1:end-2)) / (2*hx);
    gx(:, 1) = 0; gx(:, end) = 0;
    
    % Y-direction
    gy(2:end-1, :) = (m(3:end, :) - m(1:end-2, :)) / (2*hy);
    gy(1, :) = 0; gy(end, :) = 0;
end

% 修改 calc_energy 函数，使其接受 L 矩阵
function E = calculate_energy(m1, m2, m3, hx, hy, L, W) 
    % 注意：需要在调用处增加 L 参数
    Lm1 = reshape(L * m1(:), size(m1));
    Lm2 = reshape(L * m2(:), size(m2));
    Lm3 = reshape(L * m3(:), size(m3));
    
    total = sum(sum( W.*(m1.*Lm1 + m2.*Lm2 + m3.*Lm3) ));
    E = -0.5 * total * hx * hy;
end

% 相应地修改 energy_wrapper
function E = energy_xi_wrapper(xi, m1h, m2h, m3h, hx, hy, L, W)
    t1 = m1h + xi; t2 = m2h + xi; t3 = m3h + xi;
    nt = sqrt(t1.^2 + t2.^2 + t3.^2);
    % 传入 L
    E = calculate_energy(t1./nt, t2./nt, t3./nt, hx, hy, L, W);
end
