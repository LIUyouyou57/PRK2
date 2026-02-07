% function LL_Scheme_Type1_Corrected()
    % 不带能量修正的LM2_CN格式
    
    clear; clc; close all;
    %% ================= 1. Parameters =================
    Nx = 48; Ny = 48; h=1/Nx;
    Lx = 1; Ly = 1;
    hx = Lx/Nx; hy = Ly/Ny;
    
    dt = 2e-4; 
    T  = 0.5; 
    Nt = round(T/dt);
    
    % 若要测试 HeatFlow，请设置: gamma=1, beta=0, S=0
    gamma = 1.0;
    beta  = -1.0; % 设为 0 以匹配 HeatFlow
    S     = 0; % 设为 0 以匹配 HeatFlow

    %% ================= Grid =================
    x = linspace(-Lx/2, Lx/2, Nx+1);
    y = linspace(-Ly/2, Ly/2, Ny+1);
    [X,Y] = meshgrid(x,y);
    
    %% ================= Operators =================
    [L, ~, ~] = NeumannLaplacian2D_Matrix(Nx, Ny, hx, hy);
    N = numel(X);
    I = speye(N);
    
    % Step 1 Linear Operator Matrix (LHS)
    % Equation: (1/dt - (S+gamma)/2 * Delta) * mt_new = ...
    % 当 S=0, gamma=1 时，变为 (1/dt - 0.5*Delta)，与 HeatFlow 一致
    A_step1 = (1/dt)*I - 0.5*(S+gamma)*L;
    
    % Operator for RHS: (S+gamma)/2 * Delta
    B_step1 = 0.5*(S+gamma)*L;

    %% ================= Initial Condition =================
    % % % 保持与 HeatFlow_Solver 完全一致的初值生成
    % LL_r = sqrt(X.^2 + Y.^2);
    % AA = (1 - 2*LL_r).^4;
    % BB = AA.^2 + LL_r.^2;
    % mask = (LL_r <= 0.5);
    % 
    % BB(BB==0) = 1; % 避免除零，与 HeatFlow 一致
    % 
    % m1 = mask.*(2*X.*AA./BB); 
    % m2 = mask.*(2*Y.*AA./BB);
    % m3 = mask.*((AA.^2 - LL_r.^2)./BB) - (~mask).*(ones(size(LL_r)));
   
    m1 = 0.3 * sin(pi*X) .* sin(pi*Y);
    m2 = 0.3 * sin(3*pi*X) .* sin(pi*Y);
    m3 = 1.0 + 0.2 * cos(2*pi*X) .* cos(2*pi*Y);
    
    [m1, m2, m3] = normalize_m(m1, m2, m3);
    
    % --- Initialization of History Variables ---
    % m_n: 当前物理时间步 (n) 的值
    m1_n = m1; m2_n = m2; m3_n = m3;
    
    % mt_n: 公示中的 \tilde{m}^n
    % 修正点：初始时，\tilde{m}^n 就是物理值 m^n
    mt1_n = m1; mt2_n = m2; mt3_n = m3;
    
    % mt_nm1: 公式中的 \tilde{m}^{n-1}，用于外推
    mt1_nm1 = m1; mt2_nm1 = m2; mt3_nm1 = m3;

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
    
    fprintf('Start LL_Type1: S=%.1f, beta=%.1f, gamma=%.1f\n', S, beta, gamma);

    %% ================= Time Loop =================
    tic;
    for n = 1:Nt
        
        %% === Step 1: Gauss-Seidel Predictor ===
        % 计算外推项 (Extrapolation) \tilde{m}^\dagger
        if n == 1
            mt1_dag = mt1_n; mt2_dag = mt2_n; mt3_dag = mt3_n;
        else
            mt1_dag = 1.5*mt1_n - 0.5*mt1_nm1;
            mt2_dag = 1.5*mt2_n - 0.5*mt2_nm1;
            mt3_dag = 1.5*mt3_n - 0.5*mt3_nm1;
        end
        
        % 计算 Dagger 的拉普拉斯项
        Lap_mt1_dag = reshape(L*mt1_dag(:), size(m1));
        Lap_mt2_dag = reshape(L*mt2_dag(:), size(m2));
        Lap_mt3_dag = reshape(L*mt3_dag(:), size(m3));
        
        % RHS 公共项: m^n/dt + (S+gamma)/2 * Delta * \tilde{m}^n
        % 修正点：这里的 \tilde{m}^n 现在等价于 m^n (因为我们在循环末尾做了重置)
        Diff_mt1_n = reshape(B_step1*mt1_n(:), size(m1));
        Diff_mt2_n = reshape(B_step1*mt2_n(:), size(m2));
        Diff_mt3_n = reshape(B_step1*mt3_n(:), size(m3));

        % --- Solve \tilde{m}_1^{n+1} ---
        Cross1 = mt2_dag.*Lap_mt3_dag - mt3_dag.*Lap_mt2_dag;
        
        RHS1 = m1_n(:)/dt + Diff_mt1_n(:) ...
             + gamma * lambda(:) .* m1_n(:) ...
             - beta * Cross1(:) ...
             - S * Lap_mt1_dag(:);
         
        mt1_np1 = reshape(A_step1 \ RHS1, size(m1));
        
        % GS Half step: \tilde{m}_1^{n+1/2}
        mt1_half = 0.5 * (mt1_np1 + mt1_n);
        Lap_mt1_half = reshape(L*mt1_half(:), size(m1));
        
        % --- Solve \tilde{m}_2^{n+1} ---
        Cross2 = mt1_half.*Lap_mt3_dag - mt3_dag.*Lap_mt1_half;
        
        RHS2 = m2_n(:)/dt + Diff_mt2_n(:) ...
             + gamma * lambda(:) .* m2_n(:) ...
             + beta * Cross2(:) ...
             - S * Lap_mt2_dag(:);
        
        mt2_np1 = reshape(A_step1 \ RHS2, size(m2));
        
        mt2_half = 0.5 * (mt2_np1 + mt2_n);
        Lap_mt2_half = reshape(L*mt2_half(:), size(m2));
        
        % --- Solve \tilde{m}_3^{n+1} ---
        Cross3 = mt1_half.*Lap_mt2_half - mt2_half.*Lap_mt1_half;
        
        RHS3 = m3_n(:)/dt + Diff_mt3_n(:) ...
             + gamma * lambda(:) .* m3_n(:) ...
             - beta * Cross3(:) ...
             - S * Lap_mt3_dag(:);
         
        mt3_np1 = reshape(A_step1 \ RHS3, size(m3));

        %% === Step 2: Corrector ===
        % v = \tilde{m}^{n+1} - dt*gamma/2 * lambda^n * m^n
        factor = 0.5 * dt * gamma;
        v1 = mt1_np1 - factor * lambda .* m1_n;
        v2 = mt2_np1 - factor * lambda .* m2_n;
        v3 = mt3_np1 - factor * lambda .* m3_n;
        
        nv = sqrt(v1.^2 + v2.^2 + v3.^2);
        
        lambda_next = (1 - nv) ./ factor;
        
        mhat1 = v1 ./ nv;
        mhat2 = v2 ./ nv;
        mhat3 = v3 ./ nv;
        
        % %% === Step 3: Energy Dissipation Correction ===
        % % m_bar = (mhat + m^n) / 2
        % m_bar1 = 0.5 * (mhat1 + m1_n);
        % m_bar2 = 0.5 * (mhat2 + m2_n);
        % m_bar3 = 0.5 * (mhat3 + m3_n);
        % 
        % L_bar1 = reshape(L*m_bar1(:), size(m1));
        % L_bar2 = reshape(L*m_bar2(:), size(m2));
        % L_bar3 = reshape(L*m_bar3(:), size(m3));
        % 
        % c1 = m_bar2.*L_bar3 - m_bar3.*L_bar2;
        % c2 = m_bar3.*L_bar1 - m_bar1.*L_bar3;
        % c3 = m_bar1.*L_bar2 - m_bar2.*L_bar1;
        % 
        % Diss = sum(sum(W.*(c1.^2 + c2.^2 + c3.^2))) * hx * hy;
        % E_target = Energy(n) - dt * gamma * Diss;
        % 
        % % Secant method for xi
        % func_xi = @(xi_val) energy_xi_wrapper(xi_val, mhat1, mhat2, mhat3, hx, hy, L, W) - E_target;
        % 
        % xi0 = 0;
        % xi1 =  -dt^2;% Guess
        % f0 = func_xi(xi0);
        % f1 = func_xi(xi1);
        % 
        % iter = 0;
        % while abs(f1)>1e-12
        %     iter = iter + 1;
        %     xi_new = xi1 - f1*(xi1-xi0)/(f1-f0);
        %     xi0 = xi1; f0 = f1;
        %     xi1 = xi_new; f1 = func_xi(xi1);
        %     xi_sol = xi1;
        % end

        xi_sol = 0;
        
        % Final m^{n+1}
        temp1 = mhat1 + xi_sol;
        temp2 = mhat2 + xi_sol;
        temp3 = mhat3 + xi_sol;
        nn = sqrt(temp1.^2 + temp2.^2 + temp3.^2);
        
        m1_next = temp1 ./ nn;
        m2_next = temp2 ./ nn;
        m3_next = temp3 ./ nn;
        
        %% === Update History (关键修正位置) ===
        
        % 1. 更新 n-1 时刻的值 
        % 修正后：使用 m_n (物理值) 而非 mt_n (预测值)
        mt1_nm1 = mt1_n; mt2_nm1 = mt2_n; mt3_nm1 = mt3_n; 
        
        % 2. 更新 n 时刻的值 
        % 修正后：强制将 \tilde{m}^n 重置为物理变量 m^{n+1}
        % 这保证了 Step 1 的 RHS 总是从物理修正后的状态出发
        mt1_n = m1_next; mt2_n = m2_next; mt3_n = m3_next; 
        
        % 3. 更新物理变量 m
        m1_n = m1_next; m2_n = m2_next; m3_n = m3_next;
        
        % 4. 更新 Lambda
        lambda = lambda_next;
        
        Energy(n+1) = calculate_energy(m1_n, m2_n, m3_n, hx, hy, L, W);
        
        if mod(n, 1) == 0
            fprintf('Step %d/%d, E=%.5e, |xi|=%.2e \n', n, Nt, Energy(n+1), abs(xi_sol));
        end
    end
    toc;

    % %% ================= Plotting =================
    figure;
    plot((0:Nt)*dt, Energy, '-r', 'LineWidth', 1.5);
    xlabel('Time'); ylabel('Energy');
    title('Energy Dissipation (LL Type-I)');
    grid on;
    % 
    % figure;
    % quiver(X(1:2:end,1:2:end), Y(1:2:end,1:2:end), ...
    %        m1_n(1:2:end,1:2:end), m2_n(1:2:end,1:2:end));
    % axis equal; xlim([-0.5 0.5]); ylim([-0.5 0.5]);
    % title(['Spin Field at T = ' num2str(T)]);
% end

%% ================= Helper Functions =================

function [m1, m2, m3] = normalize_m(m1, m2, m3)
    mod_m = sqrt(m1.^2 + m2.^2 + m3.^2);
    m1 = m1 ./ mod_m;
    m2 = m2 ./ mod_m;
    m3 = m3 ./ mod_m;
end

function [L, Dx, Dy] = NeumannLaplacian2D_Matrix(Nx, Ny, hx, hy)
    e = ones(Nx+1, 1);
    Lx_1d = spdiags([e -2*e e], -1:1, Nx+1, Nx+1);
    Lx_1d(1,2) = 2; Lx_1d(Nx+1, Nx) = 2; 
    Lx_1d = Lx_1d / hx^2;

    e = ones(Ny+1, 1);
    Ly_1d = spdiags([e -2*e e], -1:1, Ny+1, Ny+1);
    Ly_1d(1,2) = 2; Ly_1d(Ny+1, Ny) = 2; 
    Ly_1d = Ly_1d / hy^2;

    Ix = speye(Nx+1);
    Iy = speye(Ny+1);

    L = kron(Iy, Lx_1d) + kron(Ly_1d, Ix);
    Dx = []; Dy = [];
end

function E = calculate_energy(m1, m2, m3, hx, hy, L, W) 
    % 注意：需要在调用处增加 L 参数
    Lm1 = reshape(L * m1(:), size(m1));
    Lm2 = reshape(L * m2(:), size(m2));
    Lm3 = reshape(L * m3(:), size(m3));
    
    total = sum(sum( W.*(m1.*Lm1 + m2.*Lm2 + m3.*Lm3) ));
    E = -0.5 * total * hx * hy;
end

