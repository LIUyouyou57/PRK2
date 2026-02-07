% 基于论文 Length Preserving Numerical Schemes for Landau-Lifshitz Equation（arXiv:2206.02882） 
% 中对方程(5.4)的二阶length and energy preserving数值格式(5.8)-(5.11)的代码
% 计算区域为[-1/2,1/2]^2
% 边界条件为零纽曼边界条件

clear; clc; 

%% ================= Parameters =================
Nx = 24; Ny = 24;
Lx = 1; Ly = 1;
hx = Lx/Nx; hy = Ly/Ny;

dt = 1e-4;
T  = 0.1;
Nt = round(T/dt);

gamma = 1.0;
beta  = -1.0;
S     = 0;

%% ================= Grid =================
x = linspace(-0.5,0.5,Nx+1);
y = linspace(-0.5,0.5,Ny+1);
[X,Y] = meshgrid(x,y);

%% ================= Operators =================
L = NeumannLaplacian2D(Nx,Ny,hx,hy);
N = numel(X);
I = speye(N);
A = (1/dt)*I - (S+gamma)*L;

%% ================= Initial condition |m|=1 =================
% m1 = sin(2*pi*X).*cos(2*pi*Y);
% m2 = -cos(2*pi*X).*sin(2*pi*Y);
% m3 = sqrt(max(1-m1.^2-m2.^2,0));

LL = sqrt(X.^2 + Y.^2);
AA = (1 - 2*LL).^4;
BB = AA.^2 + LL.^2;
m1 = (LL<=0.5).*(2*X.*AA./BB);
m2 = (LL<=0.5).*(2*Y.*AA./BB);
m3 = (LL<=0.5).*((AA.^2 - LL.^2)./BB) - (LL>0.5).*(ones(size(LL)));

[m1,m2,m3] = normalize_m(m1,m2,m3);

m1_old = m1; m2_old = m2; m3_old = m3;

% 1. 计算正确的初始 lambda = |grad m|^2
% 注意：这里需要计算中心差分梯度
[dm1x, dm1y] = numerical_grad(m1, hx, hy); % 使用之前定义的子函数
[dm2x, dm2y] = numerical_grad(m2, hx, hy);
[dm3x, dm3y] = numerical_grad(m3, hx, hy);
lambda = (dm1x.^2 + dm1y.^2 + dm2x.^2 + dm2y.^2 + dm3x.^2 + dm3y.^2);

Energy = zeros(Nt+1,1);
xi_hist = zeros(Nt,1);

Energy(1) = discrete_energy(m1,m2,m3,hx,hy);

%% ================= Time loop =================
tic;
for n = 1:Nt
    fprintf('Step %d / %d\n',n,Nt);

    %% ===== Step 1: Gauss–Seidel predictor =====
    m1_dag = 1.5*m1 - 0.5*m1_old;
    m2_dag = 1.5*m2 - 0.5*m2_old;
    m3_dag = 1.5*m3 - 0.5*m3_old;

    Lap_m1_dag = reshape(L*m1_dag(:),size(m1));
    Lap_m2_dag = reshape(L*m2_dag(:),size(m2));
    Lap_m3_dag = reshape(L*m3_dag(:),size(m3));

    % ---- tilde m1 ----
    RHS1 = m1(:)/dt ...
         + gamma*lambda(:).*m1(:) ...
         - beta*( m2_dag(:).*Lap_m3_dag(:) ...
                - m3_dag(:).*Lap_m2_dag(:) ) ...
         - S*Lap_m1_dag(:);
    m1_tilde = reshape(A\RHS1,size(m1));

    % ---- tilde m2 ----
    m1_half = 0.5*(m1_tilde + m1);
    Lap_m1_half = reshape(L*m1_half(:),size(m1));

    RHS2 = m2(:)/dt ...
         + gamma*lambda(:).*m2(:) ...
         + beta*( m1_half(:).*Lap_m3_dag(:) ...
                - m3_dag(:).*Lap_m1_half(:) ) ...
         - S*Lap_m2_dag(:);
    m2_tilde = reshape(A\RHS2,size(m2));

    % ---- tilde m3 ----
    m2_half = 0.5*(m2_tilde + m2);
    Lap_m2_half = reshape(L*m2_half(:),size(m2));

    RHS3 = m3(:)/dt ...
         + gamma*lambda(:).*m3(:) ...
         - beta*( m1_half(:).*Lap_m2_half(:) ...
                - m2_half(:).*Lap_m1_half(:) ) ...
         - S*Lap_m3_dag(:);
    m3_tilde = reshape(A\RHS3,size(m3));

    %% ----- Step 2: correct lambda and m -----
    v1 = m1_tilde - 0.5*gamma*dt.*lambda.*m1;
    v2 = m2_tilde - 0.5*gamma*dt.*lambda.*m2;
    v3 = m3_tilde - 0.5*gamma*dt.*lambda.*m3;

    nv = sqrt(v1.^2+v2.^2+v3.^2);

    lambda = (1-nv)./(0.5*gamma*dt);

    mhat1 = v1./nv; mhat2 = v2./nv; mhat3 = v3./nv;

    % %% ===== Step 3: solve scalar xi^{n+1} =====
    % % dissipation term
    % Lap_1 = reshape(L*(mhat1(:)+m1(:))/2,size(mhat1));
    % Lap_2 = reshape(L*(mhat2(:)+m2(:))/2,size(mhat2));
    % Lap_3 = reshape(L*(mhat3(:)+m3(:))/2,size(mhat3));
    % 
    % G1 = 1/2*(mhat2+m2).*Lap_3 - 1/2*(mhat3+m3).*Lap_2;
    % G2 = 1/2*(mhat3+m3).*Lap_1 - 1/2*(mhat1+m1).*Lap_3;
    % G3 = 1/2*(mhat1+m1).*Lap_2 - 1/2*(mhat2+m2).*Lap_1;
    % 
    % Diss = sum(G1(:).^2 + G2(:).^2 + G3(:).^2) * hx * hy;
    % 
    % Eold = Energy(n);
     Energy_xi = @(xi) energy_with_xi(mhat1,mhat2,mhat3,xi,hx,hy);
    % F = @(xi) Energy_xi(xi) - Eold + dt*gamma*Diss;
    % 
    % % secant iteration
    % xi0 = -dt^2;
    % xi1 = 0;
    % F0 = F(xi0);
    % F1 = F(xi1);
    % 
    % for k = 1:20
    %     xi2 = xi1 - F1*(xi1-xi0)/(F1-F0);
    %     if abs(xi2-xi1) < 1e-12
    %         break;
    %     end
    %     xi0 = xi1; F0 = F1;
    %     xi1 = xi2; F1 = F(xi1);
    % end
    % xi = xi2;
    % xi_hist(n) = abs(xi);

    xi = 0;

    % final update
    m1_old = m1; m2_old = m2; m3_old = m3;

    m1 = mhat1 + xi;
    m2 = mhat2 + xi;
    m3 = mhat3 + xi;
    [m1,m2,m3] = normalize_m(m1,m2,m3);

    Energy(n+1) = Energy_xi(xi);

    fprintf('   E = %.6e, |xi| = %.3e\n',Energy(n+1),abs(xi));
end
toc;

% %% ================= Plots =================
figure;
plot((0:Nt)*dt,Energy,'LineWidth',1.8);
xlabel('t'); ylabel('E(t)');
title('Energy dissipation');
grid on;
% 
% figure;
% semilogy((1:Nt)*dt,xi_hist,'LineWidth',1.8);
% xlabel('t'); ylabel('|\xi^{n+1}|');
% title('Scalar correction');
% grid on;
% 
% NematicPlot(m1,m2,m3, X, Y, 1)
% 
% toc;

function [mx, my] = numerical_grad(m, hx, hy)
    [Ny, Nx] = size(m); % 注意 Matlab size 是 (行, 列) -> (y, x)
    
    mx = zeros(Ny, Nx);
    my = zeros(Ny, Nx);
    
    %% === 1. X 方向梯度 (mx) ===
    % 1.1 内部点：标准中心差分 (i+1 - i-1) / 2h
    mx(:, 2:end-1) = (m(:, 3:end) - m(:, 1:end-2)) / (2*hx);
    
    % 1.2 左右边界：零纽曼条件 => 法向导数为 0
    mx(:, 1)   = 0; 
    mx(:, end) = 0;
    
    %% === 2. Y 方向梯度 (my) ===
    % 2.1 内部点：标准中心差分 (j+1 - j-1) / 2h
    my(2:end-1, :) = (m(3:end, :) - m(1:end-2, :)) / (2*hy);
    
    % 2.2 上下边界：零纽曼条件 => 法向导数为 0
    my(1, :)   = 0;
    my(end, :) = 0;
    
    %% === 3. 角点处理 (可选，但推荐) ===
    % 对于正方形区域的 Neumann 边界，角点处的两个导数通常都被视为 0
    % 上述代码已经自动涵盖了角点（因为 mx 在角点属边界，my 在角点也属边界）
    % 例如 mx(1,1) 被 1.2 置为 0， my(1,1) 被 2.2 置为 0。
end
