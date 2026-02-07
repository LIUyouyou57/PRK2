clear; clc;

%% ========== parameters ==========
Nx = 64; Ny = 64;
hx = 1/Nx; hy = 1/Ny;
dt = 2e-4; T = 0.01;
Nt = round(T/dt);

beta  = 0;
gamma = 1.0;
alpha2 = 3/2;

%% ========== grid ==========
x = linspace(-0.5,0.5,Nx+1);
y = linspace(-0.5,0.5,Ny+1);
[X,Y] = meshgrid(x,y);
sz = size(X);
N  = numel(X);

%% ========== Neumann Laplacian ==========
L = NeumannLaplacian2D(Nx,Ny,hx,hy);
A = (alpha2/dt)*speye(N) - gamma*L;

%% ========== initial data ==========
LL = sqrt(X.^2 + Y.^2);
AA = (1 - 2*LL).^4;
BB = AA.^2 + LL.^2;
m1 = (LL<=0.5).*(2*X.*AA./BB);
m2 = (LL<=0.5).*(2*Y.*AA./BB);
m3 = (LL<=0.5).*((AA.^2 - LL.^2)./BB) - (LL>0.5).*(ones(size(LL)));

[m1,m2,m3] = normalize_m(m1,m2,m3);

% 1. 计算正确的初始 lambda = |grad m|^2
% 注意：这里需要计算中心差分梯度
[dm1x, dm1y] = numerical_grad(m1, hx, hy); % 使用之前定义的子函数
[dm2x, dm2y] = numerical_grad(m2, hx, hy);
[dm3x, dm3y] = numerical_grad(m3, hx, hy);
lambda = (dm1x.^2 + dm1y.^2 + dm2x.^2 + dm2y.^2 + dm3x.^2 + dm3y.^2);

% one-step history
m1_old = m1; m2_old = m2; m3_old = m3;

Energy = zeros(Nt+1,1); xi_hist = zeros(Nt,1);
Energy(1) = discrete_energy(m1,m2,m3,hx,hy);

%% ========== time loop ==========
for n = 1:Nt
    fprintf('BDF2 step %d / %d\n',n,Nt);

    %% ---- Laplacian of m^n ----
    tmp = m1; Lap1 = reshape(L*tmp(:),sz);
    tmp = m2; Lap2 = reshape(L*tmp(:),sz);
    tmp = m3; Lap3 = reshape(L*tmp(:),sz);

    %% ---- m x Delta m ----
    C1 = m2.*Lap3 - m3.*Lap2;
    C2 = m3.*Lap1 - m1.*Lap3;
    C3 = m1.*Lap2 - m2.*Lap1;

    %% ---- A2(m^n) ----
    A2m1 = 2*m1 - 0.5*m1_old;
    A2m2 = 2*m2 - 0.5*m2_old;
    A2m3 = 2*m3 - 0.5*m3_old;

    %% ================= Step 1: Predictor =================
    RHS1 = A2m1/dt + gamma*lambda.*m1 - beta*C1;
    RHS2 = A2m2/dt + gamma*lambda.*m2 - beta*C2;
    RHS3 = A2m3/dt + gamma*lambda.*m3 - beta*C3;

    m1t = reshape(A\RHS1(:),sz);
    m2t = reshape(A\RHS2(:),sz);
    m3t = reshape(A\RHS3(:),sz);

    %% ================= Step 2: Corrector =================
    v1 = alpha2*m1t - gamma*dt*lambda.*m1;
    v2 = alpha2*m2t - gamma*dt*lambda.*m2;
    v3 = alpha2*m3t - gamma*dt*lambda.*m3;

    nv = sqrt(v1.^2 + v2.^2 + v3.^2);

    lambda = (alpha2 - nv)/(gamma*dt);

    mh1 = v1./nv;
    mh2 = v2./nv;
    mh3 = v3./nv;

    %% ===== Step 3: solve scalar xi^{n+1} =====
    % dissipation term
    Lap_1 = reshape(L*(mh1(:)+m1(:))/2,size(mh1));
    Lap_2 = reshape(L*(mh2(:)+m2(:))/2,size(mh2));
    Lap_3 = reshape(L*(mh3(:)+m3(:))/2,size(mh3));

    G1 = 1/2*(mh2+m2).*Lap_3 - 1/2*(mh3+m3).*Lap_2;
    G2 = 1/2*(mh3+m3).*Lap_1 - 1/2*(mh1+m1).*Lap_3;
    G3 = 1/2*(mh1+m1).*Lap_2 - 1/2*(mh2+m2).*Lap_1;

    Diss = sum(G1(:).^2 + G2(:).^2 + G3(:).^2) * hx * hy;

    Eold = Energy(n);
    Energy_xi = @(xi) energy_with_xi(mh1,mh2,mh3,xi,hx,hy);
    F = @(xi) Energy_xi(xi) - Eold + dt*gamma*Diss;

    % secant iteration
    xi0 = -dt^2;
    xi1 = 0;
    F0 = F(xi0);
    F1 = F(xi1);

    for k = 1:20
        xi2 = xi1 - F1*(xi1-xi0)/(F1-F0);
        if abs(xi2-xi1) < 1e-12
            break;
        end
        xi0 = xi1; F0 = F1;
        xi1 = xi2; F1 = F(xi1);
    end
    xi = xi2;
    xi_hist(n) = abs(xi);

    %% ---- update ----
    m1_old = m1; m2_old = m2; m3_old = m3;
    m1 = mh1 + xi;
    m2 = mh2 + xi;
    m3 = mh3 + xi;
    [m1,m2,m3] = normalize_m(m1,m2,m3);

    Energy(n+1) = Energy_xi(xi);
    fprintf('   E = %.6e, |xi| = %.3e\n',Energy(n+1),abs(xi));
end

%% ========== plot ==========
figure;
plot((0:Nt)*dt,Energy,'LineWidth',1.8);
xlabel('t'); ylabel('Energy'); grid on;

figure;
semilogy((1:Nt)*dt,xi_hist,'LineWidth',1.8);
xlabel('t'); ylabel('|\xi^{n+1}|');
title('Scalar correction');
grid on;



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