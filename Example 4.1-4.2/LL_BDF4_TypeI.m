clear; clc;

%% ================= parameters =================
Nx = 48; Ny = 48;
hx = 1/Nx; hy = 1/Ny;
dt = 1e-6; T = 0.008;
Nt = round(T/dt);

beta  = -1;
gamma = 1.0;

alpha4 = 25/12;

%% ================= grid =================
x = linspace(-0.5,0.5,Nx+1);
y = linspace(-0.5,0.5,Ny+1);
[X,Y] = meshgrid(x,y);

sz = size(X);
N = numel(X);

%% ================= Laplacian (Neumann) =================
L = NeumannLaplacian2D(Nx,Ny,hx,hy);
A = (alpha4/dt)*speye(N) - gamma*L;

%% ================= initial data =================
    % m1 = 0.3 * sin(pi*X) .* sin(pi*Y);
    % m2 = 0.3 * sin(3*pi*X) .* sin(pi*Y);
    % m3 = 1.0 + 0.2 * cos(2*pi*X) .* cos(2*pi*Y);

% % % 之前的涡旋初值
LL = sqrt(X.^2 + Y.^2);
AA = (1 - 2*LL).^4;
BB = AA.^2 + LL.^2;
mask = (LL <= 0.5);

% 避免除零
BB(BB==0) = 1; 

m1 = mask.*(2*X.*AA./BB);
m2 = mask.*(2*Y.*AA./BB);
m3 = mask.*((AA.^2 - LL.^2)./BB) - (~mask).*(ones(size(LL)));

[m1,m2,m3] = normalize_m(m1,m2,m3);

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

% history (BDF4 needs 4 levels)
m1_hist = zeros([sz,4]);
m2_hist = zeros([sz,4]);
m3_hist = zeros([sz,4]);
lam_hist = zeros([sz,3]);

m1_hist(:,:,1) = m1;
m2_hist(:,:,1) = m2;
m3_hist(:,:,1) = m3;
lam_hist(:,:,1) = lambda;

Energy = zeros(Nt+1,1);
Energy(1) = discrete_energy(m1,m2,m3,hx,hy,L,W);

for k = 1:3
    [m1_hist(:,:,k+1), ...
     m2_hist(:,:,k+1), ...
     m3_hist(:,:,k+1), ...
     lam_hist(:,:,k+1)] = ...
        LL_second_order_step( ...
            m1_hist(:,:,k), ...
            m2_hist(:,:,k), ...
            m3_hist(:,:,k), ...
            lam_hist(:,:,k), ...
            dt, L, beta, gamma);
    Energy(k+1) = discrete_energy(m1_hist(:,:,k+1),m2_hist(:,:,k+1),m3_hist(:,:,k+1),hx,hy,L,W);

end


%% ================= BDF4 time stepping =================
for n = 4:Nt
    fprintf('BDF4 step %d / %d\n',n,Nt);

    %% ----- A4(m^n) -----
    A4m1 = 4*m1_hist(:,:,4) - 3*m1_hist(:,:,3) ...
         + 4/3*m1_hist(:,:,2) - 1/4*m1_hist(:,:,1);
    A4m2 = 4*m2_hist(:,:,4) - 3*m2_hist(:,:,3) ...
         + 4/3*m2_hist(:,:,2) - 1/4*m2_hist(:,:,1);
    A4m3 = 4*m3_hist(:,:,4) - 3*m3_hist(:,:,3) ...
         + 4/3*m3_hist(:,:,2) - 1/4*m3_hist(:,:,1);

    %% ----- B3(lambda m) -----
    B3lm1 = 3*lam_hist(:,:,3).*m1_hist(:,:,4) ...
          - 3*lam_hist(:,:,2).*m1_hist(:,:,3) ...
          +   lam_hist(:,:,1).*m1_hist(:,:,2);
    B3lm2 = 3*lam_hist(:,:,3).*m2_hist(:,:,4) ...
          - 3*lam_hist(:,:,2).*m2_hist(:,:,3) ...
          +   lam_hist(:,:,1).*m2_hist(:,:,2);
    B3lm3 = 3*lam_hist(:,:,3).*m3_hist(:,:,4) ...
          - 3*lam_hist(:,:,2).*m3_hist(:,:,3) ...
          +   lam_hist(:,:,1).*m3_hist(:,:,2);

    %% ----- Laplacian of m^n, m^{n-1}, m^{n-2} -----
    tmp = m1_hist(:,:,4); Lap1n = reshape(L*tmp(:),sz);
    tmp = m2_hist(:,:,4); Lap2n = reshape(L*tmp(:),sz);
    tmp = m3_hist(:,:,4); Lap3n = reshape(L*tmp(:),sz);

    tmp = m1_hist(:,:,3); Lap1n1 = reshape(L*tmp(:),sz);
    tmp = m2_hist(:,:,3); Lap2n1 = reshape(L*tmp(:),sz);
    tmp = m3_hist(:,:,3); Lap3n1 = reshape(L*tmp(:),sz);

    tmp = m1_hist(:,:,2); Lap1n2 = reshape(L*tmp(:),sz);
    tmp = m2_hist(:,:,2); Lap2n2 = reshape(L*tmp(:),sz);
    tmp = m3_hist(:,:,2); Lap3n2 = reshape(L*tmp(:),sz);

    %% ----- B3(m × Δm) -----
    C1 = 3*(m2_hist(:,:,4).*Lap3n - m3_hist(:,:,4).*Lap2n) ...
       - 3*(m2_hist(:,:,3).*Lap3n1 - m3_hist(:,:,3).*Lap2n1) ...
       +   (m2_hist(:,:,2).*Lap3n2 - m3_hist(:,:,2).*Lap2n2);

    C2 = 3*(m3_hist(:,:,4).*Lap1n - m1_hist(:,:,4).*Lap3n) ...
       - 3*(m3_hist(:,:,3).*Lap1n1 - m1_hist(:,:,3).*Lap3n1) ...
       +   (m3_hist(:,:,2).*Lap1n2 - m1_hist(:,:,2).*Lap3n2);

    C3 = 3*(m1_hist(:,:,4).*Lap2n - m2_hist(:,:,4).*Lap1n) ...
       - 3*(m1_hist(:,:,3).*Lap2n1 - m2_hist(:,:,3).*Lap1n1) ...
       +   (m1_hist(:,:,2).*Lap2n2 - m2_hist(:,:,2).*Lap1n2);

    %% ----- Step 1: predictor -----
    RHS1 = A4m1/dt + gamma*B3lm1 - beta*C1;
    RHS2 = A4m2/dt + gamma*B3lm2 - beta*C2;
    RHS3 = A4m3/dt + gamma*B3lm3 - beta*C3;

    m1t = reshape(A\RHS1(:),sz);
    m2t = reshape(A\RHS2(:),sz);
    m3t = reshape(A\RHS3(:),sz);

    %% ----- Step 2: corrector (exact) -----
    v1 = alpha4*m1t - gamma*dt*B3lm1;
    v2 = alpha4*m2t - gamma*dt*B3lm2;
    v3 = alpha4*m3t - gamma*dt*B3lm3;

    nv = sqrt(v1.^2 + v2.^2 + v3.^2);

    lambda_new = (alpha4 - nv)/(gamma*dt);

    m1_new = v1./nv;
    m2_new = v2./nv;
    m3_new = v3./nv;

    %% ----- update history -----
    m1_hist = cat(3,m1_hist(:,:,2:4),m1_new);
    m2_hist = cat(3,m2_hist(:,:,2:4),m2_new);
    m3_hist = cat(3,m3_hist(:,:,2:4),m3_new);
    lam_hist = cat(3,lam_hist(:,:,2:3),lambda_new);

    Energy(n+1) = discrete_energy(m1_new,m2_new,m3_new,hx,hy,L,W);
    fprintf('   Energy = %.6e\n',Energy(n+1));
end

m1_ex = m1_new; m2_ex = m2_new; m3_ex = m3_new;

%% ================= 5. 绘图 =================
figure;
subplot(1,2,1);
plot(0:10*dt:T, Energy(1:10:end), 'LineWidth', 1.5);
xlabel('Time'); ylabel('Energy'); title('Dirichlet Energy Decay');
grid on;

subplot(1,2,2);
step = 2;
quiver(X(1:step:end, 1:step:end), Y(1:step:end, 1:step:end), ...
       m1_new(1:step:end, 1:step:end), m2_new(1:step:end, 1:step:end));
axis equal; xlim([-0.5 0.5]); ylim([-0.5 0.5]);
title(['Solution at T=' num2str(T)]);
drawnow;