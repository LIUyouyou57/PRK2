% Example 4.2: plot energy of SIP1; PRK2; scheme (4.4);
clear;  clc; 

h = 1/48;
N = 1/h + 1;
x = -1/2: h: 1/2;  y = -1/2: h: 1/2;
[X, Y] = meshgrid(x, y);

beta = 1;  alpha = 1;

LL_r = sqrt(X.^2 + Y.^2);
AA = (1 - 2*LL_r).^4;
BB = AA.^2 + LL_r.^2;
mask = (LL_r <= 0.5);

BB(BB==0) = 1; 

u1 = mask.*(2*X.*AA./BB); 
u2 = mask.*(2*Y.*AA./BB);
u3 = mask.*((AA.^2 - LL_r.^2)./BB) - (~mask).*(ones(size(LL_r)));

U_in = [u1(:); u2(:); u3(:)];
% 初值投影到球面
abs_U_in =sqrt( U_in(1:end/3).^2 + U_in(end/3+1:2*end/3).^2 + U_in(2*end/3+1:end).^2 );
min(abs_U_in)
U_in = U_in./[abs_U_in; abs_U_in; abs_U_in];
E_in = Energy(U_in,h);

% RK2 系数 d11=1
a11 = 1;  a21 = 0;  a22= 1/2;  b1=1/2;  b2=1/2;  d21=-1;  d22=2;
% PRK2_2p5 系数
a11_hat = 1;  a21_hat = -1/2;  a22_hat= 1;  b1_hat=0;  b2_hat=1;  g11=1; g21=1/2;  g22=1/2;

A = laplace_neumann_2D(N, h);
% vec = [490, 1000];
% vec = [100, 200, 400, 500, 800];
%% Energy decay
T = 1;  t = 0;  dt = 1e-4;  n = fix(T/dt);  U = U_in;  %NematicPlot(U, X, Y, 0); close all;
% 一阶方法能量
E1 = zeros(1, n); tic;
for i = 1:n
    t = t + dt
    U0 = U;
    [PD] = PMatrix_D(U0, A, alpha, beta);
    M = speye(3*N^2) - dt*PD;
    U = M\U0;
    abs_U =sqrt( U(1:end/3).^2 + U(end/3+1:2*end/3).^2 + U(2*end/3+1:end).^2 );
    U = U./[abs_U; abs_U; abs_U];
    length_U1(i) = max(abs_U(:));
    E1(i) = Energy(U,h);  
end
toc;

% RK2能量
t = 0;  U = U_in;
E2 = zeros(1, n);
tic;
for i = 1:n
    t = t + dt
    U0 = U;
    % PD0 = P(u0)*Laplace
    [PD0] = PMatrix_D(U0, A, alpha, beta);
    % M = I - dt*a11*P(u0)*Laplace
    M = speye(3*N^2) - dt*a11*PD0;
    U1 = M\U0;
    % PD1 = P(u1)*Laplace
    [PD1] = PMatrix_D(U1, A, alpha, beta);
    % M = I - dt*22*d22*P(u0)*Laplace
    M = speye(3*N^2) - dt*a22*d22*PD1;
    U2 = M\(U0 + dt*a22*d21*PD1*U1);
    U = U0+dt*( b1*(PD0*U1) + b2*(d21*PD1*U1 + d22*PD1*U2) );
    abs_U =sqrt( U(1:end/3).^2 + U(end/3+1:2*end/3).^2 + U(2*end/3+1:end).^2 );
    length_U2(i) = max(abs_U(:));
    U = U./[abs_U; abs_U; abs_U];
    E2(i) = Energy(U,h);
end

%other RK2-scheme(4.2)能量
t = 0;  U = U_in;
E3 = zeros(1, n);
tic;
for i = 1:n
        t = t + dt;
        U0 = U;
        [PD0] = PMatrix_D(U0, A, alpha, beta);
        M = speye(3*N^2) - dt*a11_hat*g11*PD0;
        U1 = M\U0;
        [PD1] = PMatrix_D(U1, A, alpha, beta);
        M = speye(3*N^2) - dt*a22_hat*(g21*PD0+g22*PD1);
        U2 = M\(U0 + dt*a21_hat*g11*PD0*U1);
        U = U0+dt*( b1_hat*g11*(PD0*U1) + b2_hat*(g21*PD0 + g22*PD1)*U2 );
        abs_U =sqrt( U(1:end/3).^2 + U(end/3+1:2*end/3).^2 + U(2*end/3+1:end).^2 );
        length_U3(i) = max(abs_U(:));
        U = U./[abs_U; abs_U; abs_U];
        E3(i) = Energy(U,h);
end

toc;

% Energy
plot(0:dt:T, [E_in, E1]);
hold on
plot(0:dt:T, [E_in, E2]);
plot(0:dt:T, [E_in, E3]);
%length of [u1, u2, u3]
figure;
plot(0:dt:T, [1, length_U1]);
hold on
plot(0:dt:T, [1, length_U2]);
plot(0:dt:T, [1, length_U3]);
