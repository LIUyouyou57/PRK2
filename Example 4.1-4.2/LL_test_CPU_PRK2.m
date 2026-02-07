%% test order of PRK2 scheme to solve LL eq
clear;  clc;

data = struct();

h = 1/48;
N = 1/h + 1;
x = -1/2: h: 1/2;  y = -1/2: h: 1/2;
[X, Y] = meshgrid(x, y);

beta = 1;  alpha = 1;

% u1 = 0.3 * sin(pi*X) .* sin(pi*Y);
% u2 = 0.3 * sin(3*pi*X) .* sin(pi*Y);
% u3 = 1.0 + 0.2 * cos(2*pi*X) .* cos(2*pi*Y);

LL_r = sqrt(X.^2 + Y.^2);
AA = (1 - 2*LL_r).^4;
BB = AA.^2 + LL_r.^2;
mask = (LL_r <= 0.5);

BB(BB==0) = 1; % 避免除零，与 HeatFlow 一致

u1 = mask.*(2*X.*AA./BB); 
u2 = mask.*(2*Y.*AA./BB);
u3 = mask.*((AA.^2 - LL_r.^2)./BB) - (~mask).*(ones(size(LL_r)));

U_in = [u1(:); u2(:); u3(:)];
% 初值投影到球面
abs_U_in =sqrt( U_in(1:end/3).^2 + U_in(end/3+1:2*end/3).^2 + U_in(2*end/3+1:end).^2 );
min(abs_U_in)
U_in = U_in./[abs_U_in; abs_U_in; abs_U_in];
E_in = Energy(U_in);

% RK2 系数 d11=1
a11 = 1;  a21 = 0;  a22= 1/2;  b1=1/2;  b2=1/2;  d21=-1;  d22=2;

A = laplace_neumann_2D(N, h);

T = 0.02;  t = 0;  dt = 0.1e-5;  n = fix(T/dt);  U = U_in;

load("LL_BDF4_N_48_T_0p02.mat")
% load("LL_PRK2_N_128_T_0p01024.mat")


U_ex = [m1_ex(:); m2_ex(:); m3_ex(:)];
% Step = 25*[8 16 32 64 128 256 512];
Error = zeros(1,7);

for j = 1:7
    t = 0; dt = 2e-3/2^j; n = round(T/dt); U = U_in; tic;
    for i = 1:n
        t = t + dt;
        U0 = U;
        [PD0] = PMatrix_D(U0, A, alpha, beta);
        M = speye(3*N^2) - dt*a11*PD0;
        U1 = M\U0;
        [PD1] = PMatrix_D(U1, A, alpha, beta);
        M = speye(3*N^2) - dt*a22*d22*PD1;
        U2 = M\(U0 + dt*a22*d21*PD1*U1);
        U = U0+dt*( b1*(PD0*U1) + b2*(d21*PD1*U1 + d22*PD1*U2) );
        abs_U =sqrt( U(1:end/3).^2 + U(end/3+1:2*end/3).^2 + U(2*end/3+1:end).^2 );
        U = U./[abs_U; abs_U; abs_U];
        data(j).U = U;
    end
    wall_times(j) = toc;
    Error(j) = sqrt(sum(h^2*(U-U_ex).^2))
end
    order = log2(Error(1:end-1)./Error(2:end))

    figure('Color', 'w', 'Name', 'Work-Precision Diagram');
    loglog(wall_times, Error, '-o', 'LineWidth', 1.5, 'DisplayName', 'PRK2');


