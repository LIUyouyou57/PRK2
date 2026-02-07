%% test order of SIP1 scheme to solve LL eq
clear;  clc;

h = 1/64;
N = 1/h + 1;
x = -1/2: h: 1/2;  y = -1/2: h: 1/2;
[X, Y] = meshgrid(x, y);

beta = 1;  alpha = 1;

u1 = 0.3 * sin(pi*X) .* sin(pi*Y);
u2 = 0.3 * sin(3*pi*X) .* sin(pi*Y);
u3 = 1.0 + 0.2 * cos(2*pi*X) .* cos(2*pi*Y);

U_in = [u1(:); u2(:); u3(:)];
% 初值投影到球面
abs_U_in =sqrt( U_in(1:end/3).^2 + U_in(end/3+1:2*end/3).^2 + U_in(2*end/3+1:end).^2 );
min(abs_U_in)
U_in = U_in./[abs_U_in; abs_U_in; abs_U_in];
E_in = Energy(U_in);

% RK2 系数 d11=1
a11 = 1;  a21 = 0;  a22= 1/2;  b1=1/2;  b2=1/2;  d21=-1;  d22=2;

A = laplace_neumann_2D(N, h);

T = 0.01024;  t = 0;  dt = 0.1e-5;  n = fix(T/dt);  U = U_in;

load("LL_BDF4_N_64_T_0p01024.mat")

U_ex = [m1_ex(:); m2_ex(:); m3_ex(:)];
Step = 4*[8 16 32 64 128 256 512];
Error = zeros(1,7);

for j = 1:7
    t = 0; n = Step(j);dt = T/n;U = U_in;
    for i = 1:n
        t = t + dt;
        U0 = U;
        [PD] = PMatrix_D(U0, A, alpha, beta);
        M = speye(3*N^2) - dt*PD;
        U = M\U0;
        abs_U =sqrt( U(1:end/3).^2 + U(end/3+1:2*end/3).^2 + U(2*end/3+1:end).^2 );
        U = U./[abs_U; abs_U; abs_U];
    end
    Error(j) = sqrt(sum(h^2*(U-U_ex).^2))
end
order = log2(Error(1:end-1)./Error(2:end))

