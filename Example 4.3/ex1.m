clear;  clc;

h = 1/24;
N = 1/h - 1;
x = h: h: 1-h;  y = x;  z = x;
[X, Y, Z] = meshgrid(x, y, z);

u1 = sin(2*pi*X) + 2 + 0.5*sin(6*pi*Y) + 0.2*sin(4*pi*Z);
u2 = cos(2*pi*X) - 2 + 0.5*cos(6*pi*Y) + 0.2*cos(4*pi*Z);
u3 = sin(2*pi*X) + cos(6*pi*Y) + cos(4*pi*Z);
% u1 = rand(N, N, N);
% u2 = rand(N, N, N);
% u3 = rand(N, N, N);

U_in = [u1(:); u2(:); u3(:)];
% 初值投影到球面
abs_U_in =sqrt( U_in(1:N^3).^2 + U_in(N^3+1:2*N^3).^2 + U_in(2*N^3+1:end).^2 );
min(abs_U_in)
U_in = U_in./[abs_U_in; abs_U_in; abs_U_in];

% RK2 系数 d11=1
a11 = 1;  a21 = 0;  a22= 1/2;  b1=1/2;  b2=1/2;  d21=-1;  d22=2;

A = laplace_dirichlet_3D(N, h);
g = Dirchlet_bc(N);
E_in = Energy(U_in, g);  % 初始能量

vec = [15, 50 600];
% vec = [1:1:30, 40:10:100 100:100:1000];
%% Energy decay
T = 1;  t = 0;  dt = 1e-3;  n = fix(T/dt);  U = U_in;  NematicPlot(U, g, X, 0); close all;
% % 一阶方法能量
% E1 = zeros(1, n);  
% for i = 1:n
%     t = t + dt
%     U0 = U;
%     [PD, b_c] = PMatrix_D(U0, A, g, h);
%     M = speye(3*N^3) - dt*PD;
%     U = M\(U0 + dt*b_c);
%     abs_U =sqrt( U(1:end/3).^2 + U(end/3+1:2*end/3).^2 + U(2*end/3+1:end).^2 );
%     U = U./[abs_U; abs_U; abs_U];
%     E1(i) = Energy(U, g);
% end
% NematicPlot(U, X, Y, Z)

% RK2能量
t = 0;  U = U_in;
E2 = zeros(1, n);  
for i = 1:n
    t = t + dt
    U0 = U;
    [PD0, b_c0] = PMatrix_D(U0, A, g, h);
    % M = I - dt*a11*P(u0)*Laplace
    M = speye(3*N^3) - dt*a11*PD0;
    U1 = M\(U0 + dt*a11*b_c0);
    % PD1 = P(u1)*Laplace
    [PD1, b_c1] = PMatrix_D(U1, A, g, h);
    % M = I - dt*22*d22*P(u0)*Laplace
    M = speye(3*N^3) - dt*a22*d22*PD1;
    U2 = M\(U0 + dt*a22*d21*(PD1*U1+b_c1) + dt*a22*d22*b_c1 );
    U = U0+dt*( b1*(PD0*U1+b_c0) + b2*(d21*(PD1*U1+b_c1) + d22*(PD1*U2+b_c1)) );
    abs_U =sqrt( U(1:end/3).^2 + U(end/3+1:2*end/3).^2 + U(2*end/3+1:end).^2 );
    U = U./[abs_U; abs_U; abs_U];
    E2(i) = Energy(U, g);
    if any(vec == i)
        NematicPlot(U, g, X, i)
        close all;
    else
    end    
end
NematicPlot(U, g, X, i);

plot(0:dt:T, [E_in, E1]);
hold on
plot(0:dt:T, [E_in, E2]);

% %% Second Order Test 
% T = 0.05;  t = 0;  dt = 0.5e-4;  n = fix(T/dt);  U = U_in;
% % Reference solution
% for i = 1:n
%     t = t + dt
%     U0 = U;
%     % PD0 = P(u0)*Laplace
%     [PD0, b_c0] = PMatrix_D(U0, A, g, h);
%     % M = I - dt*a11*P(u0)*Laplace
%     M = speye(3*N^3) - dt*a11*PD0;
%     U1 = M\(U0 + dt*a11*b_c0);
%     % PD1 = P(u1)*Laplace
%     [PD1, b_c1] = PMatrix_D(U1, A, g, h);
%     % M = I - dt*22*d22*P(u0)*Laplace
%     M = speye(3*N^3) - dt*a22*d22*PD1;
%     U2 = M\(U0 + dt*a22*d21*(PD1*U1+b_c1) + dt*a22*d22*b_c1 );
%     U = U0+dt*( b1*(PD0*U1+b_c0) + b2*(d21*(PD1*U1+b_c1) + d22*(PD1*U2+b_c1)) );
%     abs_U =sqrt( U(1:end/3).^2 + U(end/3+1:2*end/3).^2 + U(2*end/3+1:end).^2 );
%     U = U./[abs_U; abs_U; abs_U];
% end
% U_ex = U;
% NematicPlot(U, X, Y, Z)
% 
% 
% Step = [10 20 40 80 160 320];
% % Step = [5 10 20 40 80 160];
% Error = zeros(1,6);
% for j = 1:6
%     t = 0; n = Step(j);dt = T/n;U = U_in;
%     for i = 1:n
%         t = t + dt
%         U0 = U;
%     % PD0 = P(u0)*Laplace
%     [PD0, b_c0] = PMatrix_D(U0, A, g, h);
%     % M = I - dt*a11*P(u0)*Laplace
%     M = speye(3*N^3) - dt*a11*PD0;
%     U1 = M\(U0 + dt*a11*b_c0);
%     % PD1 = P(u1)*Laplace
%     [PD1, b_c1] = PMatrix_D(U1, A, g, h);
%     % M = I - dt*22*d22*P(u0)*Laplace
%     M = speye(3*N^3) - dt*a22*d22*PD1;
%     U2 = M\(U0 + dt*a22*d21*(PD1*U1+b_c1) + dt*a22*d22*b_c1 );
%     U = U0+dt*( b1*(PD0*U1+b_c0) + b2*(d21*(PD1*U1+b_c1) + d22*(PD1*U2+b_c1)) );
%     abs_U =sqrt( U(1:end/3).^2 + U(end/3+1:2*end/3).^2 + U(2*end/3+1:end).^2 );
%     U = U./[abs_U; abs_U; abs_U];
%     end
%     Error(j) = sqrt(sum(h^3*(U-U_ex).^2))
% end
% order = log2(Error(1:end-1)./Error(2:end))

