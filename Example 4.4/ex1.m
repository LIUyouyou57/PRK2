% clear;

h = 1/24;
N = 1/h + 1;
x = 0: h: 1;  y = x;  z = h: h: 1-h;
[X, Y, Z] = meshgrid(x, y, z);

% u1 = sin(2*pi*X) + 2 + 0.5*sin(6*pi*Y) + 0.2*sin(4*pi*Z);
% u2 = cos(2*pi*X) - 2 + 0.5*cos(6*pi*Y) + 0.2*cos(4*pi*Z);
% u1 = (Z<0.5).*ones(N, N, N-2) + (Z>=0.5).*zeros(N, N, N-2);
% u2 = (Z<0.5).*zeros(N, N, N-2) + (Z>=0.5).*ones(N, N, N-2);
% u3 = zeros(N, N, N-2);
u1 = 2*rand(N,N,N-2)-1; u2 = 2*rand(N,N,N-2)-1;  u3 = 2*rand(N,N,N-2)-1;

U_in = [u1(:); u2(:); u3(:)];
% 初值投影到球面
abs_U_in =sqrt( U_in(1:end/3).^2 + U_in(end/3+1:2*end/3).^2 + U_in(2*end/3+1:end).^2 );
min(abs_U_in)
U_in = U_in./[abs_U_in; abs_U_in; abs_U_in];

% RK2 系数 d11=1
a11 = 1;  a21 = 0;  a22= 1/2;  b1=1/2;  b2=1/2;  d21=-1;  d22=2;

A = laplace_neumann_dirichlet_3D(N, h);
g = Dirchlet_bc(N);

%% Energy decay
% T = 0.05;  t = 0;  dt = 1e-3;  n = fix(T/dt);  U = U_in;
% 一阶方法能量
% E = zeros(1,n);
% for i = 1:n
%     t = t + dt
%     U0 = U;
%     [PD, b_c] = PMatrix_D(U0, A, g);
%     M = speye(3*N^2) - dt*PD;
%     U = M\(U0 + dt*b_c);
%     abs_U =sqrt( U(1:N^2).^2 + U(N^2+1:2*N^2).^2 + U(2*N^2+1:end).^2 );
%     U = U./[abs_U; abs_U; abs_U];
%     u10 = reshape(U(1:N^2), [N,N]);
%     u20 = reshape(U(N^2+1:2*N^2), [N,N]);
%     u30 = reshape(U(2*N^2+1:end), [N,N]);
%     u1x = (u10 - [u10(end,:); u10(1:end-1,:)]);?????????
%     u1y = (u10 - [u10(:,end) u10(:,1:end-1)]);
%     u2x = (u20 - [u20(end,:); u20(1:end-1,:)]);
%     u2y = (u20 - [u20(:,end) u20(:,1:end-1)]);
%     E(i) = 0.5*sum(sum(u1x.^2 + u2x.^2 + u1y.^2 + u2y.^2));
%     % E0 = 0.5*U0'*(-Dh)*U0*h^2
%     % E1 = 0.5*U'*(-Dh)*U*h^2
%     % U2 = U/sqrt(sum(U.^2)*h^2);
%     % E2 = 0.5*U2'*(-Dh)*U2*h^2
% end
% RK2能量
% for i = 1:n
%     t = t + dt
%     U0 = U;
%     PD0 = PMatrix(U0,N,h);
%     M = spdiags(ones(2*N^2,1),0,2*N^2,2*N^2) - dt*a11*PD0;
%     U1 = M\U0;
%     PD1 = PMatrix(U1,N,h);
%     M = spdiags(ones(2*N^2,1),0,2*N^2,2*N^2) - dt*a22*d22*PD1;
%     U2 = M\(U0+dt*a22*d21*PD1*U1);
%     U = U0+dt*(b1*PD0*U1+b2*PD1*(d21*U1+d22*U2));
%     abs_U =sqrt( U(1:N^2).^2 + U(N^2+1:end).^2);
%     U = U./[abs_U;abs_U];
%     u10 = reshape(U(1:end/2),[N,N]);
%     u20 = reshape(U(end/2+1:end),[N,N]);
%     u1x = (u10 - [u10(end,:);u10(1:end-1,:)]);
%     u1y = (u10 - [u10(:,end) u10(:,1:end-1)]);
%     u2x = (u20 - [u20(end,:);u20(1:end-1,:)]);
%     u2y = (u20 - [u20(:,end) u20(:,1:end-1)]);
%     E(i) = 0.5*sum(sum(u1x.^2+u2x.^2+u1y.^2+u2y.^2));
% end
% vec = [1:1:30, 40:10:100 100:100:1000];
vec = [10 30 100];

%% Second Order Test
T = 1;  t = 0;  dt = 5e-3;  n = fix(T/dt);  U = U_in;  NematicPlot(U, X, 0); close all;
% Reference solution
for i = 1:n
    t = t + dt
    U0 = U;
    % PD0 = P(u0)*Laplace
    [PD0, b_c0] = PMatrix_D(U0, A, g, h);
    % M = I - dt*a11*P(u0)*Laplace
    M = speye(3*N^2*(N-2)) - dt*a11*PD0;
    U1 = M\(U0 + dt*a11*b_c0);
    % PD1 = P(u1)*Laplace
    [PD1, b_c1] = PMatrix_D(U1, A, g, h);
    % M = I - dt*22*d22*P(u0)*Laplace
    M = speye(3*N^2*(N-2)) - dt*a22*d22*PD1;
    U2 = M\(U0 + dt*a22*d21*(PD1*U1+b_c1) + dt*a22*d22*b_c1 );
    U = U0+dt*( b1*(PD0*U1+b_c0) + b2*(d21*(PD1*U1+b_c1) + d22*(PD1*U2+b_c1)) );
    abs_U =sqrt( U(1:end/3).^2 + U(end/3+1:2*end/3).^2 + U(2*end/3+1:end).^2 );
    U = U./[abs_U; abs_U; abs_U];

    if any(vec == i)
        NematicPlot(U, X, i)
        close all;
    else
    end
end
U_ex = U;
NematicPlot(U, X, i)

% Step = [100 200 400 800 1600 3200];
% % Step = [10 20 40 80 160 320];
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
%     M = speye(3*N^2*(N-2)) - dt*a11*PD0;
%     U1 = M\(U0 + dt*a11*b_c0);
%     % PD1 = P(u1)*Laplace
%     [PD1, b_c1] = PMatrix_D(U1, A, g, h);
%     % M = I - dt*22*d22*P(u0)*Laplace
%     M = speye(3*N^2*(N-2)) - dt*a22*d22*PD1;
%     U2 = M\(U0 + dt*a22*d21*(PD1*U1+b_c1) + dt*a22*d22*b_c1 );
%     U = U0+dt*( b1*(PD0*U1+b_c0) + b2*(d21*(PD1*U1+b_c1) + d22*(PD1*U2+b_c1)) );
%     abs_U =sqrt( U(1:end/3).^2 + U(end/3+1:2*end/3).^2 + U(2*end/3+1:end).^2 );
%     U = U./[abs_U; abs_U; abs_U];
%     end
%     Error(j) = sqrt(sum(h^3*(U-U_ex).^2))
% end
% order = log2(Error(1:end-1)./Error(2:end))

