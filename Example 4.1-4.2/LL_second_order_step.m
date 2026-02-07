function [m1_new,m2_new,m3_new,lambda_new] = ...
    LL_second_order_step(m1,m2,m3,lambda,dt,L,beta,gamma)
% -------------------------------------------------
% One step second-order Type-I scheme
% For bootstrapping BDF4
%
% m_t + beta m x Delta m = gamma(Delta m + lambda m)
% |m| = 1
% -------------------------------------------------

sz = size(m1);
N  = numel(m1);

%% ---------- Laplacian ----------
Lap1 = reshape(L*m1(:),sz);
Lap2 = reshape(L*m2(:),sz);
Lap3 = reshape(L*m3(:),sz);

%% ---------- cross product m x Delta m ----------
C1 = m2.*Lap3 - m3.*Lap2;
C2 = m3.*Lap1 - m1.*Lap3;
C3 = m1.*Lap2 - m2.*Lap1;

%% ---------- Step 1: predictor (linear Helmholtz) ----------
A = speye(N)/dt - gamma*L;

RHS1 = m1/dt + gamma*lambda.*m1 - beta*C1;
RHS2 = m2/dt + gamma*lambda.*m2 - beta*C2;
RHS3 = m3/dt + gamma*lambda.*m3 - beta*C3;

m1t = reshape(A\RHS1(:),sz);
m2t = reshape(A\RHS2(:),sz);
m3t = reshape(A\RHS3(:),sz);

%% ---------- Step 2: corrector (exact formula) ----------
v1 = m1t - 0.5*gamma*dt*lambda.*m1;
v2 = m2t - 0.5*gamma*dt*lambda.*m2;
v3 = m3t - 0.5*gamma*dt*lambda.*m3;

nv = sqrt(v1.^2 + v2.^2 + v3.^2);

m1_new = v1./nv;
m2_new = v2./nv;
m3_new = v3./nv;

lambda_new = (1 - nv)/(0.5*gamma*dt);

end
