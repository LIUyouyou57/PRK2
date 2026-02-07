function [PD] = PMatrix_D(u, A, alpha, beta)
% PD={alpha*(I-u*u'/|u|^2) + beta*u/|u|×} Laplace
% u = (u1; u2; u3)

K = length(u);
u1 = u(1:end/3);  u2 = u(end/3+1:end/3*2);  u3 = u(end/3*2+1:end);
abs_u =sqrt( u1.^2 + u2.^2 + u3.^2 ); % |u|;
u = u./[abs_u; abs_u; abs_u]; % u/|u|;
u1 = u(1:end/3);  u2 = u(end/3+1:end/3*2);  u3 = u(end/3*2+1:end);

u11 = spdiags(u1.^2, 0, K/3, K/3);
u12 = spdiags(u1.*u2, 0, K/3, K/3);
u13 = spdiags(u1.*u3, 0, K/3, K/3);
u22 = spdiags(u2.^2, 0, K/3, K/3);
u23 = spdiags(u2.*u3, 0, K/3, K/3);
u33 = spdiags(u3.^2, 0, K/3, K/3);

uu1 = spdiags(u1, 0, K/3, K/3);
uu2 = spdiags(u2, 0, K/3, K/3);
uu3 = spdiags(u3, 0, K/3, K/3);


zero = sparse(zeros(size(A)));
Dh = [A zero zero; zero A zero; zero zero A];
P = alpha*( spdiags(ones(K,1),0,K,K) - [u11 u12 u13; u12 u22 u23; u13 u23 u33] ) + beta*( [zero -uu3 uu2; uu3 zero -uu1; -uu2 uu1 zero] );
PD = P*Dh;

end