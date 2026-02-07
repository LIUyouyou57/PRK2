function [PD, b_c] = PMatrix_D(u, A, g, h)
% PD=(I-u*u'/|u|^2)Laplace
% b_c = (I-u*u'/|u|^2)Laplace项中边界条件多出的值
% g_bottom = [g{1,1}, g{1,2}, g{1,3}];
% g_top = [g{2,1}, g{2,2}, g{2,3}];

[N, ~] = size(g{1});
% z-Dirichlet B.C.
b_c1 = zeros(N, N, N-2);  b_c2 = zeros(N, N, N-2);  b_c3 = zeros(N, N, N-2);  
b_c1(:,:,end) = b_c1(:,:,end) + g{2,1}(:,:,end);  b_c1(:,:,1) = b_c1(:,:,1) + g{1,1}(:,:,1);
b_c2(:,:,end) = b_c2(:,:,end) + g{2,2}(:,:,end);  b_c2(:,:,1) = b_c2(:,:,1) + g{1,2}(:,:,1);
b_c3(:,:,end) = b_c3(:,:,end) + g{2,3}(:,:,end);  b_c3(:,:,1) = b_c3(:,:,1) + g{1,3}(:,:,1);
b_c = [b_c1(:); b_c2(:); b_c3(:)]; % Laplace项中边界条件多出的值
b_c = 1/h^2*sparse(b_c);

% u = (u1; u2; u3)
K = length(u);
abs_u =sqrt( u(1 : end/3).^2 + u(end/3+1 : end/3*2).^2 + u(end/3*2+1 : end).^2 ); % |u|;
u = u./[abs_u; abs_u; abs_u]; % u/|u|;

u11 = spdiags(u(1:end/3).^2, 0, K/3, K/3);
u12 = spdiags(u(1:end/3).*u(end/3+1:end/3*2), 0, K/3, K/3);
u13 = spdiags(u(1:end/3).*u(end/3*2+1:end), 0, K/3, K/3);
u22 = spdiags(u(end/3+1:end/3*2).*u(end/3+1:end/3*2), 0, K/3, K/3);
u23 = spdiags(u(end/3+1:end/3*2).*u(end/3*2+1:end), 0, K/3, K/3);
u33 = spdiags(u(end/3*2+1:end).*u(end/3*2+1:end), 0, K/3, K/3);

zero = sparse(zeros(size(A)));
Dh = [A zero zero; zero A zero; zero zero A];
P = spdiags(ones(K,1),0,K,K) - [u11 u12 u13; u12 u22 u23; u13 u23 u33];
PD = P*Dh;
b_c = P*b_c;

end