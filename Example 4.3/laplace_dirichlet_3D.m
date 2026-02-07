function A = laplace_dirichlet_3D(N, h)
    % N: 每个方向内部点的数量
    % 返回A: 三维拉普拉斯矩阵，满足Dirichlet边界条件
    
    % 构造一维拉普拉斯矩阵
    e = ones(N, 1);
    T = spdiags([e, -2*e, e], -1:1, N, N);

    % 构造三维拉普拉斯矩阵
    I = speye(N);
    A = kron(kron(I, I), T) + kron(kron(I, T), I) + kron(kron(T, I), I);
    A = 1/h^2*sparse(A);
end
