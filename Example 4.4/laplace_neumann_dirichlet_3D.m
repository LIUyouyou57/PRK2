% x,y-Neumann; z-Dirihlet
function A = laplace_neumann_dirichlet_3D(N, h)
    % N: 每行/列的网格点数（包含边界）
    % 返回二维拉普拉斯矩阵 A，满足诺伊曼边界条件

    % 构造一维拉普拉斯矩阵 T
    e = ones(N, 1);
    T1 = spdiags([e, -2*e, e], -1:1, N, N);

    % 修改边界条件: 左右边界满足 Neumann 条件
    T1(1, 1:2) = [-2, 2];      % 左边界
    T1(end, end-1:end) = [2, -2]; % 右边界

    % N-2: 内部点的数量（每行/列）
    % 输出A: (N-2)^2 x (N-2)^2 的二维拉普拉斯矩阵 (狄利克雷条件)

    % 一维拉普拉斯矩阵 (内部点)
    e = ones(N-2, 1);
    T2 = spdiags([e, -2*e, e], -1:1, N-2, N-2);

    % 构造三维拉普拉斯矩阵
    I1 = speye(N);  I2 = speye(N-2); 

    A = kron(I2, kron(I1, T1)) + kron(I2, kron(T1, I1)) + kron(T2, kron(I1, I1));
    A = 1/h^2*sparse(A);
end
