% x, y-Neumann
function A = laplace_neumann_2D(N, h)
    % N: 每行/列的网格点数（包含边界）
    % 返回二维拉普拉斯矩阵 A，满足诺伊曼边界条件

    % 构造一维拉普拉斯矩阵 T
    e = ones(N, 1);
    T1 = spdiags([e, -2*e, e], -1:1, N, N);

    % 修改边界条件: 左右边界满足 Neumann 条件
    T1(1, 1:2) = [-2, 2];      % 左边界
    T1(end, end-1:end) = [2, -2]; % 右边界

    % 构造二维拉普拉斯矩阵
    I1 = speye(N);  
    A = kron(I1, T1) + kron(T1, I1);
    A = 1/h^2*A;
end
