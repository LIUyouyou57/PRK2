function L = NeumannLaplacian2D(Nx,Ny,hx,hy)
    Nx1 = Nx + 1;
    Ny1 = Ny + 1;
    
    % x方向
    ex = ones(Nx1, 1);
    Tx = spdiags([ex, -2*ex, ex], -1:1, Nx1, Nx1) / hx^2;
    Tx(1, 1:2) = [-2, 2] / hx^2;         % 左边界
    Tx(end, end-1:end) = [2, -2] / hx^2; % 右边界
    
    % y方向
    ey = ones(Ny1, 1);
    Ty = spdiags([ey, -2*ey, ey], -1:1, Ny1, Ny1) / hy^2;
    Ty(1, 1:2) = [-2, 2] / hy^2;         % 下边界
    Ty(end, end-1:end) = [2, -2] / hy^2; % 上边界
    
    % 二维拉普拉斯
    L = kron(speye(Ny1), Tx) + kron(Ty, speye(Nx1));
end