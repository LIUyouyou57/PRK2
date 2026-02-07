function g = Dirchlet_bc(N)
g = cell(6,3);
% g_left = [g{1,1}, g{1,2}, g{1,3}];
% g_right = [g{2,1}, g{2,2}, g{2,3}];
% g_bottom = [g{3,1},...;
% g_top = [g{4,1},...;
% g_front = [g{5,1},...;
% g_back = [g{6,1},...;

% % g{i} = [u1, u2, u3];
% g{1} = [ones(N,N), zeros(N,N), zeros(N,N)]; 
% g{2} = g{1};
% g{3} = g{1};
% g{4} = g{1};
% g{5} = g{1};
% g{6} = g{1};

% ex1
h = 1/(N+1);
x = 0: h: 1;  y = x;  z = x;
[X, Y, Z] = meshgrid(x, y, z);
X = X - 0.5;  Y = Y - 0.5;  Z = Z - 0.5;
temp = sqrt( X.^2 + Y.^2 + Z.^2 );
X = X./temp;  Y = Y./temp;  Z = Z./temp;  
g{1,1} = X(2:end-1, 1, 2:end-1);  g{1,2} = Y(2:end-1, 1, 2:end-1);  g{1,3} = Z(2:end-1, 1, 2:end-1);
g{2,1} = X(2:end-1, end, 2:end-1);  g{2,2} = Y(2:end-1, end, 2:end-1);  g{2,3} = Z(2:end-1, end, 2:end-1);
g{3,1} = X(2:end-1, 2:end-1, 1);  g{3,2} = Y(2:end-1, 2:end-1, 1);  g{3,3} = Z(2:end-1, 2:end-1, 1);
g{4,1} = X(2:end-1, 2:end-1, end);  g{4,2} = Y(2:end-1, 2:end-1, end);  g{4,3} = Z(2:end-1, 2:end-1, end);
g{5,1} = X(end, 2:end-1, 2:end-1);  g{5,2} = Y(end, 2:end-1, 2:end-1);  g{5,3} = Z(end, 2:end-1, 2:end-1);
g{6,1} = X(1, 2:end-1, 2:end-1);  g{6,2} = Y(1, 2:end-1, 2:end-1);  g{6,3} = Z(1, 2:end-1, 2:end-1);

end