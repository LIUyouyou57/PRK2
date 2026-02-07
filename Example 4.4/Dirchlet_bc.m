function g = Dirchlet_bc(N)
g = cell(2,3);

% g_bottom = [g{1,1}, g{1,2}, g{1,3}];
% g_top = [g{2,1}, g{2,2}, g{2,3}];

% 测阶
% g{1} = [ones(N,N), zeros(N,N), zeros(N,N)]; 
% g{2} = [zeros(N,N), ones(N,N), zeros(N,N)]; 
g{1,1} = zeros(N, N, N);  g{1,2} = zeros(N, N, N);  g{1,3} = zeros(N, N, N);
g{2,1} = zeros(N, N, N);  g{2,2} = zeros(N, N, N);  g{2,3} = zeros(N, N, N);
g{1,1}(:,:,1) = ones(N,N,1);
g{2,2}(:,:,end) = ones(N,N,1);


% % ex1
% g{1} = [ones(N,N), zeros(N,N), zeros(N,N)]; 
% g{2} = [zeros(N,N), ones(N,N), zeros(N,N)]; 

end