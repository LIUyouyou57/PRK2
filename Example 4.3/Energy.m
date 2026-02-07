function  result = Energy(U, g)
% g_left = [g{1,1}, g{1,2}, g{1,3}];
% g_right = [g{2,1}, g{2,2}, g{2,3}];
% g_bottom = [g{3,1},...;
% g_top = [g{4,1},...;
% g_front = [g{5,1},...;
% g_back = [g{6,1},...;

N = round( (length(U)/3)^(1/3) );

u10 = reshape(U(1:end/3), [N, N, N]);
u20 = reshape(U(end/3+1:2*end/3), [N, N, N]);
u30 = reshape(U(2*end/3+1:end), [N, N, N]);

u1x = cat(1, u10, g{5,1}) - cat(1, g{6,1}, u10);
u1y = cat(2, u10, g{2,1}) - cat(2, g{1,1}, u10);
u1z = cat(3, u10, g{4,1}) - cat(3, g{3,1}, u10);
u2x = cat(1, u20, g{5,2}) - cat(1, g{6,2}, u20);
u2y = cat(2, u20, g{2,2}) - cat(2, g{1,2}, u20);
u2z = cat(3, u20, g{4,2}) - cat(3, g{3,2}, u20);
u3x = cat(1, u30, g{5,3}) - cat(1, g{6,3}, u30);
u3y = cat(2, u30, g{2,3}) - cat(2, g{1,3}, u30);
u3z = cat(3, u30, g{4,3}) - cat(3, g{3,3}, u30);

result = 0.5*( sum(sum(sum(u1x.^2 + u2x.^2 + u3x.^2))) + sum(sum(sum(u1y.^2 + u2y.^2 + u3y.^2))) + sum(sum(sum(u1z.^2 + u2z.^2 + u3z.^2))) );

end
