function  result = Energy(U,h)

N = sqrt(length(U)/3);

u10 = reshape(U(1:end/3), [N,N]);
u20 = reshape(U(end/3+1:2*end/3), [N,N]);
u30 = reshape(U(2*end/3+1:end), [N,N]);

u1x = u10(:, 2:end) - u10(:, 1:end-1);
u1y = u10(2:end, :) - u10(1:end-1, :);
u2x = u20(:, 2:end) - u20(:, 1:end-1);
u2y = u20(2:end, :) - u20(1:end-1, :);
u3x = u30(:, 2:end) - u30(:, 1:end-1);
u3y = u30(2:end, :) - u30(1:end-1, :);
result = 0.5*( sum(sum(u1x.^2 + u2x.^2 + u3x.^2)) + sum(sum(u1y.^2 + u2y.^2 + u3y.^2)) );

% [u1x,u1y] = gradient(u10,h,h);
% [u2x,u2y] = gradient(u20,h,h);
% [u3x,u3y] = gradient(u30,h,h);
% result = 0.5*sum(u1x(:).^2 + u2x(:).^2 + u3x(:).^2 + u1y(:).^2 + u2y(:).^2 + u3y(:).^2)*h^2 ;

end
