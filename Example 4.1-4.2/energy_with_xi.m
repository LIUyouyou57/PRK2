function E = energy_with_xi(m1,m2,m3,xi,hx,hy)
u1 = m1 + xi; u2 = m2 + xi; u3 = m3 + xi;
n = sqrt(u1.^2+u2.^2+u3.^2);
u1=u1./n; u2=u2./n; u3=u3./n;
u1x = u1(:, 2:end) - u1(:, 1:end-1);
u1y = u1(2:end, :) - u1(1:end-1, :);
u2x = u2(:, 2:end) - u2(:, 1:end-1);
u2y = u2(2:end, :) - u2(1:end-1, :);
u3x = u3(:, 2:end) - u3(:, 1:end-1);
u3y = u3(2:end, :) - u3(1:end-1, :);
E = 0.5*( sum(sum(u1x.^2 + u2x.^2 + u3x.^2)) + sum(sum(u1y.^2 + u2y.^2 + u3y.^2)) );
end