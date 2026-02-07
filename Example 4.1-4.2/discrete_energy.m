function E = discrete_energy(m1,m2,m3,hx,hy,L,W)
% m1x = m1(:, 2:end) - m1(:, 1:end-1);
% m1y = m1(2:end, :) - m1(1:end-1, :);
% m2x = m2(:, 2:end) - m2(:, 1:end-1);
% m2y = m2(2:end, :) - m2(1:end-1, :);
% m3x = m3(:, 2:end) - m3(:, 1:end-1);
% m3y = m3(2:end, :) - m3(1:end-1, :);
% E = 0.5*( sum(sum(m1x.^2 + m2x.^2 + m3x.^2)) + sum(sum(m1y.^2 + m2y.^2 + m3y.^2)) );
    Lm1 = reshape(L * m1(:), size(m1));
    Lm2 = reshape(L * m2(:), size(m2));
    Lm3 = reshape(L * m3(:), size(m3));
    
    total = sum(sum( W.*(m1.*Lm1 + m2.*Lm2 + m3.*Lm3) ));
    E = -0.5 * total * hx * hy;
end