function [m1,m2,m3] = normalize_m(m1,m2,m3)
n = sqrt(m1.^2+m2.^2+m3.^2);
m1 = m1./n; m2 = m2./n; m3 = m3./n;
end