h = 1/64; 
error = (m1_new-m1).^2 + (m2_new-m2).^2 + (m3_new-m3).^2;
error1 = h^2*sqrt( sum(error(:)) )
error2 = ( norm(m1(:)-m1_new(:),inf)+norm(m2(:)-m2_new(:),inf)+norm(m3(:)-m3_new(:),inf) )/3

% h = 1/64; 
% error = (m1_new-m1_n).^2 + (m2_new-m2_n).^2 + (m3_new-m3_n).^2;
% error1 = h^2*sqrt( sum(error(:)) )
% error2 = ( norm(m1_n(:)-m1_new(:),inf)+norm(m2_n(:)-m2_new(:),inf)+norm(m3_n(:)-m3_new(:),inf) )/3