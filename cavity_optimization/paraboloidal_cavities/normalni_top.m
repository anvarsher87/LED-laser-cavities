function n=normalni_top(kordinata, yuza_tur)

p_r=kordinata;
s_a=yuza_tur;
% coefs; check Granino A. Kom and Theresa M. Kom  MATHEMATICAL HANDBOOK
% Eq. 3.5.15

A=s_a(1,1)*p_r(1)+s_a(1,2)*p_r(2)+s_a(1,3)*p_r(3)+s_a(1,4);
B=s_a(2,1)*p_r(1)+s_a(2,2)*p_r(2)+s_a(2,3)*p_r(3)+s_a(2,4);
C=s_a(3,1)*p_r(1)+s_a(3,2)*p_r(2)+s_a(3,3)*p_r(3)+s_a(3,4);

% vectoring modulini topish

m=sqrt(A^2+B^2+C^2);

cos_x=A/m;
cos_y=B/m;
cos_z=C/m;

n=[cos_x,cos_y,cos_z]; % normal of surface

end