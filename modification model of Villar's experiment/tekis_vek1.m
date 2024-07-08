function [masofa,kordinat]=tekis_vek(x0,y0,z0,cos_f,a,b,c,d)
masofa=(d-a*x0-b*y0-c*z0)/(a*cos_f(1)+b*cos_f(2)+c*cos_f(3));
%masofa=abs(mas);
x=x0+masofa*cos_f(1);
y=y0+masofa*cos_f(2);
z=z0+masofa*cos_f(3);
kordinat=[x y z];
