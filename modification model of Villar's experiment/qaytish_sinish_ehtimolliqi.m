
function R_ref=qaytish_sinish_ehtimolliqi(yonalish,nv,n1,n2)

% alfa=pi-acos(dot(yonalish,nv)/(norm(yonalish)*norm(nv)))
% 
% beta=asin(n1*si
% n(alfa)/n2)
% 
% R_ref=0.5*(0*((sin(alfa-beta))^2/(sin(alfa+beta))^2 )+((tan(alfa-beta))^2/(tan(alfa+beta))^2));
a=yonalish; n=nv;   
     
     alfa=abs(acos((dot(a,n))));
     
%      c=0;
     if alfa>pi/2; alfa=pi-alfa; end 
     
     beta=abs(asin(n1/n2*sin(alfa)));
     
     if beta>pi/2; beta=pi-beta; end
     
R_ref=0.5*((sin(alfa-beta)/sin(alfa+beta))^2+(tan(alfa-beta)/tan(alfa+beta))^2);

%pause(1)


