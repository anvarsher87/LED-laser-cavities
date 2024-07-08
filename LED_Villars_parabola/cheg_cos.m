% biror yuzadan qaytgan nuring yo'nalish cosinusini topadi

function [b,c]=cheg_cos(a,v,n1,n2)

%n=-n;
% a=tushuvchi_nur_cosinusi; v=normal;  n1=sindrish_korsatkich_tushuvchi; n2=sindrish_korsatgich_sinuvchi; 

%  alfa=abs(acos(dot(a,n));
dt=dot(a,v); 

d=1-n1^2*(1-dt^2)/(n2^2);

if d>0
    %(n2<n1)&&(abs(acos(dot(a,n))) > asin(n2/n1))
    
     b=a-2*v*dot(a,v);
     d=sqrt(d);  
%      c=(n1*(a-v*dt))/n2+v*d;   
     c=(n1*(a-v*dt))/n2+v*d;
else

b=a-2*v*dot(a,v);  
c=[0 0 0];

end
    
end