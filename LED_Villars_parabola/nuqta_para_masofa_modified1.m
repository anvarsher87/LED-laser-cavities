% joy 5 uchun tepasi
function d=nuqta_para_masofa_modified1(r,dir,para)

p=para(1); 
xc=para(2);
yc=para(3);
zc=para(4);
surish=para(5);
r_i=r;

a=p*dir(1)^2; b=p*2*dir(1)*(r_i(1)+surish)+dir(3); c=p*(r_i(1)+surish)^2+r_i(3)-zc;

D=b^2-4*a*c;

if D<0
    d=Inf;
    1;
elseif D==0
    d=-b/a/2;
    2;
else
    3;
    if a==0
        d=-c/b;4;
    else
        dr(1)=(-b+sqrt(b^2-4*a*c))/(2*a); 
        if dr(1)<0.000001; dr(1)=Inf; end
        dr(2)=(-b-sqrt(b^2-4*a*c))/(2*a) ;
        if dr(2)<0.000001; dr(2)=Inf; end
   
        d=min(dr); 
        dmax = max(dr);
        
    end
%     dr
%     555
    r_f=r_i+d*dir;
    
    if (r_f(1)<0) && (dmax < Inf)
        d=dmax;
        55551;
    elseif (r_f(1)<0) && (dmax == Inf)
        55552;
        d=Inf;
    end
end
end