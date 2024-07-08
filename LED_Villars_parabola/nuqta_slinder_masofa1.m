% slindiracha bogan masofani topish

function out1=nuqta_slinder_masofa1(kordinata,yonalish, radius_aktiv_muhit, xc, yc)

in1 = kordinata;

in2 = yonalish;

% a=in2(1)^2+in2(2)^2;
% b=2*(in1(1)-xc)*in2(1)+2*(in1(2)-yc)*in2(2);
% c=(in1(1)-xc)^2+(in1(2)-yc)^2-radius_aktiv_muhit^2;
a=in2(1)^2+in2(3)^2;
b=2*(in1(1)-xc)*in2(1)+2*(in1(3)-yc)*in2(3);
c=(in1(1)-xc)^2+(in1(3)-yc)^2-radius_aktiv_muhit^2;

disk=b^2-4*a*c;

if disk <= 0
    
    out1=Inf;
    
else
    
    out11(1) = (-b+sqrt(disk))/(2*a);
    out11(2) = (-b-sqrt(disk))/(2*a);
    
    if out11(1)<1e-9; out11(1)=0; end; if out11(2)<1e-9; out11(2)=0; end;
%     out11
    if out11(1)>0 && out11(2)>0
        out1=min(out11);
    elseif out11(1)<=0 && out11(2)<=0
        out1=Inf;
    else
        out1=max(out11);
    end
end
