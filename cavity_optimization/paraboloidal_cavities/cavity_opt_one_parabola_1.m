clear; clc; close all;
% global D L dd l r_muhit A C sa teflon_kaytarish led_kaytarish nuqta_fotonlar_soni_max
% D=50;
% L=25;
% dd=5;
% l=8;
% r_muhit=2.5;
% %%--------------------------------crystal----------------------------------
% A=dlmread('15_12_2021_BLUE_LED_ems.txt');
% C=dlmread('Ce_Nd_YAG_abs.txt');
% %%------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------
% sa=0.9*pi/2;
% teflon_kaytarish=0.9;  %blue uchun
% led_kaytarish=0.75;    %blue uchun
% %%--------------sanash_uchun----------------------------------------------------
% nuqta_fotonlar_soni_max=1000;
% N=46; 
% aa=9; %parabola shohlari extremumlar
% surish=12.5-aa; %surish uchun
% % p=0.12; %parabola parametri %0.12 dan 0.3 gacha
% zc=p*aa^2; % surish uchun
% xx2=linspace(0,aa+surish,20); xp(1,:)=xx2; xp(2,:)=xx2;
% zz2=-p*(xx2-surish).^2+zc; zp(1,:)=zz2; zp(2,:)=zz2;
aa=12.5;
aa1=(aa^2)/2;
% p=linspace(0.03, 0.15, 25);
p=linspace(3, 11, 41);
% kk=1/(2*p);
% vv=p()
% vbnm
% % zz=floor(68.75*p);
zck = zeros(size(p));
h =linspace(6, 20, 15); 
% for index=1:41
%     kfk=1:(2*p);
% %     index
% end
% sxdcfvgbh
tic
for i=1:length(p)
      
      zck(i)=aa1/p(i);
      
    for k=1:length(h)
        
        
        if (zck(i)-h(k))<1
%            (zck(i)-h(k))
           eff(i,k)=0;
        else
           eff(i,k) = invers_eff_single_parabola(p(i), h(k));
%            eff(i, k)=effektivlik_topish_one_parabola(p(i), h(k));
        end
        
    end
    %          toc
    %       tgyubhjknl
    i
    %     pause(0.2)
    
end
view(0, 90)

[Hh, Pp]=meshgrid(h, p);

surf(Hh,Pp, eff)
colormap jet