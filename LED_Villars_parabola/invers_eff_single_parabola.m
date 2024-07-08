function eff = invers_eff_single_parabola(p, h)
% clear; clc; close all;
% % all in mm
% % kk=linspace(0.03, 0.15, 13); 
% % hh=linspace(3.5, 5.5, 3);
% % for k=1:13
% xlabel('x')
% ylabel('y')
% zlabel('z')
% % % view(40, 30)
% % view(0, 0)
% % hold on
kk=1/(2*p);
% p=12; % 3 dan 13 gacha
% k=0.03; 0.03 dan 0.15 gacha

% h=20;
D=25;
aa=12.5; %parabola shohlari extremumlar
% surish=12.5-aa; %surish uchun
surish=0;
% p=0.2; %parabola parametri
% p=0.04;
zc=kk*aa^2; % surish uchun
%---------------------parabola _1------------------------------------------
xx1=linspace(-aa,aa,20); xp(1,:)=xx1; xp(2,:)=xx1;
zz1=-kk*(xx1).^2+zc; zp(1,:)=zz1; zp(2,:)=zz1;
mm=2; n=19;
y=(0:mm-1)'/(mm-1) * ones(1,n+1);
para_coefs_1=zeros(4,4); para_coefs_1(1,1)=kk; para_coefs_1(1,4)=-kk*surish; para_coefs_1(4,1)=-kk*surish; para_coefs_1(4,4)=-zc+kk*surish^2; para_coefs_1(3,4)=1/2; para_coefs_1(4,3)=1/2;
% surf(zp,25*y,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.1')
%----------------------------rasm yon devor (qopqoq va o'rtakash)----------
x=[xx1];
z=[zz1];
y=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
%--------------------------------------------------------------------------
x1=[xx1];
z1=[zz1];
y1=[D D D D D D D D D D D D D D D D D D D D];
% patch(z1, y1,  x1,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
%----------------------LED-------------------------------------------------
a=0; b=0; c=1; d=0;
yy=linspace(0,D,2); xx=linspace(-12.5,12.5,2);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z=(d-a*x-b*y)/c;
% surf(z, y, x, 'EdgeColor', 'k', 'FaceColor', 'w', 'FaceAlpha', '1'); % Plot the surface
p1=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tan11=[-1 0 0]; tan12=[0 -1 0];
tek_1=[a b c d];
%----------------------------- LED rasmi ----------------------------------
yy=linspace(-2,D+2,50); xx=linspace(-12.5-2,12.5+2,50);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
% surf(z+0.1, y, x,  'EdgeColor', 'None', 'FaceColor', 'b', 'FaceAlpha', '0.1') % Plot the surface
% % surf(z+D/2-2,  y+h+LED_rod_masofa, x,'EdgeColor', 'None', 'FaceColor', 'w', 'FaceAlpha', '1') % Plot the surface
%------------------------------LED rasmi (orqa devor)----------------------
yy=linspace(-2,D+2,2); xx=linspace(-12.5-2,12.5+2,2);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
% surf(z-3, y, x,  'EdgeColor', 'None', 'FaceColor', 'c', 'FaceAlpha', '0.1') % Plot the surface
%------------------------------LED rasmi (yon devor)-----------------------
r_k=20.5; theta=45; %for rotation plane
[xk, yk, zk]=cylinder(r_k, 4);
X2 = xk*cosd(theta) - yk*sind(theta);
Y2 = xk*sind(theta) + yk*cosd(theta);
% surf(zk*3-3, Y2+12.5, X2,'EdgeColor', 'k', 'FaceColor',  'c', 'FaceAlpha', '0.05')%[0.1 0.155 0.198]
%--------------------qopqoq------------------------------------------------
a=0; b=-1; c=0; d=0;
xx=linspace(-12.5,12.5,50); zz=linspace(0,zc,50);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-c*z-a*x)/b;
% % surf(z, y, x, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.3'); % Plot the surface
p2=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_2=[a b c d];
% xdfcgvbh
%--------------------mavhum devor------------------------------------------
a=0; b=1; c=0; d=25;
xx=linspace(-12.5,12.5,50); zz=linspace(0,zc,50);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-c*z-a*x)/b;
% surf(z, y, x, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.3'); % Plot the surface
p3=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_3=[a b c d];
%---------------------------rod--------------------------------------------
r_muhit=2.5;
l=h-r_muhit; %Led dan muhit o'rtasigacha masofa
[xs,ys,zs]=cylinder(r_muhit,20); %surf(xs+l, 25*zs, ys, 'EdgeColor', 'None', 'FaceColor', 'm', 'FaceAlpha', '0.3')
% plotCircle3D_1([0, D+0.1, l],[0,1,0], r_muhit)
% plotCircle3D_1([0, 0-0.1, l],[0,1,0], r_muhit)
%-------------------------crystal------------------------------------------
A=dlmread('15_12_2021_BLUE_LED_ems.txt');
C=dlmread('Ce_Nd_YAG_abs.txt');
%--------------------------------------------------------------------------
% led_array_parabola
sa=0.9*pi/2;%pi/3;
teflon_kaytarish=0.9;  %blue uchun
led_kaytarish=0.75;    %blue uchun
LEDdan_qaytdi=0;
% lost_at_led=0;
% lost_at_teflon=0;
% chizish=0;
% used_f1=0;
nuqta_fotonlar_soni_max=100000;
% axis equal
sp=zeros(1,2500);
s=zeros(1,nuqta_fotonlar_soni_max);
E_total=0;
E_det1=0; % 
dd=h-2*r_muhit;

for ti=1: nuqta_fotonlar_soni_max
    tasodifiy_son=rand;
    for i=1:1000
        if tasodifiy_son<=A(i)
            ftu=i;  % ftu - fotoning to'lqin uzunligi deganini bildiradi
            break;
        end
    end
    sp(ftu)=sp(ftu)+1;
    s(1,ti)=ftu;
end

% x_grid_number  = 49;        x_grid_size = 2*r_muhit/x_grid_number;
% z_grid_number  = 49;        z_grid_size = 2*r_muhit/z_grid_number;
% y_grid_number  = 499;       y_grid_size = D/y_grid_number;
% 
% E_xyz = zeros(x_grid_number+1,y_grid_number+1,z_grid_number+1);
% sexdrcftvgbhj

for f=1:nuqta_fotonlar_soni_max
    ftu=s(1,f);
    E_total=E_total+(1/ftu);
    %-------------------------------------absorbtion---------------
    abs_coef=C(ftu); abs_length=absorption_length(abs_coef); abs_leng_m=abs_length;
    if (abs_leng_m>1900);  yutilmaydigan_f=yutilmaydigan_f+1; continue; end
    %--------------------------------------------------------------
    %------------- boshlang'ich kordinata va yonalish--------------
 
    r_i=led_random_square_1(0, 12.5, 0, 12.5, 12.5);
      yonalish = lambert_shape_4 (p1, tan11, tan12, sa);
    
    %--------------------------------------------------------------
    
    foton_in_system=1;
    while (foton_in_system == 1)
        %                 1
        %---------------------LEDs-------------------------
        dm(1)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_1);  tk1(1,:)=[p1(1) p1(2) p1(3)];    %past_1 Led turgan joy
        %-------------------------qopqoq-------------------
        dm(2)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish, tek_2);  tk1(2,:)=-[p2(1) p2(2) p2(3)];
        %--------------------------o'rta-------------------
        dm(3)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_3);  tk1(3,:)=-[p3(1) p3(2) p3(3)];
        %---------------------Laser rod--------------------
        dm(4)=nuqta_slinder_masofa1(r_i, yonalish, r_muhit, 0, l); tk1(4,:)= [1/r_muhit 1/r_muhit 0];              %aktiv muhit;
        %---------------------parabolalar------------------
        dm(5)=nuqta_para_masofa1(r_i,yonalish,[kk 0 0 zc]);   tk1(5,:)= [0 0 0 ];   %tepa_1
        %                         dm(6)=nuqta_para_masofa_modified2(r_i,yonalish,[p 0 0 zc surish]);   tk1(6,:)= [0 0 0];  % yon_chap_1
        %--------------------------------------------------
        for i=1:length(dm); if dm(i)<0.001; dm(i)=Inf; end; end
        %                         2
        [db,joy]=min(dm);  normal=tk1(joy,:); r_f=r_i + yonalish*db;
        
        if (db==Inf)||(db==0); foton_in_system=0; break; end
        
%                                  if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','b','LineWidth',0.5); end
%                                  pause
        
        r_i=r_f;
        
        %                             joy
        if (joy==1)
            %                            11
            r_i_n=r_f; r_f_n=r_i_n + normal*1;
%             if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',1.5); end
            if led_kaytarish > rand
                %                                                                                                  112
                yonalish=yonalish-2*normal*dot(yonalish,normal);
%                 LEDdan_qaytdi=LEDdan_qaytdi+1;
            else
                %                                                                                                  113
                foton_in_system=0; %lost_at_led=lost_at_led+1;
            end
        end
        
        if (joy>=2)&&(joy<=3)
            %                             12
            %                          dfghgnxcvbnm
            %                          r_f
            r_i_n=r_f; r_f_n=r_i_n + normal*1;
%             if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',1.5); end
            %                                 joy
            normal;
            
            %                                                                 121
            %
            if teflon_kaytarish > rand
                %                                                                        122
                yonalish=yonalish-2*normal*dot(yonalish,normal);
            else
                %                                                                                                              123
                foton_in_system=0; %lost_at_teflon=lost_at_teflon+1;
            end
            %                                 sdfgbhnjmk
            
        end
        if (joy==4)
            
            %                                                                                      14
            nor=[(r_f(1))/r_muhit 0 (r_f(3)-l)/r_muhit ];
            %                                 nor=-[(r_f(2)-h)/r_muhit  0 (r_f(1))/r_muhit];
            if (rand()>(qaytish_sinish_ehtimolliqi(yonalish,nor,1,1.8)))
                %                                                                                            140
                %                                 s_n=nor;
                %                                 r_i_n=r_f; r_f_n=r_i_n + s_n*0.5;
                %                                 if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',1.5); end
                foton_muhitda=1;
                [~,yonalish]=cheg_cos(yonalish, -nor,1,1.8);
                shs=0;
                while foton_muhitda==1
                    %                                                                                                             1401
                    mm=nuqta_slinder_masofa1(r_i, yonalish, r_muhit, 0, l);
                    if (mm==Inf)||(mm==0); foton_in_system=0; foton_muhitda=0; break; end
                    if (r_f(2)<0)||(r_f(2)>D)                %esdan chiqmasin
                        %
                        %                                                                                                                         143
                        foton_in_system=0;  foton_muhitda=0;
                        break
                    end
                    
                    if (r_i(2)<0)||(r_i(2)>D)                %esdan chiqmasin
                        %                                                                                                                    144
                        foton_in_system=0; foton_muhitda=0;
                        break
                    end
                    if mm>abs_leng_m
                        mm=abs_leng_m;
                        r_f=r_i+mm*yonalish;
                        if (r_f(2)<0)||(r_f(2)>D)
                            foton_muhitda=0; foton_in_system=0;
                        else
                        %                                         scatter3(r_f(3),r_f(2),r_f(1), 15,'fill')
                        %                                             14011
%                         used_f1=used_f1+1;
                        %                                              13111
                        E_det1=E_det1+(1/ftu);
                        foton_muhitda=0; foton_in_system=0;
                        r_f;
                     
%                         x_index = round((r_f(1)+ r_muhit)/x_grid_size)+1;
%                         z_index = round((r_f(3)-dd)/z_grid_size)+1;
%                         y_index = round((r_f(2))/y_grid_size)+1;
%                         E_xyz(x_index,y_index,z_index)=E_xyz(x_index,y_index,z_index)+1;
%                                  if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','r','LineWidth',0.1); end
                        % %        if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','r','LineWidth',0.1); gs=gs+1; end
                        end
                    else
                        %                                                                                                                         14012
                        abs_leng_m=abs_leng_m-mm;
                        r_f=r_i+mm*yonalish;
%                                                               if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','y','LineWidth',0.1); end
                        % %                                         if chizish==1; is(gs)=plot3([r_i(1),r_f(1)],[r_i(2),r_f(2)],[r_i(3),r_f(3)],'Color','r','LineWidth',0.1); gs=gs+1; end
                        nor=-[(r_f(1))/r_muhit 0 (r_f(3)-l)/r_muhit ];
                        %                                    nor=-[(r_f(2)-h)/r_muhit  0 (r_f(1))/r_muhit];
                        %                                         r_i_n=r_f; r_f_n=r_i_n + nor*1;
                        % %                                         if chizish==1; plot3([r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],[r_i_n(3),r_f_n(3)],'Color','y','LineWidth',1); end
                        if  ((acos(dot(-yonalish,nor)/(norm(yonalish)*norm(nor)))) < asin(1/1.8))&&(rand()>(qaytish_sinish_ehtimolliqi(yonalish,nor,1.8,1)))
                            [~,yonalish]=cheg_cos(yonalish,-nor,1.8,1);
                            %                                             r_i_n=r_f; r_f_n=r_i_n + yonalish*2;
                            foton_muhitda=0;
                            %                                               19013
                            %                                             if chizish==1; plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','r','LineWidth',1.5);end
                            
                        else
                            %                                                                                                                                     14014
                            yonalish=yonalish-2*nor*dot(yonalish,nor);
                        end
                        r_i=r_f;
                    end
                    
                    
                    shs=shs+1; if shs>5; foton_muhitda=0; foton_in_system=0; break; end
                end
            else
                %                                                                                                 141
                yonalish=yonalish-2*nor*dot(yonalish,nor);
            end
        end
        if (joy==5)
            %                                                              15
            normal=-normalni_top(r_f, para_coefs_1);
            %                             r_i_n=r_f; r_f_n=r_i_n + normal*1;
            %                             if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',1.5); end
            if teflon_kaytarish > rand
                %                                                                        151
                yonalish=yonalish-2*normal*dot(yonalish,normal);
            else
                %                                                                                                              152
                foton_in_system=0; % lost_at_teflon=lost_at_teflon+1;
            end
        end
        
        
        
    end
end
% end
% used_f1/(nuqta_fotonlar_soni_max*100)
eff=E_det1/E_total;
% % close all
% % end
% axis equal
% % axis off