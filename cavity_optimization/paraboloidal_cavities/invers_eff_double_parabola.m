function eff = invers_eff_double_parabola(p, h)
% clear; close all;
%Villars_parabola_3 / Villars_parabola_working dan olindi
% all in mm
% xlabel('x')
% ylabel('y')
% zlabel('z')
% % view(40, 30)
% view(0, 0)
% hold on
KK=1/(2*p);
D=25;
aa=9.5; %parabola shohlari extremumlar
surish=12.5-aa; %surish uchun
% p=0.12; %parabola parametri %0.12 dan 0.3 gacha
zc=KK*aa^2; % surish uchun
%---------------------parabola _1------------------------------------------
xx1=linspace(-aa-surish,0,20); xp(1,:)=xx1; xp(2,:)=xx1;
zz1=-KK*(xx1+surish).^2+zc; zp(1,:)=zz1; zp(2,:)=zz1;
mm=2; n=19;
y=(0:mm-1)'/(mm-1) * ones(1,n+1);
para_coefs_1=zeros(4,4); para_coefs_1(1,1)=KK; para_coefs_1(1,4)=-KK*surish; para_coefs_1(4,1)=-KK*surish; para_coefs_1(4,4)=-zc+KK*surish^2; para_coefs_1(3,4)=1/2; para_coefs_1(4,3)=1/2;
% surf(zp,25*y,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.1')
%-------------------parabola_2---------------------------------------------
xx2=linspace(0,aa+surish,20); xp(1,:)=xx2; xp(2,:)=xx2;
zz2=-KK*(xx2-surish).^2+zc; zp(1,:)=zz2; zp(2,:)=zz2;
mm=2; n=19;
% y=(0:mm-1)'/(mm-1) * ones(1,n+1);
para_coefs_2=zeros(4,4); para_coefs_2(1,1)=KK; para_coefs_2(1,4)=KK*surish; para_coefs_2(4,1)=KK*surish; para_coefs_2(4,4)=-zc+KK*surish^2; para_coefs_2(3,4)=1/2; para_coefs_2(4,3)=1/2; %para_coefs(4,4)=1/2;
% surf(zp,25*y,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.1')
%----------------------------rasm yon devor (qopqoq va o'rtakash)----------
x=[xx1 xx2];
z=[zz1 zz2];
y=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
%----------------------------------------------------------------------------
x=[xx1 xx2];
z=[zz1 zz2];
y=[D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D];
% patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
%----------------------LED-------------------------------------------------
a=0; b=0; c=1; d=0;
yy=linspace(0,D,2); xx=linspace(-12.5,12.5,2);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z=(d-a*x-b*y)/c;
% surf(z, y, x, 'EdgeColor', 'k', 'FaceColor', 'w', 'FaceAlpha', '1'); % Plot the surface
p1=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tan11=[-1 0 0]; tan12=[0 -1 0];
tek_1=[a b c d];
%----------------------------- LED rasmi -------------------------------------------------------
yy=linspace(-2,D+2,50); xx=linspace(-12.5-2,12.5+2,50);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
% surf(z+0.1, y, x,  'EdgeColor', 'None', 'FaceColor', 'b', 'FaceAlpha', '0.1') % Plot the surface
%------------------------------LED rasmi (orqa devor)--------------------------------------------
yy=linspace(-2,D+2,2); xx=linspace(-12.5-2,12.5+2,2);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
% surf(z-3, y, x,  'EdgeColor', 'None', 'FaceColor', 'c', 'FaceAlpha', '0.1') % Plot the surface
%------------------------------LED rasmi (yon devor)---------------------------------------------
r_k=20.5; theta=45; %for rotation plane
[xk, yk, zk]=cylinder(r_k, 4);
X2 = xk*cosd(theta) - yk*sind(theta);
Y2 = xk*sind(theta) + yk*cosd(theta);
% surf(zk*3-3, Y2+12.5, X2,'EdgeColor', 'k', 'FaceColor',  'c', 'FaceAlpha', '0.05')%[0.1 0.155 0.198]
%--------------------qopqoq------------------------------------------------
a=0; b=-1; c=0; d=0;
% xx=linspace(-12.5,12.5,50); zz=linspace(0,zc,50);
% [x, z] = meshgrid(xx,zz); % Generate x and y data
% y=(d-c*z-a*x)/b;
% % surf(z, y, x, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.3'); % Plot the surface
% p2=[a b c d];
p2=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_2=[a b c d];
% xdfcgvbh
%--------------------mavhum devor------------------------------------------------
a=0; b=1; c=0; d=25;
% xx=linspace(-12.5,12.5,50); zz=linspace(0,zc,50);
% [x, z] = meshgrid(xx,zz); % Generate x and y data
% y=(d-c*z-a*x)/b;
% % surf(z, y, x, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.3'); % Plot the surface
p3=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_3=[a b c d];
%---------------------------rod--------------------------------------------
r_muhit=2.5; %h=6;
l=h-r_muhit; % l 3.5 dan boshlab o'zgaradi
dd=h-2*r_muhit;
% if (zz2(1)-h)<1
%         eff=0;
%     continue
% end
% l=4.5; %Led dan muhit o'rtasigacha masofa %3.5 dan zz2(1)-3.5
[xs,ys,zs]=cylinder(r_muhit,20); %surf(xs+l, 25*zs, ys, 'EdgeColor', 'None', 'FaceColor', 'm', 'FaceAlpha', '0.1')
% plotCircle3D_1([0, D+0.1, l],[0,1,0], r_muhit)
% plotCircle3D_1([0, 0-0.1, l],[0,1,0], r_muhit)
%-------------------------crystal------------------------------------------
A=dlmread('15_12_2021_BLUE_LED_ems.txt');
C=dlmread('Ce_Nd_YAG_abs.txt');
%--------------------------------------------------------------------------
% led_array_parabola
sa=0.9*pi/2;
teflon_kaytarish=0.9;  %blue uchun
led_kaytarish=0.75;    %blue uchun
LEDdan_qaytdi=0;
lost_at_led=0;
lost_at_teflon=0;
chizish=0;
used_f1=0;
nuqta_fotonlar_soni_max=100;
axis equal
sp=zeros(1,2500);
s=zeros(1,nuqta_fotonlar_soni_max);
E_total=0;
E_det1=0; %




% x_grid_number  = 49;        x_grid_size = 2*r_muhit/x_grid_number;
% z_grid_number  = 49;        z_grid_size = 2*r_muhit/z_grid_number;
% y_grid_number  = 499;       y_grid_size = D/y_grid_number;
% 
% E_xyz = zeros(x_grid_number+1,y_grid_number+1,z_grid_number+1);



for ti=1: nuqta_fotonlar_soni_max
    tasodifiy_son=rand;
    for i=1:2500
        if tasodifiy_son<=A(i)
            ftu=i;  % ftu - fotoning to'lqin uzunligi deganini bildiradi
            break;
        end
    end
    sp(ftu)=sp(ftu)+1;
    s(1,ti)=ftu;
end


for f=1:nuqta_fotonlar_soni_max
    ftu=s(1,f);
    E_total=E_total+(1/ftu);
    %-------------------------------------absorbtion---------------------------
    abs_coef=C(ftu); abs_length=absorption_length(abs_coef); abs_leng_m=abs_length;
    if (abs_leng_m>1900);  yutilmaydigan_f=yutilmaydigan_f+1; continue; end
    %--------------------------------------------------------------------------
    %------------- boshlang'ich kordinata va yonalish--------------
    
    r_i=led_random_square_1(0, 12.5, 0, 12.5, 12.5);
   
    yonalish = lambert_shape_4 (p1, tan11, tan12, sa);
    
    %         r_f = r_i + yonalish*10;
    %         if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','b','LineWidth',0.5); end
    %--------------------------------------------------------------
    
    foton_in_system=1;
    while (foton_in_system == 1)
        %---------------------LEDs--------------------------------
        dm(1)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish, tek_1);  tk1(1,:)=[p1(1) p1(2) p1(3)];    %past_1 Led turgan joy
        %-------------------------qopqoq-------------------------------
        dm(2)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish, tek_2);  tk1(2,:)=-[p2(1) p2(2) p2(3)];
        %--------------------------o'rta---------------------------------
        dm(3)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish, tek_3);  tk1(3,:)=-[p3(1) p3(2) p3(3)];
        %---------------------Laser rod----------------------------
        dm(4)=nuqta_slinder_masofa1(r_i, yonalish, r_muhit, 0, l); tk1(4,:)= [1/r_muhit 1/r_muhit 0];              %aktiv muhit;
        %---------------------parabolalar-------------------------
        dm(5)=nuqta_para_masofa_modified1(r_i,yonalish,[KK 0 0 zc -surish]);   tk1(5,:)= [0 0 0 ];   %tepa_1
        dm(6)=nuqta_para_masofa_modified2(r_i,yonalish,[KK 0 0 zc surish]);   tk1(6,:)= [0 0 0];  % yon_chap_1
        %----------------------------------------------------------
        for i=1:length(dm); if dm(i)<0.001; dm(i)=Inf; end; end
        %                         2
        [db,joy]=min(dm);  normal=tk1(joy,:); r_f=r_i + yonalish*db;
        
        if (db==Inf)||(db==0); foton_in_system=0; break; end
        
%         if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','b','LineWidth',0.5); end
        %                          pause
        
        r_i=r_f;
        
        %                             joy
        if (joy==1)
            %                            11
            %                                r_i_n=r_f; r_f_n=r_i_n + normal*1;
            %                                  if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',1.5); end
            if led_kaytarish > rand
                %                                                                                                  112
                yonalish=yonalish-2*normal*dot(yonalish,normal);
                LEDdan_qaytdi=LEDdan_qaytdi+1;
            else
                %                                                                                                  113
                foton_in_system=0; lost_at_led=lost_at_led+1;
            end
        end
        
        if (joy>=2)&&(joy<=3)
            %                             12
            %                          dfghgnxcvbnm
            %                          r_f
            %                                r_i_n=r_f; r_f_n=r_i_n + normal*1;
            %                             if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',1.5); end
            %                                 joy
            %                             normal;
            
            %                                                                 121
            %
            if teflon_kaytarish > rand
                %                                                                        122
                yonalish=yonalish-2*normal*dot(yonalish,normal);
            else
                %                                                                                                              123
                foton_in_system=0; lost_at_teflon=lost_at_teflon+1;
            end
            %                                 sdfgbhnjmk
            
        end
        if (joy==4)
            
            %                                                                                      14
            nor=[(r_f(1))/r_muhit 0 (r_f(3)-l)/r_muhit ];
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
                    if (r_f(2)>D)||(r_f(2)<0)                %esdan chiqmasin
                        %
                        %                                                                                                                         143
                        foton_in_system=0; foton_muhitda=0;
                        break
                    end
                    
                    if (r_i(2)>D)||(r_i(2)<0)                %esdan chiqmasin
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
                            used_f1=used_f1+1;
                            %                                              13111
                            E_det1=E_det1+(1/ftu);
                            foton_muhitda=0; foton_in_system=0;
                            
%                             x_index = round((r_f(1)+ r_muhit)/x_grid_size)+1;
%                             z_index = round((r_f(3)-dd)/z_grid_size)+1;
%                             y_index = round((r_f(2))/y_grid_size)+1;
%                             E_xyz(x_index,y_index,z_index)=E_xyz(x_index,y_index,z_index)+1;

                            if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','r','LineWidth',0.1); end
                            % %                                                 if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','r','LineWidth',0.1); gs=gs+1; end
                        end
                    else
                        %                                                                                                                         14012
                        abs_leng_m=abs_leng_m-mm;
                        r_f=r_i+mm*yonalish;
                        if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','r','LineWidth',0.1); end
                        %  if chizish==1; is(gs)=plot3([r_i(1),r_f(1)],[r_i(2),r_f(2)],[r_i(3),r_f(3)],'Color','r','LineWidth',0.1); gs=gs+1; end
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
                    %                                     if r_i(3) > D/2
                    % %                                         led_number=2;
                    % %                                         r_i=r_f;
                    % %                                          192
                    %
                    %                                         foton_muhitda=0; foton_in_system=0;
                    %                                         break
                    %
                    %                                     end
                    
                    if (r_f(2)>D)||(r_f(2)<0)                %esdan chiqmasin
                        foton_muhitda=0;
                        %                                       193
                        foton_in_system=0;
                        
                        break
                    end
                    
                    if (r_i(2)>D)||(r_i(2)<0)                %esdan chiqmasin
                        %                                                                      144
                        
                        foton_muhitda=0; foton_in_system=0;
                        break
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
                foton_in_system=0; lost_at_teflon=lost_at_teflon+1;
            end
        end
        if (joy==6)
            %                                                                          16
            %                                                                          wsefrtgyhuj
            normal=-normalni_top(r_f, para_coefs_2);
            %                             r_i_n=r_f; r_f_n=r_i_n + normal*1;
            %                             if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',1.5); end
            if teflon_kaytarish > rand
                %                                                                        161
                yonalish=yonalish-2*normal*dot(yonalish,normal);
            else
                %                                                                                                              162
                foton_in_system=0; lost_at_teflon=lost_at_teflon+1;
            end
        end
        
        
    end
end
% end
% used_f1/(nuqta_fotonlar_soni_max)
eff=E_det1/E_total;
% hold on
% view(26, 25)