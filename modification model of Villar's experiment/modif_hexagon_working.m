clear; clc; close all;
%  hh=linspace(6, 100,20); 
% alfa1=linspace(0, (pi/2)*0.8, 20); 
%   for i=1:20  
xlabel('x')
ylabel('y')
zlabel('z')
D=50;
L=25;

%  alfa=alfa1(i);
    % o dan ~90 gradus
alfa=32*pi/180; 
% h=hh(i); %in mm
h=10;
dd=5;
l=6; if l<=dd; xato; end
r_muhit=2.5; %in mm
H=h+l;
ly=l*tan(alfa);
LED_1=H-2*(l-r_muhit);
LED_2=H;
lyy=((L/2+r_muhit)*(LED_2-LED_1-l))/h;
hold on;
% axis off
%%------------------------------------------------rod----------------------
[xs, ys, zs]=cylinder( r_muhit,40);  surf(D/2*zs, h+r_muhit+xs, ly+L/2+ys, 'EdgeColor', 'None', 'FaceColor', 'm', 'FaceAlpha', '0.3')
plotCircle3D_5([h+r_muhit, ly+L/2, D/2+0.1],[0,0,1], r_muhit)
plotCircle3D_5([h+r_muhit, ly+L/2, 0-0.1],[0,0,1], r_muhit)
%%----------------------------plane_1---------------------------------------
% LED_1=H-2*(l-r_muhit);
%%------------------------past_1 yani LED 1 turgan joy---------------------
a=1; b=0; c=0; d=LED_1;
yy=linspace(ly, ly+L, 50); zz=linspace(0, D/2, 50);
[y, z] = meshgrid(yy,zz); % Generate x and y data
x=(d-b*y-c*z)/a;
tashqi_plot(1)=surf(z, x, y, 'EdgeColor', 'None', 'FaceColor', 'w', 'FaceAlpha', '1'); % Plot the surface
p1=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tan11=[0 0 -1]; tan12=[0 -1 0];
tek_1=[a b c d];
%%----------------------------------LED_1 rasmi------------------------------
a=1; b=0; c=0; d=LED_1-3;
yy=linspace(ly-2,ly+L+2,2); zz=linspace(0-2, D/2+2,2);
[y, z] = meshgrid(yy,zz); % Generate x and y data
x=(d-b*y-c*z)/a;
surf(z, x, y, 'EdgeColor', 'None', 'FaceColor', 'c', 'FaceAlpha', '0.1'); % Plot the surface
surf(z, x+2.9, y, 'EdgeColor', 'None', 'FaceColor', 'b', 'FaceAlpha', '0.1'); % Plot the surface
%%-------------------------------------------------------------------------
r_k=20.5;
theta=45;
[xk, yk, zk]=cylinder(r_k, 4);
X2 = xk*cosd(theta) - yk*sind(theta);
Y2 = xk*sind(theta) + yk*cosd(theta);
% surf( X2+ly+12.5, zk*3+LED_1-3, Y2+12.5,  'EdgeColor', 'k', 'FaceColor',  'c', 'FaceAlpha', '0.1')%[0.153 0.255 0.204]%[0 0.4470 0.7410]
surf( X2+ly+10, zk*3+LED_1-3, Y2+14.7,  'EdgeColor', 'k', 'FaceColor',  'c', 'FaceAlpha', '0.1')%[0.153 0.255 0.204]%[0 0.4470 0.7410]
%%----------------------tepa_1 yani kichik asos----------------------------
a=1; b=0; c=0; d=LED_1+H;
yy=linspace(ly+L/2-dd/2,ly+L/2+dd/2,2); zz=linspace(0,D/2,2);
[y, z] = meshgrid(yy,zz); % Generate x and y data
x=(d-b*y-c*z)/a;
tashqi_plot(2)=surf(z, x, y,  'EdgeColor', 'k', 'FaceColor', [0.72 0.80 0.99], 'FaceAlpha', '0.2', 'LineStyle', '--'); % Plot the surface
p2=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_2=[a b c d];
%%------------------------ yon_o'ng_1a---------------------------------------
k=-tan(alfa); 
bd=-k*(h+2*r_muhit);
a=-k; b=1; c=0; d=bd;
xx=linspace(LED_1,LED_1+l,2); zz=linspace(0,D/2,2);
[x, z] = meshgrid(xx,zz); 
y=(d-a*x-c*z)/b;
tashqi_plot(3)=surf(z, x, y,  'EdgeColor', 'k', 'FaceColor', [0.52 0.61 0.82], 'FaceAlpha', '0.2', 'LineStyle', ':'); % Plot the surface
p3=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_3=[a b c d];
%%------------------------ yon_chap_1a---------------------------------------
k=-tan(alfa);
bd=L+2*ly+k*(h+2*r_muhit);
a=k; b=1; c=0; d=bd;
xx=linspace(LED_1, LED_1+l,2); zz=linspace(0,D/2,2);
[x, z] = meshgrid(xx,zz);
y=(d-a*x-c*z)/b;
tashqi_plot(4)=surf(z, x, y,  'EdgeColor', 'k', 'FaceColor', [0.52 0.61 0.82], 'FaceAlpha', '0.2', 'LineStyle', ':'); % Plot the surface
p4=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_4=[a b c d];
%%------------------------ yon_o'ng_1b---------------------------------------
k=(ly+(L-dd)/2)/h;
bd=-k*(h+2*r_muhit);
a=-k; b=1; c=0; d=bd;
xx=linspace(LED_1+l,LED_1+H,50); zz=linspace(0,D/2,50);
[x, z] = meshgrid(xx,zz); 
y=(d-a*x-c*z)/b;
tashqi_plot(6)=surf(z, x, y,  'EdgeColor', 'None', 'FaceColor', [0.39 0.47 0.58], 'FaceAlpha', '0.2', 'LineStyle', '--'); % Plot the surface
p5=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_5=[a b c d];
%%------------------------ yon_chap_1b---------------------------------------
k=-((dd-L)/2-ly)/h;
% bd=(L+dd)/2+2*ly+k*(LED_1+H);
bd=2*ly+L+k*(h+2*r_muhit);
a=k; b=1; c=0; d=bd;
xx=linspace(LED_1+l,LED_1+H,100); zz=linspace(0,D/2,50);
[x, z] = meshgrid(xx,zz);
y=(d-a*x-c*z)/b;
tashqi_plot(7)=surf(z, x, y, 'EdgeColor', 'None', 'FaceColor', [0.39 0.47 0.58], 'FaceAlpha', '0.2', 'LineStyle', '--'); % Plot the surface
p6=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_6=[a b c d];
%%-----------------------------qopqoq_1------------------------------------------
a=0; b=0; c=-1; d=0;
xx=linspace(LED_1,LED_1+H,50); yy=linspace(0,2*ly+L,50);
[x, y] = meshgrid(xx,yy); % Generate x and y data
z=(d-a*x-b*y)/c;
% surf(x, y, z, 'EdgeColor', 'None', 'FaceColor', 'r', 'FaceAlpha', '0.1') % Plot the surface
 tashqi_plot(8)=patch([0 0 0 0 0 0],[LED_1 LED_1 LED_1+l LED_1+H LED_1+H LED_1+l], [ly L+ly L+2*ly  L/2+ly+dd/2 L/2+ly-dd/2 0], 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
 p7=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_7=[a b c d];
 tashqi_plot(9)=patch([D/2 D/2 D/2 D/2 D/2 D/2],[LED_1 LED_1 LED_1+l LED_1+H LED_1+H LED_1+l], [ly L+ly L+2*ly  L/2+ly+dd/2 L/2+ly-dd/2 0], 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
 %%------------------------mavhum_o'rta_devor-------------------------------------
% LED_2=H;
x_medium=(LED_1+LED_2)/2;
y_m_1=-tan(alfa)*(h-x_medium);
y_m_2=tan(alfa)*(H-x_medium)+L+ly;
a=0; b=0; c=1; d=D/2;
xx=linspace(LED_1,LED_2,50); yy=linspace(y_m_1, y_m_2,50);
[x, y] = meshgrid(xx,yy); % Generate x and y data
z=(d-a*x-b*y)/c;
%  surf(x, y, z, 'EdgeColor', 'None', 'FaceColor', 'r', 'FaceAlpha', '0.1') % Plot the surface
% tashqi_plot(9)=patch([D/2 D/2 D/2 D/2 D/2 D/2],[LED_1 LED_1 x_medium LED_2 LED_2 x_medium], [ly L+ly y_m_2 L+ly ly y_m_1],  'w', 'FaceAlpha', '0.001');
p8=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_8=[a b c d];
%%--------------------------------crystal----------------------------------
A=dlmread('15_12_2021_BLUE_LED_ems.txt');
C=dlmread('Ce_Nd_YAG_abs.txt');

%%------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
sa=0.9*pi/2;
teflon_kaytarish=0.9;  %blue uchun
led_kaytarish=0.75;    %blue uchun
%%--------------sanash_uchun----------------------------------------------------
LEDdan_qaytdi=0;
lost_at_led=0;
loss_at_teflon=0;
lost_orta_devor=0;
used_f1=0;
used_f2=0;
foton_muhitga_bordi1=0;
foton_muhitga_bordi2=0;
yutilmaydigan_f=0;
%%---------------------LED_ni chaqirish------------------------------------
nuqta_fotonlar_soni_max=100;
axis equal

view(90,0)
sp=zeros(1,2500);
s=zeros(1,nuqta_fotonlar_soni_max);
E_total=0;
E_det1=0;
E_det2=0;
E_det_muhit1=0;
E_det_muhit2=0;
chizish=0;

foton_soni=0;

% bjn
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

DD=D/2;
x_grid_number  = 49;        x_grid_size = 2*r_muhit/x_grid_number;
y_grid_number  = 49;        y_grid_size = 2*r_muhit/y_grid_number;
z_grid_number  = 199;       z_grid_size = DD/z_grid_number;

E_xyz = zeros(x_grid_number+1, y_grid_number+1, z_grid_number+1);

    
    for f=1:nuqta_fotonlar_soni_max
        if rand > 0.999995; clc; disp("Algaritim bajarilish prosesi " + round(100 * f/nuqta_fotonlar_soni_max) + " %"); pause(0.1);  end  % bi xator loop nera galganini go'rsatadi
        foton_soni=foton_soni+1;
        ftu=s(1,f);
        E_total=E_total+(1/ftu);
        %-------------------------------------absorbtion---------------------------
        abs_coef=C(ftu); abs_length=absorption_length(abs_coef); abs_leng_m=abs_length;
        if (abs_leng_m>1900);  yutilmaydigan_f=yutilmaydigan_f+1; continue; end
        %---------------------------------------------------------------------------
        
        %------------- boshlang'ich kordinata va yonalish--------------

          r_i=led_random_square_2(LED_1, 12.5, 12.5, 12.5, 12.5, ly);
          yonalish = lambert_shape_4 (p1, tan11, tan12, sa);
%           r_f = r_i + yonalish*10;
%           if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','b','LineWidth',0.5); end
        
        %--------------------------------------------------------------
        foton_in_system=1;
        while (foton_in_system == 1)
            %%---------------------------------------case 1 --------------
            dm(1)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_1);  tk1(1,:)=[p1(1) p1(2) p1(3)];       %past_1 Led turgan joy
            %---------------------Cavity_walls-------------------------
            dm(2)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_2);  tk1(2,:)= -[p2(1) p2(2) p2(3)];     %tepa_1
            dm(3)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_3);  tk1(3,:)= [p3(1) p3(2) p3(3)];      %yon_o'ng_1a
            dm(4)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_4);  tk1(4,:)=-[p4(1) p4(2) p4(3)];      %yon_chap_1a
            dm(5)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_5);  tk1(5,:)= [p5(1) p5(2) p5(3)];      %yon_o'ng_1b
            dm(6)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_6);  tk1(6,:)= -[p6(1) p6(2) p6(3)];     %yon_chap_1b
            dm(7)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_7);  tk1(7,:)= -[p7(1) p7(2) p7(3)];     %qopqoq_1
            %-----------------------O'rtakash--------------------------
            dm(8)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_8);  tk1(8,:)= -[p8(1) p8(2) p8(3)];     %mavhum o'rta devor
            %---------------------Laser rod----------------------------
            dm(9)=nuqta_slinder_masofa(r_i, yonalish, r_muhit, h+r_muhit, L/2+ly); tk1(9,:)= [1/r_muhit 1/r_muhit 0];  %aktiv muhit;
            %----------------------------------------------------------
            for i=1:length(dm); if dm(i)<0.00001; dm(i)=Inf; end; end
            [db,joy]=min(dm); normal=tk1(joy,:);
            r_f=r_i + yonalish*db;
            
            if (db==Inf)||(db==0); foton_in_system=0; break; end
            %                          if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','b','LineWidth',0.1); gs=gs+1; end
%                                      if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','b','LineWidth',0.5); end
            r_i=r_f;
            if (joy==1)
                %                            11
                %                                  r_i_n=r_f; r_f_n=r_i_n + normal*1;
                %                                  if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','c','LineWidth',0.8); end
                if led_kaytarish > rand
                    %                                                                                                  112
                    yonalish=yonalish-2*normal*dot(yonalish,normal);
                    LEDdan_qaytdi=LEDdan_qaytdi+1;
                else
                    %                                                                                                  113
                    foton_in_system=0; lost_at_led=lost_at_led+1;
                end
            end
            
            
            if (joy>=2)&&(joy<=8)
                
                normal;
                
                if teflon_kaytarish > rand
                    %                                                                        122
                    yonalish=yonalish-2*normal*dot(yonalish,normal);
                else
                    %                                                                                                              123
                    foton_in_system=0; loss_at_teflon=loss_at_teflon+1;
                end
                %                             end
            end
            
            if (joy==9)
                %                                                                                      19
                foton_muhitga_bordi1=foton_muhitga_bordi1+1;
                E_det_muhit1=E_det_muhit1+(1/ftu);
                nor=[(r_f(1)-(h+r_muhit))/r_muhit (r_f(2)-(L/2+ly))/r_muhit 0];
                if (rand()>(qaytish_sinish_ehtimolliqi(yonalish,nor,1,1.8)))
                    %                                                                                            190
                    s_n=nor;
                    %                                 r_i_n=r_f; r_f_n=r_i_n + s_n*1;
                    %                                 if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','c','LineWidth',0.8); end
                    foton_muhitda=1;
                    [~,yonalish]=cheg_cos(yonalish, -nor,1,1.8);
                    shs=0;
                    while foton_muhitda==1
                        %                                                                                                             1901
                        mm=nuqta_slinder_masofa(r_i, yonalish, r_muhit, h+r_muhit, L/2+ly);
                        if (mm==Inf)||(mm==0); foton_in_system=0; foton_muhitda=0; break; end
                        if (r_f(3)>D/2)||(r_f(3)<0)                %esdan chiqmasin
                            %                                                                                                                         193
                            foton_in_system=0;
                            foton_muhitda=0;
                            break
                        end
                        
                        if (r_i(3)>D/2)||(r_i(3)<0)                %esdan chiqmasin
                            %                                                                                                                          194
                            foton_in_system=0; foton_muhitda=0;
                            break
                        end
                        if mm>abs_leng_m
                            mm=abs_leng_m;
                            r_f=r_i+mm*yonalish;
                            if (r_f(3)<0)||(r_f(3)>D/2)
                                foton_muhitda=0; foton_in_system=0;
                            else
                            %                                         scatter3(r_f(3),r_f(1),r_f(2), 15,'fill')
                            
                            x_index = round((r_f(1)-h)/x_grid_size)+1;
                            z_index = round((r_f(3))/z_grid_size)+1;
                            y_index = round((r_f(2)-(ly+L/2-r_muhit))/y_grid_size)+1;
                            E_xyz(x_index,y_index,z_index)=E_xyz(x_index,y_index,z_index)+1;
                            %
                            %                                                                                                                         13111
                            foton_muhitda=0; foton_in_system=0; used_f1=used_f1+1; E_det1=E_det1+(1/ftu);
%                                                                           if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',0.1); end
                            %                                               if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',0.1); gs=gs+1; end
                            end
                        else
                            %                                                                                                                         19012
                            abs_leng_m=abs_leng_m-mm;
                            r_f=r_i+mm*yonalish;
%                                                                   if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',0.1); end
                            %                                       if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',0.1); gs=gs+1; end
                            nor=-[(r_f(1)-(h+r_muhit))/r_muhit (r_f(2)-(L/2+ly))/r_muhit 0];
                            %                                         r_i_n=r_f; r_f_n=r_i_n + nor*1;
                            %                                         if chizish==1; plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','y','LineWidth',1); end
                            if  ((acos(dot(-yonalish,nor)/(norm(yonalish)*norm(nor)))) < asin(1/1.8))&&(rand()>(qaytish_sinish_ehtimolliqi(yonalish,nor,1.8,1)))
                                [~,yonalish]=cheg_cos(yonalish,-nor,1.8,1);
                                %                                             r_i_n=r_f; r_f_n=r_i_n + yonalish*2;
                                foton_muhitda=0;
                                %                                               19013
                                %                                             if chizish==1; plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','r','LineWidth',1.5);end
                                
                            else
                                %                                                                                                                                     19014
                                yonalish=yonalish-2*nor*dot(yonalish,nor);
                            end
                            r_i=r_f;
                        end
                        
                        if (r_f(3)>D/2)||(r_f(3)<0)                %esdan chiqmasin
                            %                                                                                                                         193
                            foton_in_system=0;
                            foton_muhitda=0;
                            break
                        end
                        
                        if (r_i(3)>D/2)||(r_i(3)<0)                %esdan chiqmasin
                            %                                                                                                                          194
                            foton_in_system=0; foton_muhitda=0;
                            break
                        end
                        shs=shs+1; if shs>5; foton_muhitda=0; foton_in_system=0; break; end
                    end
                else
                    %                                                                                                 191
                    yonalish=yonalish-2*nor*dot(yonalish,nor);
                end
            end
            %                     end
            %
            %                 end
            
            
            if (r_f(3)>D/2)||(r_f(3)<0)                %esdan chiqmasin
                %                                                                                                                         193
                foton_in_system=0;
                foton_muhitda=0;
                break
            end
            
            if (r_i(3)>D/2)||(r_i(3)<0)                %esdan chiqmasin
                %                                                                                                                          194
                foton_in_system=0; foton_muhitda=0;
                break
            end
        end
        
        
    end
    %               if gs>8; continue; end
    %             pause(3)
% end
%         delete(is)
%         pause(1)
% end
% E_det1
% E_det2
E_det=E_det1+E_det2;
E_det_muhit=E_det_muhit1+E_det_muhit2;
eff=E_det1/E_total
E_det_muhit/E_total;
% view (-120,60)
% view (-205, 17)
axis equal
grid off
% axis off
zlim([-3 4*L])
ylim([-3 2*H])
xlim([-1 D])
% pause(3)
% close all
%  end





