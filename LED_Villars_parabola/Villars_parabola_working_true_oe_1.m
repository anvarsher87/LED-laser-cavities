% function eff = effektivlik_topish_double_parabola(p, h)
clear; clc; close all;
%Villars_parabola_2 dan olindi
% all in cm
xlabel('x')
ylabel('y')
zlabel('z')
% % view(40, 30)
% view(0, 0)
p=3.2;
KK=1/(2*p);
hold on
D=25;

aa=9.5; %parabola shohlari extremumlar
surish=12.5-aa; %surish uchun
% KK=0.16; %parabola parametri %0.12 dan 0.3 gacha
zc1=KK*aa^2; % surish uchun
%---------------------parabola _1------------------------------------------
xx11=linspace(-aa-surish,0,20); xp(1,:)=xx11; xp(2,:)=xx11;
zz11=-KK*(xx11+surish).^2+zc1; zp(1,:)=zz11; zp(2,:)=zz11;
LL=max(zz11);
mm=2; n=19;
y=(0:mm-1)'/(mm-1) * ones(1,n+1);
para_coefs_11=zeros(4,4); para_coefs_11(1,1)=KK; para_coefs_11(1,4)=-KK*surish; para_coefs_11(4,1)=-KK*surish; para_coefs_11(4,4)=-zc1+KK*surish^2; para_coefs_11(3,4)=1/2; para_coefs_11(4,3)=1/2;
% surf(zp,D*y,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.1')
surf(zp,(D-6)*y+3,xp,'EdgeColor', 'None', 'FaceColor', 'b', 'FaceAlpha', '0.2')
surf(zp+0.01,3*y,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.4')
surf(zp+0.01,3*y+D-3,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.4')
% plot3(zp(20),[0,D],0,'Color','k','LineWidth',1);
%-------------------parabola_2---------------------------------------------
xx12=linspace(0,aa+surish,20); xp(1,:)=xx12; xp(2,:)=xx12;
zz12=-KK*(xx12-surish).^2+zc1; zp(1,:)=zz12; zp(2,:)=zz12;
mm=2; n=19;
% y=(0:mm-1)'/(mm-1) * ones(1,n+1);
para_coefs_12=zeros(4,4); para_coefs_12(1,1)=KK; para_coefs_12(1,4)=KK*surish; para_coefs_12(4,1)=KK*surish; para_coefs_12(4,4)=-zc1+KK*surish^2; para_coefs_12(3,4)=1/2; para_coefs_12(4,3)=1/2; %para_coefs(4,4)=1/2;
% surf(zp,D*y,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.1')
surf(zp,(D-6)*y+3,xp,'EdgeColor', 'None', 'FaceColor', 'b', 'FaceAlpha', '0.3')
surf(zp+0.001,3*y,xp,'EdgeColor', 'None', 'FaceColor', 'k	', 'FaceAlpha', '0.4')
surf(zp+0.001,3*y+D-3,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.4')
%----------------------------rasm yon devor (qopqoq va o'rtakash)----------
x=[xx11 xx12];
z=[zz11 zz12];
y=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% % patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
% patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineWidth', 1.5);
y1=y+2; y1_1=y+3;
y2=y+D-2; y1_2=y+D-3;
y3=y+D;
% y1=[3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];
% patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.6','LineWidth',1.5);
patch(z, y3,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.6','LineWidth',1.5);
% plot3(z, y1,  x,  'k', 'LineWidth',1.5, 'LineStyle', '--');
plot3(z, y1_2-0.1,  x,  'k', 'LineWidth',1.5, 'LineStyle', '-');
plot3(z, y1_1-0.1,  x,  'k', 'LineWidth',1.5,'LineStyle', '-');
% plot3(z, y2-0.1,  x,  'k', 'LineWidth',1.5, 'LineStyle', '--');
%----------------------------------------------------------------------------
x=[xx11 xx12];
z=[zz11 zz12];
y=[D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D]-0.01;
% patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.6','LineWidth', 1.5);
%----------------------------------------------------------------------------
% %--------------------------orta_chiziq ------------------------------------
% Define the constants for x, y, and D
x = 0; % constant x coordinate
y = 12.6953; % constant y coordinate
% D = 10; % maximum z coordinate

% Create z coordinates from 0 to D
z = linspace(0, D, 100); % 100 points from 0 to D

% Plot the 3D line
plot3(y * ones(size(z)), z, x * ones(size(z)), 'k', 'LineWidth', 1.5);
%----------------------LED_1-------------------------------------------------
a=0; b=0; c=1; d=0;
yy=linspace(0,D,2); xx=linspace(-12.5,12.5,2);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z=(d-a*x-b*y)/c;
surf(z, y, x, 'EdgeColor', 'k', 'FaceColor', 'w', 'FaceAlpha', '1', 'LineWidth', 1.5); % Plot the surface
p11=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tan11=[1 0 0]; tan12=[0 1 0];
tek_11=[a b c d];
%----------------------------- LED rasmi -------------------------------------------------------
yy=linspace(-2,D+2,50); xx=linspace(-12.5-2, 12.5+2,50);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
surf(z+0.01, y, x,  'EdgeColor', 'None', 'FaceColor', 'c', 'FaceAlpha', '0.1') % Plot the surface
% surf(z+D/2-2,  y+h+LED_rod_masofa, x,'EdgeColor', 'None', 'FaceColor', 'w', 'FaceAlpha', '1') % Plot the surface
%------------------------------LED rasmi (orqa devor)--------------------------------------------
yy=linspace(-2,D+2,2); xx=linspace(-12.5-2,12.5+2,2);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
surf(z-3, y, x,  'EdgeColor', 'None', 'FaceColor', 'c', 'FaceAlpha', '0.4') % Plot the surface
%------------------------------LED rasmi (yon devor)-----------------------
r_k=20.5; theta=45; %for rotation plane
[xk, yk, zk]=cylinder(r_k, 4);
X2 = xk*cosd(theta) - yk*sind(theta);
Y2 = xk*sind(theta) + yk*cosd(theta);
surf(zk*3-3, Y2+12.5, X2,'EdgeColor', 'k', 'FaceColor',  'c', 'FaceAlpha', '0.1', 'LineWidth', 1.5)%[0.1 0.155 0.198]
%--------------------------------------------------------------------------
%--------------------qopqoq_1----------------------------------------------
a=0; b=-1; c=0; d=0;
% xx=linspace(-12.5,12.5,50); zz=linspace(0,zc,50);
% [x, z] = meshgrid(xx,zz); % Generate x and y data
% y=(d-c*z-a*x)/b;
% % surf(z, y, x, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.3'); % Plot the surface
% p2=[a b c d];
p12=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_12=[a b c d];
% xdfcgvbh
%--------------------mavhum devor------------------------------------------
a=0; b=1; c=0; d=D;
% xx=linspace(-12.5,12.5,50); zz=linspace(0,zc,50);
% [x, z] = meshgrid(xx,zz); % Generate x and y data
% y=(d-c*z-a*x)/b;
% % surf(z, y, x, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.3'); % Plot the surface
p13=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_13=[a b c d];
%---------------------------rod--------------------------------------------
r_muhit=2.5; h=11;
l=h-r_muhit; % l 3.5 dan boshlab o'zgaradi
dd=h-2*r_muhit;
% if (zz2(1)-h)<1
%         eff=0;
%     continue
% end
% l=4.5; %Led dan muhit o'rtasigacha masofa %3.5 dan zz2(1)-3.5
[xs,ys,zs]=cylinder(r_muhit,20); surf(xs+l, 2*D*1.2*zs-5, ys, 'EdgeColor', 'None', 'FaceColor', 'm', 'FaceAlpha', '0.2')
plotCircle3D_1([0, 2*D+5+0.1, l],[0,1,0], r_muhit)
plotCircle3D_1([0, 0-5-0.1, l],[0,1,0], r_muhit)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% %--------------------------orta_chiziq ------------------------------------
% Define the constants for x, y, and D
x = 0; % constant x coordinate
y = 4.33; % constant y coordinate
% D = 10; % maximum z coordinate

% Create z coordinates from 0 to D
z = linspace(D, 2*D, 100); % 100 points from 0 to D

% Plot the 3D line
plot3(y * ones(size(z)), z, x * ones(size(z)), 'k', 'LineWidth', 1.5);
%----------------------LED_2-----------------------------------------------
L=50; %D ning davomi
LED_z=h+dd;
a=0; b=0; c=1; d=LED_z;
yy=linspace(D,L,2); xx=linspace(-12.5,12.5,2);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z=(d-a*x-b*y)/c;
surf(z, y, x, 'EdgeColor', 'k', 'FaceColor', 'w', 'FaceAlpha', '1', 'LineWidth', 1.5); % Plot the surface
p21=-[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tan21=[0 1 0]; tan22=[1 0 0];
tek_21=[a b c d];
%----------------------------- LED rasmi ----------------------------------
yy=linspace(D-2,L+2,50); xx=linspace(-12.5-2, 12.5+2,50);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
surf(z+LED_z-0.01, y, x,  'EdgeColor', 'None', 'FaceColor', 'c', 'FaceAlpha', '0.1') % Plot the surface
% surf(z+D/2-2,  y+h+LED_rod_masofa, x,'EdgeColor', 'None', 'FaceColor', 'w', 'FaceAlpha', '1') % Plot the surface
%------------------------------LED rasmi (orqa devor)----------------------
yy=linspace(D-2,L+2,2); xx=linspace(-12.5-2,12.5+2,2);
[y, x] = meshgrid(yy,xx); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
surf(z+LED_z+3, y, x,  'EdgeColor', 'None', 'FaceColor', 'c', 'FaceAlpha', '0.4') % Plot the surface
%------------------------------LED rasmi (yon devor)---------------------------------------------
r_k=20.5; theta=45; %for rotation plane
[xk, yk, zk]=cylinder(r_k, 4);
X2 = xk*cosd(theta) - yk*sind(theta);
Y2 = xk*sind(theta) + yk*cosd(theta);
surf(zk*3+LED_z+0.01, Y2+D+12.5, X2,'EdgeColor', 'k', 'FaceColor',  'c', 'FaceAlpha', '0.1', 'LineWidth', 1.5)%[0.1 0.155 0.198]
%-----------------------------------------------------------------------------------
%----------------------------CASE_2-------------------------------------------------
%-----------------------------------------------------------------------------------
aa=9.5; %parabola shohlari extremumlar
surish=12.5-aa; %surish uchun
% KK=0.16; %parabola parametri %0.12 dan 0.3 gacha
dif_LL=2*(LL-h+r_muhit);
zc2=KK*aa^2-dif_LL; % surish uchun
%---------------------parabola _1-LED_2  uchun--------------------------------------
xx21=linspace(-aa-surish,0,20); xp(1,:)=xx21; xp(2,:)=xx21;
zz21=KK*(xx11+surish).^2+zc2 ;
zp(1,:)=zz21; zp(2,:)=zz21;
mm=2; n=19;
y=(0:mm-1)'/(mm-1) * ones(1,n+1);
para_coefs_21=zeros(4,4); para_coefs_21(1,1)=-KK; para_coefs_21(1,4)=KK*surish; para_coefs_21(4,1)=KK*surish; para_coefs_21(4,4)=zc2-KK*surish^2; para_coefs_21(3,4)=1/2; para_coefs_21(4,3)=1/2;
% surf(zp,D*y+D,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.1')
surf(zp,(D-3)*y+D,xp,'EdgeColor', 'None', 'FaceColor', 'b', 'FaceAlpha', '0.2', 'LineWidth',1.5)
surf(zp+0.01,(3)*y+2*D-3,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5)
surf(zp+0.01,(3)*y+D,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5)
%-------------------parabola_2---------------------------------------------
xx22=linspace(0,aa+surish,20); xp(1,:)=xx22; xp(2,:)=xx22;
zz22=KK*(xx22-surish).^2+zc2; 
zp(1,:)=zz22; zp(2,:)=zz22;
mm=2; n=19;
% y=(0:mm-1)'/(mm-1) * ones(1,n+1);
para_coefs_22=zeros(4,4); para_coefs_22(1,1)=-KK; para_coefs_22(1,4)=-KK*surish; para_coefs_22(4,1)=-KK*surish; para_coefs_22(4,4)=zc2-KK*surish^2; para_coefs_22(3,4)=1/2; para_coefs_22(4,3)=1/2; %para_coefs(4,4)=1/2;
% surf(zp,D*y+D,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.1')
surf(zp,(D-6)*y+D+3,xp,'EdgeColor', 'None', 'FaceColor', 'b', 'FaceAlpha', '0.2', 'LineWidth',1.5)
surf(zp+0.01,(3)*y+2*D-3,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5)
surf(zp+0.01,(3)*y+D,xp,'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5)

%----------------------------rasm yon devor (qopqoq va o'rtakash)----------
x=[xx21 xx22];
z=[zz21 zz22];

y=[D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D D];
% % patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
% patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineWidth', 1.5);

y11= y +2;
y2_1 = y +3;
% patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5);
% patch(z, y11,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05');
% plot3(z, y11,  x,  'k', 'LineStyle', '--','LineWidth',1.5);
plot3(z, y2_1,  x,  'k', 'LineStyle', '-', 'LineWidth',1.5);
% plot3(z, y2,  x,  'k', 'LineStyle', '--');


%----------------------------------------------------------------------------
x=[xx21 xx22];
z=[zz21 zz22];
y=[L L L L L L L L L L L L L L L L L L L L L L L L L L L L L L L L L L L L L L L L];
% % patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
% patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineWidth', 1.5);
ld=L-2;
y1=[ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld ld];
y2_2=y1-1;
% patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.05', 'LineStyle', '--');
patch(z, y,  x,  's', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5);
% plot3(z, y1,  x,  'k', 'LineStyle', '--', 'LineWidth',1.5);
plot3(z, y2_2,  x,  'k', 'LineStyle', '-', 'LineWidth',1.5);

%--------------------------------------------------------------------------
%--------------------qopqoq_2----------------------------------------------
a=0; b=1; c=0; d=L;
% xx=linspace(-12.5,12.5,50); zz=linspace(0,zc,50);
% [x, z] = meshgrid(xx,zz); % Generate x and y data
% y=(d-c*z-a*x)/b;
% % surf(z, y, x, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.3'); % Plot the surface
% p2=[a b c d];
p22=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_22=[a b c d];
% xdfcgvbh
%--------------------mavhum devor------------------------------------------
a=0; b=1; c=0; d=D;
% xx=linspace(-12.5,12.5,50); zz=linspace(0,zc,50);
% [x, z] = meshgrid(xx,zz); % Generate x and y data
% y=(d-c*z-a*x)/b;
% % surf(z, y, x, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.3'); % Plot the surface
p23=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_23=[a b c d];
%-------------------------crystal------------------------------------------
A=dlmread('15_12_2021_BLUE_LED_ems.txt');
C=dlmread('Ce_Nd_YAG_abs.txt');
%--------------------------------------------------------------------------
led_array_parabola_1
sa=pi/3;%0.9*pi/2;
teflon_kaytarish=0.9;  %blue uchun
led_kaytarish=0.75;    %blue uchun
LEDdan_qaytdi=0;
lost_at_led=0;
lost_at_teflon=0;
chizish=1;
used_f1=0;
nuqta_fotonlar_soni_max=100;
axis equal
hold on
sp=zeros(1,2500);
s=zeros(1,nuqta_fotonlar_soni_max);
E_total=0;
E_det1=0; %

% ctfvgybhnkml


x_grid_number  = 49;        x_grid_size = 2*r_muhit/x_grid_number;
z_grid_number  = 49;        z_grid_size = 2*r_muhit/z_grid_number;
y_grid_number  = 199;       y_grid_size = L/y_grid_number;

E_xyz = zeros(x_grid_number+1,z_grid_number+1,y_grid_number+1);

% dfcgvhb

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

% cfgbh
for led_number0=1:2
    %     led_point
    for f=1:nuqta_fotonlar_soni_max
        if rand > 0.999995; clc; disp("Algaritim bajarilish prosesi " + round(100 * f/nuqta_fotonlar_soni_max) + " %"); pause(0.1);  end  % bi xator loop nera galganini go'rsatadi
        ftu=s(1,f);
        E_total=E_total+(1/ftu);
        %-------------------------------------absorbtion-----------------------
        abs_coef=C(ftu); abs_length=absorption_length(abs_coef); abs_leng_m=abs_length;
        if (abs_leng_m>1900);  yutilmaydigan_f=yutilmaydigan_f+1; continue; end
        %----------------------------------------------------------------------
        %------------- boshlang'ich kordinata va yonalish----------------------
        if led_number0==1
            r_i=led_random_square_1(0, 12.5, 0, 12.5, 12.5);
%             scatter3 (r_i(3), r_i(2), r_i(1),  15,'fill')
            yonalish = lambert_shape_4 (p11, tan11, tan12, sa);
%                 r_f=r_i + yonalish*10;
%                 if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','b','LineWidth',0.5); end
            led_number=led_number0;
        end
        
        if led_number0==2
            r_i=led_random_square_1(0, D+12.5, h+dd, 12.5, 12.5);
%             scatter3 (r_i(3), r_i(2), r_i(1),  15,'fill')
            yonalish = lambert_shape_4 (p21, tan21, tan22, sa);
%                  r_f=r_i + yonalish*10;
%                 if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','b','LineWidth',0.5); end
            led_number=led_number0;
        end
        
        %------------------------------------------------------------------
        
        foton_in_system=1;
        while (foton_in_system == 1)
            %--------------------case_1------------------------------------
            if (led_number==1)&&(foton_in_system == 1)
                
                while (led_number==1)&&(foton_in_system ==1)
                    %---------------------LEDs-----------------------------
                    dm(1)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish, tek_11);  tk1(1,:)=[p11(1) p11(2) p11(3)];    %past_1 Led turgan joy
                    %-------------------------qopqoq-------------------------------
                    dm(2)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish, tek_12);  tk1(2,:)=-[p12(1) p12(2) p12(3)];
                    %--------------------------o'rta---------------------------------
                    dm(3)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish, tek_13);  tk1(3,:)=-[p13(1) p13(2) p13(3)];
                    %---------------------Laser rod----------------------------
                    dm(4)=nuqta_slinder_masofa1(r_i, yonalish, r_muhit, 0, l); tk1(4,:)= [1/r_muhit 1/r_muhit 0];              %aktiv muhit;
                    %---------------------parabolalar-------------------------
                    dm(5)=nuqta_para_masofa_modified1(r_i,yonalish,[KK 0 0 zc1 -surish]);   tk1(5,:)= [0 0 0 ];   %tepa_1
                    dm(6)=nuqta_para_masofa_modified2(r_i,yonalish,[KK 0 0 zc1 surish]);   tk1(6,:)= [0 0 0];  % yon_chap_1
                    %------------------------------------------------------
                    for i=1:length(dm); if dm(i)<0.001; dm(i)=Inf; end; end
                    
                    %                         2
                    [db,joy]=min(dm);  normal=tk1(joy,:); r_f=r_i + yonalish*db;
                    
                    if (db==Inf)||(db==0); foton_in_system=0; break; end
                    
                    if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','b','LineWidth',1); end
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
                                if (r_f(2)>D)||(r_f(2)<0)                %esdan chiqmasin
                                    %
                                    %                                                                                                                         143
                                    foton_in_system=0;
                                    foton_muhitda=0;
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
%                                                                                 scatter3(r_f(3),r_f(2),r_f(1), 15,'fill')
                                        %                                             14011
                                        used_f1=used_f1+1;
                                        %                                              13111
                                        E_det1=E_det1+(1/ftu);
                                        foton_muhitda=0; foton_in_system=0;
                                        
                                        x_index = round((r_f(1)+ r_muhit)/x_grid_size)+1;
                                        z_index = round((r_f(3)-dd)/z_grid_size)+1;
                                        y_index = round((r_f(2))/y_grid_size)+1;
                                        E_xyz(x_index,z_index,y_index)=E_xyz(x_index,z_index,y_index)+1;
                                        r_f;
                                        if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','r','LineWidth',1); end
                                        % %                                                 if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','r','LineWidth',0.1); gs=gs+1; end
                                    end
                                else
                                    %                                                                                                                         14012
                                    abs_leng_m=abs_leng_m-mm;
                                    r_f=r_i+mm*yonalish;
                                    if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','r','LineWidth',1); end
                                    % %                                         if chizish==1; is(gs)=plot3([r_i(1),r_f(1)],[r_i(2),r_f(2)],[r_i(3),r_f(3)],'Color','r','LineWidth',0.1); gs=gs+1; end
                                    nor=-[(r_f(1))/r_muhit 0 (r_f(3)-l)/r_muhit ];
                                    %                                    nor=-[(r_f(2)-h)/r_muhit  0 (r_f(1))/r_muhit];
                                                                            r_i_n=r_f; r_f_n=r_i_n + nor*1;
                                    %                                         if chizish==1; plot3([r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],[r_i_n(3),r_f_n(3)],'Color','y','LineWidth',1); end
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
                        normal=-normalni_top(r_f, para_coefs_11);
                                                    r_i_n=r_f; r_f_n=r_i_n + normal*1;
                                                    if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',0.5); end
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
                        normal=-normalni_top(r_f, para_coefs_12);
                                                    r_i_n=r_f; r_f_n=r_i_n + normal*1;
                                                    if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',0.5); end
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
            %--------------------------------case_2------------------------
            if (led_number==2) && (foton_in_system == 1)
                
                while (led_number==2)&&(foton_in_system ==1)
                    %---------------------LEDs-----------------------------
                    dm(1)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish, tek_21);  tk1(1,:)=[p21(1) p21(2) p21(3)];    %past_1 Led turgan joy
                    %-------------------------qopqoq-------------------------------
                    dm(2)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish, tek_22);  tk1(2,:)=-[p22(1) p22(2) p22(3)];
                    %--------------------------o'rta---------------------------------
                    dm(3)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish, tek_23);  tk1(3,:)=[p23(1) p23(2) p23(3)];
                    %---------------------Laser rod----------------------------
                    dm(4)=nuqta_slinder_masofa1(r_i, yonalish, r_muhit, 0, l); tk1(4,:)= [1/r_muhit 1/r_muhit 0];              %aktiv muhit;
                    %---------------------parabolalar-------------------------
                    dm(5)=nuqta_para_masofa_modified1(r_i,yonalish,[-KK 0 0 zc2 -surish]);   tk1(5,:)= [0 0 0 ];   %tepa_1
                    dm(6)=nuqta_para_masofa_modified2(r_i,yonalish,[-KK 0 0 zc2 surish]);   tk1(6,:)= [0 0 0];  % yon_chap_1
                    %------------------------------------------------------
                    for i=1:length(dm); if dm(i)<0.0001; dm(i)=Inf; end; end
%                     dm
                    %                         2
                    [db,joy]=min(dm);  normal=tk1(joy,:); r_f=r_i + yonalish*db;
                    
                    if (db==Inf)||(db==0); foton_in_system=0; break; end
                    
                    if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','b','LineWidth',1); end
                    %                          pause
                    
                    r_i=r_f;
                    
%                                                 joy
                    if (joy==1)
                        %                            11
                                                       r_i_n=r_f; r_f_n=r_i_n + normal*1;
                                                         if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',0.5); end
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
                                                       r_i_n=r_f; r_f_n=r_i_n + normal*1;
                                                    if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',0.5); end
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
                        %                                 nor=-[(r_f(2)-h)/r_muhit  0 (r_f(1))/r_muhit];
                        if (rand()>(qaytish_sinish_ehtimolliqi(yonalish,nor,1,1.8)))
                            %                                                                                            140
                                                            s_n=nor;
                                                            r_i_n=r_f; r_f_n=r_i_n + s_n*0.5;
                                                            if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',0.5); end
                            foton_muhitda=1;
                            [~,yonalish]=cheg_cos(yonalish, -nor,1,1.8);
                            shs=0;
                            while foton_muhitda==1
                                %                                                                                                             1401
                                mm=nuqta_slinder_masofa1(r_i, yonalish, r_muhit, 0, l);
                                if (mm==Inf)||(mm==0); foton_in_system=0; foton_muhitda=0; break; end
                                if (r_f(2)>L)||(r_f(2)<D)                %esdan chiqmasin
                                    %
                                    %                                                                                                                         143
                                    foton_in_system=0;
                                    foton_muhitda=0;
                                    break
                                end
                                
                                if (r_i(2)>L)||(r_i(2)<D)                %esdan chiqmasin
                                    %                                                                                                                    144
                                    foton_in_system=0; foton_muhitda=0;
                                    break
                                end
                                if mm>abs_leng_m
                                    mm=abs_leng_m;
                                    r_f=r_i+mm*yonalish;
                                    if (r_f(2)<D)||(r_f(2)>L)
                                        foton_muhitda=0; foton_in_system=0;
                                    else
%                                                                                 scatter3(r_f(3),r_f(2),r_f(1), 15,'fill')
                                        %                                             14011
                                        used_f1=used_f1+1;
                                        %                                              13111
                                        E_det1=E_det1+(1/ftu);
                                        foton_muhitda=0; foton_in_system=0;
                                        
                                        x_index = round((r_f(1)+ r_muhit)/x_grid_size)+1;
                                        z_index = round((r_f(3)-dd)/z_grid_size)+1;
                                        y_index = round((r_f(2))/y_grid_size)+1;
                                        E_xyz(x_index,z_index,y_index)=E_xyz(x_index,z_index,y_index)+1;
                                        r_f;
                                        if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','r','LineWidth',1); end
                                        % %                                                 if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','r','LineWidth',0.1); gs=gs+1; end
                                    end
                                else
                                    %                                                                                                                         14012
                                    abs_leng_m=abs_leng_m-mm;
                                    r_f=r_i+mm*yonalish;
                                    if chizish==1; plot3([r_i(3),r_f(3)],[r_i(2),r_f(2)],[r_i(1),r_f(1)],'Color','r','LineWidth',1); end
                                    % %                                         if chizish==1; is(gs)=plot3([r_i(1),r_f(1)],[r_i(2),r_f(2)],[r_i(3),r_f(3)],'Color','r','LineWidth',0.1); gs=gs+1; end
                                    nor=-[(r_f(1))/r_muhit 0 (r_f(3)-l)/r_muhit ];
                                    
                                    %                                    nor=-[(r_f(2)-h)/r_muhit  0 (r_f(1))/r_muhit];
                                                                            r_i_n=r_f; r_f_n=r_i_n + nor*0.1;
                                                                            if chizish==1; plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','y','LineWidth',1); end
                                    if  ((acos(dot(-yonalish,nor)/(norm(yonalish)*norm(nor)))) < asin(1/1.8))&&(rand()>(qaytish_sinish_ehtimolliqi(yonalish,nor,1.8,1)))
                                        [~,yonalish]=cheg_cos(yonalish,-nor,1.8,1);
                                                                                    r_i_n=r_f; r_f_n=r_i_n + yonalish*2;
                                        foton_muhitda=0;
                                        %                                               19013
%                                                                                     if chizish==1; plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','r','LineWidth',1.5);end
                                        
                                    else
                                        %                                                                                                                                     14014
                                        yonalish=yonalish-2*nor*dot(yonalish,nor);
                                    end
                                    r_i=r_f;
                                end
                               
                                
                                if (r_f(2)>L)||(r_f(2)<D)                %esdan chiqmasin
                                    foton_muhitda=0;
                                    %                                       193
                                    foton_in_system=0;
                                    
                                    break
                                end
                                
                                if (r_i(2)>L)||(r_i(2)<D)                %esdan chiqmasin
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
                        normal=normalni_top(r_f, para_coefs_21);
                                                    r_i_n=r_f; r_f_n=r_i_n + normal*1;
                                                    if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',0.5); end
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
                        normal=normalni_top(r_f, para_coefs_22);
                                                    r_i_n=r_f; r_f_n=r_i_n + normal*1;
                                                    if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(2),r_f_n(2)],[r_i_n(1),r_f_n(1)],'Color','c','LineWidth',0.5); end
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
        end
    end
end

used_f1/(2*nuqta_fotonlar_soni_max)
eff=E_det1/E_total
hold on
view(60, 30)
axis off
xlim([-5 22])
ylim([-8 58])
zlim([-15 15])