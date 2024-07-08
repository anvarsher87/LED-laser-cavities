clear; clc; close all;
% ll=linspace(6, 20, 14); 
%  for i=1:14  
xlabel('x')
ylabel('y')
zlabel('z')
D=50;
L=25; %in mm
% h=hh(i); %in mm
h=10;
% h=10; %in mm
r_muhit=2.5; %in mm
dd=5;
% l=ll(i);
l=6;
H=h+l;
LED_1=H-2*(l-r_muhit);
LED_2=H;
ly=(L-dd)/2;
hold on;
%------------------------------LED_1 rasmi (orqa devor)--------------------------------------------
a=1; b=0; c=0; d=LED_1;
yy=linspace(-2,L+2,2); zz=linspace(-2,D/2+2,2);
[y, z] = meshgrid(yy,zz); % Generate x and y data
x=(d-b*y-c*z)/a;
surf(z, x-0.1, y, 'EdgeColor', 'None', 'FaceColor', 'c', 'FaceAlpha', '0.1','LineWidth',1.5); % Plot the surface
surf(z, x-3, y, 'EdgeColor', 'None', 'FaceColor', 'c', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface
%------------------------------LED_2 rasmi (orqa devor)--------------------------------------------
a=1; b=0; c=0; d=LED_2;
yy=linspace(-2,L+2,2); zz=linspace(-2,D/2+2,2);
[y, z] = meshgrid(yy,zz); % Generate x and y data
x=(d-b*y-c*z)/a;
surf(z+D/2, x+0.1, y, 'EdgeColor', 'None', 'FaceColor', 'c', 'FaceAlpha', '0.1', 'LineWidth',1.5); % Plot the surface
surf(z+D/2, x+3, y, 'EdgeColor', 'None', 'FaceColor', 'c', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface
%------------------------------LED rasmi (yon devor)---------------------------------------------
r_k=20.5; theta=45; %for rotation plane
[xk, yk, zk]=cylinder(r_k, 4);
X2 = xk*cosd(theta) - yk*sind(theta);
Y2 = xk*sind(theta) + yk*cosd(theta);
surf(Y2+D/4, zk*3+LED_1-3, X2+L/2,'EdgeColor', 'k', 'FaceColor',  'c', 'FaceAlpha', '0.1','LineWidth',1.5)
surf(Y2+D-L/2, zk*3+H, X2+L/2,'EdgeColor', 'k', 'FaceColor',  'c', 'FaceAlpha', '0.1', 'LineWidth',1.5)
%%----------------------------plane_1---------------------------------------
%%------------------------past_1 yani LED 1 turgan joy---------------------
a=1; b=0; c=0; d=LED_1;
yy=linspace(0, L, 2); zz=linspace(0, D/2, 2);
[y, z] = meshgrid(yy,zz); % Generate x and y data
x=(d-b*y-c*z)/a;
tashqi_plot(1) = surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'w', 'FaceAlpha', '1', 'LineWidth',1.5); % Plot the surface
p1=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tan11=[0 0 -1]; tan12=[0 -1 0];
tek_1=[a b c d];
%%------------------------ tepa_1------------------------------------------
a=1; b=0; c=0; d=LED_1+H;
yy=linspace(ly,ly+dd,2); zz=linspace(0,D/2,2);
[y, z] = meshgrid(yy,zz); % Generate x and y data
x=(d-b*y-c*z)/a;
% tashqi_plot(2) = surf(z,x, y,  'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.1', 'LineWidth',1.5); % Plot the surface
p3=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_3=[a b c d];
%-------------------chizish uchun------------------------------------------
yy=linspace(ly,ly+dd,2); zz1=linspace(3,D/2-3,2);
[y, z] = meshgrid(yy,zz1); % Generate x and y data
surf(z,x, y,  'EdgeColor', 'k', 'FaceColor', 'b', 'FaceAlpha', '0.2', 'LineWidth',1.5); % Plot the surface

yy=linspace(ly,ly+dd,2); zz2=linspace(0,3,2);
[y, z] = meshgrid(yy,zz2); % Generate x and y data
surf(z,x, y,  'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface


yy=linspace(ly,ly+dd,2); zz2=linspace(D/2-3,D/2,2);
[y, z] = meshgrid(yy,zz2); % Generate x and y data
surf(z,x, y,  'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface

%%----------------yon_chap_1-----------------------------------------------
k=ly/H;
bd=-k*LED_1;
a=-k; b=1; c=0; d=bd;
xx=linspace(LED_1,LED_1+H,2); zz=linspace(0,D/2,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
% tashqi_plot(3)=surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.1', 'LineWidth',1.5); % Plot the surface
p4=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_4=[a b c d];

%-----------------chizish uchun--------------------------------------------
xx=linspace(LED_1,LED_1+H,2); zz=linspace(3,D/2-3,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'b', 'FaceAlpha', '0.2', 'LineWidth',1.5); % Plot the surface


xx=linspace(LED_1,LED_1+H,2); zz=linspace(0,3,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface


xx=linspace(LED_1,LED_1+H,2); zz=linspace(D/2-3,  D/2, 2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface


%%----------------yon_o'ng_1-----------------------------------------------
k=ly/H;
bd=L+k*LED_1;
a=k; b=1; c=0; d=bd;
xx=linspace(LED_1,LED_1+H,2); zz=linspace(0,D/2,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
% tashqi_plot(4) = surf( z, x, y, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.1', 'LineWidth',1.5); % Plot the surface
p5=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_5=[a b c d];
%-----------------chizish uchun--------------------------------------------
xx=linspace(LED_1,LED_1+H,2); zz=linspace(3,D/2-3,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'b', 'FaceAlpha', '0.2', 'LineWidth',1.5); % Plot the surface


xx=linspace(LED_1,LED_1+H,2); zz=linspace(0,3,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface


xx=linspace(LED_1,LED_1+H,2); zz=linspace(D/2-3,  D/2, 2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface




%%-----------------------qopqoq_1---------------------------------------------
a=0; b=0; c=-1; d=0;
xx=linspace(LED_1,LED_1+H,50); yy=linspace(0,L,50);
[x, y] = meshgrid(xx,yy); % Generate x and y data
z=(d-a*x-b*y)/c;
%surf(x, y, z, 'EdgeColor', 'None', 'FaceColor', 'r', 'FaceAlpha', '0.1') % Plot the surface
% patch([0 0 7 7], [0 25 16 9], [0 0 0 0], 's', 'FaceAlpha', '0.1')
tashqi_plot(5) = patch([0 0 0 0], [LED_1 LED_1 LED_1+H LED_1+H], [0 L L-ly ly],  's', 'FaceColor', 'k', 'FaceAlpha', '0.6', 'LineWidth',1.5); 
p6=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_6=[a b c d];
%%-------------------------mavhum o'rta devor_1-----------------------------
a=0; b=0; c=1; d=D/2;
xx=linspace(LED_1,LED_1+H,50); yy=linspace(0,L,50);
[x, y] = meshgrid(xx,yy); % Generate x and y data
z=(d-a*x-b*y)/c;
% surf(x, y, z, 'EdgeColor', 'None', 'FaceColor', 'r', 'FaceAlpha', '0.1') % Plot the surface
%patch([0 0 3.5 7, 7, 3.5], [0 25 20.5 25 0 4.5], [25 25 25 25 25 25], 's', 'FaceAlpha', '0.1')
tashqi_plot(6) = patch([D/2 D/2 D/2 D/2], [LED_1 LED_1 LED_1+H LED_1+H], [0 L L-ly ly],  's', 'FaceColor', 'k', 'FaceAlpha', '0.6', 'LineWidth',1.5);
p7=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];

tek_7=[a b c d];



%%----------------------------plane_2---------------------------------------
%%------------------------past_2----------------------------------------
a=-1; b=0; c=0; d=0;
yy=linspace(ly,ly+dd,2); zz=linspace(D/2,D,2);
[y, z] = meshgrid(yy,zz); % Generate x and y data
x=(d-b*y-c*z)/a;
% tashqi_plot(7) =surf(z, x, y,  'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', '0.1',  'LineWidth',1.5); % Plot the surface
p8=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_8=[a b c d];
%-----------------chizish uchun--------------------------------------------
yy=linspace(ly,ly+dd, 2); zz=linspace(D/2+3,D-3, 2);
[y, z] = meshgrid(yy,zz); % Generate x and y data
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'b', 'FaceAlpha', '0.2', 'LineWidth',1.5); % Plot the surface


yy=linspace(ly,ly+dd,2); zz=linspace(D/2,D/2+3,2);
[y, z] = meshgrid(yy,zz); % Generate x and y data
x=(d-b*y-c*z)/a;
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface


yy=linspace(ly,ly+dd,2); zz=linspace(D-3,  D, 2);
[y, z] = meshgrid(yy,zz); % Generate x and y data
x=(d-b*y-c*z)/a;
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface



%%------------------------ tepa_2 yani LED 2 turgan joy----------------------
a=1; b=0; c=0; d=LED_2;
yy=linspace(0,L,50); zz=linspace(D/2,D,50);
[y, z] = meshgrid(yy,zz); % Generate x and y data
x=(d-b*y-c*z)/a;
tashqi_plot(8) = surf(z, x, y,  'EdgeColor', 'None', 'FaceColor', 'w', 'FaceAlpha', '1', 'LineWidth',1.5); % Plot the surface
p2=-[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tan21=[0 0 1]; tan22=[0 -1 0];
tek_2=[a b c d];
%%----------------yon_chap_2-----------------------------------------------
k=-ly/H;
bd=ly;
a=-k; b=1; c=0; d=bd;
xx=linspace(0,LED_2,2); zz=linspace(D/2,D,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
% tashqi_plot(9) = surf(z,x, y,  'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.1',  'LineWidth',1.5); % Plot the surface
p9=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_9=[a b c d];
%---------------chizish uchun ---------------------------------------------
xx=linspace(0,LED_2,2); zz=linspace(D/2+3,D-3,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z,x, y,  'EdgeColor', 'k', 'FaceColor', 'b', 'FaceAlpha', '0.2',  'LineWidth',1.5); % Plot the surface


xx=linspace(0,LED_2,2); zz=linspace(D/2,D/2+3,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z,x, y,  'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4',  'LineWidth',1.5); % Plot the surface

xx=linspace(0,LED_2,2); zz=linspace(D-3,D,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z,x, y,  'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4',  'LineWidth',1.5); % Plot the surface

%%----------------yon_o'ng_2-----------------------------------------------
k=ly/H;
bd=-(L-k*H);
a=-k; b=1; c=0; d=-bd;
xx=linspace(LED_2,0,2); zz=linspace(D,D/2,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
% tashqi_plot(10) = surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.1', 'LineWidth',1.5); % Plot the surface
p10=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_10=[a b c d];
%-------------------chizish uchun------------------------------------------
xx=linspace(LED_2,0,2); zz=linspace(D-3,D/2+3,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'b', 'FaceAlpha', '0.2', 'LineWidth',1.5); % Plot the surface

xx=linspace(LED_2,0,2); zz=linspace(D,D-3,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface

xx=linspace(LED_2,0,2); zz=linspace(D/2+3,D/2,2);
[x, z] = meshgrid(xx,zz); % Generate x and y data
y=(d-a*x-c*z)/b;
surf(z, x, y, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', '0.4', 'LineWidth',1.5); % Plot the surface


%%-----------------------qopqoq_2---------------------------------------------
a=0; b=0; c=1; d=D;
xx=linspace(0,h,50); yy=linspace(0,L,50);
[x, y] = meshgrid(xx,yy); % Generate x and y data
z=(d-a*x-b*y)/c;
% surf(x, y, z, 'EdgeColor', 'None', 'FaceColor', 'r', 'FaceAlpha', '0.1') % Plot the surface
% patch([0 0 7 7], [9 16 25 0], [50 50 50 50], 's', 'FaceAlpha', '0.1')
tashqi_plot(11) = patch([D D D D], [LED_2 LED_2 0 0], [0 L L-ly ly], 's', 'FaceColor', 'k', 'FaceAlpha', '0.6', 'LineWidth',1.5);
p11=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_11=[a b c d];
%%-------------------------mavhum o'rta devor_2-----------------------------
a=0; b=0; c=1; d=D/2;
xx=linspace(0,h,50); yy=linspace(0,L,50);
[x, y] = meshgrid(xx,yy); % Generate x and y data
z=(d-a*x-b*y)/c;
% surf(x, y, z, 'EdgeColor', 'None', 'FaceColor', 'r', 'FaceAlpha', '0.1') % Plot the surface
%patch([0 0 3.5 7, 7, 3.5], [0 25 20.5 25 0 4.5], [25 25 25 25 25 25], 's', 'FaceAlpha', '0.1')
tashqi_plot(12) = patch([D/2 D/2 D/2 D/2], [LED_2 LED_2 0 0], [0 L L-ly ly], 's', 'FaceColor', 'k', 'FaceAlpha', '0.6',  'LineWidth',1.5);
p12=[a/(sqrt(a^2+b^2+c^2)) b/(sqrt(a^2+b^2+c^2)) c/(sqrt(a^2+b^2+c^2))];
tek_12=[a b c d];
%%--------------------rod-------------------------------------------------------
[xs, ys, zs]=cylinder(r_muhit,40); surf( D*1.2*zs-5, h+r_muhit+xs, L/2+ys,'EdgeColor', 'None', 'FaceColor', 'm', 'FaceAlpha', '0.2')
% plotCircle3D_2([ L/2, h+r_muhit,D+0.01, ],[0,0,1], r_muhit-0.1)
plotCircle3D_2([ L/2, h+r_muhit,D+0.01+5, ],[0,0,1], r_muhit-0.1)
% plotCircle3D_2([L/2, h+r_muhit, 0-0.01],[0,0,1], r_muhit-0.1)
plotCircle3D_2([L/2, h+r_muhit, 0-0.01-5],[0,0,1], r_muhit-0.1)
hold on;

%%--------------------------------crystal----------------------------------
A=dlmread('15_12_2021_BLUE_LED_ems.txt');
% A=dlmread('blue_led.txt');
% % B=dlmread('Nd_YAG_abs.txt'); % yordamichi malumot, lazer qristalining yutish
C=dlmread('Ce_Nd_YAG_abs.txt');
%%------------------------------------------------------------------------------
sa=0.9*pi/2;
teflon_kaytarish=0.9;  %blue uchun
led_kaytarish=0.75;    %blue uchun
LEDdan_qaytdi=0;
lost_at_led=0;
loss_at_teflon=0;
lost_orta_devor=0;
led_array
chizish=1;
used_f1=0; used_f2=0;
nuqta_fotonlar_soni_max=100;
axis equal
% % axis off
xlim([-10 58])
zlim([-10 30])
ylim([-10 LED_2+10])
%  view (57,17)
view(-150,30)
axis off
% view(90,0)
sp=zeros(1,2500);
s=zeros(1,nuqta_fotonlar_soni_max);
E_total=0;
E_det1=0;  E_det2=0;

% bnm

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

% y{10 15}; z{0 25}; x{h h+5}
DD=D/2;
x_grid_number  = 49;        x_grid_size = 2*r_muhit/x_grid_number;
y_grid_number  = 49;        y_grid_size = 2*r_muhit/y_grid_number;
z_grid_number  = 199;       z_grid_size = D/z_grid_number;

E_xyz = zeros(x_grid_number+1, y_grid_number+1, z_grid_number+1);

        tic
        for led_number0=1:2
            for f=1:nuqta_fotonlar_soni_max
                if rand > 0.999995; clc; disp("Algaritim bajarilish prosesi " + round(100 * f/nuqta_fotonlar_soni_max) + " %"); pause(0.1);  end  % bi xator loop nera galganini go'rsatadi
                
                ftu=s(1,f);
                E_total=E_total+(1/ftu);
                %-------------------------------------absorbtion---------------------------
                abs_coef=C(ftu); abs_length=absorption_length(abs_coef); abs_leng_m=abs_length;
                if (abs_leng_m>1900);  yutilmaydigan_f=yutilmaydigan_f+1; continue; end
                
                %--- -----------------------------------------------------------------------
                
                %------------- boshlang'ich kordinata va yonalish--------------
                if led_number0==1
                %     r_i=reshape(led_nuqta_kor(led_point, :),[1,3]);
                r_i=led_random_square_1(LED_1, 12.5, 12.5, 12.5, 12.5);
                %        scatter3 (r_i(3), r_i(1), r_i(2),  15,'fill')
                yonalish = lambert_shape_4 (p1, tan11, tan12, sa);
%                     r_f = r_i + yonalish*30;
%                     if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','b','LineWidth',0.5); end
                 led_number=led_number0;
                end
                
                 if led_number0==2
                %     r_i=reshape(led_nuqta_kor(led_point, :),[1,3]);
                r_i=led_random_square_1(LED_2, 12.5, D/2+12.5, 12.5, 12.5);
%                        scatter3 (r_i(3), r_i(1), r_i(2),  15,'fill')
                yonalish = lambert_shape_4 (p2, tan21, tan22, sa);
%                     r_f = r_i + yonalish*30;
%                     if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','b','LineWidth',0.5); end
                 led_number=led_number0;
                end
                %--------------------------------------------------------------
                
                %             abs_leng_m=10*abs_length;
                
                
                foton_in_system=1;
                while (foton_in_system == 1)
                    %--------------------case_1------------------------------------
                    if (led_number==1)&&(foton_in_system == 1)
                        
                        while (led_number==1)&&(foton_in_system ==1)
                            
                            %---------------------LEDs--------------------------------
                            dm(1)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_1);  tk1(1,:)=[p1(1) p1(2) p1(3)];    %past_1 Led turgan joy
                            %---------------------Laser rod----------------------------
                            dm(2)=nuqta_slinder_masofa(r_i, yonalish, r_muhit, h+r_muhit, L/2); tk1(2,:)= [1/r_muhit 1/r_muhit 0];              %aktiv muhit;
                            %---------------------Cavity_walls-------------------------
                            dm(3)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_3);  tk1(3,:)= -[p3(1) p3(2) p3(3)];     %tepa_1
                            dm(4)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_4);  tk1(4,:)= [p4(1) p4(2) p4(3)];    % yon_chap_1
                            dm(5)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_5);  tk1(5,:)=- [p5(1) p5(2) p5(3)];   %yon_o'ng_1
                            dm(6)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_6);  tk1(6,:)=- [p6(1) p6(2) p6(3)];  %qopqoq_1
                            %-----------------------O'rtakash--------------------------
                            dm(7)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_7);  tk1(7,:)= -[p7(1) p7(2) p7(3)];  %mavhum o'rta devor
                            %----------------------------------------------------------
                            for i=1:length(dm); if dm(i)<0.000001; dm(i)=Inf; end; end
                            
                            [db,joy]=min(dm); normal=tk1(joy,:); r_f=r_i + yonalish*db;
                            
                            if (db==Inf)||(db==0); foton_in_system=0; break; end
                            
                            %                          if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','b','LineWidth',0.1); gs=gs+1; end
                                                     if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','b','LineWidth',1); end
                            r_i=r_f;
                            if (joy==1)
                                     r_i_n=r_f; r_f_n=r_i_n + normal*1; if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','c','LineWidth',0.5); end
                                if led_kaytarish > rand
                                    %                                                                                                  112
                                    yonalish=yonalish-2*normal*dot(yonalish,normal);
                                    LEDdan_qaytdi=LEDdan_qaytdi+1;
                                else
                                    %                                                                                                  113
                                    foton_in_system=0; lost_at_led=lost_at_led+1;
                                end
                            end
                            
                            if (joy==2)
                                %                                                                                      13
                                nor=[(r_f(1)-(h+r_muhit))/r_muhit (r_f(2)-L/2)/r_muhit 0];
                                if (rand()>(qaytish_sinish_ehtimolliqi(yonalish,nor,1,1.8)))
                                    %                                                                                            131
                                                                    s_n=nor;
                                                                    r_i_n=r_f; r_f_n=r_i_n + s_n*1;
                                     if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','c','LineWidth',0.5); end
                                    foton_muhitda=1;
                                    [~,yonalish]=cheg_cos(yonalish, -nor,1,1.8);
                                    shs=0;
                                    while foton_muhitda==1
                                        %                                                                                                             1311
                                        mm=nuqta_slinder_masofa(r_i, yonalish, r_muhit, h+r_muhit, L/2);
                                        if (mm==Inf)||(mm==0); foton_in_system=0; foton_muhitda=0; break; end
                                        if (r_f(3)>D/2)||(r_f(3)<0)                %esdan chiqmasin
                                            %                                                                                                                         1312
                                            foton_in_system=0;
                                            foton_muhitda=0;
                                            break
                                        end
                                        
                                        if (r_i(3)>D/2)||(r_i(3)<0)                %esdan chiqmasin
                                            %                                                                                                                          2312
                                            foton_in_system=0; foton_muhitda=0;
                                            break
                                        end
                                        if mm>abs_leng_m
                                            mm=abs_leng_m;
                                            r_f=r_i+mm*yonalish;
                                            if (r_f(3)<0)||(r_f(3)>D/2)
                                                foton_muhitda=0; foton_in_system=0;
                                            else
                                                
%                                                 scatter3(r_f(3),r_f(1),r_f(2), 30,'fill')
                                               
                                                x_index = round((r_f(1)-h)/x_grid_size)+1;
                                                z_index = round((r_f(3))/z_grid_size)+1;
                                                y_index = round((r_f(2)-(L/2-r_muhit))/y_grid_size)+1;
                                                E_xyz(x_index,y_index,z_index)=E_xyz(x_index,y_index,z_index)+1;
                                                
                                                foton_muhitda=0; foton_in_system=0; used_f1=used_f1+1; E_det1=E_det1+(1/ftu);
                                            end
                                            %                                                                                                                         13111
                                            
                                                                                          if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',1); end
                                            % %                                                 if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',0.1); gs=gs+1; end
                                        else
                                            %                                                                                                                         13112
                                            abs_leng_m=abs_leng_m-mm;
                                            r_f=r_i+mm*yonalish;
                                                                                  if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',1); end
                                            % %                                         if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',0.1); gs=gs+1; end
                                            nor=-[(r_f(1)-(h+r_muhit))/r_muhit (r_f(2)-L/2)/r_muhit 0];
                                            %                                         r_i_n=r_f; r_f_n=r_i_n + nor*1;
                                            %                                         if chizish==1; plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','y','LineWidth',0.8); end
                                            if  ((acos(dot(-yonalish,nor)/(norm(yonalish)*norm(nor)))) < asin(1/1.8))&&(rand()>(qaytish_sinish_ehtimolliqi(yonalish,nor,1.8,1)))
                                                [~,yonalish]=cheg_cos(yonalish,-nor,1.8,1);
                                                %                                             r_i_n=r_f; r_f_n=r_i_n + yonalish*2;
                                                foton_muhitda=0;
                                                %                                               131121
                                                %                                             if chizish==1;plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','r','LineWidth',1.5);end
                                            else
                                                %                                                                                                                                     131122
                                                yonalish=yonalish-2*nor*dot(yonalish,nor);
                                            end
                                            r_i=r_f;
                                        end
                                        
                                        if (r_f(3)>D/2)||(r_f(3)<0)                %esdan chiqmasin
                                            %                                                                                                                         1312
                                            foton_in_system=0;
                                            foton_muhitda=0;
                                            break
                                        end
                                        
                                        if (r_i(3)>D/2)||(r_i(3)<0)                %esdan chiqmasin
                                            %                                                                                                                          2312
                                            foton_in_system=0; foton_muhitda=0;
                                            break
                                        end
                                        shs=shs+1; if shs>10; foton_muhitda=0; foton_in_system=0; break; end
                                    end
                                else
                                    %                                                                                                 132
                                    yonalish=yonalish-2*nor*dot(yonalish,nor);
                                end
                            end
                            
                            if (joy>=3)&&(joy<=7)
                                                            r_i_n=r_f; r_f_n=r_i_n + normal*1; if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','c','LineWidth',0.5); end
                                %                                                         14
                                %                             joy
                                normal;
                                %                             if joy==7
                                %                                                                                                  18
                                %                                 if r_f(1)<h/2
                                %                                     kichik_asos_yoni = r_f(1);
                                % %                                                                         181
                                %                                 else
                                %                                     kichik_asos_yoni = h/2-(r_f(1)-h/2);
                                % %                                                                         182
                                %                                 end
                                %
                                %                                 if abs(r_f(2)-L/2)<(dd+2*kichik_asos_yoni/tan(burchak))/2
                                %                                     led_number=2;
                                %                                     r_i=r_f;
                                % %                                                                         183
                                %                                 else
                                %                                     yonalish=yonalish-2*normal*dot(yonalish,normal);
                                % %                                                                         184
                                %                                 end
                                %
                                %                             else
                                %                                 joy
                                %                                                                 141
                                %
                                if teflon_kaytarish > rand
                                    %                                                                         1411
                                    yonalish=yonalish-2*normal*dot(yonalish,normal);
                                else
                                    %                                                                                                              1412
                                    foton_in_system=0; loss_at_teflon=loss_at_teflon+1;
                                end
                                %                             end
                            end
                            
                        end
                    end
                    
                    
% % %                     --------------------------- case 2------------------------
                                    if (led_number==2) && (foton_in_system == 1)
                    
                                        while (led_number==2)&&(foton_in_system ==1)
                    %                         200
                                            %---------------------LEDs--------------------------------
                                            dm(1)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_2);  tk2(1,:)=-[p2(1) p2(2) p2(3)];     %tepa_2 Led turgan joy
                                            %---------------------Laser rod----------------------------
                                            dm(2)=nuqta_slinder_masofa(r_i, yonalish, r_muhit, h+r_muhit, L/2); tk2(2,:)=[1/r_muhit 1/r_muhit 0];              %aktiv muhit;
                                            %---------------------Cavity_walls-------------------------
                                            dm(3)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_8);  tk2(3,:)=-[p8(1) p8(2) p8(3)];    %past_2
                                            dm(4)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_9);  tk2(4,:)= [p9(1) p9(2) p9(3)];    % yon_chap_2
                                            dm(5)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_10);  tk2(5,:)= -[p10(1) p10(2) p10(3)];   %yon_o'ng_2
                                            dm(6)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_11);  tk2(6,:)=- [p11(1) p11(2) p11(3)];  %qopqoq_2
                                            %-----------------------O'rtakash--------------------------
                                            dm(7)=tekis_vek(r_i(1),r_i(2),r_i(3),yonalish,tek_12);  tk2(7,:)= -[p12(1) p12(2) p12(3)];  %mavhum o'rta devor
                                            %----------------------------------------------------------
                                            for i=1:length(dm); if dm(i)<0.000001; dm(i)=Inf; end; end
                    
                    
                                            [db,joy]=min(dm); normal=tk2(joy,:);
                                            if db == Inf; foton_in_system=0; break; end
                                           
                                            r_f=r_i + yonalish*db;
                    %                          if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','b','LineWidth',0.1); gs=gs+1; end
                                            if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','b','LineWidth',1); end
                                              r_i=r_f;
                                            if (joy==1)
                                                r_i_n=r_f; r_f_n=r_i_n + normal*1; if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','c','LineWidth',0.5); end
                    %                                                                                      211
                                                if led_kaytarish > rand
                    %                                                                                                   212
                                                    yonalish=yonalish-2*normal*dot(yonalish,normal);
                                                    LEDdan_qaytdi=LEDdan_qaytdi+1;
                                                else
                    %                                                                                                 213
                                                    foton_in_system=0;lost_at_led=lost_at_led+1;
                                                end
                                            end
                                            if (joy==2)
                    %                                                                                      23
                                                nor=[(r_f(1)-(h+r_muhit))/r_muhit (r_f(2)-L/2)/r_muhit 0];
                                                if (rand()>(qaytish_sinish_ehtimolliqi(yonalish, nor,1,1.8)))
                    %                                                                                                 231
                                                    s_n=nor;
                                                    r_i_n=r_f; r_f_n=r_i_n + s_n*1;
                                                    if chizish==1; plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','c','LineWidth',0.5); end
                                                        foton_muhitda=1;
                                                    [~,yonalish]=cheg_cos(yonalish, -nor,1,1.8);
                    %                                 r_i_n=r_f; r_f_n=r_i_n + yonalish*1;
                    %                                 plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','k','LineWidth',1.5);
                                                    shs2=0;
                                                    while foton_muhitda==1
                    %                                                                                                              2311
                                                        mm=nuqta_slinder_masofa(r_i, yonalish, r_muhit, h+r_muhit, L/2);
                                                        if (mm==Inf)||(mm==0); foton_in_system=0; foton_muhitda=0; break; end
                                                        if mm==Inf; xatolik_yuz_berdi; end
                                                        yonalish;
                                                        abs_leng_m;
                                                        
                                                        if (r_f(3)>D)||(r_f(3)<D/2)                %esdan chiqmasin
                                                            %                                                                                                                         1312
                                                            foton_in_system=0;
                                                            foton_muhitda=0;
                                                            break
                                                        end
                                                        
                                                        if (r_i(3)>D)||(r_i(3)<D/2)                %esdan chiqmasin
                                                            %                                                                                                                          2312
                                                            foton_in_system=0; foton_muhitda=0;
                                                            break
                                                        end
                                                        
                                                        if mm>abs_leng_m
                                                            mm=abs_leng_m;
                                                            r_f=r_i+mm*yonalish;
                                                            if (r_f(3)>D)||(r_f(3)<D/2)                %esdan chiqmasin
                                                            %                                                                                                                         1312
                                                            foton_in_system=0;
                                                            foton_muhitda=0;
                                                            else
%                                                             scatter3(r_f(3),r_f(1),r_f(2),30,'fill')
                                                            
                                                            x_index = round((r_f(1)-h)/x_grid_size)+1;
                                                            z_index = round((r_f(3))/z_grid_size)+1;
                                                            y_index = round((r_f(2)-(L/2-r_muhit))/y_grid_size)+1;
                                                            E_xyz(x_index,y_index,z_index)=E_xyz(x_index,y_index,z_index)+1;
                    
                                                            foton_muhitda=0; foton_in_system=0; used_f2=used_f2+1; E_det2=E_det2+(1/ftu);
                    %                                                                                                                          23111
                    
                                                            if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',1); end
                    %                                          if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',0.1); gs=gs+1; end
                                                            end
                                                        else
                    %                                                                                                                          23112
                                                            abs_leng_m=abs_leng_m-mm;
                                                            r_f=r_i+mm*yonalish;
                                                            if chizish==1; plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',1); end
                    %                                        if chizish==1; is(gs)=plot3([r_i(3),r_f(3)],[r_i(1),r_f(1)],[r_i(2),r_f(2)],'Color','r','LineWidth',0.1); gs=gs+1; end
                                                           nor=-[(r_f(1)-(h+r_muhit))/r_muhit (r_f(2)-L/2)/r_muhit 0];
                    %                                         r_i_n=r_f; r_f_n=r_i_n + nor*1;
                                                            if chizish==1; plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','y','LineWidth',0.5); end
                                                            if  ((acos(dot(-yonalish,nor)/(norm(yonalish)*norm(nor)))) < asin(1/1.8))&&(rand()>(qaytish_sinish_ehtimolliqi(yonalish,nor,1.8,1)))
                                                                [~,yonalish]=cheg_cos(-yonalish,-nor,1.8,1);
                    %                                             r_i_n=r_f; r_f_n=r_i_n + yonalish*2;
                    %                                             231121
%                                                                 if chizish==1; plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],'[r_i_n(2),r_f_n(2)],Color','k','LineWidth',1.5); end
                                                                foton_muhitda=0;
                                                            else
                    %                                                                                                                                      231122
                                                                yonalish=yonalish-2*nor*dot(yonalish,nor);
                                                            end
                                                            r_i=r_f;
                                                        end
                                                       
                                                        if (r_f(3)>D)||(r_f(3)<D/2)                %esdan chiqmasin
                    %                                                                                                                          2312
                                                            foton_in_system=0; foton_muhitda=0;
                                                            break
                                                        end
                    
                                                        if (r_i(3)>D)||(r_i(3)<D/2)                %esdan chiqmasin
                    %                                                                                                                          2312
                                                            foton_in_system=0; foton_muhitda=0;
                                                            break
                                                        end
                                                        shs2=shs2+1; if shs2>10; foton_muhitda=0; foton_in_system=0; break; end
                                                    end
                                                else
                    %                                                                                                  232
                                                    yonalish=yonalish-2*nor*dot(yonalish,nor);
                                                end
                                            end
                                                if (joy>=3)&&(joy<=7)
                    %                             r_i_n=r_f; r_f_n=r_i_n + normal*1; if chizish==1;  plot3([r_i_n(3),r_f_n(3)],[r_i_n(1),r_f_n(1)],[r_i_n(2),r_f_n(2)],'Color','r','LineWidth',1.5); end
                                                %                             joy
                    %                                                         29
                    %                             if joy==7
                    % %                                                                 2213
                    %                                 if r_f(1)<h/2
                    %                                     kichik_asos_yoni = r_f(1);
                    % %                                                                         22131
                    %                                 else
                    %                                     kichik_asos_yoni = h/2-(r_f(1)-h/2);
                    % %                                                                         22132
                    %                                 end
                    %
                    %                                 if abs(r_f(2)-L/2)<(dd+2*kichik_asos_yoni/tan(burchak))/2
                    %                                     led_number=1;
                    %                                     r_i=r_f;
                    % %                                                                         22133
                    %                                 else
                    %                                     yonalish=yonalish-2*normal*dot(yonalish,normal);
                    % %                                                                         22134
                    %                                 end
                    %                             else
                    %                                                                    290
                                                    %                                    joy
                                                    if teflon_kaytarish > rand
                    %                                                                                                              291
                                                        yonalish=yonalish-2*normal*dot(yonalish,normal);
                                                    else
                    %                                                                                                             292
                                                        foton_in_system=0;  loss_at_teflon=loss_at_teflon+1;
                                                    end
                                                end
                    
                    %                         end
                    
                                        end
                                    end
                    
                    
%                                 end
%                     %               if gs>8; continue; end
%                     %             pause(10)
%                             end
                    %         delete(is)
                    %         pause(1)
                    
                    
                end
%                     pause(0.001)
            end
        end
   toc             
       
%     end
% end
% pause(3)
used_f1
used_f2
eff_u=used_f1/(nuqta_fotonlar_soni_max*2)
eff=(E_det1+E_det2)/E_total
% close all
%  end
% axis off