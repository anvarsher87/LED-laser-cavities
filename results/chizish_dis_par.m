clear; clc; close all;

% load('dist_ce_for_tem_2.mat')



% load('modeldan_garaki_10mln_single_50_50_200.mat')
load('modeldan_garaki_1mln_s_50_50_200_tek.mat')


ram=2.5;
L=50;
r_muhit = ram;

%-----------------------------------


x_grid_number  = 50;        x_grid_size = 2*r_muhit/x_grid_number;
y_grid_number  = 50;        y_grid_size = 2*r_muhit/y_grid_number;
z_grid_number= 200;         z_grid_size  = L/z_grid_number;

ss0=size(E_xyz);
E_xyz=smooth3(E_xyz);
clear X Y Z
ui=1;
% hold on
for ix=1:ss0(1)

    for iy=1:ss0(2)
        
        for iz=1:ss0(3)
                     

                
                if (((ix-1)*x_grid_size-ram)^2+((iy-1)*y_grid_size-ram)^2 <= ram^2)
%                     if (ix<(ss0(1)/2))&&(iz<=(ss0(3)/2))
%                         continue
%                     else
                    
                    
                    xc=(ix-1)*x_grid_size-ram;
                    yc=(iy-1)*y_grid_size-ram;
                    zc=(iz-1)*z_grid_size-0.5*L;
                    %                 xc
                    
                    %                 Rho  = hypot(xc,zc);
                    %                 Teta = atan2(zc,xc); if Teta<0; Teta=Teta+2*pi; end
                    %
                    %                 %----------------------------------------------------------
                    %                 rho_index = round(Rho/rho_grid_size)+1;
                    %                 teta_index = round((Teta)/teta_grid_size)+1;
                    %                 z_s_index = round((yc+0.5*L)/z_grid_size)+1;
                    %                 E_slindr(rho_index,teta_index,z_s_index)=E_slindr(rho_index,teta_index,z_s_index)+1;
                    %                 %----------------------------------------------------------
                    
                    
                    X(ui)=xc;
                    Y(ui)=yc;
                    Z(ui)=zc;
                    %                 S21A(ui)=mode_dist_x_y_z(ix,iy,iz);
                    %                 S21A(ui)=u(ix,iy,iz);
                    
                    S21A(ui)=E_xyz(ix,iy,iz);
                    %                 C=E_xyz(ix,iy,iz);
                    %
                    ui=ui+1;
                    %                 if E_xyz(ix,iy,iz)>0 && rand>0.9999
                    %                     %                     scatter3(xc,yc,zc)
                    %                     C
                    %                     scatter3(xc,yc,zc, 20, C, 'filled','MarkerEdgeColor','None')
                    %                 end
                    
                    
                    
%                     end
               
                end
            
            
        end
    end
    
end

figure
scatter3( Z,  Y, X,  20, S21A, 'filled','MarkerEdgeColor','None')

% scatter3(Z,Y,X)
xlabel('x'); ylabel('y'); zlabel('z');
% axis off
% grid minor
axis equal
% xlim([-2.5 2.5])
zlim([-2.5 0])

colormap jet
view(-154, 29)
adasdad
% %----------------------------------------

% 
% load('modeldan_garaki_10_5_nd_yag.mat')
% 
% ss=size(E_rtz);
% E_xyz=smooth3(E_rtz);
% 
% x=linspace(-ram,ram,ss(1));
% y=linspace(-ram,ram,ss(2));
% z=linspace(-L/2,L/2,ss(3));
% 
% ro=linspace(0,ram,25); dr=ro(2)-ro(1);
% dro   = dr;
%  
% teta0 =linspace(0,2*pi,72);
% dteta = teta0(2)-teta0(1);
% dzs   = z(2)-z(1);
% 
% 
% for i=1:length(x)
%     for j=1:length(y)
%         for k=1:length(z)
%             
%             xd=x(i); yd=y(j); zd=z(k);
%             rr=(xd^2+yd^2)^(0.5);
%             if rr <= ram
%                 %-----------------------------------
%                 if (xd>0)&&(yd>0)
%                     teta=180*atan(yd/xd)/pi;
%                 elseif (xd==0)&&(yd>0)
%                     teta=90;
%                 elseif (xd==0)&&(yd<0)
%                     teta=270;
%                 elseif (xd>0)&&(yd==0)
%                     teta=0;
%                 elseif (xd<0)&&(yd==0)
%                     teta=180;
%                 elseif (xd<0)&&(yd>0)
%                     teta=180+180*atan(yd/xd)/pi;
%                 elseif (xd<0)&&(yd<0)
%                     teta=180+180*atan(yd/xd)/pi;
%                 elseif (xd>0)&&(yd<0)
%                     teta=2*180+180*atan(yd/xd)/pi;
%                 elseif xd==0 && yd==0
%                     teta=0;
%                 else
%                     disp("Xatolik");
%                 end
%                 teta;
%                 teta=pi*teta/180;
%                 %-----------------------------------
%                 teta_index = round(teta/dteta)+1;
%                 ro_index   = round(rr/dro)+1;
%                 zs_index   = k;
%                 %---------------------------------------
%                 
% %                 if (E_xyz(i,j,k)==0)&&(i~=1)&&(j~=1)&&(k~=1)&&(i~=length(x))&&(j~=length(y))&&(k~=length(z))
% % %                     hggh
% %                     E_xyz(i,j,k) = ((E_xyz(i-1,j,k)+E_xyz(i+1,j,k))/2 + ... 
% %                         (E_xyz(i,j-1,k)+E_xyz(i,j+1,k))/2 + ...
% %                         (E_xyz(i,j,k-1)+E_xyz(i,j,k+1))/2)/3;
% %                     
% %                 end
%                 
%                 
%                 E_ro_teta_z(ro_index,teta_index,zs_index) = E_xyz(i,j,k);
%             end
%         end
%     end
% end
% 
% 
% 
% % figure
% Q=E_ro_teta_z;
% S2p=Q;
% 
% ss=size(Q);
% 
% %----------------------------------------------------
% %filling defective locations
% for ii=1:ss(1)
%     for jj=1:ss(2)
%         for kk=1:ss(3)
%             
%  
%             
%             if (Q(ii,jj,kk) == 0)&& (ii < ss(1)-5)
%                 
%                 in=1;
%                 while (Q(ii,jj,kk) == 0)&&( in<5)
%                     
%                     Q(ii,jj,kk)=Q(ii+in,jj,kk);
%                     in=in+1;
%                 end
%                 
% 
%                 
%             end
%             
%         end
%     end
% end

%----------------------------------------------------



ui=1;

clear X Y Z
load('modeldan_garaki_10_5_nd_yag_3.mat')

ss=size(E_rtz);
Q=smooth3(E_rtz);

ro=linspace(0, ram, ss(1)); dr=ro(2)-ro(1);
dro   = dr;
z=linspace(-L/2,L/2,ss(3)); 
teta0 =linspace(0, 2*pi, ss(2));
% dteta = teta0(2)-teta0(1);
% dzs   = z(2)-z(1);

for ix=1:ss(1)
    
    for iy=1:ss(2)
        
        for iz=1:ss(3)
            
            X(ui) = ro(ix)*cos(teta0(iy));
            Y(ui) = ro(ix)*sin(teta0(iy));
            Z(ui) = z(iz);
            S31A(ui)=Q(ix,iy,iz);
            ui=ui+1;
            
        end
    end
    
end

scatter3(Z,Y,X, 20, S31A, 'filled','MarkerEdgeColor','None')

xlabel('x'); ylabel('y'); zlabel('z');
axis equal; 
% xlim([0 5.5])
% zlim([-0.25 0 ])

%-----------------------------------------

