clear; close all; 
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% load('modeldan_garaki_10_6_nd_yag.mat')





clear X Y Z
ui=1; %loop=1;
load('modeldan_garaki_trapetsiya_10mln_50_50_200_oe_check.mat')
ad=E_xyz;
S2=smooth3(ad);
ss0=size(E_xyz);

ram=2.5;
L=50;


x_grid_number  = 50;        x_grid_size    = 2*ram/x_grid_number;
z_grid_number  = 200;       z_grid_size  = L/z_grid_number;
y_grid_number  = 50;        y_grid_size    = 2*ram/y_grid_number;

subplot(4,2,1);
for ix=1:ss0(1)
    
    
    for iy=1:ss0(2)
        
        for iz=1:ss0(3)
            
            if (ix>ss0(1)/2)&&(iz<ss0(3)/2.1)%&&(iy<ss0(2)/2)
                     continue  
            else
            
            X(ui)=ix-0;
            Z(ui)=iz-0;
            Y(ui)=iy;
            if (X(ui)-25)^2+(Y(ui)-25)^2 > 25^2
                
                D21A(ui)=0;
                
            else
                D21A(ui)=S2(ix,iy,iz);
                ui=ui+1;
            end
            end
%             S21A(ui)=S2(ix,iy,iz);
%                 ui=ui+1;
           

%            

        end
    end
end
hold on
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% scatter3((Z-500)*0.1,(Y-32)*0.1+3,(X-32)*0.1+3, 15, S21A, 'filled','MarkerEdgeColor','None');
% scatter3((Z-500)*0.05,(Y-25)*0.05,(X-25)*0.05, 15, S21A, 'filled','MarkerEdgeColor','None');
% scatter3((Z-0)*0.05+80, (Y-0)*0.05+2.5, X*0.05+2.5, 15, S21A, 'filled','MarkerEdgeColor','None');
% scatter3(X*0.05-0.1, Y*0.05-0.1, Z*0.05-5, 15, S21A, 'filled','MarkerEdgeColor','None');
scatter3( (Z)*z_grid_size, (Y)*y_grid_size, (X)*x_grid_size, 15, D21A, 'filled','MarkerEdgeColor','None');
% scatter3( (Z)*z_d_grid_size, (Y)*y_grid_size, (X)*x_grid_size, 15, S21A, 'filled','MarkerEdgeColor','None');
% scatter3( (Z), (Y-1), (X-1), 15, S21A, 'filled','MarkerEdgeColor','None');
caxis([0 60])
colormap jet

axis equal
hold on
axis off






% % grid on
% % axis off
% % view(90,200)
% % view(-108,13)
view(-112,11)
xlim([0,50])
ylim([0, 5])
zlim([0, 5])
% axis on

% XZCsfz
% erdctfvgybhujnk
% 
subplot(4,2,2);
load('modeldan_garaki_trapetsiya_10mln_50_50_200_oe_check.mat')
ad=E_xyz;
S2=smooth3(ad);
integral_array_1= squeeze(sum(S2, 3));
% integral_array_1=rescale(integral_array_1);
surf(integral_array_1, 'EdgeColor', 'None');
axis off
axis equal
view(-90, 90)
% colorbar
caxis([0*10^3 1*10^4])
% caxis([0.1*10^4 5*10^4])

% Plot the surface


% dfghbjnkm
clear X Y Z
ui=1; loop=1;
load('modeldan_garaki_hexagon_10mln_50_50_200_oe_check.mat')
bd=E_xyz;
C2=smooth3(bd);
ss1=size(E_xyz);


ram=2.5;
L=50;


x_grid_number  = 50;        x_grid_size    = 2*ram/x_grid_number;
z_grid_number  = 200;       z_grid_size  = L/z_grid_number;
y_grid_number  = 50;        y_grid_size    = 2*ram/y_grid_number;
subplot(4,2,3);
for ix=1:ss1(1)
    
    
    for iy=1:ss1(2)
        
        for iz=1:ss1(3)
            if (ix>ss1(1)/2)&&(iz<ss1(3)/2.1)%&&(iy<ss1(2)/2)
                     continue  
            else
            
            X(ui)=ix-0;
            Z(ui)=iz-0;
            Y(ui)=iy;
            if (X(ui)-25)^2+(Y(ui)-25)^2 > 25^2
                
                C21A(ui)=0;
                
            else
                C21A(ui)=C2(ix,iy,iz);
                ui=ui+1;
            end
            end
%             S21A(ui)=S2(ix,iy,iz);
%                 ui=ui+1;
           

%            

        end
    end
end


scatter3( (Z)*z_grid_size, (Y)*y_grid_size, (X)*x_grid_size, 15, C21A, 'filled','MarkerEdgeColor','None');
axis equal
axis off
caxis([0 60])
colormap jet

view(-112,11)
xlim([0,50])
ylim([0, 5])
zlim([0, 5])


subplot(4,2,4);
load('modeldan_garaki_hexagon_10mln_50_50_200_oe_check.mat')
bd=E_xyz;
C2=smooth3(bd);
integral_array_1= squeeze(sum(C2, 3));
% integral_array_1=rescale(integral_array_1);
surf(integral_array_1, 'EdgeColor', 'None');
axis off
axis equal
view(-90, 90)
% colorbar
% caxis([0.1*10^4 5*10^4])
caxis([0*10^3 1*10^4])



clear X Y Z
ui=1; %loop=1;
load('modeldan_garaki_10mln_single_50_50_200_oe_check.mat')
ad=E_xyz;
S2=smooth3(ad);
ss0=size(E_xyz);

ram=2.5;
L=50;


x_grid_number  = 50;        x_grid_size    = 2*ram/x_grid_number;
y_grid_number  = 50;        y_grid_size    = 2*ram/y_grid_number;
z_grid_number  = 200;       z_grid_size  = L/z_grid_number;

subplot(4,2,5);
for ix=1:ss0(1)
    
    
    for iy=1:ss0(2)
        
        for iz=1:ss0(3)
            
            if (ix>ss0(1)/2)&&(iz<ss0(3)/2.1)%&&(iy<ss0(2)/2)
                     continue  
            else
            
            X(ui)=ix-0;
            Z(ui)=iz-0;
            Y(ui)=iy;
            if (X(ui)-25)^2+(Y(ui)-25)^2 > 25^2
                
                D21A(ui)=0;
                
            else
                D21A(ui)=S2(ix,iy,iz);
                ui=ui+1;
            end
            end
%             S21A(ui)=S2(ix,iy,iz);
%                 ui=ui+1;
           

%            

        end
    end
end
hold on
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% scatter3((Z-500)*0.1,(Y-32)*0.1+3,(X-32)*0.1+3, 15, S21A, 'filled','MarkerEdgeColor','None');
% scatter3((Z-500)*0.05,(Y-25)*0.05,(X-25)*0.05, 15, S21A, 'filled','MarkerEdgeColor','None');
% scatter3((Z-0)*0.05+80, (Y-0)*0.05+2.5, X*0.05+2.5, 15, S21A, 'filled','MarkerEdgeColor','None');
% scatter3(X*0.05-0.1, Y*0.05-0.1, Z*0.05-5, 15, S21A, 'filled','MarkerEdgeColor','None');
scatter3( (Z)*z_grid_size, (Y)*y_grid_size, (X)*x_grid_size, 15, D21A, 'filled','MarkerEdgeColor','None');
% scatter3( (Z)*z_d_grid_size, (Y)*y_grid_size, (X)*x_grid_size, 15, S21A, 'filled','MarkerEdgeColor','None');
% scatter3( (Z), (Y-1), (X-1), 15, S21A, 'filled','MarkerEdgeColor','None');
% caxis([0 300])
caxis([0 60])
colormap jet
axis equal
axis off

view(-112,11)
xlim([0,50])
ylim([0, 5])
zlim([0, 5])




% 
subplot(4,2,6);
load('modeldan_garaki_10mln_single_50_50_200_oe_check.mat')
ad=E_xyz;
S2=smooth3(ad);
integral_array_1= squeeze(sum(S2, 3));
% integral_array_1=rescale(integral_array_1);
surf(integral_array_1, 'EdgeColor', 'None');
axis off
axis equal
view(-90, 90)
% colorbar
% caxis([0*10^4 4*10^4])
% caxis([0.1*10^4 5*10^4])
caxis([0*10^3 1*10^4])

% Plot the surface


clear X Y Z
ui=1; loop=1;
load('modeldan_garaki_10mln_double_50_50_200_oe_check.mat')
bd=E_xyz;
C2=smooth3(bd);
ss1=size(E_xyz);


ram=2.5;
L=50;


x_grid_number  = 50;        x_grid_size    = 2*ram/x_grid_number;
y_grid_number  = 50;        y_grid_size    = 2*ram/y_grid_number;
z_grid_number  = 200;       z_grid_size  = L/z_grid_number;

subplot(4,2,7);
for ix=1:ss1(1)
    
    
    for iy=1:ss1(2)
        
        for iz=1:ss1(3)
            if (ix>ss1(1)/2)&&(iz<ss1(3)/2.1)%&&(iy<ss1(2)/2)
                     continue  
            else
            
            X(ui)=ix-0;
            Z(ui)=iz-0;
            Y(ui)=iy;
            if (X(ui)-25)^2+(Y(ui)-25)^2 > 25^2
                
                C21A(ui)=0;
                
            else
                C21A(ui)=C2(ix,iy,iz);
                ui=ui+1;
            end
            end
%             S21A(ui)=S2(ix,iy,iz);
%                 ui=ui+1;
           

%            

        end
    end
end


scatter3( (Z)*z_grid_size, (Y)*y_grid_size, (X)*x_grid_size, 15, C21A, 'filled','MarkerEdgeColor','None');
axis equal
axis off
% caxis([0 300])
caxis([0 60])
colormap jet
% view(-79,6)

view(-112,11)
xlim([0,50])
ylim([0,5])
zlim([0,5])


subplot(4,2,8);
load('modeldan_garaki_10mln_double_50_50_200_oe_check.mat')
bd=E_xyz;
C2=smooth3(bd);
integral_array_1= squeeze(sum(C2, 3));
% integral_array_1=rescale(integral_array_1);
surf(integral_array_1, 'EdgeColor', 'None');
axis off
axis equal
view(-90, 90)
% colorbar
caxis([0*10^3 1*10^4])


