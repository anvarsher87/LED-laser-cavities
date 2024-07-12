function plotCircle3D_2(center,normal,radius)
% center is an array [x,y,z]
% normal is an array [nx,ny,nz]
% radius is a value
% E.g., plotCircle3D([0,0,0],[-1,0,0],1)
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
%plot3(points(1,:),points(2,:),points(3,:),'b-','linewidth',0.002);
fill3(points(3,:), points(2,:), points(1,:), 'k','LineWidth',1.5,'EdgeColor', 'k','FaceColor', 'm','FaceAlpha', '0.4');
