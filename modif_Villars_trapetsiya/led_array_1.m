% D=50;
% h1=0; %pastdagi led 
% L1=L/10;
% D1=D/20;

chizish=1;
%  for led_number=1:2
led_number=1;
dy=2;
dz=2;
birinchi_led_loc = [0 3.5 3.5]; 
%--------------------- bar 1 p--------------------------
 led_nuqta=1;
%1 dan 100
% hold on; 
xi=LED_1;
for index_z=1:10  
    for index_y=1:10
        yi=3.5 + dy*(index_y-1);
        zi=3.5 + dz * (index_z-1);
         led_nuqta_kor(led_number, led_nuqta,:) = [xi, yi, zi];
         led_nuqta=led_nuqta+1;
        if chizish ==1
        scatter3(zi, xi,yi, 10, 's',  'filled', 'MarkerEdgeColor','None',...
        'MarkerFaceColor','b')
        end
    end
end


% %---------------------------bar 2 t--------------------------
% led_nuqta 101 dan 200; 

hold on;
led_number=2;
led_nuqta=1;
xi=LED_2; %h
for index_z=1:10  
    for index_y=1:10
        yi=3.5 + dy*(index_y-1);
        zi=25+3.5 + dz * (index_z-1);
         led_nuqta_kor(led_number, led_nuqta,:) = [xi, yi, zi];
         led_nuqta=led_nuqta+1; 
        if chizish ==1
%         scatter3(zi,xi,yi, 10, 's',  'filled', 'MarkerEdgeColor','None',...
%         'MarkerFaceColor','k')
        end
    end
end
%    led_joy=F(led_number, led_nuqta,:);

% end
% xlabel('x')
% ylabel('y')
% zlabel('z')