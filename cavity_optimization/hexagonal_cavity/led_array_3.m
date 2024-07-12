% function led_joy=led_array_1 (LED_1_x, LED_2_x, y_boyicha_surish)

chizish=1;
%  for led_number=1:2
led_number=1;
dy=2;
dz=2;
% birinchi_led_loc = [LED_1_x y_boyicha_surish+3.5 3.5]; 
birinchi_led_loc = [LED_1 ly+3.5 3.5]; 
%--------------------- bar 1 p--------------------------
 led_nuqta=1;
%1 dan 100
% hold on; 
xi=LED_1;
for index_z=1:10  
    for index_y=1:10
        yi=ly+3.5 + dy*(index_y-1);
        zi=3.5 + dz * (index_z-1);
         led_nuqta_kor(led_nuqta,:) = [xi, yi, zi];
         led_nuqta=led_nuqta+1;
%         if chizish ==1
        scatter3(zi,xi,yi, 20, 'o',  'filled', 'MarkerEdgeColor','None',...
        'MarkerFaceColor','b')
%         end
    end
end


%--------------------- bar 1 p--------------------------
 led_nuqta=1;
%1 dan 100
% hold on; 
xi=LED_2;
for index_z=1:10  
    for index_y=1:10
        yi=ly+3.5 + dy*(index_y-1);
        zi=D/2+3.5 + dz * (index_z-1);
         led_nuqta_kor(led_nuqta,:) = [xi, yi, zi];
         led_nuqta=led_nuqta+1;
%         if chizish ==1
        scatter3(zi,xi,yi, 20, 'o',  'filled', 'MarkerEdgeColor','None',...
        'MarkerFaceColor','b')
%         end
    end
end

% end
% xlabel('x')
% ylabel('y')
% zlabel('z')
% end