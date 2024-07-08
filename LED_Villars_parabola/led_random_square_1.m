function [r_initial]=led_random_square_1(xc, yc1, zc, led_eni, led_boyi)
%r_i=led_random_for_square_1(x_center, LED_1y, LED_1z, led_eni, led_eni)led_eni=12.5;

%         led_ortasi=L(leds_soni,:);
%         xi=led_ortasi(1);
%         yi=led_ortasi(2);
%         zi=led_ortasi(3);


aa=rand*led_eni;
bb=rand*led_boyi;
% r=sqrt(rand);
r_ix=2*aa+xc-led_eni;
r_iy=2*bb+yc1-led_boyi;
% r_iy=2*b+yc-boyi;

r_iz=zc+0; 


r_initial=[r_ix, r_iy, r_iz];
