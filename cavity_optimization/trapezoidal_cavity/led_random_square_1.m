function [r_initial]=led_random_square_1(xc, yc, zc, led_eni, led_boyi)
%r_i=led_random_for_square_1(x_center, LED_1y, LED_1z, led_eni, led_eni)led_eni=12.5;


aa=rand*led_eni;
bb=rand*led_boyi;
% r=sqrt(rand);
r_iz=2*aa+zc-led_eni;
r_iy=2*bb+yc-led_boyi;
% r_iy=2*b+yc-boyi;

r_ix=xc+0; 


r_initial=[r_ix, r_iy, r_iz];
