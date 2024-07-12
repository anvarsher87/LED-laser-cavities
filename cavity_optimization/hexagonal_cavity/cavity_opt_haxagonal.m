clear all; clc; close all;
global D L dd l r_muhit A C sa teflon_kaytarish led_kaytarish nuqta_fotonlar_soni_max
D=50;
L=25;
dd=5;
l=6;
r_muhit=2.5;
%%--------------------------------crystal----------------------------------
A=dlmread('15_12_2021_BLUE_LED_ems.txt');
C=dlmread('Ce_Nd_YAG_abs.txt');
%%------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
sa=0.9*pi/2;
teflon_kaytarish=0.9;  %blue uchun
led_kaytarish=0.75;    %blue uchun
%%--------------sanash_uchun----------------------------------------------------
nuqta_fotonlar_soni_max=10000;
% N=46; 
alfa=linspace(0, pi/4, 46);
h =linspace(10, 60, 51); 
tic
for i=1:46
    
    for k=1:51
        
      eff(i,k) = effektivlik_hexagon_working(alfa(i), h(k)); 

      
    end
         
%       tgyubhjknl
    i
%     pause(0.2)
   
end
toc
% [Aa, Hh]=meshgrid(h, alfa);
% 
% surf(Aa,Hh, eff)
[Hh, Aa]=meshgrid(h, alfa*180/pi);

surf(Hh,Aa, eff)
view(0, 90)
