clear all; clc; close all;
global D L dd r_muhit A C sa teflon_kaytarish led_kaytarish nuqta_fotonlar_soni_max
D=50;
L=25;
dd=5;
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
nuqta_fotonlar_soni_max=1;
% N=46; 
l=linspace(6, 20, 30);
h =linspace(10, 50, 41); 
tic
for i=1:30
    
    for k=1:41
        
      eff(i,k) = effektivlik_trapetsiya_working_cm(l(i), h(k)); 

      
    end
         
%       tgyubhjknl
    i
    pause(0.2)
   
end
toc
[Hh, Ll]=meshgrid(h, l);

surf(Hh,Ll, eff)


