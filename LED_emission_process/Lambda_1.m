clear; clc; close all;
A=zeros(2500,1);
b=A;
a=A;
b_sumulated=A;

 A=dlmread('BLUE_LED_emission.txt');
 i=1:2500;
pause(0.1);
hold off
for i=2:2400   
a(i)=A(i-1); 
b(i)=a(i)-a(i-1);
end;
i=360:560;
figure(1)
 plot(i,a(i),'b');
 grid on;
 figure(2)
 plot(i,1000*b(i)/31.89,'b','LineWidth',2);
%  sdfcgv
 grid on;
 hold on
  for k=1:1000
     k
     %%%%%%%%%% opredelenie lyambda %%%%%%%%%%%
  qq=rand; 
  lyamb=2500;
for i=1:1000
    if qq<=a(i)      
        lyamb=i;
        break;
    end
end
b_sumulated(lyamb)=b_sumulated(lyamb)+1;
  end
  

  sum=0;
for i=1:1000
  sum=sum+b_sumulated(i);
end

 for i=1:1000
   b_sumulated(i)=b_sumulated(i)/sum;  
 end
 hold on 
%   i=400:500;
% %  plot(i,1000*b_sumulated(i),'m.','LineWidth',1,...
% %                 'MarkerEdgeColor','m',...
% %                 'MarkerFaceColor','m',...
% %                 'MarkerSize',15);
%  grid on;
 
 nuqta_fotonlar_soni_max=1000;
% axis equal
sp=zeros(1,2500);
s=zeros(1,nuqta_fotonlar_soni_max);


for ti=2: nuqta_fotonlar_soni_max
    
    tasodifiy_son=rand;
    for i=2:2500
        if tasodifiy_son<=A(i-1)
            ftu=i;  % ftu - fotoning to'lqin uzunligi deganini bildiradi
            break;
        end
    end
    sp(ftu)=sp(ftu)+1;
    s(1,ti)=ftu;
end


plot(sp/max(sp), 'g','LineWidth',1)
grid on
xlim([360 560])
