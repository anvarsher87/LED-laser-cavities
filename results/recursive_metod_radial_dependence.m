% % program opower;

% clear ; clc; close all;

%experimental_points


%--------------------Ce: Nd: YAG-------------------------------------------------

% load('qayta_trap_2.5_optical_engineering.mat')
% load('qayta_hexagon_2.5_optical_engineering.mat')
% load('qayta_s_par_optical_engineering.mat')
load('qayta_d_par_optical_engineering.mat')


%--------------------------------------------------------------------------
% load('qayta_d_par_rad_1.5.mat')
% load('qayta_d_par_rad_2.mat')
% 
% 
% load('qayta_d_par_rad_3.mat')
% load('qayta_d_par_rad_3.5.mat')
pump_efficiency=eff
% pump_efficiency=0.702
% fgvhbj
% zdxcfvgbh
quantum_yeald_ce = 0.8;
% eff
% pump_efficiency=eff; %for Nd_KGW (3 %)
pump_efficiency =  quantum_yeald_ce*pump_efficiency; 
internal_loss_c= 0.003; % for Ce:Nd:YAG
refractive_index=1.8; 
% ----------------describing pumping rate-----------------------------------
Max_power=300;      % in W
plank=6.62607004e-34;   % metre^2 kg / s
ave_wavelength_pump=4.6e-7;     % in metres
laser_wavelength=1.064e-6;  % in metres
% --------------------------------------------------------------------------
laser_level_lifetime=230.0e-06;  tau=laser_level_lifetime; % in sec
Nt=1.38+20; sigma_em=2.8e-19;



% load('modeldan_olingan_yutlish_taqsimoti.mat')
size_of_array = size(E_xyz);

% T=zeros(size_of_array)+300; % for Nd:YAG
T=zeros(size_of_array)+300; % for Ce:Nd:YAG


%----------------describing an active medium-------------------------------
Length_of_crystal=5;    % in sm
radius_of_crystal=0.25;    % in sm
     % of the active medium    
%----------------describing a cavity parameters----------------------------
length_of_cavity = 12;    % in sm, this is free space inside the cavity 
ref_coef_of_exit_mirror=0.94;  % reflection villars 0.96
tran_coef_of_exit_mirror1=0.002; % transmission
tran_coef_of_exit_mirror2=1-ref_coef_of_exit_mirror; % transmission
%----------------------Griding---------------------------------------------

x=linspace(-radius_of_crystal,radius_of_crystal, size_of_array(1)); y=x;  % this devision is comes from T calcualtion
z=linspace(0,Length_of_crystal, size_of_array(3));

[X,Y,Z] = meshgrid(x,y,z);


kristal_turi='bir narsa';



% flow_rate=100; 
% heatable_c = 0.15; % for Ce: 0.125; for Nd: 0.055; 




% alfa=1;  % ayutilish koefisinti
% str = {'Flow Rate: 50','Pumping Efficiency: 0.15','Internal loss: 0.003', 'Heatable Part: 0.12'};

                % in sec
% total_numer_of_pump_photons=Pump_input*time/(plank*3e+8/ave_wavelength_pump);
% pump_size=0.5*(2*radius_of_crystal); % all in cm and focal spot is equal to diameter of the rod
% wp=pump_size; 

% Pump dis and its normalization
% 
% for k=1:length(z)
%     for i=1:length(x)
%         for j=1:length(y)
%             if x(i)^2+y(j)^2>radius_of_crystal^2
%                 pump_dist_x_y_z(i,j,k)=0;
%             else
%                 pump_dist_x_y_z(i,j,k)=exp(-2.*(x(i).^2+y(j).^2)./wp).*exp(-alfa.*z(k));
%             end
%         end
%     end
% end
pump_dist_x_y_z = E_xyz;
pump_dist_x_y=trapz(z,pump_dist_x_y_z,3);
sum_pump = trapz(y,trapz(x,pump_dist_x_y,2));
pump_dist_x_y_z_normalized=pump_dist_x_y_z./sum_pump;

%------------------------------------------------------------------------------
% Ve  - effective volume of the mode in the cavity
w0 = 0.99*radius_of_crystal;
ll=linspace(0,length_of_cavity, size_of_array(3));
for k=1:length(ll)
    for i=1:length(x)
        for j=1:length(y)
            if x(i)^2+y(j)^2>w0^2
                mode_dist_x_y_z(i,j,k)=0;
            else
                mode_dist_x_y_z(i,j,k)=exp(-2.*(0*x(i).^2+0*y(j).^2)./w0).*(0*ll(k)+1);
            end
            
        end
    end
end
mode_dist_x_y=trapz(ll,mode_dist_x_y_z,3);
sum_mode_cavity = trapz(y,trapz(x,mode_dist_x_y,2));

%------------------------------------------------------------------------------
% Ve  - effective volume of the mode inside the rod
for k=1:length(z)
    for i=1:length(x)
        for j=1:length(y)
            if x(i)^2+y(j)^2>w0^2
                mode_dist_x_y_z(i,j,k)=0;
            else
                mode_dist_x_y_z(i,j,k)=exp(-2.*(0*x(i).^2+0*y(j).^2)./w0).*(0*z(k)+1);
            end
        end
    end
end
mode_dist_x_y=trapz(z,mode_dist_x_y_z,3);
sum_mode_rod = refractive_index*trapz(y,trapz(x,mode_dist_x_y,2));
%--------------------------------------------------------------------------

Ve=sum_mode_rod+sum_mode_cavity; % Effective volume 

% --------------------------------------------------------------------------
loss_at_mirror1=-log(1-tran_coef_of_exit_mirror1);
loss_at_mirror2=-log(1-tran_coef_of_exit_mirror2);
a=0.002; Li=internal_loss_c*Length_of_crystal;
internal_loss=-log(1-a)-log(1-Li);   % Zvelto 7.2.7
loss=internal_loss+(loss_at_mirror1+loss_at_mirror2)/2;
% --------------------------------------------------------------------------
Optical_length = (refractive_index)*Length_of_crystal+1*length_of_cavity;
photon_lifetime=Optical_length/loss/3e+10; tauc=photon_lifetime;


P_input=linspace(0, Max_power,50);
tic;
for step=1:length(P_input)
    total_numer_of_pump_photons=P_input(step)*1/(6.62607004e-34*3e+8/ave_wavelength_pump);
    cavity_photons_0=total_numer_of_pump_photons*pump_efficiency*tauc;
    Rate_pump =  pump_efficiency * total_numer_of_pump_photons.*pump_dist_x_y_z./sum_pump;
    
    %T=temperature_distribution(X,Y,Z, radius_of_crystal, P_input(step)*pump_efficiency,pump_size,mh,alfa,K,Tc, h);
    P_input_c=P_input(step);
    
%     Teplorod_end_pumping_P_variable1;
%     
%     aylantirgich;
%     
%     T=T_x_y_z;
%     Tmax=max(max(max(T)))
%     hgfhfh
    step
    %     f=ab_cross_section(T)./sm_cross_section(T);
    f=0; %sigma_absorption_thermal_population (T)/sigma_em;
    
    lup=0;
    
    while lup < 50 %cavity_photons-cavity_photons0
        
        
        under_integral_x_y_z = (Rate_pump.*tau.*(1+f)-f.*Nt)./((1+f).*tau+Ve./(cavity_photons_0.*3e+10.*(2.8e-19).*mode_dist_x_y_z));
        
        under_integral_x_y_z(isnan(under_integral_x_y_z))=0;
        under_integral_x_y=trapz(z,under_integral_x_y_z,3);
        cavity_photons = photon_lifetime*trapz(y,trapz(x,under_integral_x_y,2));
        lup=lup+1;
        cavity_photons_0=cavity_photons;
        %         pause
    end
    
    if cavity_photons<0
        cavity_photons=0;
    end
    
    out_cavity_photons(step)=cavity_photons;
    
    
end



Pout=(loss_at_mirror2*3e+10/2/Optical_length) * out_cavity_photons  * (plank*3e+8/1.064e-6);

% if kristal_turi=="Nd:YAG"
%    
%     plot(P_input,Pout,'-.m','LineWidth', 2) 
%     
% elseif kristal_turi=="Ce:Nd:YAG"
%   
%     plot(P_input,Pout,'MarkerFaceColor','#EDB120','LineWidth', 2) 
% end

% plot(P_input,Pout,'MarkerFaceColor','#EDB120','LineWidth', 2) 
plot(P_input,Pout,'LineWidth', 2) 
xlabel('The pump power [W]','FontSize',14,'FontWeight','bold');
ylabel('The output power [W]','FontSize',14,'FontWeight','bold');
hold on;
% plot(Power_input,P_out_total,'ro', 'LineWidth', 2)
% xlabel('Solar input power [W]','FontSize',18,'FontWeight','bold');
% ylabel('Laser output power [W]','FontSize',18,'FontWeight','bold');
grid on;
% box(axes, 'on');
% Set the remaining axes properties
% set(axes,'FontSize',14,'FontWeight','bold');
% gridn;
% title( o'Ce:Nd:YAG')
% dim = [0.2 0.5 0.3 0.3];

%  annotation('textbox',...
%     [0.151041670819124 0.74832216967177 0.13854166251421 0.135570568364432],...
%     'String',{'The rod material is', kristal_turi, 'Flow Rate =' flow_rate,'cm^3/s','Pumping Efficiency =' 100*pump_efficiency, '%' },...
%     'LineStyle',':',...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'FitBoxToText','on',...
%     'BackgroundColor',[0.9 0.9 0.9]);
sound(sin(1:100));














