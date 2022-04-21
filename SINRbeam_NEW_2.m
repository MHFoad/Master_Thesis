clear all; 
% close all; 
clc;
rng(2755);
%constant / parameters 
f = 3500e6; %frequency
c= (3e8)/f; %wavelength
%Noise = 12; %in db
pi=3.1416;
PT= 30 ;%dBm transmit power of base station
nBS= 6;
nUS= 50;
XBS=100*rand(2,nBS);
XUS=100*rand(2,nUS);
n=4; %path losss component

antenna_orientation_deg = wrapTo180(360*rand(1,nBS));
N_array = 8;

%theta_deg =thetabus_deg;
beam_angle_deg = 0;

% distance between the base station

for i=1:nBS
   for j=1:nBS
    Dbs(i,j)=sqrt(sum((XBS(:,i)-XBS(:,j)).^2));
    thetabs_deg =atan2((XBS(:,i)),(XBS(:,j)));
   end
end
ad= mean(Dbs(:));
% distance between the users
for ii=1:nUS
   for jj=1:nUS
    Dus(ii,jj)=sqrt(sum((XUS(:,ii)-XUS(:,jj)).^2));
    thetaus_deg =atan2((XUS(:,ii)),(XUS(:,jj)));
   end
end

% distance between base station and the users
for k=1:nBS
   for m=1:nUS
    Dbu(k,m)=sqrt(sum((XBS(:,k)-XUS(:,m)).^2));
    % theta= atan2(( YXBS-YXUS),(XXBS-XXUS))
%   thetabus_deg(k,m) =atan2((XBS(2,k)-XUS(2,m)),(XBS(1,k)-XUS(1,m)));
  theta_rad(k,m) =atan2((XBS(2,k)-XUS(2,m)),(XBS(1,k)-XUS(1,m)));
  %theta_deg(k,m) =thetabus_deg; 
%   theta_rad(k,m) = thetabus_deg(k,m)/180*pi;
  beam_angle_rad = theta_rad(k,m) ;%beam_angle_deg/180*pi;
antenna_orientation_rad = antenna_orientation_deg(k)/180*pi;
antenna_power(k,m) = antenna_gain(theta_rad(k,m),beam_angle_rad,antenna_orientation_rad,N_array);
%  antenna_power_lin(k,m)=10.^(antenna_power(k,m)/10);    
PL(k,m)=32.4+10*n*log10(Dbu(k,m)/1000)+20*log10(f/1e6);%in db Free space loss
    PRx(k,m)= PT-PL(k,m)+antenna_power(k,m);
     PRx_lin(k,m)= 10.^( PRx(k,m)/10);
%      PRx_total_lin(k,m)= PRx_lin(k,m).*antenna_power_lin(k,m);
%    PRx(k,m)=10*log10(PRx_total_lin(k,m));
    
   end
end
ad= mean(Dbs(:)); PRx_lin(k,m);
PLc=32.4+10*n*log10(ad/1000)+20*log10(f/1e6);%in db Free space loss
 PRxc= PT-PLc;
 SNR=0;% considering at the cell edge
    % N=PRxc-SNR;
      N=-31;
    N_lin= 10.^(N/10);
    
 
% % Define which BSs are serving which users
[~,US_serving_BS_index] =  min(Dbu,[],1);
    

 % define which users are under each BS
for kk=1:nBS 
    US_ind_per_BS{kk} = find(US_serving_BS_index == kk);
end
 
 num_of_users_per_BS = cellfun(@length,US_ind_per_BS);
 if any(num_of_users_per_BS<1)
     error('At least one user per BS is needed. Try running again...')
 end
 
 num_of_combinations = prod(num_of_users_per_BS(num_of_users_per_BS>0));
 
 combination_cell = getCombinations(US_ind_per_BS);

 SNR_served = NaN(num_of_combinations,nBS);
 SINR_served = SNR_served;
 interference_served = SNR_served;
 US_ind_for_each_combination = SNR_served;
 for combination_ind = 1:num_of_combinations
     

     
      %served US for this combination index
     served_US = cellfun(@(x) x(combination_ind),combination_cell);
     
     US_ind_for_each_combination(combination_ind,:) = served_US;
     
     for kk = 1:nBS
         
         SNR_served(combination_ind,kk) = PRx(kk,served_US(kk)) - N;
         
         % interfering BS ind
         interfering_BS_ind=[ 1:(kk-1),(kk+1):nBS];
         P_interference_lin_tot = 0;
         for ii = interfering_BS_ind
             beam_angle_rad = theta_rad(ii,served_US(ii));
             
             antenna_orientation_rad = antenna_orientation_deg(ii)/180*pi;
             ant_gain = antenna_gain(theta_rad(ii,served_US(kk)),beam_angle_rad,antenna_orientation_rad,N_array);
 
            P_interference = PT-PL(ii,served_US(kk))+ ant_gain;
            P_interference_lin = 10.^(P_interference/10);
            P_interference_lin_tot = P_interference_lin_tot + P_interference_lin;
         end
         P_interference_tot = 10*log10(P_interference_lin_tot);
         interference_served(combination_ind,kk) = P_interference_tot;
         SINR_served(combination_ind,kk) = 10*log10(PRx_lin(kk,served_US(kk))/(N_lin+P_interference_lin_tot));      
     end
 end

average_SNR_per_combination = 10*log10(mean(10.^(SNR_served/10),2));
average_interference_per_combination = 10*log10(mean(10.^(interference_served/10),2));
average_SINR_per_combination = 10*log10(mean(10.^(SINR_served/10),2));
%h = histogram(x,'Normalization','probability')

figure
histogram(average_SNR_per_combination,10,'Normalization','probability')
hold on
histogram(average_SINR_per_combination(:),10,'Normalization','probability')
hold off
legend('SNR','SINR')
xlabel('dB')
ylabel('probability')

% figure
% histogram(average_SNR_per_combination,20)
% hold on
% histogram(average_interference_per_combination(:),20)
% histogram(average_SINR_per_combination(:),20)
% hold off
% legend('SNR','Interf','SINR')

return
