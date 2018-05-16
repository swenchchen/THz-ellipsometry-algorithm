function output=polarizer_combine(time,data1,data2,mag_gradient,phase_shift,phase_gradient,method)
maxF=1/(time(2)-time(1));
data1_fd=rfft(data1);
data2_fd=rfft(data2);
freq=linspace(0,maxF/2,numel(data1_fd));
high_freq_region=find(freq>3);

% interpolation model
load('pol_cal_data.mat');
mag_itp=interp1(freq_use,pol_mag,freq,'spline');
phase_itp=interp1(freq_use,pol_phase,freq,'spline');
if method==0
    pol_effect=mag_itp.*exp(1i*phase_itp);
elseif method==1
% linear model
    pol_effect=(mag_gradient*freq).*exp(1i*(phase_shift*pi/180+phase_gradient*pi*freq/180));
end

pol_effect(high_freq_region)=0;
combine_fd=data1_fd+data2_fd.*pol_effect;
output=irfft(combine_fd,numel(time));