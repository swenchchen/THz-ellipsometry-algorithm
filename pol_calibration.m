function [p_output,s_output]=pol_calibration(p_input,s_input,mag_gradient,phase_shift,phase_gradient,method)
%%
p_fd=rfft(p_input(2,:));
s_fd=rfft(s_input(2,:));
maxF=1/(p_input(1,2)-p_input(1,1));
freq=linspace(0,maxF/2,numel(p_fd));
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

p_fd_correct=(p_fd-s_fd.*pol_effect);
s_fd_correct=(s_fd-p_fd.*pol_effect);

% p_fd_correct=(p_fd-s_fd.*pol_effect)./(1-pol_effect.^2);
% s_fd_correct=(s_fd-p_fd.*pol_effect)./(1-pol_effect.^2);
%%

p_output(1,:)=p_input(1,:);
s_output(1,:)=s_input(1,:);
p_output(2,:)=irfft(p_fd_correct,numel(p_input(2,:)));
s_output(2,:)=irfft(s_fd_correct,numel(p_input(2,:)));