%% Code explanation

% This code is used for processing THz ellipsometry data using the
% algorithm avaiable in the paper "Robust and Accurate Terahertz
% Time-Domain Spectroscopic Ellipsometry " by X. Chen et. al.

% The code is written for data acquired from Menlo Systems K15 THz-TDS system

% Some sub-functions are also attached for necessary data processing
%% Fill the experimental parameters
clear all
theta_i=67.1; % incident angle
P1=45;  % direction of P1
P2_align1=16;  % direction of P2
%% Read, interpolarization and HP filter
% read data from files with interpolation in the time-domain (TD)
s_load=read_step_data(0,'cond_side_s.txt',2.5,80);
p_load=read_step_data(s_load.data_itp(1,:),'cond_side_p.txt',2.5,0);
align_load=read_step_data(s_load.data_itp(1,:),'cond_side_pos16.txt',2.5,0); % align data was measured at the beta direction

% Filter the low frequency offset before applying window function
s_itp=filter_slow(s_load.data_itp(1,:),s_load.data_itp(2,:),0.1,0.04);
p_itp=filter_slow(p_load.data_itp(1,:),p_load.data_itp(2,:),0.1,0.04);
align_itp=filter_slow(align_load.data_itp(1,:),align_load.data_itp(2,:),0.1,0.04);

%% align sample data
% find the center point for generating a time-domain window function
cent_point=round(mean([find(s_itp(2,1:1800)==max(s_itp(2,1:1800))),find(s_itp(2,1:1800)==min(s_itp(2,1:1800)))]));
% window TD data and do fft
s_process=selected_win_fft(s_itp,cent_point,1280,200,0);
p_process=selected_win_fft(p_itp,cent_point,1280,200,0);
align_process=selected_win_fft(align_itp,cent_point,1280,200,0);
% extract data for applying the algorithm
s_itp_use=s_process.td;
p_itp_use=p_process.td;
align_itp_use=align_process.td;

phase_step=0.5*pi/180; % degree at 1THz
time_shift_range=0.15; %+-range in ps
phase_shift_range=2*pi*time_shift_range;
s_fd=rfft(s_itp_use(2,:));
p_fd=rfft(p_itp_use(2,:));
maxF=1/(s_itp_use(1,2)-s_itp_use(1,1));
freq=linspace(0,maxF/2,numel(s_fd));
align_fd=rfft(align_itp_use(2,:));
shift_array=-phase_shift_range:phase_step:phase_shift_range;
j=1;
clear p_shift_td align_shift_td align_shift_fd
for phase_shift=shift_array
    p_shift_fd=p_fd.*exp(1i*freq*phase_shift);
    p_shift_td(j,:)=irfft(p_shift_fd,numel(p_itp_use(1,:)));
    align_shift_fd(j,:)=align_fd.*exp(1i*freq*phase_shift);
    align_shift_td(j,:)=irfft(align_shift_fd(j,:),numel(align_itp_use(1,:)));
    current_filter=LP_filter(s_itp_use(1,:),align_shift_td(j,:),1.2,0.5);
    align_shift_td_filter(j,:)=current_filter(2,:);
    j=j+1;
end

time_res=s_itp_use(1,2)-s_itp_use(1,1);
compare_region=(cent_point-round(5/time_res)):(cent_point+round(5/time_res))-80;
fd_region=find(freq>0.2 & freq<1);

P3_error=0; % correct the orientation of P3 if necessary
p_shift_td=p_shift_td./cos((45+P3_error)*pi/180)*cos(45*pi/180);
s_itp_use=s_itp_use./cos((45-P3_error)*pi/180)*cos(45*pi/180);
align_shift_fd=align_shift_fd./cos((45-P2_align1+P3_error)*pi/180)*cos((45-P2_align1)*pi/180);

clear vs error_vs_2shift
for m=1:numel(p_shift_td(:,1))
    [p_correct_itp,s_correct_itp]=pol_calibration([p_itp_use(1,:);p_shift_td(m,:)],s_itp_use,0.0521,90.19,9.675,0);
    vector_sum1=sin(P2_align1*pi/180)*s_correct_itp(2,:)+cos(P2_align1*pi/180)*p_correct_itp(2,:); %vector sum at align direction
    vector_sum1=vector_sum1/cos(45*pi/180)*cos((45-P2_align1)*pi/180); %calibrate according to the align direction
    vector_sum2=cos(P2_align1*pi/180)*s_correct_itp(2,:)-sin(P2_align1*pi/180)*p_correct_itp(2,:);%vector sum at align perpendicular direction
    vector_sum2=vector_sum2/cos(45*pi/180)*cos((45+P2_align1)*pi/180); %calibration
    vector_sum=polarizer_combine(p_itp_use(1,:),vector_sum1,vector_sum2,0.0521,90.19,9.675,0); %polarizer calibration
    vector_sum_filter=LP_filter(s_itp_use(1,:),vector_sum,1.5,0.5);
    vector_sum_fd=rfft(vector_sum);
    vs(m,:)=vector_sum_filter(2,:);
    for n=1:numel(align_shift_td_filter(:,1))
        error_vs_2shift(m,n)=mean(abs(vector_sum_fd(fd_region)-align_shift_fd(n,fd_region)));
    end
end

[m_min,n_min]=find(error_vs_2shift==min(min(error_vs_2shift)));
% check the error matrix
figure
imagesc(error_vs_2shift)
% check the best fitting result
figure
plot(vs(m_min,:))
hold on
plot(align_shift_td_filter(n_min,:),'r--')

% output the corresponding calibrated p and s components
[p_cal,s_cal]=pol_calibration([p_itp_use(1,:);p_shift_td(m_min,:)],s_itp_use,0.0521,90.19,9.675,0);