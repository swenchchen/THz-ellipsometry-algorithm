function output=LP_filter(time,data,filter_bd,width)
%%
signal_length=numel(data);
deltaT=time(2)-time(1);
MaxFreq=1/deltaT;
freq = linspace(0,MaxFreq,signal_length);
if filter_bd~=0
    LP_filter=exp(-((freq-filter_bd).^2)./width^2);
    f_ind=find(freq<filter_bd);
    LP_filter(f_ind)=1;

    data_fd=rfft(data);
    filter_use=LP_filter(1:numel(data_fd));
    data_fd_filter=data_fd.*filter_use;
    data_td_filter=irfft(data_fd_filter,signal_length);
else
    data_td_filter=data;
end

output(1,:)=time;
output(2,:)=data_td_filter;