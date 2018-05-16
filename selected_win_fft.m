function output=selected_win_fft(data,cent_point,win_length,winPara,zerofill)

% win_length must be even!!

time=data(1,:);
signal=data(2,:);

maxF=1/(time(2)-time(1));
if winPara~=0;
    if win_length/2<cent_point
        signal=data(2,(cent_point-win_length/2):(cent_point+win_length/2-1));
    else
        signal=[linspace(0,0,win_length/2-cent_point),data(2,1:(cent_point+win_length/2))];
    end
    window_ftn=chebwin(win_length,winPara);
    % make passing window more transparent and blocking part more opaque
    window_ftn=0.5*sin(window_ftn*pi-pi/2)+0.5;
    signal_win=window_ftn'.*signal;
    if win_length/2<cent_point
        signal_out=[linspace(0,0,cent_point-win_length/2-1),signal_win,linspace(0,0,numel(data(2,:))-cent_point-win_length/2+1)];
    else
        signal_out=[signal_win(win_length/2-cent_point+1:end),linspace(0,0,numel(data(2,:))-cent_point-win_length/2)];
    end
    time_out=time;
    if zerofill>numel(signal_out)
        signal_out=[signal_out,linspace(0,0,zerofill-numel(signal_out))];
        time_out=linspace(time(1),(time(end)-time(1))/(numel(time)-1)*(numel(signal_out)-1),numel(signal_out));
    end
    
    signal_fd=rfft(signal_out);
    
elseif winPara==0
%     if zerofill<numel(signal)
%         time_out=time(1:zerofill);
%         signal_out=signal(1:zerofill);
%         signal_fd=rfft(signal_out);
%     else
        signal_out=signal;
        time_out=time;
        signal_fd=rfft(signal_out);
%     end
end

freq=linspace(0,maxF/2,numel(signal_fd));
output.td(1,:)=time_out;
output.td(2,:)=signal_out;
output.fd(1,:)=freq;
output.fd(2,:)=signal_fd;




