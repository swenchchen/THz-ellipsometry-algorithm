function output=read_step_data(time_axis_define,file_name,time_res_multiple,time_end)


data_read=load(file_name);
data_read=data_read';

data_number=round(time_res_multiple*numel(data_read(1,:)));
if round(data_number/2)==data_number/2
    data_number=data_number;
else
    data_number=data_number+1;
end
    
if numel(time_axis_define)>1
    time_axis=time_axis_define;
elseif time_end~=0;
    time_axis=linspace(0,time_end,data_number);
else
    time_axis=linspace(0,data_read(1,end),data_number);
end

data_itp(1,:)=time_axis;
data_itp(2,:)=interp1(data_read(1,:),data_read(2,:),data_itp(1,:),'spline');

nan_ind=find(isnan(data_itp(2,:)));
data_itp(2,nan_ind)=0;

output.data_read=data_read;
output.data_itp=data_itp;