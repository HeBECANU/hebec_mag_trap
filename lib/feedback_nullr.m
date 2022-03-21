function [btrap,nuller_status]=feedback_nullr(btrap,nullr_opt)
%
if ~isfield(nullr_opt,'do_opt') || nullr_opt.do_opt
    options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-6,'MaxFunEvals',1e3);
    nuller_curr_opt=fminsearch(@(x) evaluate_trail_nuller(btrap,nullr_opt,x),...
        nullr_opt.current_guess,options);
    nullr_opt.current=nuller_curr_opt;
    btrap.nullr.current=nuller_curr_opt;
else
    nullr_opt.current=nullr_opt.current_guess;
    btrap.nullr.current=nullr_opt.current_guess;
end

[~,nuller_status]=evaluate_trail_nuller(btrap,nullr_opt,nullr_opt.current);

btrap.nullr.sensor=nullr_opt.sensor;
btrap=build_nuller(btrap,nullr_opt);



end


function [rms_err,details]=evaluate_trail_nuller(btrap,nullr_opt,curr_mat)

%convert the set point from voltage to feild
set_pt_v=col_vec(nullr_opt.set_pt_v);
mag_to_voltage_conversion=nullr_opt.b_to_volt_conversion;

nullr_opt.current=curr_mat;
trial_btrap=build_nuller(btrap,nullr_opt);  % build nuller

scal_prop_opt.type='b_comp';
scal_prop_opt.btrap=trial_btrap;
pos_tmp=arrayfun(@(x) x.pos,nullr_opt.sensor,'UniformOutput',0);
scal_prop_opt.xyz_list=cat(1,pos_tmp{:});
dirn_tmp=arrayfun(@(x) x.dirn,nullr_opt.sensor,'UniformOutput',0);
scal_prop_opt.component_vec=cat(1,dirn_tmp{:});
scal_res=compute_scalar_property(scal_prop_opt);

b_val_sensors=scal_res.val;
b_val_chans=b_val_sensors;
mean_ch=0.5*(b_val_chans(nullr_opt.avg_sensors(1))...
    -b_val_chans(nullr_opt.avg_sensors(2)));
b_val_chans(nullr_opt.avg_sensors)=[];
b_val_chans=cat(1,b_val_chans,mean_ch);
chan_voltages=b_val_chans*mag_to_voltage_conversion;
rms_err=norm(chan_voltages-set_pt_v);
if nargout>1
    details=[];
    details.chan_voltages=chan_voltages;
    details.b_val_sensors=b_val_sensors;
    details.rmse=rms_err;
    details.currents=curr_mat;
end 
end