function btrap=feedback_nullr(btrap,nullr_opt)
%

options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-6,'MaxFunEvals',1e3);
nuller_curr_opt=fminsearch(@(x) evaluate_trail_nuller(btrap,nullr_opt,x),...
    nullr_opt.current_guess,options);
nullr_opt.current=nuller_curr_opt;
btrap.nullr.current=nuller_curr_opt;
btrap.nullr.sensor=nullr_opt.sensor;
btrap=build_nuller(btrap,nullr_opt);

end


function rms_err=evaluate_trail_nuller(btrap,nullr_opt,curr_mat)
nullr_opt.current=curr_mat;

trial_btrap=build_nuller(btrap,nullr_opt);  % build nuller

scal_prop_opt.type='b_comp';
scal_prop_opt.btrap=trial_btrap;
pos_tmp=arrayfun(@(x) x.pos,nullr_opt.sensor,'UniformOutput',0);
scal_prop_opt.xyz_list=cat(1,pos_tmp{:});
dirn_tmp=arrayfun(@(x) x.dirn,nullr_opt.sensor,'UniformOutput',0);
scal_prop_opt.component_vec=cat(1,dirn_tmp{:});
scal_res=compute_scalar_property(scal_prop_opt);
mean_ch=0.5*(scal_res.val(nullr_opt.avg_sensors(1))...
    -scal_res.val(nullr_opt.avg_sensors(2)));
comb_ch=scal_res.val';
comb_ch(nullr_opt.avg_sensors)=[];
comb_ch=[comb_ch,mean_ch];
rms_err=norm(comb_ch);
end