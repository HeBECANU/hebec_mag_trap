%Simple optimisesation of angle uniformity along the z-axis( i.e. axis of
%the coils) with x- field bias and shunt size
options = optimset('PlotFcns',@optimplotfval);
x0 = [0.2,0,3.4,5]; %what are the params
x = fminsearch(@cost_calc,x0,options);
%%
[theta,omega,btrap] = theta_calc(x);
max_theta=(theta(2)-theta(1)).*180/pi;
fprintf('\n max angle change %u \n \n',max_theta)
%%
x0 = [0.2,0,3.4,2.0];
[theta,omega,btrap]=theta_calc(x0)
%%
function Cost=cost_calc(params)
[theta,omega,btrap] = theta_calc(params);
if or((btrap(1,2).power_bias>40.6),(btrap(1,1).power_quad>27.9))
    burn_bias = 1e15;
else
    burn_bias = 1;
end
Cost = abs((theta(2)-theta(1)).*180/pi)*burn_bias*(exp(2*pi*400/omega)+1)*1/params(4);
end

function [theta,omega,btrap]=theta_calc(params);
[anal_out,btrap]=main(params(1:3));
omega = anal_out.trap_freq(3)*2*pi;
trap_cent=anal_out.trap_cent;
osc_amp = params(4)*1e-3/0.417/omega;
%%
pts_list = [trap_cent; trap_cent+[0 0 1].*osc_amp];
[Bmag,Bvec]=trap_eval(btrap,pts_list);
%%
k = [1,0,0]; %light field vector
k = k./norm(k);
theta = acos(Bvec*k'./Bmag);
end