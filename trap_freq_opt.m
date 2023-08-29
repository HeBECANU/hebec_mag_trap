


%% mag trap
cost_func = @(V) prod(mag_trap_freq(V(2),V(1)).^[1/3,1/3,-2/3]);
mag_trap_freq(3.0,4)
cost_func([3.0,0.4])

options = optimset('PlotFcns',@optimplotfval,'TolX',1e-9);
x0=[3.0,0.2];
[x_opt]=fminsearch(@(x) cost_func(x),x0,options);

% mag_trap_freq(6.228007325020534e-02,4.912118337797798e+00)

function trap_freq = mag_trap_freq(Vsh,Vq)

hebec_constants
trap_config=[];
trap_config.dlen_num=1e-4;
trap_config.dlen_num=1e-4; % the finite linear segments that are used in the helix code

trap_config.Bext=1e-4*[0,0,0];     % external bias field [T] (uniform)


trap_config.v_quad=Vq; %3.4 used in 'normal trap'
trap_config.v_shunt=Vsh;  

%ML extreme trap
% trap_config.v_quad=5.7; %start of evap
% trap_config.v_shunt=0;

%ML damping trap
% trap_config.v_quad=0.25; %3.4 used in 'normal trap'
% trap_config.v_shunt=0.75; 
%------------- END USER Config-------------------------


% construct the biquick trap
btrap=[];
btrap=biquic_trap_loops(btrap,trap_config);  % build biquic

solve_trapdepth=0;
anal_out=[];
verbose = 0;
if Vsh<0 || Vq<0.5 || Vsh>4 || Vq>5
    trap_freq = [nan nan nan];
else
anal_out=trap_characterise(anal_out,btrap,[-4e-3,0,0],solve_trapdepth,verbose);
trap_freq = anal_out.trap_freq;
end
end