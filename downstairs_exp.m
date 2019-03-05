% function [anal_out,btrap]=main(params)
% quad_val=params(3);
% shunt_val=params(1);
% x_bias=params(2);
%% BiQUIC trap magnetic field calculator
% D K Shin
% refer to https://doi.org/10.1103/PhysRevA.35.1535 for expressions used
% B M Henson

%creates a structure bfield that specifies the currents

%features currently implemented
%find trap depth 3 diffrent methods
    %2d thresholding/region detection
    %3d thresholding/region detection
    %saddle point finder using opt of derivative landscape
%plots in 2,3d
%finds trap frequency in cartesian axis


% Known Bugs/Errors
% trap freq not reproduced
%trap ratio not reproduced


% TODO:
%break out potential from B feild for future optical & grav potential
%add in gravitational,optical potentials for hybrid trapping
%reproduce numbers from https://www.sciencedirect.com/science/article/pii/S0030401806009680?via%3Dihub
%find principle axes for trap
%finding net trap osc period with amplitude
    %use integration of potential landscape


% DONE:
%code runs on arrays
%check that the Energy of a He* atom in a B feild is E=2*ub*B  
%looks to be the case http://iopscience.iop.org.virtual.anu.edu.au/article/10.1088/1464-4266/5/2/360/pdf
%close all

%%% CORD system
% x +ve towards LVIS
% Y -ve towards desk/ source
% Y axis along ZS

t_start=tic;
% ------------------START USER Config--------------

solve_trapchar=1;
solve_trapdepth=0;
solve_stpt=0;

plot2d_opts.do=true;
plot2d_opts.type='b_scal'; %['b_scal','b_vec_comp','b_angle','hess_l1norm']
plot2d_opts.hess_delt=1e-7;
plot2d_opts.plot_cen=[0.0,0.0,0e-3];
plot2d_opts.rot=pi/2 * [0, 0, 0];
plot2d_opts.range=[[-1,1];
                   [-1,1]]*10e-3;
plot2d_opts.nsamp=[1,1]*60; 
plot2d_opts.zero_on_cen=false;   
plot2d_opts.vec_plot.do=true;
plot2d_opts.vec_plot.nsamp=[1,1]*10;


plot3d_opts.do=true;
plot3d_opts.show_coils=true;
plot3d_opts.type='b_scal'; %['b_scal','b_vec_comp','b_angle','hess_l1norm']
plot3d_opts.hess_delt=1e-7;
plot3d_opts.cen=[0,0,3];
plot3d_opts.range=[[-1,1];
                   [-1,1];
                   [-1,1]]*200*1e-3;
plot3d_opts.nsamp=[1,1,1]*10;           
plot3d_opts.zero_on_cen=false;    


%% mag trap
trap_config.v_quad=3.4; %3.4 used in 'normal trap' 14.178 A
trap_config.v_shunt=0.2;%0.2 used in TO 

trap_config.Bext=1e-4*[0,0,0];     % external bias field [T] (uniform assumption)

% trap_config.v_quad=3.4; %3.4 used in 'normal trap'
% trap_config.v_shunt=0.0;  

%ML extreme trap
% trap_config.v_quad=5.7; %start of evap
% trap_config.v_shunt=0;

%ML damping trap
%trap_config.v_quad=0.25; %3.4 used in 'normal trap'
%trap_config.v_shunt=0.75; 
%------------- END USER Config-------------------------

%add all subfolders to the path
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));

hebec_constants


if exist('btrap','var') && isfield(btrap,'nullr') && isfield(btrap.nullr,'current')
    curr_guess=btrap.nullr.current;
else
    curr_guess=[[1,1];[1,1];[1,1]];
end

btrap=[];

%btrap=biquic_trap_loops(btrap,trap_config);  % build biquic
btrap=biquic_trap_helix(btrap,trap_config);  % build biquic


nullr_opt=[];
nullr_opt.radius=500e-3;
nullr_opt.inset=90e-3;
nullr_opt.do_opt=true;
nullr_opt.sensor=[];


%%% CORD system
% x +ve towards LVIS
% Y -ve towards desk/ source
% Y axis along ZS
% YELLOW SENSOR
nullr_opt.sensor(1).pos=[90,-30,0]*1e-3;
nullr_opt.sensor(1).dirn=[-1,0,0];
%GREEN SENSOR
% x 94±3mm away from LVIS
nullr_opt.sensor(2).pos=[-94,30,55]*1e-3;
nullr_opt.sensor(2).dirn=[1,0,0];
%BLUE SENSOR
% x 51±5mm away from LVIS
nullr_opt.sensor(3).pos=[-51,55,0]*1e-3;
nullr_opt.sensor(3).dirn=[0,-1,0];
% WHITE SENSOR 
% x 43mm towards LVIS
% displacement along Y
nullr_opt.sensor(4).pos=[43,-55,0]*1e-3; 
nullr_opt.sensor(4).dirn=[0,1,0];
%RED SENSOR
% x 37mm towards LVIS
% z 45 up
nullr_opt.sensor(5).pos=[55,55,45]*1e-3;
nullr_opt.sensor(5).dirn=[0,0,-1];
%unknown colour
% 96-148
nullr_opt.sensor(6).pos=[-45,-55,-52]*1e-3;
nullr_opt.sensor(6).dirn=[0,0,1];

nullr_opt.avg_sensors=[2,3];
nullr_opt.current_guess=curr_guess;
btrap=feedback_nullr(btrap,nullr_opt);



anal_out=[];
if solve_trapchar>0
    % evaluate trap center, 1D trap potential, trap frequencies
    % NOTE: there are multiple points of potential minimia
    verbose=0;
    anal_out=trap_characterise(anal_out,btrap,[-5e-3,0,0],solve_trapdepth,verbose);
end
%%% 2D grid trap B-field
if plot2d_opts.do
    plot2d_opts.btrap=btrap;
    visualise_2d(plot2d_opts)
end

%%% 3D grid trap B-field
if plot3d_opts.do
    plot3d_opts.btrap=btrap;
    visualise_3d(plot3d_opts)
end


if solve_stpt>0
    st_pts=stationary_points(btrap,anal_out.trap_cent,solve_stpt);
end




%% End of code
t_end=toc(t_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_end);
disp('===================ALL TASKS COMPLETED===================');
%end