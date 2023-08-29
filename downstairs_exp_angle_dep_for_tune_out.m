
%add all subfolders to the path
addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path
set_up_project_path

hebec_constants

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


%% mag trap
trap_config=[];
trap_config.dlen_num=1e-3;

trap_config.v_quad=3.4; %3.4 used in 'normal trap' 14.178 A
trap_config.v_shunt=-0.1;%0.2 used in TO 

trap_config.dlen_num=0.01e-4; % the finite linear segments that are used in the helix code

trap_config.Bext=1e-4*[0,0,0];     % external bias field [T] (uniform)



% no field helix (wrong wires)
% {-3.232367,-0.076193,-0.002435} mm
% vec angle 6.174
% 1 gauss y
% vec angle  4.0723
% {-3.278436,-0.389554,-0.004517} mm


% trap_config.v_quad=2.5; %3.4 used in 'normal trap'
% trap_config.v_shunt=0.7;  

%trap_config.v_quad=3.4;%0.7; %2.5;%3.4 used in 'normal trap'
%trap_config.v_shunt=0.2;%0.85; %4;%0.75 

%ML extreme trap
% trap_config.v_quad=5.7; %start of evap
% trap_config.v_shunt=0;

%ML damping trap
%trap_config.v_quad=3.4; %3.4 used in 'normal trap'
%trap_config.v_shunt=0.0; 
%------------- END USER Config-------------------------


% construct the biquick trap
btrap=[];
btrap=biquic_trap_loops(btrap,trap_config);  % build biquic
%btrap=biquic_trap_helix(btrap,trap_config);  % build biquic trap


plot2d_opts.do=true;
plot2d_opts.type='b_scal'; %['b_scal','b_vec_comp','b_angle','hess_l1norm']
plot2d_opts.cen=[0,0,0]; %anal_out.trap_cen.pos;

plot2d_opts.rot=pi/2 * [0, 0, 0];
plot2d_opts.range=[[-1,1];
                   [-0.1,0.1]]*10e-3;
plot2d_opts.nsamp=[1,1]*50; 
plot2d_opts.zero_on_cen=false;  %does not do anything  




% if plot2d_opts.do
%     plot2d_opts.btrap=btrap;
%     visualise_2d(plot2d_opts)
%     axis equal
% end

disp("---------------------start--------------------")
%%
anal_out=[];
if solve_trapchar>0
    % evaluate trap center, 1D trap potential, trap frequencies
    % NOTE: there are multiple points of potential minimia
    verbose=1;
    anal_out=trap_characterise(anal_out,btrap,[-4e-3,0,0],solve_trapdepth,verbose);
    fprintf('trap position (x,y,z)=%s\n',sprintf('%.2g,',anal_out.trap_cen.pos))
end

disp("---------------------end of trap characterise--------------------")


%% add the nuller sensors
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
%BLACK SENSOR
% 96-148
nullr_opt.sensor(6).pos=[-45,-55,-52]*1e-3;
nullr_opt.sensor(6).dirn=[0,0,1];
nullr_opt.avg_sensors=[3,4];

nullr_opt.b_to_volt_conversion=1e4*10; %1e4 G/T  * 10v/gauss


if exist('btrap','var') && exist('nuller_curr_tmp','var')  
    curr_guess=nuller_curr_tmp;
    %curr_guess=[[1,1];[1,1];[1,1]]*0.1;
else
    curr_guess=[[-4,4];[4,4];[4,4]]*10;
end

nullr_opt.current_guess=curr_guess;
%[yellow (1), green (0) , red (2), black (3), (blue, white)]
% for the TO experiment the main data was taken with both averaged sensors (B,W) at 2.72v
nullr_opt.set_pt_v=[-3.54,3.54,-1.27,0,1.8];%[0.0,0,0,0,20.72];%%[-3.54,3.54,0,0,2.3];%[2.8,-2.8,2.8,-2.8,0.72];%[3.54,-3.54,3.54,-3.54,0.72];%
% apply nuller feedback, the value of the optimzation function is the rms voltage error 
[btrap,nuller_status]=feedback_nullr(btrap,nullr_opt);
nuller_curr_tmp=nuller_status.currents;


anal_out=trap_characterise(anal_out,btrap,[-4e-3,0,0],solve_trapdepth,verbose);

%


%%
calc_props=[];
calc_props.type='b_angle'; %['b_scal','b_vec_comp','b_angle','hess_l1norm','b_angle','b_polar_angle']
calc_props.ref_vec=[1,0,0];
calc_props.convert_to_deg=1;
calc_props.xyz_list= anal_out.trap_cen.pos;
calc_props.btrap=btrap;
scal_out=compute_scalar_property(calc_props);
vec_angle=scal_out.val;

%%
% calc_props=[];
% calc_props.type='b_polar_angle'; %['b_scal','b_vec_comp','b_angle','hess_l1norm','b_angle','b_polar_angle']
% calc_props.ref_vec_pointing=[1,0,0];
% calc_props.ref_vec_angle=[0,1,0];
% calc_props.convert_to_deg=1;
% calc_props.xyz_list= anal_out.trap_cen.pos;
% calc_props.btrap=btrap;
% scal_out=compute_scalar_property(calc_props);
% vec_polar_angle=scal_out.val




%% 2D contor/surface plot of the trap B-field

plot2d_opts.do=true;
plot2d_opts.type='b_scal'; %['b_scal','b_vec_comp','b_angle','hess_l1norm']
plot2d_opts.hess_delt=1e-7;
plot2d_opts.cen=anal_out.trap_cen.pos;

plot2d_opts.rot=pi/2 * [0, 0, 0];
plot2d_opts.range=[[-1,1];
                   [-0.1,0.1]]*20e-3;
plot2d_opts.nsamp=[1,1]*30; 
plot2d_opts.zero_on_cen=false;  %does not do anything  




% if plot2d_opts.do
%     plot2d_opts.btrap=btrap;
%     visualise_2d(plot2d_opts)
% end

%%

plot3d_opts.do=true;
plot3d_opts.show_coils=true;
plot3d_opts.type='b_scal'; %['b_scal','b_vec_comp','b_angle','hess_l1norm']
%plot3d_opts.hess_delt=1e-7;
plot3d_opts.cen=[10,-10,7]*1e-3;
%plot3d_opts.cen=[0,0,0]
plot3d_opts.range=[[-1,1];
                   [-1,1];
                   [-1,1]]*2*1e-3;
plot3d_opts.nsamp=[1,1,1]*30;           
plot3d_opts.zero_on_cen=false;    

if plot3d_opts.do
    plot3d_opts.btrap=btrap;
    visualise_3d(plot3d_opts)
end




%% 3D plot of the trap B-field
plot3d_opts.do=true;
plot3d_opts.show_coils=true;
plot3d_opts.type='b_scal'; %['b_scal','b_vec_comp','b_angle','hess_l1norm']
%plot3d_opts.hess_delt=1e-7;
plot3d_opts.cen=anal_out.trap_cen.pos
%plot3d_opts.cen=[0,0,0]
plot3d_opts.range=[[-1,1];
                   [-1,1];
                   [-1,1]]*30*1e-3;
plot3d_opts.nsamp=[1,1,1]*5;           
plot3d_opts.zero_on_cen=false;    

if plot3d_opts.do
    plot3d_opts.btrap=btrap;
    visualise_3d(plot3d_opts)
end

return 
%%

if solve_stpt>0
    st_pts=stationary_points(btrap,anal_out.trap_cent,solve_stpt);
end


%% narrow plots

plot2d_opts.do=true;
plot2d_opts.type='b_scal'; %['b_scal','b_vec_comp','b_angle','hess_l1norm']
plot2d_opts.hess_delt=1e-7;
plot2d_opts.plot_cen= anal_out.trap_cen.pos;
plot2d_opts.rot=pi/2 * [0, 0, 0];
plot2d_opts.range=[[-1,1];
                   [-1,1]]*2*1e-3;
plot2d_opts.nsamp=[1,1]*60; 
plot2d_opts.zero_on_cen=false;   
plot2d_opts.vec_plot.do=true;
plot2d_opts.vec_plot.nsamp=[1,1]*10;


%%% 2D grid trap B-field
if plot2d_opts.do
    plot2d_opts.btrap=btrap;
    visualise_2d(plot2d_opts)
end


%% plot the field angle
plot1d_opts.do=true;
plot1d_opts.type='b_angle'; %['b_scal','b_vec_comp','b_angle','hess_l1norm','b_angle','b_polar_angle']
plot1d_opts.ref_vec=[1,0,0]; %relativeto the x axis
plot1d_opts.unwrap_phase=true;
plot1d_opts.convert_to_deg=true;

plot1d_opts.plot_cen= anal_out.trap_cen.pos;
plot1d_opts.rot=pi/2 * [0, 0, 1];
%plot1d_opts.range=[-1,1]*10e-3/420; % approx scale of oscillation
plot1d_opts.range=[-1,1]*9e-6; % approx scale of beam
plot1d_opts.nsamp=[1,1]*60; 


%%% 2D grid trap B-field
if plot2d_opts.do
    plot1d_opts.btrap=btrap;
    visualise_1d(plot1d_opts)
end

%% polar angle plot
plot1d_opts.do=true;
plot1d_opts.type='b_polar_angle'; %['b_scal','b_vec_comp','b_angle','hess_l1norm','b_angle','b_polar_angle']
plot1d_opts.ref_vec_pointing=[1,0,0];
plot1d_opts.ref_vec_angle=[0,1,0];
plot1d_opts.unwrap_phase=true;
plot1d_opts.convert_to_deg=true;

plot1d_opts.plot_cen= anal_out.trap_cen.pos;
plot1d_opts.rot=pi/2 * [0, 0, 1];
%plot1d_opts.range=[-1,1]*10e-3/420; % approx scale of oscillation
plot1d_opts.range=[-1,1]*9e-6; % approx scale of beam
plot1d_opts.nsamp=[1,1]*60; 

%%% 2D grid trap B-field
if plot2d_opts.do
    plot1d_opts.btrap=btrap;
    vis_dat=visualise_1d(plot1d_opts)
end



%%
a=[];
a.xyz_list=[-3.711805903984313e-03,6.997183069798872e-06,-2.258008604850403e-04];%anal_out.trap_cen.pos;%+[0,10,0]*1e-6
% a.type='b_polar_angle'
% a.ref_vec_pointing=[1,0,0];
% a.ref_vec_angle=[0,1,0];
a.type='b_angle'
a.ref_vec=[1,0,1]./sqrt(2);%[1,0,0];
a.convert_to_deg=1;

a.btrap=btrap;
b=compute_scalar_property(a)
b.val



%% visualize the vector field
plot_opts=[]
plot_opts.cen= [-3.711805903984313e-03,6.997183069798872e-06,-2.258008604850403e-04];%anal_out.trap_cen.pos;
plot_opts.rot=pi/2 * [0, 0, 0];
plot_opts.range=[[-1,1];
                   [-1,1];
                   [-1,1]]*20*1e-6;
plot_opts.nsamp=[1,1,1]*5; 
plot_opts.zero_on_cen=true;   

plot_opts.btrap=btrap;
visualise_vector_3d(plot_opts)




%% End of code
t_end=toc(t_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_end);
disp('===================ALL TASKS COMPLETED===================');
%end