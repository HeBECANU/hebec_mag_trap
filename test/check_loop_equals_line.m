solve_trapchar=0;
solve_trapdepth=0;
solve_stpt=0;
plot2d_opts.do=true;
plot2d_opts.type='b_scal'; %['b_scal','b_vec_comp','b_angle','hess_l1norm']
plot2d_opts.hess_delt=1e-7;
plot2d_opts.plot_cen=[0.1,0.0,0];
plot2d_opts.rot=pi/2 * [0, 1, 0];
plot2d_opts.range=[[-2,2];
                   [-2,2]];
plot2d_opts.nsamp=[1,1]*100; 
plot2d_opts.zero_on_cen=false;   
plot2d_opts.vec_plot.do=true;
plot2d_opts.vec_plot.nsamp=[1,1]*15;


plot3d_opts.do=false;
plot3d_opts.type='b_scal'; %['b_scal','b_vec_comp','b_angle','hess_l1norm']
plot3d_opts.hess_delt=1e-7;
plot3d_opts.cen=[0,0,0.5];
plot3d_opts.range=[[-2,2];
                   [-2,2];
                   [-2,2]];
plot3d_opts.nsamp=[30,30,30];           
plot3d_opts.zero_on_cen=false;    

%------------- END USER Config-------------------------
t_start=tic;
%add all subfolders to the path
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));

hebec_constants
% Build BiQUIC trap
%single loop
% btrap=[];
% simple_loop.type='loop';
% simple_loop.param.radius=1;
% simple_loop.param.current=1;
% simple_loop.param.position=[0.5,0.2,0];
% simple_loop.param.rot=pi/2*[0.0,1,0];
% btrap.b_src=[simple_loop];

%quic trap
btrap=[];
loopa=[];

loopa.type='loop';
loopa.param.radius=1;
loopa.param.current=1;
loopa.param.position=[0.0,0,0];
loopa.param.rot=pi/2*[0,1,0];
loopb=loopa;
loopc=loopa;
loopb.param.current=1;
loopb.param.position=[0.1,0,0];

% linea.type='line';
% linea.param.length=0.2;
% linea.param.current=1;
% linea.param.position=[0,0.0,-0.1];
% linea.param.rot=pi/2*[0,0,0];



%loop_from_lines(1,100)
%loopb
btrap.b_src=loopa;
b_single_loop=trap_eval(btrap,[0.1,0,0]);
btrap.b_src=[loopa,loop_from_lines(1,50000,1)];
b_both_loops=trap_eval(btrap,[0.1,0,0]);

fprintf('cancels down to %g\n',b_both_loops/(2*b_single_loop))

anal_out=[];
if solve_trapchar>0
    % evaluate trap center, 1D trap potential, trap frequencies
    % NOTE: there are multiple points of potential minimia
    verbose=0;
    anal_out=trap_characterise(anal_out,btrap,[5e-3,0,0],solve_trapdepth,verbose);
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