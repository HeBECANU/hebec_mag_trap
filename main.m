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
%move st_point calculator to trap charaterize
%change btrap from a vector of structures to a structure with cells
%break out potential from B feild for future optical & grav potential
%add in gravitational,optical potentials for hybrid trapping
%reproduce numbers from https://www.sciencedirect.com/science/article/pii/S0030401806009680?via%3Dihub
%currently trap freq way higher than reported there 
%factor out plot code
    %add options for various 2d slices

%find principle axes for trap
%more flexibility of orentation of mag cois
%finding net trap osc period with amplitude
    %use integration of potential landscape
%should estimate ideal numerical derivative step size with https://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating_point_arithmetic    



% DONE:
%code runs on arrays
%check that the Energy of a He* atom in a B feild is E=2*ub*B  
%looks to be the case http://iopscience.iop.org.virtual.anu.edu.au/article/10.1088/1464-4266/5/2/360/pdf
%close all
t_start=tic;
% ------------------START USER Config--------------
%%% Flags
verbose=1;          % graphical output

solve_trapchar=1;   % characterise trap params including freq and center
plot_2D=0;
plot_3D=0;         % solve full 3D vector B-field (takes a while)
solve_trapdepth=0;  %VERY SLOW and a can require fiddling to get working ok
solve_stpt=0;
solve_hessian=0;

%% mag trap
trap_config.v_quad=6.2; %3.4 used in 'normal trap' 14.178 A
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
% Build BiQUIC trap
btrap=biquic_trap([],trap_config);  % build biquic
anal_out=[];
if solve_trapchar>0
    % evaluate trap center, 1D trap potential, trap frequencies
    % NOTE: there are multiple points of potential minimia
    anal_out=trap_characterise(anal_out,btrap,[5e-3,0,0],solve_trapdepth,verbose);
end
%%% 2D grid trap B-field
if plot_2D
    visualise_2d(btrap,anal_out.trap_cent)
end

%%% 3D grid trap B-field
if plot_3D
    visualise_3d(btrap,anal_out.trap_cent)
end


if solve_hessian>0
    fprintf('starting hessian plot \n')
    %plot_range=[[-2,12];[-8.5,8.5]]*1e-3; %x,z
     plot_range=[[-2,12];[-8.5,8.5]]*1e-3; %x,z
    cluster_range=1e-5;
    zero_feild_thresh=1e-7; %1uT is when we say it is a zero crossing minima
    %2d derivative grid
    ngrid=100;          % 50 - med; 300 - very fine;
    delt=1e-8;
    clip=1e3;
    
    
    % grid in trap centered ref frame
    plot_range=[plot_range(1,:)+trap_cent(1);plot_range(2,:)+trap_cent(3)];%center on trap cent
    xyz_grid=[];
    [xyz_grid(:,:,:,1),xyz_grid(:,:,:,3)]=...
    meshgrid(linspace(plot_range(1,1),plot_range(1,2),ngrid),...
             linspace(plot_range(2,1),plot_range(2,2),ngrid));    % meshgrid

    xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
    H_list=num_hessian(@(x) trap_eval(btrap,x),xyz_list,delt);
    det_list=zeros(size(H_list,3),1);
    for i=1:size(H_list,3) %find the determinacy of each of these hessians
        %det_list(i)=range(eig(H_list(:,:,i)));
        det_list(i)=sum(sum(abs(H_list(:,:,i))));
    end
    %det_list(det_list>clip)=clip;
    det_list=reshape(det_list,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
    figure(22)
    clf;
    h=surf((xyz_grid(:,:,:,1)-trap_cent(1))*1e3,(xyz_grid(:,:,:,3)-trap_cent(3))*1e3,det_list,'facealpha',0.9);
    colormap(viridis())
    set(h,'LineStyle','none')
    box on;
    xlabel('X (mm)');
    ylabel('Z (mm)');
    title('Hessian')
    zlabel('RSS Grad. (T/m) ');
    set(gcf,'Color',[1 1 1]);
    pause(0.001);

    figure(24)
    clf
    contour((xyz_grid(:,:,:,1)-trap_cent(1))*1e3,(xyz_grid(:,:,:,3)-trap_cent(3))*1e3,det_list,300)
    set(gcf,'Color',[1 1 1]);
    title('Hessian')
    xlabel('X (mm)');
    ylabel('Z (mm)');
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