solve_trapchar=0;

plot_2d_opt.cen=[0,0,0.5];
plot_2d_opt.rot=pi/2 * [0, 0, 0];
plot_2d_opt.range=[[-2,2];
                   [-2,2]];
plot_2d_opt.zero_on_cen=false;              
%------------- END USER Config-------------------------

%add all subfolders to the path
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));

hebec_constants
% Build BiQUIC trap
btrap=[];
simple_loop.type='loop';
simple_loop.param.radius=1;
simple_loop.param.current=1;
simple_loop.param.position=[0,0,0];
simple_loop.param.rot=pi/2*[1,0,0];
btrap.b_src=[simple_loop];

anal_out=[];
if solve_trapchar>0
    % evaluate trap center, 1D trap potential, trap frequencies
    % NOTE: there are multiple points of potential minimia
    anal_out=trap_characterise(anal_out,btrap,[5e-3,0,0],solve_trapdepth,verbose);
end
%%% 2D grid trap B-field
if plot_2D
    visualise_2d(btrap,plot_2d_opt)
end

%%% 3D grid trap B-field
if plot_3D
    visualise_3d(btrap,[0,0,0])
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