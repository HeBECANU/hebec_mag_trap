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
%reproduce numbers from https://www.sciencedirect.com/science/article/pii/S0030401806009680?via%3Dihub
%Reduce function calls to trap_eval in trap_characterise? 
%currently trap freq way higher than reported there 
%factor out plot code
    %add options for various 2d slices
%add in gravitational,optical potentials for hybrid trapping
%find principle axes for trap
%more flexibility of orentation of mag cois
%finding net trap osc period with amplitude
    %use integration of potential landscape
%should estimate ideal numerical derivative step size with https://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating_point_arithmetic    
global const  


% DONE:
%code runs on arrays
%check that the Energy of a He* atom in a B feild is E=2*ub*B  
%looks to be the case http://iopscience.iop.org.virtual.anu.edu.au/article/10.1088/1464-4266/5/2/360/pdf
close all
t_start=tic;
% ------------------START USER Config--------------
%%% Flags
verbose=1;          % graphical output

solve_trapchar=1;   % characterise trap params including freq and center
solve_2D=1;
solve_3D=0;         % solve full 3D vector B-field (takes a while)
solve_trapdepth=0;  %VERY SLOW and a can require fiddling to get working ok
solve_stpt=1;
solve_hessian=0;

%%% mag trap
v_quad=7.2; %3.4 used in 'normal trap'
v_shunt=0;  
Bext=1e-4*[0.0,0,0];     % external bias field [T] (uniform assumption)

%------------- END USER Config-------------------------

%%Add dependencies
addpath('Colormaps') 
constants
% Build BiQUIC trap
btrap=biquic_trap(v_quad,v_shunt,Bext);  % build biquic

if solve_trapchar>0
    % evaluate trap center, 1D trap potential, trap frequencies
    % NOTE: there are multiple points of potential minimia
    [f0,trap_cent,B_cent]=trap_characterise(btrap,5e-3,solve_trapdepth,verbose);
end

%%% 2D grid trap B-field
if solve_2D
    trap_profile_2d(btrap,trap_cent)
end

%%% 3D grid trap B-field
if solve_3D
    visualise_3D(btrap,trap_cent)
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
    stationary_points(btrap,trap_cent,solve_stpt)
%     fprintf('starting derivative plot \n')
%     %plot_range=[[-2,12];[-8.5,8.5]]*1e-3; %x,z
%      plot_range=[[-2,12];[-8.5,8.5]]*1e-3; %x,z
%     min_cluster_range=1e-5;
%     st_cluster_range=1e-5; %saddle pt cluster range
%     zero_feild_thresh=1e-7; %1uT is when we say it is a zero crossing minima
%     %2d derivative grid
%     ngrid=100;          % 50 - med; 300 - very fine;
%     delt=1e-10;
% 
%     
%     
%     % grid in trap centered ref frame
%     plot_range=[plot_range(1,:)+trap_cent(1);plot_range(2,:)+trap_cent(3)];%center on trap cent
%     xyz_grid=[];
%     [xyz_grid(:,:,:,1),xyz_grid(:,:,:,3)]=...
%     meshgrid(linspace(plot_range(1,1),plot_range(1,2),ngrid),...
%              linspace(plot_range(2,1),plot_range(2,2),ngrid));    % meshgrid
% 
%     xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
%     Bgrad_list=num_grad(@(x) trap_eval(btrap,x),xyz_list,delt);
%     
%     Bgrad_mag=sum(abs(Bgrad_list),2);
%     Bgrad_grid=reshape(Bgrad_mag,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
%     figure(12)
%     clf;
%     h=surf((xyz_grid(:,:,:,1)-trap_cent(1))*1e3,(xyz_grid(:,:,:,3)-trap_cent(3))*1e3,Bgrad_grid,'facealpha',0.9);
%     colormap(viridis())
%     set(h,'LineStyle','none')
%     box on;
%     xlabel('X (mm)');
%     ylabel('Z (mm)');
%     title('Jacobian')
%     zlabel('RSS Grad. (T/m) ');
%     set(gcf,'Color',[1 1 1]);
%     pause(0.001);
% 
%     figure(14)
%     clf
%     contour((xyz_grid(:,:,:,1)-trap_cent(1))*1e3,(xyz_grid(:,:,:,3)-trap_cent(3))*1e3,Bgrad_grid,300)
%     set(gcf,'Color',[1 1 1]);
%     title('Jacobian')
%     xlabel('X (mm)');
%     ylabel('Z (mm)');
%     pause(1e-3)
%     if solve_stpt>1
%         fprintf('starting optimization on potential landscape \n')
%         %find the minimum points of the potential landscape
%         min_pts=[];
%         num_pts=30; 
%         %rng(123456) %for protoryping
%         x0=[plot_range(1,1)+rand(1,num_pts)'*(plot_range(1,2)-plot_range(1,1)),plot_range(2,1)+rand(1,num_pts)'*(plot_range(2,2)-plot_range(2,1))];
%         for n=1:num_pts
%             options=optimset('MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-10);
%             min_pts(n,:)=fminsearch(@(x) trap_eval(btrap,[[x(1),0,x(2)]]),x0(n,:),options);
%         end
%         %filter the points that have stayed in the range
%         mask=min_pts(:,1)>plot_range(1,1) & min_pts(:,1)<plot_range(1,2) & min_pts(:,2)>plot_range(2,1) & min_pts(:,2)<plot_range(2,2);
%         st_pts=min_pts(mask,:);
%         fprintf('surviving points %i out of %i \n',size(st_pts,1),num_pts)
%         %custer the points
%         % now cluster points
%         %break this clustering out as a function once 3d
%         min_pt_clusters={};
%         clust_idx=1;
%         while size(min_pts,1)>0
%             mask=sum((min_pts(1,:)-min_pts).^2,2)<min_cluster_range; %are the points within some threshold distance
%             min_pt_clusters{clust_idx}=min_pts(mask,:);
%             min_pts=min_pts(mask~=1,:);
%             clust_idx=clust_idx+1;
%         end
%         min_pts=[];
%         for n=1:size(min_pt_clusters,2) 
%             min_pts(n,:)=mean(min_pt_clusters{n},1);
%         end
%         min_pts(:,3)=trap_eval(btrap,[min_pts(:,1),zeros(size(min_pts,1),1),min_pts(:,2)]); 
%         min_pts(:,4)=sum(abs(num_grad(@(x)trap_eval(btrap,x),[min_pts(:,1),zeros(size(min_pts,1),1),min_pts(:,2)],delt)),2);
%         fprintf('minima clustered to %i points, %i zero crossing \n',size(min_pts,1),sum(min_pts(:,3)<zero_feild_thresh))
% 
%         figure(14)
%         hold on
%         scatter3((min_pts(:,1)-trap_cent(1))*1e3,(min_pts(:,2)-trap_cent(3))*1e3, min_pts(:,4),20,'r','filled')
%         pause(1e-3)
%         hold off       
%               
%     end
%     
%     if solve_stpt>2
%         fprintf('starting optimization on derivative landscape \n')        
%         st_pts=[];
%         num_pts=30;
%         %rng(123456) %for protoryping
%         x0=[plot_range(1,1)+rand(1,num_pts)'*(plot_range(1,2)-plot_range(1,1)),plot_range(2,1)+rand(1,num_pts)'*(plot_range(2,2)-plot_range(2,1))];
%         for n=1:num_pts
%             options=optimset('MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',1e-5);
%             st_pts(n,:)=fminsearch(@(y) sum(abs(num_grad(@(x) trap_eval(btrap,x),[[y(1),0,y(2)]],delt)),2),x0(n,:),options);
%             %PSoptions = optimoptions(@patternsearch,'Display','iter','TolFun',1e-10,'TolMesh',1e-10);
%             %[st_pts(n,1:2),st_pts(n,3)] = patternsearch(@(y)  num_grad(@(x) trap_eval(btrap,x),y,delt),x0(n,:),[],[],[],[],plot_range(:,1),plot_range(:,2),PSoptions);
%         end
%         mask=plot_range(1,1)<st_pts(:,1) & st_pts(:,1)<plot_range(1,2) & plot_range(2,1)<st_pts(:,2) & st_pts(:,2)<plot_range(2,2);
%         fprintf('removing %i points that are out of ROI \n',sum(mask))
%         st_pts=st_pts(mask,:);
%         st_pts_raw=st_pts;
%         
%         %remove points that are close to the previously found minima
%         mask=false(1,size(st_pts,1));
%         for n=1:size(st_pts,1)
%             dist=sqrt(sum((st_pts(n,1:2)-min_pts(:,1:2)).^2,2));
%             mask(n)=min(dist)<st_cluster_range;
%         end
%         fprintf('removing %i points that are close to minima \n',sum(mask))
%         st_pts=st_pts(~mask,:);
%         
%         H_list=num_hessian(@(x) trap_eval(btrap,x),[st_pts(:,1),zeros(size(st_pts,1),1),st_pts(:,2)],delt);
%         mask=false(1,size(st_pts,1));
%         for i=1:size(H_list,3) %find the determinacy of each of these hessians
%             eh=eig(H_list(:,:,i));
%             mask(i)=sum(eh>1)>0 & sum(eh<0)>0; %check that contains +ve and -ve
%         end
%         clear('eh')
%         fprintf('removing %i points  that do not pass hessian test \n',sum(~mask))
%         st_pts=st_pts(mask,:);
%         
%         fprintf('surviving points %i out of %i \n',size(st_pts,1),num_pts)
%         
%         st_pt_clusters={};
%         clust_idx=1;
%         mask=false(1,size(st_pts,1));
%         while size(st_pts,1)>0
%             mask=sum((st_pts(1,:)-st_pts).^2,2)<st_cluster_range; %are the points within some threshold distance
%             st_pt_clusters{clust_idx}=st_pts(mask,:);
%             st_pts=st_pts(mask~=1,:);
%             clust_idx=clust_idx+1;
%         end
%         st_pts=[];
%         for n=1:size(st_pt_clusters,2) 
%             st_pts(n,:)=mean(st_pt_clusters{n},1);
%         end
%           
%         st_pts(:,3)=trap_eval(btrap,[st_pts(:,1),zeros(size(st_pts,1),1),st_pts(:,2)]); 
%         st_pts(:,4)=sum(abs(num_grad(@(x)trap_eval(btrap,x),[st_pts(:,1),zeros(size(st_pts,1),1),st_pts(:,2)],delt)),2);
%         fprintf('saddles clustered to %i points, %i zero crossing \n',size(st_pts,1),sum(st_pts(:,3)<zero_feild_thresh))
% 
% 
%         for i=1: size(st_pts,1)
%             fprintf('Trap saddle found %2.3fG (%2.3f MHz)(%2.3E K)\n',st_pts(i,3)*1e4,st_pts(i,3)*const.b_freq*1e-6,(2*st_pts(i,3)*const.mub)*2/const.kb)
%         end
%         
%         figure(14)
%         hold on
%         scatter3((st_pts(:,1)-trap_cent(1))*1e3,(st_pts(:,2)-trap_cent(3))*1e3,st_pts(:,4),20,'k','filled')
%         pause(1e-3)
%         hold off
%         
%       
%     end
end




%% End of code
t_end=toc(t_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_end);
disp('===================ALL TASKS COMPLETED===================');