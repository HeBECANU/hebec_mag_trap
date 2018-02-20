%% BiQUIC trap magnetic field calculator
% D K Shin
% refer to https://doi.org/10.1103/PhysRevA.35.1535 for expressions used
% B M Henson

%creates a structure bfeild that specifies the currents


% Known Bugs/Errors
% trap freq not reproduced
%trap ratio not reprotuced


% TODO:
%reproduce numbers from https://www.sciencedirect.com/science/article/pii/S0030401806009680?via%3Dihub
%currently trap freq way higher than reported there 
%factor out plot code
    %add options for various 2d slices
%add in gravitational potential
%find trap depth
    %to start with using the connected region method
    %then using the saddle following method


% DONE:
%code runs on arrays
%check that the Energy of a He* atom in a B feild is E=2*ub*B  
%looks to be the case http://iopscience.iop.org.virtual.anu.edu.au/article/10.1088/1464-4266/5/2/360/pdf

t_start=tic;
% ------------------START USER Config--------------
%%% Flags
verbose=1;          % graphical output

solve_3D=0;         % solve full 3D vector B-field (takes a while)
solve_2D=1;
solve_trapchar=1;   % characterise trap params including freq and center
solve_trapdepth=0;  %VERY SLOW and a can require fidling to get working ok

%%% mag trap
v_quad=3.4;%2.4%3.4 %3.4 used in 'normal trap'
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
    [f0,trap_cent,B_cent]=trap_characterise(btrap,-0e-3,solve_trapdepth,verbose);
end



%%% 3D grid trap B-field
if solve_3D
    ngrid=100;          % 50 - med; 300 - very fine;
    % grid in trap centered ref frame
    xyz_grid=[];
    [xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3)]=...
    meshgrid(trap_cent(1)+linspace(-5e-3,20e-3,ngrid),...
             trap_cent(2)+linspace(-5e-3,5e-3,ngrid),...
             trap_cent(3)+linspace(-5e-3,5e-3,ngrid));    % meshgrid
    xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);

    [Bmag_list,Bxyz]=trap_eval(btrap,xyz_list);

    Bmag_grid=reshape(Bmag_list,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
    figure(2)
    set(gcf,'Color',[1 1 1]);
    clf;
    nBisosurf=10;
    Bmax=max(Bmag_list(isfinite(Bmag_list)));
    Bmin=min(Bmag_list(isfinite(Bmag_list)));
    Bisoval=logspace(log10(Bmin),log10(Bmax),nBisosurf);
    Bisoval=Bisoval(1:end);   % cull the min and max
    cc=viridis(nBisosurf);
    p={};
    pp=[];
    for ii=1:nBisosurf
        p{ii}=isosurface(1e3*xyz_grid(:,:,:,1),1e3*xyz_grid(:,:,:,2),1e3*xyz_grid(:,:,:,3),Bmag_grid,Bisoval(ii));
        pp(ii)=patch(p{ii},'FaceColor',cc(ii,:),'EdgeColor','none','FaceAlpha',0.15,...
            'DisplayName',sprintf('%0.1g',1e4*Bisoval(ii)));
    end
    box on;
    daspect([1,1,1]);
    view(3);
    camlight;
    lighting gouraud;

    xlabel('X [mm]');
    ylabel('Y [mm]');
    zlabel('Z [mm]');

    nnring=100;
    phi=linspace(0,2*pi,nnring);
    [x_ring,y_ring]=pol2cart(phi,1);
    z_ring=zeros(1,nnring);
    for ii=1:numel(btrap)
        % only plot coils
        if isequal(btrap(ii).type,'coil')
            % transform from unit ring
            R_coil_this=btrap(ii).param{1};
            pos_coil_this=btrap(ii).param{3};
            xthis=R_coil_this*x_ring+pos_coil_this(1);
            ythis=R_coil_this*y_ring+pos_coil_this(2);
            zthis=R_coil_this*z_ring+pos_coil_this(3);

            % draw this coil
            hold on;
            plot3(1e3*xthis,1e3*ythis,1e3*zthis,...
                'Color','k','LineWidth',2);
        end
    end
pause(0.01)
end







%%% 2D grid trap B-field
if solve_2D
    ngrid=100;          % 50 - med; 300 - very fine;
    % grid in trap centered ref frame
    xyz_grid=[];
    [xyz_grid(:,:,:,1),xyz_grid(:,:,:,3)]=...
    meshgrid(trap_cent(1)+linspace(-20e-3,20e-3,ngrid),...
             trap_cent(3)+linspace(-7e-3,7e-3,ngrid));    % meshgrid

    xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
    [Bmag_list,Bxyz]=trap_eval(btrap,xyz_list);
    Bmag_grid=reshape(Bmag_list,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
    figure(3)
    clf;
    h=surf((xyz_grid(:,:,:,1)-trap_cent(1))*1e3,(xyz_grid(:,:,:,3)-trap_cent(3))*1e3,Bmag_grid*1e4);
    colormap(viridis())
    set(h,'LineStyle','none')
    box on;
    xlabel('X [mm]');
    ylabel('Z [mm]');
    zlabel('B (Gauss)');
    set(gcf,'Color',[1 1 1]);
    pause(0.001);
    
    
    

    xyz_grid=zeros(ngrid,ngrid,1,3);
    [xyz_grid(:,:,:,1),xyz_grid(:,:,:,2)]=...
    meshgrid(trap_cent(1)+linspace(-20e-3,20e-3,ngrid),...
             trap_cent(3)+linspace(-5e-3,5e-3,ngrid));    % meshgrid
    xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
    [Bmag_list,Bxyz]=trap_eval(btrap,xyz_list);
    Bmag_grid=reshape(Bmag_list,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);

    figure(4)
    clf;
    h=surf((xyz_grid(:,:,:,1)-trap_cent(1))*1e3,(xyz_grid(:,:,:,2)-trap_cent(3))*1e3,Bmag_grid*1e4);
    colormap(viridis())
    set(h,'LineStyle','none')
    box on;
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('B (Gauss)');
    set(gcf,'Color',[1 1 1]);
    
end







%% Solve B-field for trap in 3D

% 
% if solve_3D>0
%     % 3D B-field
%         % calculate magnetic field (macro summary)
%     
%     %%% visualise
%     % 3D B-isosurfaces
%     if verbose>0
%         hfig_btrap=plot_B_3d(btrap,Bmag,xyz_list);
%         
%         if solve_trapchar>0
%             % display evaluated trap centre
%             scatter3(1e3*trap_cent(1),1e3*trap_cent(2),1e3*trap_cent(3),...
%                 'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none',...
%                 'SizeData',30,'DisplayName','Trap centre');
%         end
%     end
% end
% 




%% End of code
t_end=toc(t_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_end);
disp('===================ALL TASKS COMPLETED===================');