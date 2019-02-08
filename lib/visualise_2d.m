function visualise_2d(btrap,plot_2d_opt)
    trap_cen=plot_2d_opt.cen;
    plot_range=plot_2d_opt.range;

    ngrid=100;          % 50 - med; 300 - very fine;
    rot_mat=rotationVectorToMatrix(plot_2d_opt.rot);
    if plot_2d_opt.zero_on_cen
        rot_cen=trap_cen*rot_mat;
    else
        rot_cen=zeros(1,3);
    end
    
    %plot_range=[[-2,12];[-8.5,8.5]]*1e-3; %x,z
    
    % grid in trap centered ref frame
    xyz_grid=zeros(ngrid,ngrid,1,3);
    [xyz_grid(:,:,:,1),xyz_grid(:,:,:,2)]=...
    meshgrid(linspace(plot_range(1,1),plot_range(1,2),ngrid),...
             linspace(plot_range(2,1),plot_range(2,2),ngrid));    % meshgrid
    xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
    xyz_list=xyz_list+trap_cen;
    xyz_list=xyz_list*rot_mat;
    [Bmag_list,Bxyz]=trap_eval(btrap,xyz_list);
    Bmag_grid=reshape(Bmag_list,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
    figure(3)
    clf;
    h=surf((xyz_grid(:,:,:,1)-rot_cen(1))*1e3,(xyz_grid(:,:,:,2)-rot_cen(2))*1e3,Bmag_grid*1e4);
    colormap(viridis())
    set(h,'LineStyle','none')
    box on;
    xlabel('X [mm]');
    ylabel('Y [mm]');
    title('Potential')
    zlabel('B (Gauss)');
    set(gcf,'Color',[1 1 1]);
    pause(0.001);
    
    
    figure(13)
    clf
    contour((xyz_grid(:,:,:,1)-rot_cen(1))*1e3,(xyz_grid(:,:,:,2)-rot_cen(2))*1e3,Bmag_grid,300)
    set(gcf,'Color',[1 1 1]);
    title('Potential')
    xlabel('X (mm)');
    ylabel('Y (mm)');
    
%     xyz_grid=zeros(ngrid,ngrid,1,3);
%     [xyz_grid(:,:,:,1),xyz_grid(:,:,:,2)]=...
%     meshgrid(trap_cen(1)+linspace(-20e-3,20e-3,ngrid),...
%              trap_cen(3)+linspace(-5e-3,5e-3,ngrid));    % meshgrid
%     xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
%     [Bmag_list,Bxyz]=trap_eval(btrap,xyz_list);
%     Bmag_grid=reshape(Bmag_list,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
% 
%     figure(4)
%     clf;
%     h=surf((xyz_grid(:,:,:,1)-trap_cen(1))*1e3,(xyz_grid(:,:,:,2)-trap_cen(3))*1e3,Bmag_grid*1e4,'facealpha',0.9);
%     colormap(viridis())
%     set(h,'LineStyle','none')
%     box on;
%     xlabel('X (mm)');
%     ylabel('Y (mm)');
%     zlabel('B (Gauss)');
%     title('Potential')
%     set(gcf,'Color',[1 1 1]);
%     pause(0.001)
%     
%     figure(5)
%     clf
%     contour((xyz_grid(:,:,:,1)-trap_cen(1))*1e3,(xyz_grid(:,:,:,2)-trap_cen(3))*1e3,Bmag_grid,300)
%     set(gcf,'Color',[1 1 1]);
%     title('Potential')
%     xlabel('X (mm)');
%     ylabel('Y (mm)');
    
end