function  visualise_vector_3d(plot_opts)

    plot_scaling=1e3; % change position scales
    if ~isfield(plot_opts,'show_coils')
        plot_opts.show_coils=true;
    end

    plot_range=plot_opts.range;
    
    if plot_opts.zero_on_cen
        axes_offset=plot_opts.cen;
        plot_range=plot_range+repmat(plot_opts.cen',1,2);
    else
        plot_range=plot_range+repmat(plot_opts.cen',1,2);
        axes_offset=zeros(1,3);
    end
    % grid in trap centered ref frame
    xyz_grid=[];
    [xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3)]=...
    meshgrid(linspace(plot_range(1,1),plot_range(1,2),plot_opts.nsamp(1)),...
             linspace(plot_range(2,1),plot_range(2,2),plot_opts.nsamp(2)),...
             linspace(plot_range(3,1),plot_range(3,2),plot_opts.nsamp(3)));    % meshgrid
    xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
    plot_opts.xyz_list=xyz_list;
    [~,bvec_list]=trap_eval(plot_opts.btrap,xyz_list);
     
    
    Bvec_grid=reshape(bvec_list,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3),3]);
    
    contor_plot_handle=stfig('3d vector','add_stack',1);
    set(gcf,'Color',[1 1 1]);
    clf;
    
    quiver3(plot_scaling*(xyz_grid(:,:,:,1)-axes_offset(1)),...
            plot_scaling*(xyz_grid(:,:,:,2)-axes_offset(2)),...
            plot_scaling*(xyz_grid(:,:,:,3)-axes_offset(3)),...
            Bvec_grid(:,:,:,1),Bvec_grid(:,:,:,2),Bvec_grid(:,:,:,3))
    
    
    
    xlim((plot_range(1,:)-axes_offset(1))*plot_scaling)
    ylim((plot_range(2,:)-axes_offset(2))*plot_scaling)
    zlim((plot_range(3,:)-axes_offset(3))*plot_scaling)
    xlabel('X [mm]');
    ylabel('Y [mm]');
    zlabel('Z [mm]');


    hold off;
pause(0.01)

end

