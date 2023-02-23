function visualise_2d(plot_opts)
    %plots a scalar value in a plane
    %can rotate that plane about the center point 
    
    %plot_range=[[-2,12];[-8.5,8.5]]*1e-3; %x,z
    plot_scaling=1e3;
    
    plot_range=plot_opts.range;
    plot_cen=plot_opts.cen;
    is_rotated=~isequal(plot_opts.rot,[0,0,0]);
    rot_mat=rotationVectorToMatrix(plot_opts.rot);
   
    % grid in unrotated ref frame
    xyz_grid=zeros(plot_opts.nsamp(1),plot_opts.nsamp(2),3);
    [xyz_grid(:,:,1),xyz_grid(:,:,2)]=...
    meshgrid(linspace(plot_range(1,1),plot_range(1,2),plot_opts.nsamp(1)),...
             linspace(plot_range(2,1),plot_range(2,2),plot_opts.nsamp(2)));    % meshgrid
    xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2),size(xyz_grid,3)]);
    %optimizations could be implmented here to not shift or rotate in the null cases
    %xyz_list=xyz_list; 
    xyz_list=xyz_list*rot_mat;
    xyz_list=xyz_list+plot_cen; 
    %xyz_list=xyz_list+plot_cen.*[0,0,1];
    plot_opts.xyz_list=xyz_list;
    scal_res=compute_scalar_property(plot_opts);
    
    bvec_grid=reshape(scal_res.val,[size(xyz_grid,1),size(xyz_grid,2)]);
    stfig('2d potential surface plot','add_stack',1) ;
    clf;
    h=surf((xyz_grid(:,:,1))*plot_scaling,(xyz_grid(:,:,2))*plot_scaling,bvec_grid*1e4);
    colormap(viridis())
    set(h,'LineStyle','none')
    box on;
    if is_rotated % put a ' to indicate the roated cord frame
        xlabel("X' (mm)");
        ylabel("Y' (mm)");
        title(strcat(scal_res.plot_title,' (rotated cord.)'))
    else
        xlabel('X (mm)');
        ylabel('Y (mm)');
        title(scal_res.plot_title)
    end
    zlabel(scal_res.plot_zlabel);
    set(gcf,'Color',[1 1 1]);
    pause(0.001);
    
    
    stfig('2d potential contor plot','add_stack',1) 
    clf
    contour((xyz_grid(:,:,1)-plot_cen(1))*plot_scaling,(xyz_grid(:,:,2)-plot_cen(2))*plot_scaling,bvec_grid*1e4,300)
    set(gcf,'Color',[1 1 1]);
    colormap(viridis())
    c = colorbar;
    ylabel(c, scal_res.plot_zlabel)
<<<<<<< HEAD
    title(scal_res.plot_title)
    xlabel('X (mm)');
    ylabel('Y (mm)');
    
    if plot_opts.vec_plot.do
        hold on
        % grid in unrotated ref frame
        xyz_grid=zeros(plot_opts.vec_plot.nsamp(1),plot_opts.vec_plot.nsamp(2),3);
        [xyz_grid(:,:,1),xyz_grid(:,:,2)]=...
        meshgrid(linspace(plot_range(1,1),plot_range(1,2),plot_opts.vec_plot.nsamp(1)),...
                 linspace(plot_range(2,1),plot_range(2,2),plot_opts.vec_plot.nsamp(2)));    % meshgrid
        xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2),size(xyz_grid,3)]);
        %optimizations could be implmented here to not shift or rotate in the null cases
        xyz_list=xyz_list-plot_cen; 
        xyz_list=xyz_list*rot_mat;
        xyz_list=xyz_list+plot_cen; 

        [~,bvec_list]=trap_eval(plot_opts.btrap,xyz_list);
    
        
        bvec_list=bvec_list*rot_mat;
        bvec_list = normalize(bvec_list,2);
        bvec_grid=reshape(bvec_list,[size(xyz_grid,1),size(xyz_grid,2),3]);

        h=quiver((xyz_grid(:,:,1)-plot_cen(1))*1e3,(xyz_grid(:,:,2)-plot_cen(2))*1e3,bvec_grid(:,:,1),bvec_grid(:,:,2));
        
        hold off
=======
    if is_rotated % put a ' to indicate the roated cord frame
        xlabel("X' (mm)");
        ylabel("Y' (mm)");
        title(strcat(scal_res.plot_title,' (rotated cord.)'))
    else
        xlabel('X (mm)');
        ylabel('Y (mm)');
        title(scal_res.plot_title)
>>>>>>> fb94d8a7023f83e406b2cf3e332982b4585859fc
    end
    
 
        
       
    xlim((plot_range(1,:))*plot_scaling)
    ylim((plot_range(2,:))*plot_scaling)
    

    
end