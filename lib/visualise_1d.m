function visualise_1d(plot_opts)
    %plots a scalar value along a line
    %can rotate that line about the center point 
    
    plot_range=plot_opts.range;
    plot_cen=plot_opts.plot_cen;
    is_rotated=~isequal(plot_opts.rot,[0,0,0]);
    rot_mat=rotationVectorToMatrix(plot_opts.rot);
   
    % grid in unrotated ref frame
    xyz_list=zeros(plot_opts.nsamp(1),3);
    xyz_list(:,1)=linspace(plot_range(1,1),plot_range(1,2),plot_opts.nsamp(1));
    %optimizations could be implmented here to not shift or rotate in the null cases
    xyz_list=xyz_list*rot_mat;
    xyz_list=xyz_list+plot_cen; 
    %xyz_list=xyz_list+plot_cen.*[0,0,1];
    plot_opts.xyz_list=xyz_list;
    scal_res=compute_scalar_property(plot_opts);
    
    result_vec=scal_res.val;
    stfig('1d plot','add_stack',1) 
    clf;
    plot((xyz_list(:,1)-plot_cen(1))*1e3,result_vec,'k');
    if is_rotated % put a ' to indicate the roated cord frame
        xlabel("X' (mm)");
        title(strcat(scal_res.plot_title,' (rotated cord.)'))
    else
        xlabel('X (mm)');
    end
    ylabel(scal_res.plot_zlabel);
    set(gcf,'Color',[1 1 1]);
    pause(0.001);
    
        
    
end
