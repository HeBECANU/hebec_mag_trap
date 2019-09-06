function visualise_3d(plot_opts)
    plot_scaling=1e3;
    if ~isfield(plot_opts,'show_coils')
        plot_opts.show_coils=true;
    end

    plot_range=plot_opts.range;
    if plot_opts.zero_on_cen
        trap_cen=plot_opts.cen;
    else
        trap_cen=zeros(1,3);
    end
    % grid in trap centered ref frame
    xyz_grid=[];
    [xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3)]=...
    meshgrid(linspace(plot_range(1,1),plot_range(1,2),plot_opts.nsamp(1)),...
             linspace(plot_range(2,1),plot_range(2,2),plot_opts.nsamp(2)),...
             linspace(plot_range(3,1),plot_range(3,2),plot_opts.nsamp(3)));    % meshgrid
    xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
    xyz_list=xyz_list+trap_cen;
    plot_opts.xyz_list=xyz_list;
    scal_res=compute_scalar_property(plot_opts);
    
    Bmag_grid=reshape(scal_res.val,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
    contor_plot_handle=stfig('3d potential contors','add_stack',1);
    set(gcf,'Color',[1 1 1]);
    clf;
    nisosurf=10;
    Bmax=max(Bmag_grid(isfinite(Bmag_grid(:))));
    Bmin=min(Bmag_grid(isfinite(Bmag_grid(:))));
    Bisoval=logspace(log10(Bmin),log10(Bmax),nisosurf);
    Bisoval=Bisoval(1:end);   % cull the min and max
    cc=viridis(nisosurf);
    p={};
    pp=[];
    for ii=1:nisosurf
        p{ii}=isosurface(plot_scaling*xyz_grid(:,:,:,1),plot_scaling*xyz_grid(:,:,:,2),plot_scaling*xyz_grid(:,:,:,3),Bmag_grid,Bisoval(ii));
        pp(ii)=patch(p{ii},'FaceColor',cc(ii,:),'EdgeColor','none','FaceAlpha',0.15,...
            'DisplayName',sprintf('%0.1g',1e4*Bisoval(ii)));
    end
    box on;
    daspect([1,1,1]);
    view(3);
    camlight;
    lighting gouraud;
    xlim(plot_range(1,:)*plot_scaling)
    ylim(plot_range(2,:)*plot_scaling)
    zlim(plot_range(3,:)*plot_scaling)

    xlabel('X [mm]');
    ylabel('Y [mm]');
    zlabel('Z [mm]');

    nnring=100;
    phi=linspace(0,2*pi,nnring);
    xyz_unit_loop=zeros(numel(phi),3);
    [xyz_unit_loop(:,1),xyz_unit_loop(:,2)]=pol2cart(phi,1);
    xyz_unit_loop(:,3)=zeros(1,nnring);
    nnline=100;
    xyz_unit_line=zeros(nnline,3);
    xyz_unit_line(:,3)=linspace(0,1,nnline);
    
    stfig(contor_plot_handle);
    hold on;
    for ii=1:numel(plot_opts.btrap.b_src)
        elm_param=plot_opts.btrap.b_src(ii).param;
        % only plot coils
        if isequal(plot_opts.btrap.b_src(ii).type,'loop')
            % transform from unit ring
            rad_coil=elm_param.radius;
            pos=elm_param.position;
            xyz_pos=xyz_unit_loop*rad_coil;
            rot_vec=elm_param.rot;
            rot_mat=rotationVectorToMatrix(-rot_vec);
            xyz_pos=xyz_pos*rot_mat+pos;
            
            % draw this coil
            
            plot3(plot_scaling*xyz_pos(:,1),plot_scaling*xyz_pos(:,2),plot_scaling*xyz_pos(:,3),...
                'Color','k','LineWidth',2);
        end
        if isequal(plot_opts.btrap.b_src(ii).type,'helix')
            % transform from unit ring
            radius=elm_param.radius;
            pitch=elm_param.pitch;
            turn=[0,elm_param.turns];
            dlen=elm_param.dlen;
            pos=elm_param.position;
            rot_vec=elm_param.rot;
            rot_mat=rotationVectorToMatrix(-rot_vec);
            tlim=2*pi*turn;
            num_lines=min([1e3,ceil(range(tlim)/dlen)]);
            tvec=linspace(tlim(1),tlim(2),num_lines)';
            wire_pos=[radius*cos(tvec), radius*sin(tvec), pitch*tvec/(2*pi)];
            xyz_pos=wire_pos*rot_mat+pos;
            
            % draw this coil
            
            plot3(plot_scaling*xyz_pos(:,1),plot_scaling*xyz_pos(:,2),plot_scaling*xyz_pos(:,3),...
                'Color','k','LineWidth',2);
        end
        if isequal(plot_opts.btrap.b_src(ii).type,'line')
            % transform from unit ring
            len_line_this=plot_opts.btrap.b_src(ii).param.length;
            pos=plot_opts.btrap.b_src(ii).param.position;
            xyz_pos=xyz_unit_line*len_line_this;
            rot_mat=rotationVectorToMatrix(-plot_opts.btrap.b_src(ii).param.rot);
            xyz_pos=xyz_pos*rot_mat+pos;
            
            % draw this coil
            
            plot3(plot_scaling*xyz_pos(:,1),plot_scaling*xyz_pos(:,2),plot_scaling*xyz_pos(:,3),...
                'Color','k','LineWidth',2);
        end
    end
    sensor_size=30e-3;
    stfig(contor_plot_handle);
    if isfield(plot_opts.btrap,'nullr') && isfield(plot_opts.btrap.nullr,'sensor')
        for ii=1:numel(plot_opts.btrap.nullr.sensor)
            start_pt=plot_opts.btrap.nullr.sensor(ii).pos;
            norm_dirn=plot_opts.btrap.nullr.sensor(ii).dirn/norm(plot_opts.btrap.nullr.sensor(ii).dirn);
            end_pt=start_pt+sensor_size*norm_dirn;
            mArrow3(plot_scaling*start_pt,plot_scaling*end_pt,'color','red','stemWidth',1,'facealpha',0.5);
        end
    end
    
    xlim(plot_range(1,:)*plot_scaling)
    ylim(plot_range(2,:)*plot_scaling)
    zlim(plot_range(3,:)*plot_scaling)
    hold off;
pause(0.01)

end