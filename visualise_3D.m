function visualise_3D(btrap,trap_cent)

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