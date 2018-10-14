function H=plot_B_3d(btrap,Bmag,xyz_list)
% Plot
% H = plot_B_3d(btrap, Bmag, XYZ)
%

% config
nBisosurf=5;

%%% Parse input
% grid limits in x,y,z
xyz_lim=cellfun(@(x)[min(x(:)),max(x(:))],xyz_list,'UniformOutput',false);
xyz_lim=vertcat(xyz_lim{:});

%%% Magnetic field
% quiver plot for B field
H=figure();
% quiver3(1e3*XYZ{1},1e3*XYZ{2},1e3*XYZ{3},Bxyz{:},...
%     'Color','k','LineWidth',1.5,'Visible','off');
hold on;

%%% B field magnitude - isosurfaces (isopotentials)
Bmax=max(Bmag(isfinite(Bmag)));
Bmin=min(Bmag(isfinite(Bmag)));
% Bisoval=linspace(Bmin,Bmax,nBisosurf+2);
Bisoval=logspace(log10(Bmin),log10(Bmax),nBisosurf+2);
Bisoval=Bisoval(2:end-1);   % cull the min and max
cc=viridis(nBisosurf);
p={};
pp=[];
for ii=1:nBisosurf
    p{ii}=isosurface(1e3*xyz_list{1},1e3*xyz_list{2},1e3*xyz_list{3},Bmag,Bisoval(ii));
    pp(ii)=patch(p{ii},'FaceColor',cc(ii,:),'EdgeColor','none','FaceAlpha',0.15,...
        'DisplayName',sprintf('%0.1g',1e4*Bisoval(ii)));
end
box on;
daspect([1,1,1]);
view(3);
% xlim(1e3*[min(xyz{1}),max(xyz{1})]);
% ylim(1e3*[min(xyz{2}),max(xyz{2})]);
% zlim(1e3*[min(xyz{3}),max(xyz{3})]);
xlim(1e3*xyz_lim(1,:));
ylim(1e3*xyz_lim(2,:));
zlim(1e3*xyz_lim(3,:));
camlight;
lighting gouraud;

xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');

%%% Coils
% NOTE: coil widths in plot are not to true scale!
% unit ring
nnring=100;
phi=linspace(0,2*pi,nnring);
[x_ring,y_ring]=pol2cart(phi,1);
z_ring=zeros(1,nnring);
figure(H);
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

% legend
lgd=legend(pp(:));
title(lgd,'$B$-isosurface (G)');

end