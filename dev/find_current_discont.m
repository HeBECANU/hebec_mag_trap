%try current discontinuities finding methods

clear all
close all
% create a path of a line
syms t
len=1;
curr=1;

rot_vec=[0.0,0,0];
tlim=[0,len];
dlen=1e-1;
path_line=symfun([0,0,t],t);

%% Bfield_finite_line_analytic  vs  Bfield_path_numeric
plot_range=[[-1,1];[-1,1];[-1,1]].*0.001;
num_points=100;
xyz_grid=[];
[xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3)]=...
    meshgrid(linspace(plot_range(1,1),plot_range(1,2),num_points),...
             linspace(plot_range(2,1),plot_range(2,2),num_points),...
             linspace(plot_range(3,1),plot_range(3,2),num_points));    % meshgrid
xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);

vector_results_analytic=Bfield_finite_line_analytic(tlim,curr,rot_vec,xyz_list);

Bmag_grid_analytic=reshape(vector_results_analytic,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3),3]);
    
% find the divergence

div_grid_analytic=divergence(xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3),...
                    Bmag_grid_analytic(:,:,:,1),Bmag_grid_analytic(:,:,:,2),Bmag_grid_analytic(:,:,:,3));
div_grid_analytic=abs(div_grid_analytic);                
                
fprintf('analytic method\n')
fprintf('div mean %e\n',mean(div_grid_analytic(:)))
fprintf('div max %e \n',max(div_grid_analytic(:)))                
fprintf('div min %e \n',min(abs(div_grid_analytic(:))))

%div_grid=vecnorm(Bmag_grid,2,4);                

stfig('divergance isosurfaces')
clf
subplot(1,2,1)
nisosurf=10;
Bmax=max(Bmag_grid_analytic(isfinite(Bmag_grid_analytic(:))));
Bmin=min(Bmag_grid_analytic(isfinite(Bmag_grid_analytic(:))))+eps;
isoval=logspace(log10(Bmin),log10(Bmax),nisosurf);
isoval=isoval(1:end);   % cull the min and max
cc=viridis(nisosurf);
surface={};
surface_fig_obj=[];
for ii=1:nisosurf
    surface{ii}=isosurface(xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3),div_grid_analytic,isoval(ii));
    surface_fig_obj(ii)=patch(surface{ii},'FaceColor',cc(ii,:),'EdgeColor','none','FaceAlpha',0.15,...
        'DisplayName',sprintf('%0.1g',1e4*isoval(ii)));
end                
xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');     
xlim(plot_range(1,:))
ylim(plot_range(2,:))
zlim(plot_range(3,:))


                
%% 


vector_results_numeric=Bfield_finite_line_numeric(tlim,dlen,curr,rot_vec,xyz_list);
Bmag_grid_numeric=reshape(vector_results_numeric,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3),3]);
    
% find the divergence
div_grid_numeric=divergence(xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3),...
                    Bmag_grid_numeric(:,:,:,1),Bmag_grid_numeric(:,:,:,2),Bmag_grid_numeric(:,:,:,3));
div_grid_numeric=abs(div_grid_numeric);
%div_grid=vecnorm(Bmag_grid,2,4);                

fprintf('numeric method\n')
fprintf('div mean %e\n',mean(div_grid_numeric(:)))
fprintf('div max %e \n',max(div_grid_analytic(:)))                
fprintf('div min %e \n',min(abs(div_grid_analytic(:))))

stfig('divergance isosurfaces')
subplot(1,2,2)
nisosurf=10;
Bmax=max(div_grid_numeric(isfinite(div_grid_numeric(:))));
Bmin=min(div_grid_numeric(isfinite(div_grid_numeric(:))))+eps;
isoval=logspace(log10(Bmin),log10(Bmax),nisosurf);
isoval=isoval(1:end);   % cull the min and max
cc=viridis(nisosurf);
surface={};
surface_fig_obj=[];
for ii=1:nisosurf
    surface{ii}=isosurface(xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3),div_grid_numeric,isoval(ii));
    surface_fig_obj(ii)=patch(surface{ii},'FaceColor',cc(ii,:),'EdgeColor','none','FaceAlpha',0.15,...
        'DisplayName',sprintf('%0.1g',1e4*isoval(ii)));
end                
xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');     
xlim(plot_range(1,:))
ylim(plot_range(2,:))
zlim(plot_range(3,:))


%%

Bmag_grid_diff=frac_diff(Bmag_grid_numeric,Bmag_grid_analytic);
Bmag_grid_diff=abs(nanmean(Bmag_grid_diff,4));
              

fprintf('mag field diff\n')
fprintf('diff mean %e\n',mean(Bmag_grid_diff(:)))
fprintf('diff max %e \n',max(Bmag_grid_diff(:)))                
fprintf('diff min %e \n',min(abs(Bmag_grid_diff(:))))

stfig('mag diff')
nisosurf=20;
Bmax=max(Bmag_grid_diff(isfinite(Bmag_grid_diff(:))));
Bmin=min(Bmag_grid_diff(isfinite(Bmag_grid_diff(:))))+eps;
isoval=logspace(log10(Bmin),log10(Bmax),nisosurf);
isoval=isoval(1:end);   % cull the min and max
cc=viridis(nisosurf);
surface={};
surface_fig_obj=[];
for ii=1:nisosurf
    surface{ii}=isosurface(xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3),Bmag_grid_diff,isoval(ii));
    surface_fig_obj(ii)=patch(surface{ii},'FaceColor',cc(ii,:),'EdgeColor','none','FaceAlpha',0.15,...
        'DisplayName',sprintf('%0.1g',1e4*isoval(ii)));
end                
xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');     
xlim(plot_range(1,:))
ylim(plot_range(2,:))
zlim(plot_range(3,:))