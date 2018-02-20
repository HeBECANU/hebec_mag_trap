%Trap Depth Play
ngrid=1000;
plot3d=1;
plot2d=1;
rmax=2;
%try this ridge folowing algorithm


[x_grid,y_grid]=...
meshgrid(linspace(-rmax,rmax,ngrid),...
         linspace(-rmax,rmax,ngrid));    % meshgrid
xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
u_grid=trap(x_grid,y_grid);
%[Bmag_list,Bxyz]=trap;
%u_grid=reshape(Bmag_list,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
   
if plot3d
    figure(1)
    clf
    h=surf(x_grid,y_grid,u_grid,'facealpha',0.9);
    set(gcf,'Color',[1 1 1]);
    colormap(magma())
    set(h,'LineStyle','none')
    xlabel('x')
    ylabel('y')
    hold on
   
    scatter3(pt(1),pt(2),trap(pt(1),pt(2)),20,'k','filled')
    hold off
end

contor_start=0.8;
%if plot2d
    figure(2)
    clf;
    contour(x_grid,y_grid,u_grid,5)
    hold on
    [c,h]=contour(x_grid,y_grid,u_grid,'LevelList',[contor_start],'LineColor','k','LevelListMode','manual');
    hold off
%end

    %scatter(st_pt(1),st_pt(2),'k','filled')  
high_pot=0.8;
low_pot=0.5;
found_pot=1;
thresh=1e-6;
while found_pot
mid_pot=(high_pot+low_pot)/2;
figure(3)
clf;
bw_grid=u_grid<mid_pot;
imagesc(bw_grid)
regions=bwconncomp(bw_grid);
if regions.NumObjects>1
    low_pot=mid_pot;
else
    high_pot=mid_pot;
end
pause(0.1)
if (high_pot-low_pot)<thresh 
    found_pot=0;
end

end


function u=trap(x,y)
sigmax1=2;
sigmay1=2;
amp1=1;
sigmax2=0.5;
sigmay2=0.5;
amp2=0.5;
xcen1=0.3;
ycen1=-0.3;
u=amp1*exp(-((x.^2/(2*sigmax1^2))+(y.^2/(2*sigmay1^2))))...
    -amp2*exp(-(((x-xcen1).^2/(2*sigmax2^2))+((y-ycen1).^2/(2*sigmay2^2))));
end