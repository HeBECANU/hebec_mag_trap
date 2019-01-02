%trap angle scratch

%Have a look at a 100 um box around the trap center
% scan_range = linspace(-100*1e-6,100*1e6,25);
trap_cent=anal_out.trap_cent;
% [x y z] = meshgrid(scan_range+trap_cent(1),scan_range+trap_cent(2),scan_range+trap_cent(3));
%%
ngrid=40;          % 50 - med; 300 - very fine;
% grid in trap centered ref frame
xyz_grid=[];
xpts=trap_cent(1)+linspace(-20e-6,20e-6,ngrid);
ypts=trap_cent(2)+linspace(-20e-6,20e-6,ngrid);
zpts=trap_cent(3)+linspace(-20e-6,20e-6,ngrid);
[xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3)]=...
meshgrid(xpts,ypts,zpts);    % meshgrid
xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
[Bmag_list,Bvec_list]=trap_eval(btrap,xyz_list);
Bmag_grid=reshape(Bmag_list,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
Bx_grid=reshape(Bvec_list(:,1),[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
By_grid=reshape(Bvec_list(:,2),[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
Bz_grid=reshape(Bvec_list(:,3),[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
%%
k = [1,0,0]; %light field vector
k = k./norm(k);
theta = acos((k(1).*Bx_grid + k(2).*By_grid + k(3).*Bz_grid)./Bmag_grid);

% figure(54)
% pcolor(xpts.*1e6,ypts.*1e6,theta(:,:,end/2).*180/pi)
% xlabel('x (\mu m)')
% ylabel('y (\mu m)')
% shading interp
% colorbar
% figure(55)
% pcolor(xpts.*1e6,zpts.*1e6,squeeze(theta(:,end/2,:)).*180/pi)
% xlabel('x (\mu m)')
% ylabel('z (\mu m)')
% shading interp
% colorbar
%%
figure(56)
pcolor((xpts(:)-trap_cent(1)).*1e6,ypts(:).*1e6,theta(:,:,end/2).*180/pi)
shading interp
% hold on
% quiver(xyz_grid(:,:,end/2,1).*1e6,xyz_grid(:,:,end/2,2).*1e6,Bx_grid(:,:,end/2),By_grid(:,:,end/2))
% hold off
colormap(viridis)
colorbar
xlabel('x (\mu m)')
ylabel('y (\mu m)')
% xlim([-20 20])
% ylim([-20 20])
%caxis([0 6])
%%
figure(59)
pcolor((xpts(:)-trap_cent(1)).*1e6,zpts(:).*1e6,squeeze(theta(:,end/2,:)).*180/pi)
shading interp
% hold on
% quiver(xyz_grid(:,:,end/2,1).*1e6,xyz_grid(:,:,end/2,2).*1e6,Bx_grid(:,:,end/2),By_grid(:,:,end/2))
% hold off
colormap(viridis)
colorbar
xlabel('x (\mu m)')
ylabel('z (\mu m)')
% xlim([-20 20])
% ylim([-20 20])
%caxis([0 15])
%%
% figure(57)
% quiver3(xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3),Bx_grid,By_grid,Bz_grid)

%% 
% [bmag,bvec]=trap_eval(btrap,[3.584 0.0 0.1]*1e-3);
% [azimuth,elevation,r]=cart2sph(bvec(1),bvec(2),bvec(3));
% 
% azimuth*360/(2*pi)
% elevation*360/(2*pi)