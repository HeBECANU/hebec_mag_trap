%%% useful rotator for meshgrid
function [x,y,z]=rotmesh(Mrot,x0,y0,z0)
nn=size(x0);    % size of the meshgrid
xyz=Mrot*[x0(:)';y0(:)';z0(:)'];
x=reshape(xyz(1,:),nn);
y=reshape(xyz(2,:),nn);
z=reshape(xyz(3,:),nn);
end