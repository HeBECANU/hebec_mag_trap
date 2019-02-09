tic
ngrid=30;          % 50 - med; 300 - very fine;
% grid in trap centered ref frame
xyz=cell(3,1);
%%% MACRO
xyz{1}=linspace(-20e-3,20e-3,3);      % x-vect
xyz{2}=linspace(-15e-3,15e-3,ngrid);      % y-vect
xyz{3}=linspace(-120e-3,150e-3,ngrid);      % z-vect
% %%% MICRO
% xyz{1}=linspace(-1e-3,1e-3,ngrid);      % x-vect
% xyz{2}=linspace(-1e-3,1e-3,ngrid);      % y-vect
% xyz{3}=linspace(-1e-3,1e-3,ngrid);      % z-vect
XYZ=cell(3,1);
[XYZ{1},XYZ{2},XYZ{3}]=meshgrid(xyz{:});    % meshgrid
% permute the 3D array so that indexing goes x-y-z
%XYZ=cellfun(@(YXZ) permute(YXZ,[2,1,3]),XYZ,'UniformOutput',false);
xyz_list = [reshape(XYZ{1},[],1),reshape(XYZ{2},[],1),reshape(XYZ{3},[],1)];

r = 1e-3;
N = 100;
c = 0.1e-3;
%B_OUT=helix_coil(N,r,c,xyz_list);
t = linspace(-5e-3,5e-3)';
l = [t,zeros(length(t),1),zeros(length(t),1)+0.005];
dl = [ones(length(t),1),zeros(length(t),1),zeros(length(t),1)];
B_OUT=biot_savart(l,dl,t,xyz_list);
t = 0:0.1:2*pi*N;
x_coil = r*cos(t);
y_coil = r*sin(t);
z_coil = c*t;
%%
B_X = reshape(B_OUT(:,1),ngrid,3,ngrid);
B_Y = reshape(B_OUT(:,2),ngrid,3,ngrid);
B_Z = reshape(B_OUT(:,3),ngrid,3,ngrid);
B_mag=sqrt(B_X.^2+B_Y.^2+B_Z.^2);

B_X = B_X./sqrt(B_X.^2+B_Y.^2+B_Z.^2);
B_Y = B_Y./sqrt(B_X.^2+B_Y.^2+B_Z.^2);
B_Z = B_Z./sqrt(B_X.^2+B_Y.^2+B_Z.^2);

X=XYZ{1};
Y=XYZ{2};
Z=XYZ{3};
%%
%clf
figure
quiver3(X,Y,Z,B_X,B_Y,B_Z)
hold on
% plot3(x_coil,y_coil,z_coil)
plot3(l(:,1),l(:,2),l(:,3))
figure
pcolor(squeeze(Y(:,2,:)),squeeze(Z(:,2,:)),squeeze((B_mag(:,2,:))))
hold on
quiver(squeeze(Y(:,2,:)),squeeze(Z(:,2,:)),squeeze(B_Y(:,2,:)),squeeze(B_Z(:,2,:)),'w')
toc