function Bout=Bfield_finite_line_numeric_order_0(tlim,dlen,curr,rot_vec,xyz)
% Bout = Bfield_finite_line
% B field calculator for single line of current (located at origin, pointing Z-axis)
% Bout is a 3x1 cell-array: {Bx,By,Bz} defined at points x,y,z
% 20190210 validated to give the same answer as Bfield_path_numerical 
% uses direct integration of the Biot–Savart law which makes it fairly slow to converge
% for a line you should not use this method it is more a validation of the approach
% is slightly 

%http://web.mit.edu/viz/EM/visualizations/coursenotes/modules/guide09.pdf
% physical constants
%global const
%mu_0=const.mu0;  
%hack to get to work in parfor
mu_0=1.2566370614*10^-6;

rot_mat=rotationVectorToMatrix(rot_vec);
rev_rot_mat=rotationVectorToMatrix(-rot_vec); %reverse rotation vector
xyz=xyz*rot_mat;

Bout_trapz=xyz*nan; %initalize
tvec=linspace(tlim(1),tlim(2),ceil(range(tlim)/dlen))';
%define the parametric curve (with t is the parameter) for the current
wire_pos=tvec*[0,0,1];
%define the tangent to this curve
dLdt=tvec*0+[0,0,1];
for ii=1:size(xyz,1)
    %difference between interogation point and the wire position
    Rvec=xyz(ii,:)-wire_pos;
    %calculate the infentesimal b field vector
    dBvec=cross(dLdt,Rvec)./(vecnorm(Rvec,2,2).^3);
    %combine them all
    %Bout_trapz(ii,:)=trapz(tvec,dBvec);
    Bout_trapz(ii,:)=trapz(tvec,dBvec);
end
Bout=((mu_0*curr)/(4*pi)).*Bout_trapz;

Bout=Bout*rev_rot_mat;


end


% tic
% dL=@(t)  [0,0,t];
% %define the tangent to this curve
% dLdt=@(t) [0,0,1];
% %difference between interogation point and 
% Rvec=@(t) xyz(ii,:)-dL(t);
% dBvec=@(t) cross(dLdt(t),Rvec(t))./(norm(Rvec(t)).^3);
% 
% Bout_inbuitlt_integeral(ii,:)=integral(dBvec,tlim(1),tlim(2),'ArrayValued',1);
% if ~isequal(Bout_trapz,Bout_inbuitlt_integeral)
%     Bout_trapz-Bout_inbuitlt_integeral
%     warning('error methods do not match')
% end