function Bout=Bfield_finite_line_numeric_order_1(tlim,dlen,curr,rot_vec,xyz)
% Bout = Bfield_finite_line
% B field calculator for single line of current (located at origin, pointing Z-axis)
% Bout is a 3x1 cell-array: {Bx,By,Bz} defined at points x,y,z

% uses sum of current sticks approach
% http://web.mit.edu/6.013_book/www/chapter8/8.2.html

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
for ii=1:size(xyz,1)
    %difference between interogation point current stick start
    bvec=xyz(ii,:)-wire_pos(1:end-1,:);
    %difference between interogation point current stick end
    cvec=xyz(ii,:)-wire_pos(2:end,:);
    %vector of the differnece between start and end
    avec=wire_pos(2:end,:)-wire_pos(1:end-1,:);
    %calculate the magnetic feild from each current stick
    % from eq 22 we caulate the feild from each stick
    c_cross_a=cross(cvec,avec,2);
    Bvec=(c_cross_a./(vecnorm(c_cross_a,2,2).^2)) .*( ( dot(avec,cvec,2)./vecnorm(cvec,2,2) ) - ( dot(avec,bvec,2)./vecnorm(bvec,2,2) ) ) ;
    %combine them all
    Bout_trapz(ii,:)=sum(Bvec,1);
end
Bout=((mu_0*curr)/(4*pi)).*Bout_trapz;

Bout=Bout*rev_rot_mat;


end
