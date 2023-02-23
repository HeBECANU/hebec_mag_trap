function [Bmag,Bout]=trap_eval(trap,xyz_list)
% [Bmag,Bout]=trap_eval(trap,x,y,z)
% Evaluate B field for a trap
% Bmag is an array (same size as x,y,z)
% Bout is a 3x1 cell-array: {Bx,By,Bz} defined at points x,y,z
%


% trap configuration
ncomps=numel(trap.b_src);
Bout=repmat([0,0,0],[size(xyz_list,1),1]);
% evaluate B field from each component at the grid region
for ii=1:ncomps
    trap_elm=trap.b_src(ii); %for each element in the cell array of structs
    elm_type=trap_elm.type;
    elm_param=trap_elm.param;
    
    % evaluate B field for this object
    switch elm_type
        case 'loop'
            % get element parameters and format them for the calculator          
            rad=elm_param.radius;
            curr=elm_param.current;
            pos=elm_param.position;
            rot_vec=elm_param.rot;
            rel_xyz=xyz_list-pos;
            % call the elemnt calculator
            Bxyz_this=Bfield_loop_analytic(rad,curr,rot_vec,rel_xyz);
        case {'num line','numeric line'}   
            tlim=[0,elm_param.length];
            curr=elm_param.current;
            pos=elm_param.position;
            dlen=elm_param.dlen;
            rot_vec=elm_param.rot;
            rel_xyz=xyz_list-pos;
            Bxyz_this=Bfield_finite_line_numeric_order_1(tlim,dlen,curr,rot_vec,rel_xyz);
        case {'line','ana line','analytic line'}       
            tlim=[0,elm_param.length];
            curr=elm_param.current;
            pos=elm_param.position;
            rot_vec=elm_param.rot;
            rel_xyz=xyz_list-pos;
            Bxyz_this=Bfield_finite_line_analytic(tlim,curr,rot_vec,rel_xyz);
        case {'helix','num helix','numeric helix'} 
            rad=elm_param.radius;
            pitch=elm_param.pitch;
            turn=[0,elm_param.turns];
            dlen=elm_param.dlen;
            curr=elm_param.current;
            pos=elm_param.position;
            rot_vec=elm_param.rot;
            rel_xyz=xyz_list-pos;
            Bxyz_this=Bfield_helix_numeric(rad,pitch,turn,dlen,curr,rot_vec,rel_xyz);
        case 'uniform'
            % uniform B field
            % params: {[Bx,By,Bz]}
            Bxyz_this=elm_param;     % get the uniform B field vector
        otherwise
            error('element type of "%s" is not recognised.',elm_type);
    end
    
    % add to total B field array
    Bout=Bout+Bxyz_this;

end
% get B field vector components
Bmag=sqrt(sum(Bout.^2,2));    % absolute magnetic field strength [T]


end