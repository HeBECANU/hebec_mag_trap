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
            % get coil param: {R, I, POS, ORI}            
            R=elm_param.radius;
            I=elm_param.current;
            pos=elm_param.position;
            rot_vec=elm_param.rot;
            %=param_this{4};
            rel_xyz=xyz_list-pos;
            % call the coil calculator
            Bxyz_this=Bfield_loop(R,I,rot_vec,rel_xyz);
        case 'line'
            error('not yet implemented')
        case 'spiral'
            error('not yet implemented')
            Bxyz_this=Bfield_coil(R,I,rel_xyz);
        case 'uniform'
            % uniform B field
            % params: {[Bx,By,Bz]}
            Bxyz_this=elm_param;     % get the uniform B field vector
        otherwise
            error('<TRAP>.type of %s is not recognised.',string(typethis));
    end
    
    % add to total B field array
    Bout=Bout+Bxyz_this;

end
% get B field vector components
Bmag=sqrt(sum(Bout.^2,2));    % absolute magnetic field strength [T]


end