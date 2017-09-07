%%% Evaluate B field for a trap
function [Bmag,Bout]=trap_eval(trap,x,y,z)
% Bmag is an array (same size as x,y,z)
% Bout is a 3x1 cell-array: {Bx,By,Bz} defined at points x,y,z
%

% trap configuration
ncomps=numel(trap);

% evaluate B field from each component at the grid region
for ii=1:ncomps
    obj_this=trap(ii);
    type_this=obj_this.type;
    param_this=obj_this.param;
    
    Bxyz_this=cell(3,1);
    
    % evaluate B field for this object
    switch type_this
        case 'coil'
            % get coil param: {R, I, POS, ORI}            
            R=param_this{1};
            I=param_this{2};
            pos=param_this{3};
            ori=param_this{4};
            
            % coord translation
            xyz={x,y,z};
            xyz_tf=cell(3,1);
            for jj=1:3
                xyz_tf{jj}=xyz{jj}-pos(jj);
            end
            
            % call the coil calculator
            Bxyz_this=Bfield_coil(R,I,xyz_tf{:});
            
        case 'uniform'
            % uniform B field
            % params: {[Bx,By,Bz]}
            Buni=param_this{1};     % get the uniform B field vector
            for jj=1:3
                Bxyz_this{jj}=Buni(jj);
            end
            
        otherwise
            error('<TRAP>.type of %s is not recognised.',string(typethis));
    end
    
    % add to total B field array
    if ii==1
        Bout=Bxyz_this;
    else
        Bout=cellfun(@(x,y)x+y,Bout,Bxyz_this,'UniformOutput',false);
    end
end
% get B field vector components
Bmag=sqrt(esumsqr(Bout{:}));     % absolute magnetic field strength [T]
end