%%% Evaluate B field for a trap
function [Bxx,Byy,Bzz,Bmag]=trap_eval(trap,x,y,z)
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
            [Bxyz_this{1},Bxyz_this{2},Bxyz_this{3}]=Bfield_coil(R,I,xyz_tf{:});
            
        otherwise
            error('<TRAP>.type of %s is not recognised.',string(typethis));
    end
    
    % add to total B field array
    if ii==1
        Bxyz=Bxyz_this;
    else
        Bxyz=cellfun(@(x,y)x+y,Bxyz,Bxyz_this,'UniformOutput',false);
    end
end
% get B field vector components
Bxx=Bxyz{1};
Byy=Bxyz{2};
Bzz=Bxyz{3};
Bmag=sqrt(Bxx.^2+Byy.^2+Bzz.^2);     % absolute magnetic field strength [T]
end