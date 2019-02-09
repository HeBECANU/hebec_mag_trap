function btrap=biquic_trap(btrap,trap_config)
% btrap = biquic_trap(Iq, Ib, Bext)
% Calculates the correct elements for a biquic trap and appends them to the btrap structure.

% Input
%   btrap array of structures which sepcify the magnetic feild
%       btrap(n).type string,element specifier 'coil' or 'uniform
%       btrap(n).param cell array of coil parameters
%           for type='coil' {Radius,Current,[x,y,z],[0,0,0]} forth element reserved for angle
%           for type='uniform' {[Bx,By,Bz]}
%   trap_config.v_quad  scalar value, of the DAQ voltage for the trap current [volts]
%   trap_config.v_shunt
%   trap_config.Bext a 1x3 vector of bias magnetic field applied by the nuller (Bx,By,Bz) [T]

% Output
%   btrap array of structures which sepcify the magnetic feild. Appended to by this code.
% 


%%% Biquic geometry
% trap: 14 mm Dia x 10 turns
% Bias (Ioffe): 14 mm x 18 turns 
% pitch=500 um (wire dia)
% AH sep ~17 mm
% trap-bias sep = 18.5 mm
% Current: ~36 amp (max)
% 

%%Current
amp_per_volt_trap=4.17;
amp_per_volt_shunt=1.46;

Itrap=amp_per_volt_trap*trap_config.v_quad;
Ibias=amp_per_volt_trap*trap_config.v_quad+amp_per_volt_shunt*trap_config.v_shunt;

%v_quad,v_shunt,Bext)



%% configure coil params
pitch_coil=0.5e-3;


% coil diameters
Dtrap=13.94e-3+2*0.7e-3+pitch_coil; %
Dbias=12e-3+2*0.7e-3+pitch_coil;

% coil radii
Rtrap=Dtrap/2;
Rbias=Dbias/2;

% turns
Nturnbias=18;
Nturntrap=10;

%measurments

disp_ah=21.5e-3;      % AH separation (across the chamber)
disp_qb=19.43e-3;    % trap-bias separation (in the plane of the bi-quic)

supply_wire_length=400e-3; %supply wire length

res_trap=49.3e-3;
res_baias=71.8e-3;



%manual playing
%disp_ah=16.01e-3;      % AH separation (across the chamber)
%disp_qb=18.5e-3;    % trap-bias separation (in the plane of the bi-quic)


%disp_ah=16.75e-3;      % AH separation (across the chamber)
%disp_qb=18.5e-3;    % trap-bias separation (in the plane of the bi-quic)

% BiQUIC paper
%disp_ah=17e-3;      % AH separation (across the chamber)
%disp_qb=18.5e-3;    % trap-bias separation (in the plane of the bi-quic)
% coil diameters
%Dtrap=14e-3+pitch_coil; %
%Dbias=14e-3+pitch_coil;

% % From A.G. Manning PhD thesis
% disp_ah=9e-3;      % AH separation
% disp_qb=17e-3;    % trap-bias separation

%Itrap=Iq;
%Ibias=Ib;



fprintf('trap current %2.2f A, power %2.1f W\n',Itrap,Itrap^2*res_trap)
fprintf('bais current %2.2f A, power %2.1f W\n',Ibias,Ibias^2*res_baias)
fprintf('Total power %2.1f W\n',Itrap^2*res_trap+Ibias^2*res_baias)
btrap={};
btrap.power.trap=Itrap^2*res_trap;
btrap.power.bais=Ibias^2*res_baias;
btrap.power.total=Itrap^2*res_trap+Ibias^2*res_baias;

%%% Build trap
% traprupole - ref
trap_coil=[];
trap_coil.type='loop';
trap_coil.param.radius=Rtrap;
trap_coil.param.current=Itrap;
trap_coil.param.position=[-disp_qb/2,disp_ah/2,0];
trap_coil.param.rot=pi/2*[1,0,0]; %angle of pointing in theta,phi
% Shunt (Ioffe) - ref
bias_coil=[];
bias_coil.type='loop';
bias_coil.param.radius=Rbias;
bias_coil.param.current=Ibias;
bias_coil.param.position=[disp_qb/2,disp_ah/2,0];
bias_coil.param.rot=pi/2*[1,0,0]; %angle of pointing in theta,phi
% Trap external bias (nuller) [http://dx.doi.org/10.1063/1.2472600]
extbias.type='uniform';
extbias.param=trap_config.Bext;

btrap.b_src=[];
btrap.b_src=extbias;
%add the coils
for ii=1:Nturntrap
    trap_coil_temp=trap_coil;
    % shift by wire pitch
    trap_coil_temp.param.position=trap_coil_temp.param.position+(ii-1)*[0,pitch_coil,0];
    btrap.b_src=[btrap.b_src,trap_coil_temp];
    
    % anti-helmholtz pair - mirror symmetry around Y
    trap_coil_temp.param.current=-1*trap_coil_temp.param.current;     % flip current dir
    trap_coil_temp.param.position=[1,-1,1].*trap_coil_temp.param.position; %flip position
    btrap.b_src=[btrap.b_src,trap_coil_temp];
end
for ii=1:Nturnbias
    bias_coil_temp=bias_coil;
    % shift by wire pitch
    bias_coil_temp.param.position=bias_coil_temp.param.position+(ii-1)*[0,pitch_coil,0];
    btrap.b_src=[btrap.b_src,bias_coil_temp];
    
    % anti-helmholtz pair - mirror symmetry around Y
    bias_coil_temp.param.current=-1*bias_coil_temp.param.current;   % flip current dir
    bias_coil_temp.param.position=[1,-1,1].*bias_coil_temp.param.position;
    btrap.b_src=[btrap.b_src,bias_coil_temp];
end

%add the wires for the bias col
%TO CHECK THE CURRENTS ARE GOING IN THE RIGHT DIRECTION!
bias_supply_wires_top_pos=[];
bias_supply_wires_top_pos.type='line';
bias_supply_wires_top_pos.param.current=bias_coil.param.current;
bias_supply_wires_top_pos.param.position=bias_coil.param.position+[0,0,bias_coil.param.radius];
bias_supply_wires_top_pos.param.rot=pi/2*[1,0,0];
bias_supply_wires_top_pos.param.length=supply_wire_length;
bias_supply_wires_btm_pos=bias_supply_wires_top_pos;
bias_supply_wires_btm_pos.param.position=bias_coil.param.position+[0,(Nturnbias-1)*pitch_coil,-bias_coil.param.radius];
bias_supply_wires_btm_pos.param.current=-bias_supply_wires_top_pos.param.current;
bias_supply_wires_btm_neg=bias_supply_wires_btm_pos;
bias_supply_wires_top_neg=bias_supply_wires_top_pos;
bias_supply_wires_btm_neg.param.position=[1,-1,1].*bias_supply_wires_btm_neg.param.position;
bias_supply_wires_btm_neg.param.rot=-1*bias_supply_wires_btm_neg.param.rot;
bias_supply_wires_top_neg.param.position=[1,-1,1].*bias_supply_wires_top_neg.param.position;
bias_supply_wires_top_neg.param.rot=-1*bias_supply_wires_top_neg.param.rot;
btrap.b_src=[btrap.b_src,bias_supply_wires_top_pos,bias_supply_wires_top_neg,bias_supply_wires_btm_pos,bias_supply_wires_btm_neg];


trap_supply_wires_top_pos=[];
trap_supply_wires_top_pos.type='line';
trap_supply_wires_top_pos.param.current=trap_coil.param.current;
trap_supply_wires_top_pos.param.position=trap_coil.param.position+[0,0,trap_coil.param.radius];
trap_supply_wires_top_pos.param.rot=pi/2*[1,0,0];
trap_supply_wires_top_pos.param.length=supply_wire_length;
trap_supply_wires_btm_pos=trap_supply_wires_top_pos;
trap_supply_wires_btm_pos.param.position=trap_coil.param.position+[0,(Nturntrap-1)*pitch_coil,-trap_coil.param.radius];
trap_supply_wires_btm_pos.param.current=-trap_supply_wires_top_pos.param.current;
trap_supply_wires_btm_neg=trap_supply_wires_btm_pos;
trap_supply_wires_top_neg=trap_supply_wires_top_pos;
trap_supply_wires_btm_neg.param.position=[1,-1,1].*trap_supply_wires_btm_neg.param.position;
trap_supply_wires_btm_neg.param.rot=-1*trap_supply_wires_btm_neg.param.rot;
trap_supply_wires_top_neg.param.position=[1,-1,1].*trap_supply_wires_top_neg.param.position;
trap_supply_wires_top_neg.param.rot=-1*trap_supply_wires_top_neg.param.rot;
btrap.b_src=[btrap.b_src,trap_supply_wires_top_pos,trap_supply_wires_top_neg,trap_supply_wires_btm_pos,trap_supply_wires_btm_neg];


btrap.power={};
btrap.power.trap=Itrap^2*res_trap;
btrap.power.bias=Ibias^2*res_baias;
btrap.power.total=Itrap^2*res_trap+Ibias^2*res_baias;
end