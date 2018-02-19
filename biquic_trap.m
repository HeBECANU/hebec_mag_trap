function btrap=biquic_trap(v_quad,v_shunt,Bext)
% btrap = biquic_trap(Iq, Ib, Bext)
%
% btrap is biquic magnetic trap object created by parameters defined in the
% inputs.
%
% Iq, Ib are current [Amp]
% Bext is 1x3 array of bias magnetic field at the trap (Bx,By,Bz) [T]
%
%
%%% Biquic geometry
% Quad: 14 mm Dia x 10 turns
% Bias (Ioffe): 14 mm x 18 turns 
% pitch=500 um (wire dia)
% AH sep ~17 mm
% quad-bias sep = 18.5 mm
% Current: ~36 amp (max)
% 

%%Current

amp_per_volt_quad=4.17;
amp_per_volt_shunt=1.46;

Iquad=amp_per_volt_quad*v_quad-amp_per_volt_shunt*v_shunt;
Ibias=amp_per_volt_quad*v_quad;


%% configure coil params
pitch_coil=0.5e-3;


% coil diameters
Dquad=14e-3+pitch_coil; %
Dbias=14e-3+pitch_coil;

% coil radii
Rquad=Dquad/2;
Rbias=Dbias/2;

% turns
Nturnbias=18;
Nturnquad=10;

disp_ah=16.01e-3;      % AH separation (across the chamber)
disp_qb=18.5e-3;    % Quad-bias separation (in the plane of the bi-quic)

%disp_ah=16.75e-3;      % AH separation (across the chamber)
%disp_qb=18.5e-3;    % Quad-bias separation (in the plane of the bi-quic)

% BiQUIC paper
%disp_ah=17e-3;      % AH separation (across the chamber)
%disp_qb=18.5e-3;    % Quad-bias separation (in the plane of the bi-quic)
% % From A.G. Manning PhD thesis
% disp_ah=9e-3;      % AH separation
% disp_qb=17e-3;    % Quad-bias separation

%Iquad=Iq;
%Ibias=Ib;

% Trap external bias (nuller) [http://dx.doi.org/10.1063/1.2472600]

%%% Build trap
% Quadrupole - ref
quad_coil.type='coil';
quad_coil.param={Rquad,Iquad,[0,0,disp_ah/2],[0,0,0]};
% Shunt (Ioffe) - ref
bias_coil.type='coil';
bias_coil.param={Rbias,Ibias,[disp_qb,0,disp_ah/2],[0,0,0]};
% Bias (nuller)
extbias.type='uniform';
extbias.param={Bext};

btrap=[];
for ii=1:Nturnquad
    quad_coil_temp=quad_coil;
    % shift by wire pitch
    quad_coil_temp.param{3}=quad_coil_temp.param{3}+(ii-1)*[0,0,pitch_coil];
    btrap=[btrap,quad_coil_temp];
    
    % anti-helmholtz pair - mirror symmetry around Z
    quad_coil_temp.param{2}=-1*quad_coil_temp.param{2};     % flip current dir
    quad_coil_temp.param{3}=[1,1,-1].*quad_coil_temp.param{3};
    btrap=[btrap,quad_coil_temp];
end
for ii=1:Nturnbias
    bias_coil_temp=bias_coil;
    % shift by wire pitch
    bias_coil_temp.param{3}=bias_coil_temp.param{3}+(ii-1)*[0,0,pitch_coil];
    btrap=[btrap,bias_coil_temp];
    
    % anti-helmholtz pair - mirror symmetry around Z
    bias_coil_temp.param{2}=-1*bias_coil_temp.param{2};   % flip current dir
    bias_coil_temp.param{3}=[1,1,-1].*bias_coil_temp.param{3};
    btrap=[btrap,bias_coil_temp];
end
btrap=[btrap,extbias];

end