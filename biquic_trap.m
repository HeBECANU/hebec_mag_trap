function btrap=biquic_trap(Iq,Ish,Bbias)
% btrap = biquic_trap(Iq, Ish, Bbias)
%
% btrap is biquic magnetic trap object created by parameters defined in the
% inputs.
%
% Iq, Ish are current [Amp]
% Bbias is 1x3 array of bias magnetic field at the trap (Bx,By,Bz) [T]
%
%
%%% Biquic geometry
% Quad: 14 mm Dia x 10 turns
% Shunt (Ioffe): 14 mm x 18 turns 
% pitch=500 um (wire dia)
% AH sep ~17 mm
% quad-shunt sep = 18.5 mm
% Current: ~36 amp (max)
% 

%% configure coil params
% coil diameters
Dquad=14e-3;
Dshunt=14e-3;

% coil radii
Rquad=Dquad/2;
Rshunt=Dshunt/2;

% turns
Nturnshunt=18;
Nturnquad=10;

pitch_coil=0.5e-3;

% BiQUIC paper
disp_ah=17e-3;      % AH separation
disp_qs=18.5e-3;    % Quad-shunt separation
% % From A.G. Manning PhD thesis
% disp_ah=9e-3;      % AH separation
% disp_qs=17e-3;    % Quad-shunt separation

Iquad=Iq;
Ishunt=Ish;

% Trap bias (nuller) [http://dx.doi.org/10.1063/1.2472600]

%%% Build trap
% Quadrupole - ref
quad_coil.type='coil';
quad_coil.param={Rquad,Iquad,[0,0,disp_ah/2],[0,0,0]};
% Shunt (Ioffe) - ref
shunt_coil.type='coil';
shunt_coil.param={Rshunt,Ishunt,[disp_qs,0,disp_ah/2],[0,0,0]};
% Bias (nuller)
bias.type='uniform';
bias.param={Bbias};

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
for ii=1:Nturnshunt
    shunt_coil_temp=shunt_coil;
    % shift by wire pitch
    shunt_coil_temp.param{3}=shunt_coil_temp.param{3}+(ii-1)*[0,0,pitch_coil];
    btrap=[btrap,shunt_coil_temp];
    
    % anti-helmholtz pair - mirror symmetry around Z
    shunt_coil_temp.param{2}=-1*shunt_coil_temp.param{2};   % flip current dir
    shunt_coil_temp.param{3}=[1,1,-1].*shunt_coil_temp.param{3};
    btrap=[btrap,shunt_coil_temp];
end
btrap=[btrap,bias];

end