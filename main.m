%% BiQUIC trap magnetic field calculator
% D K Shin
% refer to https://doi.org/10.1103/PhysRevA.35.1535 for expressions used

clear all;

t_start=tic;

% TODO:
% [] Voltage to current conversion
% [] experiment params

%% Config
%%% Flags
verbose=1;          % graphical output

solve_3D=1;         % solve full 3D vector B-field (takes a while)
solve_trapchar=1;   % characterise trap params including freq and center

%%% magneitc trap
Vquad=10;
Vshunt=10;
% Vquad=linspace(0,1,10);      % control voltage
% Vshunt=linspace(0,1,10);

k_v2i=1;    % conversion factor for control V [V] to I [A] at coils
% TODO - currents in BiQUIC isn't actually constructed like this...
Iquad=k_v2i*Vquad;
Ishunt=k_v2i*Vshunt;    % TODO - this coil labelled "BIAS"
Bbias=1e-4*[0,0,0];     % external bias field [T] (uniform assumption)

%%% 3D grid trap B-field
if solve_3D>0
    ngrid=50;          % 50 - med; 300 - very fine;
    
    % grid in trap centered ref frame
    xyz=cell(3,1);
    %%% MACRO
    xyz{1}=linspace(-20e-3,20e-3,ngrid);      % x-vect
    xyz{2}=linspace(-20e-3,20e-3,ngrid);      % y-vect
    xyz{3}=linspace(-20e-3,20e-3,ngrid);      % z-vect
    
    % %%% MICRO
    % xyz{1}=linspace(-1e-3,1e-3,ngrid);      % x-vect
    % xyz{2}=linspace(-1e-3,1e-3,ngrid);      % y-vect
    % xyz{3}=linspace(-1e-3,1e-3,ngrid);      % z-vect
    
    XYZ=cell(3,1);
    [XYZ{1},XYZ{2},XYZ{3}]=meshgrid(xyz{:});    % meshgrid
    % permute the 3D array so that indexing goes x-y-z
    XYZ=cellfun(@(YXZ) permute(YXZ,[2,1,3]),XYZ,'UniformOutput',false);
end

% %% Create meshgrid of trap configs
% [IQ,IS]=meshgrid(Iquad,Ishunt);
% I_biquic=[IQ(:),IS(:)]';

%% Build BiQUIC trap
btrap=biquic_trap(Iquad,Ishunt,Bbias);  % build biquic

%% characterise trap
% initialise
f0=NaN;
trap_cent=NaN;
H_trap=NaN;

if solve_trapchar>0
    % evaluate trap center, 1D trap potential, trap frequencies
    % NOTE: there are multiple points of potential minimia
    [f0,trap_cent,H_trap]=trap_characterise(btrap,1e-3,verbose);
end

%% Solve B-field for trap in 3D
% initialise
Bmag=NaN;
Bxyz=NaN;
hfig_btrap=NaN;

if solve_3D>0
    % 3D B-field
    [Bmag,Bxyz]=trap_eval(btrap,XYZ{:});    % calculate magnetic field (macro summary)
    
    %%% visualise
    % 3D B-isosurfaces
    if verbose>0
        hfig_btrap=plot_B_3d(btrap,Bmag,XYZ);
        
        if solve_trapchar>0
            % display evaluated trap centre
            scatter3(1e3*trap_cent(1),1e3*trap_cent(2),1e3*trap_cent(3),...
                'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none',...
                'SizeData',30,'DisplayName','Trap centre');
        end
    end
end

%% End of code
t_end=toc(t_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_end);
disp('===================ALL TASKS COMPLETED===================');