function [f0,trap_cent,Bmag,XYZ,btrap]=biquic_calculator(Iquad,Ibias,Bext,x0,verbose)
% [trap_freq,trap_cent,B,XYZ,btrap] = biquic_calculator(Iquad,Ibias,Bext,x0,verbose)
%
% Function calculates BiQUIC trap properties at user defined coil current
% settings and external bias field
%
% NOTE: user should supply an estimate for trap centre along x-axis (x0)
% [m] - otherwise an incrementally search from a known solution...
%   should lie somewhere between 0mm and 10mm for practical traps
%

%% Function config
%%% 3D grid trap B-field
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


%% Build BiQUIC trap
btrap=biquic_trap(Iquad,Ibias,Bext);  % build biquic


%% characterise trap
% evaluate trap center, 1D trap potential, trap frequencies
% NOTE: there are multiple points of potential minimia
[f0,trap_cent,H_trap]=trap_characterise(btrap,x0,verbose);


%% Solve B-field for trap in 3D
% initialise
hfig_btrap=NaN;

% 3D B-field
[Bmag,Bxyz]=trap_eval(btrap,XYZ{:});    % calculate magnetic field (macro summary)

if verbose>0
    %%% visualise
    % 3D B-isosurfaces
    hfig_btrap=plot_B_3d(btrap,Bmag,XYZ);
    
    % display evaluated trap centre
    scatter3(1e3*trap_cent(1),1e3*trap_cent(2),1e3*trap_cent(3),...
        'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none',...
        'SizeData',30,'DisplayName','Trap centre');
end

end