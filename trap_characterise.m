function [f0,trap_cent,H]=trap_characterise(btrap,verbose)
% Evaluates trap frequencies and center
%
% [f0,cent]=trap_characterise(btrap,verbose)
%
% f0 is 3x1 array of trap frequencies [Hz] in x,y,z axis
% trap_cent is 1x3 array of trap center [m]
% H is a cell array - handle to figures: {trap_3d, trap_1d, trap_freq}
%
% btrap is the trap object to characterise
% verbose: set to >0 for graphics
%

%% Trap center
% trap center: point of minimum B magnitude

%%% Find trap center with current trap config
% find X to minimise $Bmag$ from trap_eval(btrap,X,0,0) - from symmetry
% initial guess param X~0.1 [mm]
% NOTE: x_cent solve in mm scale (function domain scaled to order of unity)
[x_cent,B_cent]=fminsearch(@(x) trap_eval(btrap,1e-3*x,0,0),0.1);
trap_cent=[1e-3*x_cent,0,0];    % mm-->m evaluated trap centre [m]

%%% Get B field near trap centre
% config
ngrid_trap=50;
trap_lim=[-10e-6,10e-6; -10e-6,10e-6; -10e-6,10e-6];

% build trap scale grid (~10 um each dir)
xyz_trap=cell(3,1);
for ii=1:3
    xyz_trap{ii}=trap_cent(ii)+linspace(trap_lim(ii,1),trap_lim(ii,2),ngrid_trap);
end
XYZ_trap=cell(3,1);
[XYZ_trap{1},XYZ_trap{2},XYZ_trap{3}]=meshgrid(xyz_trap{:});    % meshgrid
% permute the 3D array so that indexing goes x-y-z
XYZ_trap=cellfun(@(YXZ) permute(YXZ,[2,1,3]),XYZ_trap,'UniformOutput',false);     

% calculate magnetic field
[Bmag_trap,Bxyz_trap]=trap_eval(btrap,XYZ_trap{:});

%%% 1D Bmag profile - X,Y,Z line profile
% indices to trap centre (in grid)
[~,I0]=min(abs(Bmag_trap(:)-B_cent));
II0=zeros(1,3);        % this is ordered in YXY
[II0(1),II0(2),II0(3)]=ind2sub(size(Bmag_trap),I0);     % index to trap cent

B_trap_1d=cell(3,1);     % 1D line-profile of magnetic field potential [T]
idxcirc=[1,2,3];
Bmagcirc=Bmag_trap;      % temporary copy
for ii=1:3
    B_trap_1d{ii}=Bmagcirc(:,II0(idxcirc(2)),II0(idxcirc(3)))';     % 1xNi array 
    idxcirc=circshift(idxcirc,[-1,0]);
    Bmagcirc=permute(Bmagcirc,[2,3,1]);     % dimension circular permutation
end
clearvars Bmagcirc;

%%% Visualise trap: 3D B-isosurfaces
if verbose>0
    hfig_btrap_cent_3d=plot_B_3d(btrap,Bmag_trap,XYZ_trap);
    
    % display evaluated trap centre
    figure(hfig_btrap_cent_3d);
    hold on;
    scatter3(1e3*trap_cent(1),1e3*trap_cent(2),1e3*trap_cent(3),...
        'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none',...
        'SizeData',30,'DisplayName','Trap centre');
end

%%% Visualise B-1D: 1D B-profiles
if verbose>0
    hfig_bmag_1d=figure();
    axisstr={'X','Y','Z'};
    linestyle={'-','--',':'};
    p=[];
end
xyz_trap0=cell(3,1);
for ii=1:3
    % shift coord to approximate trap center
    xyz_trap0{ii}=xyz_trap{ii}-trap_cent(ii);       % trap center ref'd disp [m]
    if verbose>0
        hold on;
        p(ii)=plot(1e3*xyz_trap0{ii},1e4*B_trap_1d{ii},...
        'LineStyle',linestyle{ii},'LineWidth',1.5,...
        'DisplayName',axisstr{ii});
    end
end
if verbose>0
    box on;
    xlabel('Displacement [mm]');
    ylabel('$B$ [G]');
    lgd=legend(p);
    title(lgd,sprintf('(%0.2g,%0.2g,%0.2g) [mm]',1e3*trap_cent(:)));
end


%% Evaluate magnetic trap potential (1d: x-y-z)
% TODO - U_B± \propto mu_He * (±)Bmag (for weak(strong)-field seeking state)
const_prop=1;       % TODO
U_B_1d=cellfun(@(x) const_prop*x,B_trap_1d,'UniformOutput',false);


%% Trap frequency
% frequency: omega=sqrt(U''/m) [spatial 2nd order derivative]
% trap 1D: B_trap_1d - 3x1 cell array of 1D B magnitude profiles

% numerical second order differential - TODO check
d2U=cellfun(@(x)diff(x,2),U_B_1d,'UniformOutput',false);    % 2nd order DIFFERENCE in potential
dx=cellfun(@(x)diff(x),xyz_trap0,'UniformOutput',false);    % 1st ord DIFF in displacement

% evaluate trap frequency (harmonic)
m_He=1;     % TODO
f_trap=cellfun(@(D2Y,DX)(2*pi)*sqrt((D2Y./(DX(1:end-1).^2))/m_He),d2U,dx,'UniformOutput',false);     % trap frequency in [Hz]

f0=cellfun(@(x)mean(x),f_trap);     % mean trap frequency around trap centre [Hz]

%%% Visualise trap frequency - anharmonicity, etc
if verbose>0
    hfig_trap_freq=figure();
    axisstr={'X','Y','Z'};
    linestyle={'-','--',':'};
    p=[];
    for ii=1:3
        hold on;
        p(ii)=plot(1e3*xyz_trap0{ii}(2:end-1),f_trap{ii},...
            'LineStyle',linestyle{ii},'LineWidth',1.5,...
            'DisplayName',axisstr{ii});
    end
    box on;
    xlabel('Displacement [mm]');
    ylabel('$f$ [Hz]');
    lgd=legend(p);
end


%% Collate all figs
H={hfig_btrap_cent_3d, hfig_bmag_1d, hfig_trap_freq};

end