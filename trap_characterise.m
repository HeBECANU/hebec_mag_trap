function [f0,trap_cent,B_cent]=trap_characterise(btrap,x0,solve_trapdepth,verbose)
% Evaluates trap frequencies and center
%
% [f0,cent]=trap_characterise(btrap,verbose)
%
% f0 is 3x1 array of trap frequencies [Hz] in x,y,z axis
% trap_cent is 1x3 array of trap center [m]
%
% btrap is the trap object to characterise
% x0 is the initial guess for trap minimum [m]
% verbose: set to >0 for graphics
%

%% Numbers

%% Trap center
% trap center: point of minimum B magnitude

%%% Find trap center with current trap config
% find X to minimise $Bmag$ from trap_eval(btrap,X,0,0) - from symmetry
% initial guess param X~0.1 [mm]
% NOTE: x_cent solve in mm scale (function domain scaled to order of unity)
% x0=1e-3*18.5/2;     % initial estimate for trap centre (middle of two coils) 
options = optimset('TolX',1e-8);
[x_cent,B_cent]=fminsearch(@(x) trap_eval(btrap,[1e-3*x,0,0]),1e3*x0,options);
trap_cent=[1e-3*x_cent,0,0];    % mm-->m evaluated trap centre [m]
global const
 
fprintf('found trap minimum\nB=%2.3fG (%2.3fMHz) \nat {%f,%f,%f} mm\n',...
    1e4*B_cent,1e-6*B_cent*const.b_freq,trap_cent(1)*1e3,trap_cent(2)*1e3,trap_cent(3)*1e3)

if B_cent<1e-10
    warning('===============the trap minimum is close to zero=============')
end


%%% 1D Bmag profile - X,Y,Z line profile
%X slice
%select a small range
points=1000;
range=[[1,-1];[1,-1];[1,-1]]*1e-4;
labels=['x','y','z'];
figure(1)
set(gcf,'Color',[1 1 1]);
clf;
for n=1:3
xyz_points=zeros(points,3);
xyz_points(:,n)=linspace(range(n,1),range(n,2),points)';
xyz_points=trap_cent+xyz_points;
[bmag,bvec]=trap_eval(btrap,xyz_points);   
bmag=bmag-B_cent;
subplot(3,1,n)
deltx=xyz_points(:,n)-trap_cent(n);
plot(deltx,bmag) %
xlabel(labels(n))
ylabel('Bfeild')

poly=polyfit(deltx,bmag,6);
hold on
plot(deltx,polyval(poly,deltx),'r')
hold off
%hold on
dudx(n,1)=polyval(polyder(poly),0);
dudx(n,2)=polyval(polyder(polyder(poly)),0);
dudx(n,3)=polyval(polyder(polyder(polyder(poly))),0);
dudx(n,4)=polyval(polyder(polyder(polyder(polyder(poly)))),0);
dudx(n,5)=polyval(polyder(polyder(polyder(polyder(polyder(poly))))),0);
dudx(n,6)=polyval(polyder(polyder(polyder(polyder(polyder(polyder(poly)))))),0);
end
fprintf('trap curvature {%f , %f, %f} G/cm \n',dudx(1,1)*1e2,dudx(2,1)*1e2,dudx(3,1)*1e2)
fprintf('trap curvature {%f , %f, %f} G/cm^2 \n',dudx(1,2),dudx(2,2),dudx(3,2))
trap_freq=sqrt(2*const.mub*dudx(:,2)'/const.mhe)/(2*pi);
fprintf('trap freq {%f , %f, %f} \n',trap_freq(1),trap_freq(2),trap_freq(3))
fprintf('trap ratio {%f , %f}={y/x,z/x} \n',trap_freq(2)/trap_freq(1),trap_freq(3)/trap_freq(1))




if solve_trapdepth
    %find the trap depth by finding when the thresholded B feild map
    %changes the number of connected regions
    %uses the bisection method to find this potential
    %needs to have a well sampled grid (to see the other minima), a good
    %region of interest and decent gess of the inital potential
    %can do it in 2 or 3d
    %this method is a little messy but it works
    %should implement a ridge following method like
    %https://carter.princeton.edu/wp-content/uploads/sites/316/2015/08/EAC-047.pdf
    %for now just sample the potential in 3d and find the point that the
    %number of regions change
    

    high_pot=B_cent+30*1e-4;
    low_pot=B_cent+2*1e-5;
    thresh=1e-5*1e-5; %convergence threshold
    dimensions=2;
    
    
    
    
    
    
    if dimensions==2
        ngrid=1000;          % 50 - med; 300 - very fine;
        % grid in trap centered ref frame
        xyz_grid=[];
        xpts=trap_cent(1)+linspace(-10e-3,20e-3,ngrid);
        zpts=trap_cent(3)+linspace(-10e-3,10e-3,ngrid);
        [xyz_grid(:,:,:,1),xyz_grid(:,:,:,3)]=...
        meshgrid(xpts,zpts);    % meshgrid
        xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
        fprintf('trap depth calculating grid \n')
        [Bmag_list,~]=trap_eval(btrap,xyz_list);
        Bmag_grid=reshape(Bmag_list,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);

        
        figure(5)
        set(gcf,'Color',[1 1 1]);
        clf;
        bw_grid=Bmag_grid<low_pot;
        imagesc((xpts-trap_cent(1))*1e3,(zpts-trap_cent(3))*1e3,flipud(bw_grid))
        colormap(gray)
        set(gca,'YDir','normal')
        title(['Connected Regions B=',num2str(low_pot*1e4,'%2.3f'),'G'])
        xlabel('x(mm)')
        ylabel('y(mm)')
        regions=bwconncomp(bw_grid,4);
        base_regions=regions.NumObjects;
        pause(3)
        n=1;
        found_pot=1;
        while found_pot
            mid_pot=(high_pot+low_pot)/2;
            figure(6)
            set(gcf,'Color',[1 1 1]);
            clf;
            bw_grid=Bmag_grid<mid_pot; %
            imagesc((xpts-trap_cent(1))*1e3,(zpts-trap_cent(3))*1e3,flipud(bw_grid))
            xlabel('x(mm)')
            ylabel('y(mm)')
            title(['Connected Regions B=',num2str(mid_pot*1e4,'%2.3f'),'G'])
            colormap(gray)
            set(gca,'YDir','normal')
            regions=bwconncomp(bw_grid,4);
            delt_regions=regions.NumObjects-base_regions;
            if delt_regions<0
                high_pot=mid_pot;
            elseif delt_regions>0
                low_pot=mid_pot;
            elseif delt_regions==0
                low_pot=mid_pot;    
            end
            fprintf('trap depth evaluation %i at %2.3fG,delta regions %i, delta B %2.3E G \n',n,mid_pot*1e4,delt_regions,(high_pot-low_pot)*1e4)
            pause(0.5)
            if (high_pot-low_pot)<thresh 
                found_pot=0;
            end 
            n=n+1;
        end
    elseif dimensions==3
        ngrid=300;          % 50 - med; 300 - very fine;
        % grid in trap centered ref frame
        xpts=trap_cent(1)+linspace(-10e-3,20e-3,ngrid);
        ypts=trap_cent(2)+linspace(-6e-3,6e-3,ngrid);
        zpts=trap_cent(3)+linspace(-8e-3,8e-3,ngrid);
        
        xyz_grid=[];
        [xyz_grid(:,:,:,1),xyz_grid(:,:,:,2),xyz_grid(:,:,:,3)]=meshgrid(xpts,ypts,zpts);    % meshgrid
        xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
        fprintf('trap depth start calculating grid \n')
        [Bmag_list,~]=trap_eval(btrap,xyz_list);
        fprintf('trap depth done calculating grid \n')
        Bmag_grid=reshape(Bmag_list,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
        
        figure(5)
        clf;
        bw_grid=Bmag_grid<low_pot;
        isosurface(1e3*xyz_grid(:,:,:,1),1e3*xyz_grid(:,:,:,2),1e3*xyz_grid(:,:,:,3),bw_grid,[0.5]);
        colormap(gray)
        set(gcf,'Color',[1 1 1]);
        title(['Connected Regions B=',num2str(low_pot*1e4,'%2.3f'),'G'])
        xlabel('x(mm)')
        ylabel('y(mm)')
        zlabel('z(mm)')
        regions=bwconncomp(bw_grid,6);
        base_regions=regions.NumObjects;
        fprintf('trap depth base evaluation at %2.3fG, regions %i\n',low_pot*1e4,base_regions)
        pause(3)
        n=1;
        found_pot=1;
        while found_pot
            mid_pot=(high_pot+low_pot)/2;
            bw_grid=Bmag_grid<mid_pot; %
            figure(6)
            clf;
            set(gcf,'Color',[1 1 1]);
            isosurface(1e3*xyz_grid(:,:,:,1),1e3*xyz_grid(:,:,:,2),1e3*xyz_grid(:,:,:,3),bw_grid,[0.5]);
            xlabel('x(mm)')
            ylabel('y(mm)')
            zlabel('z(mm)')
            title(['Connected Regions B=',num2str(mid_pot*1e4,'%2.3f'),'G'])
            colormap(gray)
            regions=bwconncomp(bw_grid,6);
            delt_regions=regions.NumObjects-base_regions;
            if delt_regions<0
                high_pot=mid_pot;
            elseif delt_regions>0
                low_pot=mid_pot;
            elseif delt_regions==0
                low_pot=mid_pot;    
            end
            fprintf('trap depth evaluation %i at %2.3fG,delta regions %i, delta B %2.3E G \n',n,mid_pot*1e4,delt_regions,(high_pot-low_pot)*1e4)
            pause(0.01)
            if (high_pot-low_pot)<thresh 
                found_pot=0;
            end 
            n=n+1;
        end
        
        
    end
    
    trap_depth=mid_pot;
    fprintf('Trap depth found to be %2.3fG (%2.3f MHz)(%2.3EK)\n',trap_depth*1e4,trap_depth*const.b_freq*1e-6,(2*trap_depth*const.mub)*2/const.kb)
end
% Y slice
% points=100;
% y_points=trap_cent-[zeros(points,1),linspace(1E-6,-1E-6,100)', zeros(points,1)];
% [ybmag,ybvec]=trap_eval(btrap,y_points);
% ybmag=ybmag-B_cent;
% plot(y_points(:,2),k_UB*ybmag)
% 
% 
% Z slice
% points=100;
% z_points=trap_cent-[zeros(points,1), zeros(points,1),linspace(1E-6,-1E-6,100)'];
% [zbmag,zbvec]=trap_eval(btrap,z_points);
% zbmag=zbmag-B_cent;
% plot(z_points(:,3),k_UB*zbmag)
% end

%sample the trap in x,y,z




%%% Visualise trap: 3D B-isosurfaces
%if verbose>0
    % display evaluated trap centre
%    figure(hfig_btrap_cent_3d);
%    hold on;
%    scatter3(1e3*trap_cent(1),1e3*trap_cent(2),1e3*trap_cent(3),...
%        'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none',...
%        'SizeData',30,'DisplayName','Trap centre');
%end

%%% Visualise B-1D: 1D B-profiles
%if verbose>0
%    hfig_bmag_1d=figure();
%    axisstr={'X','Y','Z'};
%    linestyle={'-','--',':'};
%    p=[];
%end


% xyz_trap0=cell(3,1);
% for ii=1:3
%     shift coord to approximate trap center
%     xyz_trap0{ii}=xyz_trap{ii}-trap_cent(ii);       % trap center ref'd disp [m]
%     if verbose>0
%         hold on;
%         p(ii)=plot(1e3*xyz_trap0{ii},1e4*B_trap_1d{ii},...
%         'LineStyle',linestyle{ii},'LineWidth',1.5,...
%         'DisplayName',axisstr{ii});
%     end
% end


% if verbose>0
%     box on;
%     xlabel('Displacement [mm]');
%     ylabel('$B$ [G]');
%     lgd=legend(p);
%     title(lgd,sprintf('(%0.2g,%0.2g,%0.2g) [mm]',1e3*trap_cent(:)));
% end




% %% Trap frequency
% % frequency: omega=sqrt(U''/m) [spatial 2nd order derivative]
% % trap 1D: B_trap_1d - 3x1 cell array of 1D B magnitude profiles
% 
% % numerical second order differential - TODO check
% d2U=cellfun(@(x)diff(x,2),U_B_1d,'UniformOutput',false);    % 2nd order DIFFERENCE in potential
% dx=cellfun(@(x)diff(x),xyz_trap0,'UniformOutput',false);    % 1st ord DIFF in displacement
% 
% % evaluate trap frequency (harmonic)
% m_He=6.6465e-27;     % mass of 4-helium [kg]
% f_trap=cellfun(@(D2Y,DX)(2*pi)*sqrt((D2Y./(DX(1:end-1).^2))/m_He),d2U,dx,'UniformOutput',false);     % trap frequency in [Hz]
% 
% f0=cellfun(@(x)mean(x),f_trap);     % mean trap frequency around trap centre [Hz]
% 
% %%% Visualise trap frequency - anharmonicity, etc
% if verbose>0
%     hfig_trap_freq=figure();
%     axisstr={'X','Y','Z'};
%     linestyle={'-','--',':'};
%     p=[];
%     for ii=1:3
%         hold on;
%         p(ii)=plot(1e3*xyz_trap0{ii}(2:end-1),f_trap{ii},...
%             'LineStyle',linestyle{ii},'LineWidth',1.5,...
%             'DisplayName',axisstr{ii});
%     end
%     box on;
%     xlabel('Displacement [mm]');
%     ylabel('$f$ [Hz]');
%     lgd=legend(p);
% end
% 
f0=[trap_freq(1),trap_freq(2),trap_freq(3)];

end