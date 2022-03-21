function anal_out=trap_characterise(anal_out,btrap,x0,solve_trapdepth,verbose)
% Evaluates trap frequencies and center
% trap_characterise(btrap,verbose)
%
% btrap is the trap object to characterise
% x0 is the initial guess for trap minimum [m]
% verbose: set to >0 for graphics

global const

%% Trap center
% trap center: point of minimum B magnitude
%%% Find trap center with current trap config
% find X to minimise $Bmag$ from trap_eval(btrap,x_vec) 
options = optimset('TolX',1e-9);
anal_out.trap_cen.pos=fminsearch(@(x) trap_eval(btrap,x),x0,options);

%calculate the trap center field
[anal_out.trap_cen.b_mag,anal_out.trap_cen.b_vec]=trap_eval(btrap,anal_out.trap_cen.pos);
 
fprintf('found trap minimum\nB=%2.3fG (%2.3fMHz) \nat {%f,%f,%f} mm\n',...
    1e4*anal_out.trap_cen.b_mag,1e-6*anal_out.trap_cen.b_mag*const.b_freq,...
    anal_out.trap_cen.pos*1e3)

if anal_out.trap_cen.b_mag<1e-10 %this checks that the minimum is not too close to a zero crossing
    warning('===============the trap minimum found is close to zero=============')
end


%%% 1D Bmag profile - X,Y,Z line profile
%X slice
%select a small range
anal_out = mag_profile_1d(anal_out,btrap);


if solve_trapdepth
    %find the trap depth by finding when the thresholded B field map
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


end