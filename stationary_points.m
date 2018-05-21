function stationary_points(btrap,trap_cent,solve_stpt)

    fprintf('starting derivative plot \n')
    %plot_range=[[-2,12];[-8.5,8.5]]*1e-3; %x,z
     plot_range=[[-2,12];[-8.5,8.5]]*1e-3; %x,z
    min_cluster_range=1e-5;
    st_cluster_range=1e-5; %saddle pt cluster range
    zero_feild_thresh=1e-7; %1uT is when we say it is a zero crossing minima
    %2d derivative grid
    ngrid=100;          % 50 - med; 300 - very fine;
    delt=1e-10;

    
    
    % grid in trap centered ref frame
    plot_range=[plot_range(1,:)+trap_cent(1);plot_range(2,:)+trap_cent(3)];%center on trap cent
    xyz_grid=[];
    [xyz_grid(:,:,:,1),xyz_grid(:,:,:,3)]=...
    meshgrid(linspace(plot_range(1,1),plot_range(1,2),ngrid),...
             linspace(plot_range(2,1),plot_range(2,2),ngrid));    % meshgrid

    xyz_list=reshape(xyz_grid,[size(xyz_grid,1)*size(xyz_grid,2)*size(xyz_grid,3),3]);
    Bgrad_list=num_grad(@(x) trap_eval(btrap,x),xyz_list,delt);
    
    Bgrad_mag=sum(abs(Bgrad_list),2);
    Bgrad_grid=reshape(Bgrad_mag,[size(xyz_grid,1),size(xyz_grid,2),size(xyz_grid,3)]);
    figure(12)
    clf;
    h=surf((xyz_grid(:,:,:,1)-trap_cent(1))*1e3,(xyz_grid(:,:,:,3)-trap_cent(3))*1e3,Bgrad_grid,'facealpha',0.9);
    colormap(viridis())
    set(h,'LineStyle','none')
    box on;
    xlabel('X (mm)');
    ylabel('Z (mm)');
    title('Jacobian')
    zlabel('RSS Grad. (T/m) ');
    set(gcf,'Color',[1 1 1]);
    pause(0.001);

    figure(14)
    clf
    contour((xyz_grid(:,:,:,1)-trap_cent(1))*1e3,(xyz_grid(:,:,:,3)-trap_cent(3))*1e3,Bgrad_grid,300)
    set(gcf,'Color',[1 1 1]);
    title('Jacobian')
    xlabel('X (mm)');
    ylabel('Z (mm)');
    pause(1e-3)
    if solve_stpt>1
        fprintf('starting optimization on potential landscape \n')
        %find the minimum points of the potential landscape
        min_pts=[];
        num_pts=30; 
        %rng(123456) %for protoryping
        x0=[plot_range(1,1)+rand(1,num_pts)'*(plot_range(1,2)-plot_range(1,1)),plot_range(2,1)+rand(1,num_pts)'*(plot_range(2,2)-plot_range(2,1))];
        for n=1:num_pts
            options=optimset('MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-10);
            min_pts(n,:)=fminsearch(@(x) trap_eval(btrap,[[x(1),0,x(2)]]),x0(n,:),options);
        end
        %filter the points that have stayed in the range
        mask=min_pts(:,1)>plot_range(1,1) & min_pts(:,1)<plot_range(1,2) & min_pts(:,2)>plot_range(2,1) & min_pts(:,2)<plot_range(2,2);
        st_pts=min_pts(mask,:);
        fprintf('surviving points %i out of %i \n',size(st_pts,1),num_pts)
        %custer the points
        % now cluster points
        %break this clustering out as a function once 3d
        min_pt_clusters={};
        clust_idx=1;
        while size(min_pts,1)>0
            mask=sum((min_pts(1,:)-min_pts).^2,2)<min_cluster_range; %are the points within some threshold distance
            min_pt_clusters{clust_idx}=min_pts(mask,:);
            min_pts=min_pts(mask~=1,:);
            clust_idx=clust_idx+1;
        end
        min_pts=[];
        for n=1:size(min_pt_clusters,2) 
            min_pts(n,:)=mean(min_pt_clusters{n},1);
        end
        min_pts(:,3)=trap_eval(btrap,[min_pts(:,1),zeros(size(min_pts,1),1),min_pts(:,2)]); 
        min_pts(:,4)=sum(abs(num_grad(@(x)trap_eval(btrap,x),[min_pts(:,1),zeros(size(min_pts,1),1),min_pts(:,2)],delt)),2);
        fprintf('minima clustered to %i points, %i zero crossing \n',size(min_pts,1),sum(min_pts(:,3)<zero_feild_thresh))

        figure(14)
        hold on
        scatter3((min_pts(:,1)-trap_cent(1))*1e3,(min_pts(:,2)-trap_cent(3))*1e3, min_pts(:,4),20,'r','filled')
        pause(1e-3)
        hold off       
              
    end
    
    if solve_stpt>2
        fprintf('starting optimization on derivative landscape \n')        
        st_pts=[];
        num_pts=30;
        %rng(123456) %for protoryping
        x0=[plot_range(1,1)+rand(1,num_pts)'*(plot_range(1,2)-plot_range(1,1)),plot_range(2,1)+rand(1,num_pts)'*(plot_range(2,2)-plot_range(2,1))];
        for n=1:num_pts
            options=optimset('MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',1e-5);
            st_pts(n,:)=fminsearch(@(y) sum(abs(num_grad(@(x) trap_eval(btrap,x),[[y(1),0,y(2)]],delt)),2),x0(n,:),options);
            %PSoptions = optimoptions(@patternsearch,'Display','iter','TolFun',1e-10,'TolMesh',1e-10);
            %[st_pts(n,1:2),st_pts(n,3)] = patternsearch(@(y)  num_grad(@(x) trap_eval(btrap,x),y,delt),x0(n,:),[],[],[],[],plot_range(:,1),plot_range(:,2),PSoptions);
        end
        mask=plot_range(1,1)<st_pts(:,1) & st_pts(:,1)<plot_range(1,2) & plot_range(2,1)<st_pts(:,2) & st_pts(:,2)<plot_range(2,2);
        fprintf('removing %i points that are out of ROI \n',sum(mask))
        st_pts=st_pts(mask,:);
        st_pts_raw=st_pts;
        
        %remove points that are close to the previously found minima
        mask=false(1,size(st_pts,1));
        for n=1:size(st_pts,1)
            dist=sqrt(sum((st_pts(n,1:2)-min_pts(:,1:2)).^2,2));
            mask(n)=min(dist)<st_cluster_range;
        end
        fprintf('removing %i points that are close to minima \n',sum(mask))
        st_pts=st_pts(~mask,:);
        
        H_list=num_hessian(@(x) trap_eval(btrap,x),[st_pts(:,1),zeros(size(st_pts,1),1),st_pts(:,2)],delt);
        mask=false(1,size(st_pts,1));
        for i=1:size(H_list,3) %find the determinacy of each of these hessians
            eh=eig(H_list(:,:,i));
            mask(i)=sum(eh>1)>0 & sum(eh<0)>0; %check that contains +ve and -ve
        end
        clear('eh')
        fprintf('removing %i points  that do not pass hessian test \n',sum(~mask))
        st_pts=st_pts(mask,:);
        
        fprintf('surviving points %i out of %i \n',size(st_pts,1),num_pts)
        
        st_pt_clusters={};
        clust_idx=1;
        mask=false(1,size(st_pts,1));
        while size(st_pts,1)>0
            mask=sum((st_pts(1,:)-st_pts).^2,2)<st_cluster_range; %are the points within some threshold distance
            st_pt_clusters{clust_idx}=st_pts(mask,:);
            st_pts=st_pts(mask~=1,:);
            clust_idx=clust_idx+1;
        end
        st_pts=[];
        for n=1:size(st_pt_clusters,2) 
            st_pts(n,:)=mean(st_pt_clusters{n},1);
        end
          
        st_pts(:,3)=trap_eval(btrap,[st_pts(:,1),zeros(size(st_pts,1),1),st_pts(:,2)]); 
        st_pts(:,4)=sum(abs(num_grad(@(x)trap_eval(btrap,x),[st_pts(:,1),zeros(size(st_pts,1),1),st_pts(:,2)],delt)),2);
        fprintf('saddles clustered to %i points, %i zero crossing \n',size(st_pts,1),sum(st_pts(:,3)<zero_feild_thresh))


        for i=1: size(st_pts,1)
            fprintf('Trap saddle found %2.3fG (%2.3f MHz)(%2.3E K)\n',st_pts(i,3)*1e4,st_pts(i,3)*const.b_freq*1e-6,(2*st_pts(i,3)*const.mub)*2/const.kb)
        end
        
        figure(14)
        hold on
        scatter3((st_pts(:,1)-trap_cent(1))*1e3,(st_pts(:,2)-trap_cent(3))*1e3,st_pts(:,4),20,'k','filled')
        pause(1e-3)
        hold off
        
      
    end
    
end