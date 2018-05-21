function trap_freq = mag_profile_1d(btrap,B_cent,trap_cent)
    global const
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
        ylabel('Bfield')

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

end