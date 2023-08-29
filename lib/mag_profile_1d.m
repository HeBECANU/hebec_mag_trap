function anal_out = mag_profile_1d(anal_out,btrap,B_cent)
    global const
    points=1000;
    range=[[1,-1];[1,-1];[1,-1]]*0.5*1e-4;
    labels=['x','y','z'];
    figure(1)
    set(gcf,'Color',[1 1 1]);
    clf;
    for n=1:3
        xyz_points=zeros(points,3);
        xyz_points(:,n)=linspace(range(n,1),range(n,2),points)';
        xyz_points=anal_out.trap_cen.pos+xyz_points;
        [bmag,bvec]=trap_eval(btrap,xyz_points);   
        bmag=bmag-anal_out.trap_cen.b_mag;

        subplot(3,1,n)
        deltx=xyz_points(:,n)-anal_out.trap_cen.pos(n);
        plot(deltx,bmag) %
        xlabel(labels(n))
        ylabel('Bfield')
        
        poly=polyfit(deltx,bmag,6);
        hold on
        plot(deltx,polyval(poly,deltx),'r')
 %       hold off
        %hold on
        dudxn(n,1)=polyval(polyder(poly),0);
        dudxn(n,2)=polyval(polyder(polyder(poly)),0);
        dudxn(n,3)=polyval(polyder(polyder(polyder(poly))),0);
        dudxn(n,4)=polyval(polyder(polyder(polyder(polyder(poly)))),0);
        dudxn(n,5)=polyval(polyder(polyder(polyder(polyder(polyder(poly))))),0);
        dudxn(n,6)=polyval(polyder(polyder(polyder(polyder(polyder(polyder(poly)))))),0);
    end
    fprintf('trap grad {%f , %f, %f} G/cm \n',dudxn(1,1)*1e2,dudxn(2,1)*1e2,dudxn(3,1)*1e2)
    fprintf('trap curvature {%f , %f, %f} G/cm^2 \n',dudxn(1,2),dudxn(2,2),dudxn(3,2))

    trap_freq=sqrt(2*const.mub*dudxn(:,2)'/const.mhe)/(2*pi);
    %calculate the derivative of the trap freq with position
    anal_out.trap_freq_anh1=(trap_freq.*dudxn(:,3)')./(2*dudxn(:,2)');
    anal_out.trap_freq_anh2=trap_freq.*((dudxn(:,4)'./(2*dudxn(:,2)'))-(dudxn(:,3)'.^2/(4*dudxn(:,2)'.^2)));
    %rescale to be derivative of trap freq with energy (only makes sense
    %for the second derivative)
    anal_out.trap_freq_anh2_u=anal_out.trap_freq_anh2./(2*const.mub.*dudxn(:,2)');
    %given a velocity in the trap what it the anharmonic shift
    v_osc=0.005;
    x_osc=v_osc./trap_freq;
    anal_out.trap_freq_shift=anal_out.trap_freq_anh1.*x_osc+anal_out.trap_freq_anh2.*(x_osc.^2);
    
    fprintf('trap freq {%f , %f, %f} Hz\n',trap_freq(1),trap_freq(2),trap_freq(3))
    fprintf('trap ratio {%f , %f}={y/x,z/x} \n',trap_freq(2)/trap_freq(1),trap_freq(3)/trap_freq(1))

    anal_out.trap_freq=trap_freq;
    

end