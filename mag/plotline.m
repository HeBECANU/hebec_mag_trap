%%%% A 3-d plotting routine for the MAGNITUDE of the magnetic field

function[Bmag]=plotline(B,num,coord) 

%%%%Remeber to change back to G and cm!!!!
x=B(:,1);
y=B(:,2);
z=B(:,3);

Bmag= sqrt(sum(B(:,(4:6)).^2,2));
Bmag=Bmag/1e-4;
if coord == 1
   plot(x*100,Bmag,'g-.')
   xlabel('cm');
   ylabel('B (Gauss)');
   s=(gca);
set(s,'Fontsize',14)

   elseif coord == 2
      plot(y*100,Bmag,'b-')
      s=(gca);
set(s,'Fontsize',14)

elseif coord == 3
   plot(z*100,Bmag,'r-.')
   s=(gca);
set(s,'Fontsize',14)

elseif coord == 4
   for i = 1:length(x)
      if x(i)<0
         xy(i)=-1*sqrt(x(i)^2+y(i)^2);
      else
         xy(i)=1*sqrt(x(i)^2+y(i)^2);
      end
   end
   
   plot(xy*100,Bmag,'g:')
   s=(gca);
set(s,'Fontsize',14)

   elseif coord == 5
   for i = 1:length(z)
      if z(i)<0
         zx(i)=-1*sqrt(z(i)^2+x(i)^2);
      else
         zx(i)=1*sqrt(z(i)^2+x(i)^2);
      end
   end
   
   plot(zx*100,Bmag,'c')
s=(gca);
set(s,'Fontsize',14)

end

    xlabel('cm','fontsize',14);
   ylabel('B (Gauss)','fontsize',14);

