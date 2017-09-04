%%%% A 3-d plotting routine for the MAGNIUTDE of the magnetic field

function[]=plotmag3(B,num1,num2) 

x=B(:,1);
y=B(:,2);
z=B(:,3);

Bmag= sqrt(sum(B(:,(4:6)).^2,2));
count3=0;
count1=0;
for i = 1:num2
   count2=0;
   count1=count1+1;
   for j = 1:num1
      count2=count2+1;
         count3=count3+1;
      BB(count1,count2)=Bmag(count3); 
      yc(count1,count2)=y(count3);
      zc(count1,count2)=z(count3);
   end
   end
   max(Bmag)
mesh(zc,yc,BB)
