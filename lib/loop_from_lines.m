function lines=loop_from_lines(radius,segments,current)
dtheta=2*pi/segments;
dl=2*tan(dtheta/2)*radius;

%can adjust the circle that the line segments are inscribed in to get better convergence
%rad_correction=radius*(cos(dtheta)-1)*(7/40);
%better yet we can set the area to be equal to the circle with radius r
r_poly=sqrt(2*pi*radius^2/(segments*sin(2*pi/segments)));
radius=r_poly;


line.type='line';
line.param.length=dl;
line.param.current=current;
lines=[line,line]; %initalize the structure array
lines(segments)=line;
for ii=1:segments %do it backwards for implicit structure array prealocation
    theta=ii*dtheta;
    line.param.position=[0,sin(theta-dl/2)*radius,cos(theta-dl/2)*radius];
    line.param.rot=[theta+pi/2,0,0];
    lines(ii)=line;
end


end