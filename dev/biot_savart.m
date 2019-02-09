function B_out=biot_savart(l,dl,t,pts)
B_out = zeros(size(pts,1),3);
for ii = 1:size(pts,1)
    x = pts(ii,1);
    y = pts(ii,2);
    z = pts(ii,3);
    %simple direct biot-savart
    r_pt = repmat([x,y,z],length(t),1);
    
    r_dash = r_pt-l;
    r_dash_norm = vecnorm(r_dash')';
    B_out(ii,:) = trapz(t,cross(dl,r_dash)./(r_dash_norm.^(3)));
end
end