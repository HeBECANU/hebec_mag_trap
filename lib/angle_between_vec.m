function angle=angle_between_vec(u,b)
angle=atan2(norm(cross(u,v)),dot(u,v));
end