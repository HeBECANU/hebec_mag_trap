function btrap=square_loop(btrap,sq_opt)
%build a square loop out of line primitives
% sq_loop_opt.position=[0,30e-3,0];
% sq_loop_opt.current=1;
% sq_loop_opt.rot=pi/2*[1,0,0];
% sq_loop_opt.radius=20e-3;

if ~isfield(btrap,'b_src')
    btrap.b_src=[];
end

rot_mat_pos=rotationVectorToMatrix(sq_opt.rot);
rot_mat_neg=rotationVectorToMatrix(-sq_opt.rot);

top_wire=[];
top_wire.type='line';
top_wire.param.current=sq_opt.current;

top_wire.param.position=sq_opt.position+[1,1,0]*sq_opt.radius*rot_mat_pos;
top_wire.param.rot=rotationMatrixToVector(rot_mat_neg*rotationVectorToMatrix(pi/2*[0,1,0]));
top_wire.param.length=sq_opt.radius*2;

btm_wire=top_wire;
btm_wire.param.position=sq_opt.position+[-1,-1,0]*sq_opt.radius*rot_mat_pos;
btm_wire.param.rot=rotationMatrixToVector(rot_mat_neg*rotationVectorToMatrix(pi/2*[0,-1,0]));

r_wire=top_wire;
r_wire.param.position=sq_opt.position+[1,-1,0]*sq_opt.radius*rot_mat_pos;
r_wire.param.rot=rotationMatrixToVector(rot_mat_neg*rotationVectorToMatrix(pi/2*[1,0,0]));

l_wire=top_wire;
l_wire.param.position=sq_opt.position+[-1,1,0]*sq_opt.radius*rot_mat_pos;
l_wire.param.rot=rotationMatrixToVector(rot_mat_neg*rotationVectorToMatrix(pi/2*[-1,0,0]));

btrap.b_src=[btrap.b_src,top_wire,btm_wire,r_wire,l_wire];

end