function btrap=build_nuller(btrap,nullr_opt)

rad_inset=nullr_opt.radius-nullr_opt.inset;
%posx
loop=[];
loop.position=[rad_inset,0,0];
loop.current=nullr_opt.current(1,1);
loop.rot=pi/2*[0,1,0];
loop.radius=nullr_opt.radius;
btrap=square_loop(btrap,loop);
%negx
loop.current=nullr_opt.current(1,2);
loop.position=[-rad_inset,0,0];
btrap=square_loop(btrap,loop);
%posy
loop.rot=pi/2*[1,0,0];
loop.current=nullr_opt.current(2,1);
loop.position=[0,rad_inset,0];
btrap=square_loop(btrap,loop);
%negy
loop.current=nullr_opt.current(2,2);
loop.position=[0,-rad_inset,0];
btrap=square_loop(btrap,loop);
%posz
loop.current=nullr_opt.current(3,1);
loop.rot=pi/2*[0,0,0];
loop.position=[0,0,rad_inset];
btrap=square_loop(btrap,loop);
%negz
loop.current=nullr_opt.current(3,2);
loop.position=[0,0,-rad_inset];
btrap=square_loop(btrap,loop);



end