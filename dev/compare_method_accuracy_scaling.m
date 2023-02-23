%compare line methods numerical error convergence

rot_vec=[0.0,0,0];
len=1;
tlim=[0,len];
dlen_vec=logspace(-1,-6,20);
path_line=symfun([0,0,t],t);

xyz_list=[1,1,1]*-1e-3;

iimax=numel(dlen_vec);
b_analytic=Bfield_finite_line_analytic(tlim,curr,rot_vec,xyz_list);
line_compare=[];
line_compare.zeroth.error=nan(iimax,1);
line_compare.zeroth.time=nan(iimax,1);
line_compare.first.error=nan(iimax,1);
line_compare.first.time=nan(iimax,1);


fprintf('%02u',0)
for ii=1:iimax
    dlen=dlen_vec(ii);
    timer_obj=tic;
    b_zeroth_order=Bfield_finite_line_numeric_order_0(tlim,dlen,curr,rot_vec,xyz_list);
    line_compare.zeroth.time(ii)=toc(timer_obj);
    line_compare.zeroth.error(ii)=norm(b_zeroth_order-b_analytic)/norm(b_analytic);
    
    timer_obj=tic;
    b_first_order=Bfield_finite_line_numeric_order_1(tlim,dlen,curr,rot_vec,xyz_list);
    line_compare.first.time(ii)=toc(timer_obj);
    line_compare.first.error(ii)=norm(b_first_order-b_analytic)/norm(b_analytic);

    fprintf('\b\b%02u',ii)
end
fprintf('\n')




%%
stfig('line method error');
clf
subplot(1,2,1)
plot(dlen_vec,line_compare.zeroth.error)
hold on
plot(dlen_vec,line_compare.first.error)
hold off
legend('zeroth','first')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'Xdir', 'reverse')
xlabel('current segment length(m)')
ylabel('error')



subplot(1,2,2)
plot(dlen_vec,line_compare.zeroth.time)
hold on
plot(dlen_vec,line_compare.first.time)
hold off
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'Xdir', 'reverse')
xlabel('current segment length(m)')
ylabel('eval time(s)')



%% Compare the convergence for a loop

dlen_vec=logspace(-1,-6,20);

syms t
radius=1;
pitch=1;
path_line=symfun([radius*cos(t), radius*sin(t), 0],t);
len=2*pi;
curr=1;
dlen=1e-3;
tlim=[0,len];

%xyz_list=[1,0,0]*radius*1.001;
xyz_list=[1,0,1]*radius*1;

b_analytic=Bfield_loop_analytic(radius,curr,rot_vec,xyz_list);

iimax=numel(dlen_vec);
loop_compare=[];
loop_compare.zeroth.error=nan(iimax,1);
loop_compare.zeroth.time=nan(iimax,1);
loop_compare.first.error=nan(iimax,1);
loop_compare.first.time=nan(iimax,1);
fprintf('%02u',0)
for ii=1:iimax
    dlen=dlen_vec(ii);
    timer_obj=tic;
    b_zeroth_order=Bfield_path_numeric_order_0(path_line,tlim,dlen,curr,rot_vec,xyz_list);
    loop_compare.zeroth.time(ii)=toc(timer_obj);
    loop_compare.zeroth.error(ii)=norm(b_zeroth_order-b_analytic)/norm(b_analytic);
    
    timer_obj=tic;
    b_first_order=Bfield_path_numeric_order_1(path_line,tlim,dlen,curr,rot_vec,xyz_list);
    loop_compare.first.time(ii)=toc(timer_obj);
    loop_compare.first.error(ii)=norm(b_first_order-b_analytic)/norm(b_analytic);
    fprintf('\b\b%02u',ii)
end
fprintf('\n')




%
stfig('loop method error');
clf
subplot(1,2,1)
plot(dlen_vec,loop_compare.zeroth.error)
hold on
plot(dlen_vec,loop_compare.first.error)
hold off
legend('zeroth','first')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'Xdir', 'reverse')
xlabel('current segment length(m)')
ylabel('error')

subplot(1,2,2)
plot(dlen_vec,loop_compare.zeroth.time)
hold on
plot(dlen_vec,loop_compare.first.time)
hold off
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'Xdir', 'reverse')
xlabel('current segment length(m)')
ylabel('eval time(s)')

