%test equalities
% this script goes through and cross checks the various B feild calculators

% we will validate by samping in ±2 across all dimensions 
xyz_samp=(rand(10,3)-0.5)*4;
rot_vec=[0.1,0,0];
logic_str={'pass','!!! FAIL FAIL FAIL !!!'};
verbose=2;
fail_count=0;
frac_tolerance=1e-9;
%%
%% Bfield_finite_line_analytic  vs  Bfield_path_numeric
fprintf('Bfield_finite_line_analytic vs Bfield_path_numeric \n')
% create a path of a line
syms t
path_line=symfun([0,0,t],t);
len=1;
curr=1;
dlen=1e-6; %need lots of detail for numeric to match analytic well enough
tlim=[0,len];

tic
res_a=Bfield_finite_line_analytic(tlim,curr,rot_vec,xyz_samp);
time_a=toc;

tic
res_b=Bfield_path_numeric(path_line,tlim,dlen,curr,rot_vec,xyz_samp);
time_b=toc;

res_diff=res_a-res_b;
res_mean=mean(cat(3,res_a,res_b),3);
res_frac=res_diff./res_mean;
res_diff_max=max(res_frac(:));
res_diff_mean=mean(res_frac(:));
res_diff_norm=vecnorm(res_frac,2,2);

if verbose>1
    fprintf('INFO: eval time Bfield_finite_line_analytic  %.3g\n',time_a)
    fprintf('INFO: eval time Bfield_path_numeric          %.3g\n', time_b)
    fprintf('INFO: rel time                               %.3g\n', time_a/time_b)
end
fail_test=(abs(res_diff_max)>frac_tolerance);
if fail_test || verbose>0
    
    fprintf('TEST: %s \n',logic_str{1+(fail_test)})
end

if verbose>0 || fail_test
    fprintf('INFO: max frac diff %g\n',res_diff_max)
end

fail_count=fail_count+fail_test;
fprintf('\n\n')
%% Bfield_finite_line_numeric  vs  Bfield_path_numeric
fprintf('Bfield_finite_line_numeric vs Bfield_path_numeric \n')
% create a path of a line
syms t
path_line=symfun([0,0,t],t);
len=1;
curr=1;
dlen=1e-3;
tlim=[0,len];

tic
res_a=Bfield_finite_line_numeric(tlim,dlen,curr,rot_vec,xyz_samp);
time_a=toc;

tic
res_b=Bfield_path_numeric(path_line,tlim,dlen,curr,rot_vec,xyz_samp);
time_b=toc;

res_diff=res_a-res_b;
res_mean=mean(cat(3,res_a,res_b),3);
res_frac=res_diff./res_mean;
res_diff_max=max(res_frac(:));
res_diff_mean=mean(res_frac(:));
res_diff_norm=vecnorm(res_frac,2,2);

if verbose>1
    fprintf('INFO: eval time Bfield_finite_line_numeric   %.3g\n',time_a)
    fprintf('INFO: eval time Bfield_path_numeric          %.3g\n', time_b)
    fprintf('INFO: rel time                               %.3g\n', time_a/time_b)
end
fail_test=(abs(res_diff_max)>frac_tolerance);
if fail_test || verbose>0
    
    fprintf('TEST: %s \n',logic_str{1+(fail_test)})
end

if verbose>0 || fail_test
    fprintf('INFO: max frac diff %g\n',res_diff_max)
end

fail_count=fail_count+fail_test;
fprintf('\n\n')
%ok looks like the arb path method is the same but about 2-3x slower (asymptotically)..
%%
%% Bfield_loop_analytic  vs  Bfield_path_numeric
fprintf('Bfield_loop_analytic vs Bfield_path_numeric \n')
% create a path of a line
syms t
radius=1;
pitch=1;
path_line=symfun([radius*cos(t), radius*sin(t), 0],t);
len=2*pi;
curr=1;
dlen=1e-3;
tlim=[0,len];

tic
res_a=Bfield_loop_analytic(radius,curr,rot_vec,xyz_samp);
time_a=toc;
tic
res_b=Bfield_path_numeric(path_line,tlim,dlen,curr,rot_vec,xyz_samp);
time_b=toc;

res_diff=res_a-res_b;
res_mean=mean(cat(3,res_a,res_b),3);
res_frac=res_diff./res_mean;
res_diff_max=max(res_frac(:));
res_diff_mean=mean(res_frac(:));
res_diff_norm=vecnorm(res_frac,2,2);

if verbose>1
    fprintf('INFO: eval time Bfield_loop_analytic         %.3g\n',time_a)
    fprintf('INFO: eval time Bfield_path_numeric          %.3g\n', time_b)
    fprintf('INFO: rel time                               %.3g\n', time_a/time_b)
end
fail_test=(abs(res_diff_max)>frac_tolerance);
if fail_test || verbose>0
    
    fprintf('TEST: %s \n',logic_str{1+(fail_test)})
end

if verbose>0 || fail_test
    fprintf('INFO: max frac diff %g\n',res_diff_max)
end
fail_count=fail_count+fail_test;

fprintf('\n\n')
%% Bfield_helix_numeric  vs  Bfield_path_numeric
fprintf('Bfield_helix_numeric vs Bfield_path_numeric \n')
% create a path of a line
syms t
radius=1;
pitch=1;
path_line=symfun([radius*cos(t), radius*sin(t), pitch*t/(2*pi)],t);
len=2*pi;
curr=1;
dlen=1e-3;
tlim=[0,len];

tic
res_a=Bfield_helix_numeric(radius,pitch,tlim,dlen,curr,rot_vec,xyz_samp);
time_a=toc;
tic
res_b=Bfield_path_numeric(path_line,tlim,dlen,curr,rot_vec,xyz_samp);
time_b=toc;

res_diff=res_a-res_b;
res_mean=mean(cat(3,res_a,res_b),3);
res_frac=res_diff./res_mean;
res_diff_max=max(res_frac(:));
res_diff_mean=mean(res_frac(:));
res_diff_norm=vecnorm(res_frac,2,2);

if verbose>1
    fprintf('INFO: eval time Bfield_helix_numeric         %.3g\n',time_a)
    fprintf('INFO: eval time Bfield_path_numeric          %.3g\n', time_b)
    fprintf('INFO: rel time                               %.3g\n', time_a/time_b)
end
fail_test=(abs(res_diff_max)>frac_tolerance);
if fail_test || verbose>0
    
    fprintf('TEST: %s \n',logic_str{1+(fail_test)})
end

if verbose>0 || fail_test
    fprintf('INFO: max frac diff %g\n',res_diff_max)
end
fail_count=fail_count+fail_test;

fprintf('\n\n')
%% Bfield_finite_line_numeric  vs  Bfield_loop_analytic 
fprintf('Bfield_loop vs Bfield_finite_line_numeric\n')

% create a loop
tic
btrap=[];
loopa.type='loop';
loopa.param.radius=1;
loopa.param.current=1;
loopa.param.position=[0.0,0,0];
loopa.param.rot=pi/2*[0,1,0];
btrap.b_src=loopa;
res_a=trap_eval(btrap,xyz_samp);
time_a=toc;

%create a loop from finite lines
tic

tic
btrap=[];
btrap.b_src=loop_from_lines(1,50000,1,'ana');
res_b=trap_eval(btrap,xyz_samp);
time_b=toc;


res_diff=res_a-res_b;
res_mean=mean(cat(3,res_a,res_b),3);
res_frac=res_diff./res_mean;
res_diff_max=max(res_frac(:));
res_diff_mean=mean(res_frac(:));
res_diff_norm=vecnorm(res_frac,2,2);

if verbose>1
    fprintf('INFO: eval time Bfield_loop                    %.3g\n',time_a)
    fprintf('INFO: eval time Bfield_finite_line_numeric     %.3g\n', time_b)
    fprintf('INFO: rel time                                 %.3g\n', time_a/time_b)
end
fail_test=(abs(res_diff_max)>frac_tolerance);
if fail_test || verbose>0
    
    fprintf('TEST: %s \n',logic_str{1+(fail_test)})
end

if verbose>0 || fail_test
    fprintf('INFO: max frac diff %g\n',res_diff_max)
end
fail_count=fail_count+fail_test;
fprintf('\n\n')


if fail_count>0, error('something failed to pass the cross check'), end