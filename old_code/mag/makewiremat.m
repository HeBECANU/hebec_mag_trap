%%%% 

function[wiremat]=makewiremat(rlayers,zlayers,numofwires,howtowrap) 

wiremat=ones(rlayers,zlayers);
leftover=(rlayers*zlayers)-numofwires
if leftover<0
   error1=('You have too many wires for the number of input layers')
   pause
end
if (numofwires)<((rlayers*zlayers)-rlayers)
   error1=('You have too few wires for the number of input layers')
   pause
end

   if numofwires~= (rlayers*zlayers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Wrapping options%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if howtowrap==0
   %%%%This option is for wrapping from the outer edge first
   % so for example if you have ten total wraps and you 
   % want to use three layers in the z and 4 layers in the r direction
   % then there will be a leftover - do you want this coil
   % to be as far away from the origin or as close to the origin
   % as possible? This option is for as close.
   wiremat(((rlayers-leftover+1):rlayers),zlayers)=0;
elseif howtowrap==1
   %this option is for furthest away from the origin
   wiremat(((rlayers-leftover+1):rlayers),1)=0;
end
end