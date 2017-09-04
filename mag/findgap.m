%%%%%%% subroutine for finding the gap between each wire
function[gap]=findgap(n,w,stop,start) 

 %allow space for each wire 
 totwirethickness = n*w;
 width=stop-start;
 leftover=width-totwirethickness;
 gap=leftover/(2*n);


