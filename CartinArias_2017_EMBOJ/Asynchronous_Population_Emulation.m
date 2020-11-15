% to emulate RedN coverage in asynchronous population.
% distribution of relative cell cycle progression is based on Rodrigo's demograph of DnaN 

Ntot=8000; % total cells in the demograph
 N03=200;  % cut off of "too long" cells - wil be ignored, not essential;
 N01=6700; % approximate "average" appearance of DnaN focus in the demograph
 N02=500; % approximate "average" disappearance of DnaN focus in the demograph
 
% then assuming that before and after DnaN focus RedN distrbution is one simple triangle... 
% and during is triangle+trapezoid with zero at relative (genomically) DnaN
% position (based on our modelling, see the paper), converge according to demograph distribution.

% remove discarded cells
Ntot=Ntot-N03;
 N01=N01-N03;
 N02=N02-N03;
 
% now converge
nn=100; %bins along chromosome position 
 xx=0:1/nn:1; % chromosome coordinates
 dx=xx(2)-xx(1);

% we will get combined distribution by accumulatiung distributions from the 3 cases:

% case 1 before replication
% y=-2(x-1)
YY=-2*(xx-1)*(Ntot-N01);

% case 3 after replication
% y=-2(x-1)
YY=YY+(-2*(xx-1)*N02);

% case 2 during replication
% y=-2(x-1)+2*a for x<=x0, x0 - replisome position
% y=-2(a-x)     for x>x0
for cc=1:(N01-N02)
  x0=(cc-1)/(N01-N02-1); 
  YY=YY+(-2*(xx-1)+2*x0).*heaviside(xx-(x0-dx/2))+(2*(x0-xx)).*heaviside((x0-dx/2)-xx);
end

% now distribution is ready




