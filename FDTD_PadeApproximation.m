%===============================================================================
% function [ny]=FDTD_PadeApproximation(x,y,nx,Na,Nb)
%	Pade Approximation
%	Input:
%		x,y: The Original Data
%		nx: fine resolution of xn
%		Na: The highest order of the numerator (Best value= 2* Resonance Number)
%		Nb: The highest order of the denominator (Best value= 2* Resonance Number)
%	Output:
%		ny: the value @ nx 
%	Note: the length(x)>= Na+Nb+2
%	From Mittra, Ref: MGWL 8,415,1998
%	Edited by C. M. Lai, 2006.06.29
%===============================================================================
function [ny]=FDTD_PadeApproximation(x,y,nx,Na,Nb)
%===============================================================================
x			=x(:);
y			=y(:);

[sx,snx]	=xScaling(x,nx);
sLen		=length(sx);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tM1=-sx(:,ones(1,Na+1)).^(repmat(0:Na,sLen,1));							% -wi^j, j=0:tN
tM2= sx(:,ones(1,Nb  )).^(repmat(1:Nb,sLen,1)).*y(:,ones(1,Nb));	% P(wi)*wi^j, j=1:tN
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pCoef	=-[tM2 tM1]\y;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pBeta	=[pCoef(Nb:-1:1);1];
pAlpha=pCoef(end:-1:Nb+1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ny		=polyval(pAlpha,snx)./polyval(pBeta,snx);
%===============================================================================
%function_end
%===============================================================================
% function [xs]=xScaling(x)
%	Scaling the wi to the [-1 1] range
%===============================================================================
function [sx,snx]=xScaling(x,nx)
%-------------------------------------------------------------------------------
xMax	=max(x);
xMin	=min(x);
sx		=(2*x -(xMax+xMin))/(xMax-xMin);
snx	=(2*nx-(xMax+xMin))/(xMax-xMin);
%===============================================================================
%function_end
