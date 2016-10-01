%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function []=FDTD_11m()
%	Hard source: As refelction to the original source position, the behavior is strange.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=FDTD_11m()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KSteps	=200;
Ex			=zeros(KSteps,1);
Hy			=zeros(KSteps,1);
Kc			=floor(KSteps/2.0);
T0			=40.0;
spread	=12;
NSteps	=400;
tStep		=[	1 10 20 30 40 50 100 120 140 160 180 ...
				200 210 220 230 240 260 280 300 320 340 360 380 400];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for T=1:NSteps
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%% Pulse source %%%%%%%%%%%%%%%%
	pulse			=exp(-0.5*((T0-T)/spread)^2);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Ex(2:end)	=Ex(2:end)+0.5*(Hy(1:end-1)-Hy(2:end));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Ex(Kc)		=pulse;				% Hard Source
%	Ex(Kc)		=Ex(Kc)+pulse;		% Soft Source
%	Ex(Kc+50)	=Ex(Kc+50)+pulse;	% Two Soft Source
%	Ex(Kc-50)	=Ex(Kc-50)+pulse;	% Two Soft Source
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Hy(1:end-1)	=Hy(1:end-1)+0.5*(Ex(1:end-1)-Ex(2:end));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%% Ploting figure %%%%%%%%%%%%%%
	if length(find(T==tStep))~=0
		figure(1);
		subplot(2,1,1);plot(transpose(1:KSteps),Ex);title(num2str(T));axis([0 KSteps -1 1]);
		subplot(2,1,2);plot(transpose(1:KSteps),Hy);title(num2str(T));axis([0 KSteps -1 1]);
		pause
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
