%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=FDTD_12m()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KSteps	=200;
Ex			=zeros(KSteps,1);
Hy			=zeros(KSteps,1);
Kc			=floor(KSteps/2.0);
T0			=40.0;
spread	=12;
NSteps	=400;
ExBH		=[0 0];
ExBL		=[0 0];
Esi		=[ones(Kc,1); 4*ones(KSteps-Kc,1)];
tStep		=[	1 10 20 30 40 50 100 120 140 160 180 ...
				200 210 220 230 240 260 280 300 320 340 360 380 400];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for T=1:NSteps
	Ex(2:end)	=Ex(2:end)+0.5*(Hy(1:end-1)-Hy(2:end));
	pulse			=exp(-0.5*((T0-T)/spread)^2);
	Ex(Kc)		=Ex(Kc)+pulse;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Ex(1)		=ExBL(1);
	ExBL(1)	=ExBL(2);
	ExBL(2)	=Ex(2);

	Ex(end)	=ExBH(1);
	ExBH(1)	=ExBH(2);
	ExBH(2)	=Ex(end-1);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Hy(1:end-1)	=Hy(1:end-1)+0.5*(Ex(1:end-1)-Ex(2:end));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if length(find(T==tStep))~=0
		figure(1);
		subplot(2,1,1);plot(transpose(1:KSteps),Ex);title(num2str(T));axis([0 KSteps -1.5 1.5]);
		subplot(2,1,2);plot(transpose(1:KSteps),Hy);title(num2str(T));axis([0 KSteps -1.5 1.5]);
		pause
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
