%===============================================================================
% function [FDTDSetting]=FDTD_Setting()
%===============================================================================
function [FDTDSetting]=FDTD_JAP99_093112()
%-------------------------------------------------------------------------------
LightSpeed			=3.0e8;
kDispersionLimit	=10.0;
tSatbilityLimit	=2.0;			% Can NOT be Changed at this time. Need Modified Codes.
%-------------------------------------------------------------------------------
FDTDSetting.RootPath		=[FolderSetting('FDTDFolder')];
FDTDSetting.ProjectName	=['JAP99_093112_01'];
FDTDSetting.Notes			=[''];
%-------------------------------------------------------------------------------
FDTDSetting.nSteps		=[2^16];
FDTDSetting.Dimension	=[2];
FDTDSetting.dSize			=[10 10 10]*1e-9;	% in m unit
FDTDSetting.PMLNum		=[10];
FDTDSetting.PulseTime	=[1 1 800e12 1000];	% [Tyep Amp Freq Broaden ....] depend on Type
FDTDSetting.PulseSpace	=[ 31  43  0 0 1];	% Source Index [i j k]
FDTDSetting.zStore		=[ 43  31  0;
								   37  37  0;		% Store electric field Index [i j k]
								   35  43  0;
								   49  47  0];	
FDTDSetting.IsSoftSource=[0];
FDTDSetting.AllSpace		=[1 100;1 100;1 100];	% in Index unit
FDTDSetting.DefaultEsi	=[1.0];
FDTDSetting.dTime			=[FDTDSetting.dSize/LightSpeed/tSatbilityLimit];
FDTDSetting.MaxFreq		=LightSpeed./(FDTDSetting.dSize*kDispersionLimit);
%-------------------------------------------------------------------------------
Polygon(1).Esi				=[2.3^2];					% For ZnO
Polygon(1).Types			=['Hexagonal'];			% Triangular, rectangular, Circle, Polygon
Polygon(1).Vertices		=[50 50 0 380e-9];
% Vertices:		[Center Position(3) length(m)] for Triangular  
% Vertices:		[Center Position(3) xlength(m) ylength(m)] for Rectangular  
% Vertices:		[Center Position(3) Radius(m)] for Circle  
% Vertices:		[Position(3);Position(3);...] for Polygon. Need Closed  
%-------------------------------------------------------------------------------
FDTDSetting.Polygon		=Polygon;
%-------------------------------------------------------------------------------
FDTDSetting.PlotFFT		=[1];
FDTDSetting.PlotFFTLim	=[]*1e12;
FDTDSetting.PlotTime		=[0];
FDTDSetting.PlotEsi		=[1];
FDTDSetting.PlotPulse	=[1];
FDTDSetting.FinalField	=[1];
%-------------------------------------------------------------------------------
FDTDSetting.PadeApprox	=[];		% 0: off, 1: Method 1, 2: Method 2,...
FDTDSetting.PadeFLimit	=[];		% [Fmin Fmax divisions];
%-------------------------------------------------------------------------------
FDTDSetting.PadeApprox	=[0 16];		% [0]: off, [1 refine]: Method 1
FDTDSetting.PadeFLimit	=[]*1e12;		% [Fmin Fmax divisions];
%===============================================================================
%function_end
