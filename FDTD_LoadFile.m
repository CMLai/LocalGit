%===============================================================================
% function [tDat]=FDTD_LoadFile(tFile)
%	tFile: FileName saved by FDTD2DSolver.c
%===============================================================================
function [oDat,FDTDSetting]=FDTD_LoadFile(tFile)
%-------------------------------------------------------------------------------
FDTDSetting.dx=0;
FDTDSetting.sFreq=0;
FDTDSetting.spread=0;
FDTDSetting.nPos=0;
FDTDSetting.nsteps=0;
FDTDSetting.Esi=1;
%-------------------------------------------------------------------------------
fp=fopen(tFile,'r');tOuts=fgetl(fp);fclose(fp);
[pName,pValue]=strtok(tOuts);
%-------------------------------------------------------------------------------
if sum(isletter(pName))~=0
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	New File Version
	nsteps=0;
	fp=fopen(tFile,'r');
	while (nsteps==0)
		tOuts=fgetl(fp);
		[pName,pValue]=strtok(tOuts);
		switch (pName)
			case 'dx'
				FDTDSetting.dx=str2num(pValue);
			case 'sFerq'
				FDTDSetting.sFreq=str2num(pValue);
			case 'spread'
				FDTDSetting.spread=str2num(pValue);
			case 'nPos'
				FDTDSetting.nPos=str2num(pValue);
			case 'Esi'
				FDTDSetting.Esi=str2num(pValue);
			case 'nsteps'
				nsteps=str2num(pValue);
				FDTDSetting.nsteps=nsteps;
		end
	end
	oDat=transpose(fscanf(fp,'%e',[FDTDSetting.nPos+2 FDTDSetting.nsteps]));
	fclose(fp);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Old File Version
	oDat=load(tFile);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%===============================================================================
%function_end
