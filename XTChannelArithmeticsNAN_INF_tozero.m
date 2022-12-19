%
%
%  Channel Arithmetics Function for Imaris 8.2
%
%  Copyright Bitplane 2015
%
%
%  Installation:
%
%  - Copy this file into the XTensions folder in the Imaris installation directory
%  - You will find this function in the Image Processing menu
%
%    <CustomTools>
%      <Menu>
%        <Item name="Channel Arithmetics - NEW ratio" icon="Matlab" tooltip="Create a new channel combining the others using a regular expression.">
%          <Command>MatlabXT::XTChannelArithmeticsNAN_INF_tozero(%i)</Command>
%        </Item>
%      </Menu>
%    </CustomTools>
% 
%
%  Description:
%   
%   Create a new channel combining the others using a regular expression. 
% 

function XTChannelArithmeticsNAN_INF_tozero(aImarisApplicationID)

% get the application object
if isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
  % called from workspace
  vImarisApplication = aImarisApplicationID;
else
  % connect to Imaris interface
  javaaddpath ImarisLib.jar
  vImarisLib = ImarisLib;
  if ischar(aImarisApplicationID)
    aImarisApplicationID = round(str2double(aImarisApplicationID));
  end
  vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
end

vAnswer = inputdlg({sprintf(['Combination expression:\n\n', ...
  'Channel names: ch1, ch2, ...\nUse matlab operators, i.e. ', ...
  '+, -, .*, ./, .^, sqrt, ...\n'])}, ...
    'Channel Arithmetics', 1, {'sqrt(ch1 .* ch2)'});
if isempty(vAnswer), return, end

vDataSet = vImarisApplication.GetDataSet.Clone;

vLastC = vDataSet.GetSizeC;
vDataSet.SetSizeC(vLastC + 1);
vMin = vDataSet.GetChannelRangeMin(0);
vMax = vDataSet.GetChannelRangeMax(0);
vDataSet.SetChannelRange(vLastC, vMin, vMax);

vProgressDisplay = waitbar(0, 'Channel Arithmetics');
vProgressCount = 0;

vDataSize = [vDataSet.GetSizeX, vDataSet.GetSizeY, vDataSet.GetSizeZ];
if vDataSize(3) == 1
  vBlockSize = [1024, 1024, 1];
else
  vBlockSize = [512, 512, 32];
end
% vBlockSize = [100000, 100000, 1]; % process slice by slice (better not)
vBlockCount = ceil(vDataSize ./ vBlockSize);

vProgressTotalCount = vDataSet.GetSizeT*prod(vBlockCount);

for vTime = 1:vDataSet.GetSizeT
    for vIndexZ = 1:vBlockCount(3)
      vMinZ = (vIndexZ - 1) * vBlockSize(3);
      vSizeZ = min(vBlockSize(3), vDataSize(3) - vMinZ);
    for vIndexY = 1:vBlockCount(2)
      vMinY = (vIndexY - 1) * vBlockSize(2);
      vSizeY = min(vBlockSize(2), vDataSize(2) - vMinY);
    for vIndexX = 1:vBlockCount(1)
      vMinX = (vIndexX - 1) * vBlockSize(1);
      vSizeX = min(vBlockSize(1), vDataSize(1) - vMinX);

        for vChannel = 1:vLastC
          if strcmp(vDataSet.GetType,'eTypeUInt8')
              vData = vDataSet.GetDataSubVolumeAs1DArrayBytes(...
                vMinX, vMinY, vMinZ, vChannel-1, vTime-1, vSizeX, vSizeY, vSizeZ);
              vData = typecast(vData, 'uint8');
          elseif strcmp(vDataSet.GetType,'eTypeUInt16')
              vData = vDataSet.GetDataSubVolumeAs1DArrayShorts(...
                vMinX, vMinY, vMinZ, vChannel-1, vTime-1, vSizeX, vSizeY, vSizeZ);
              vData = typecast(vData, 'uint16');
          elseif strcmp(vDataSet.GetType,'eTypeFloat')
              vData = vDataSet.GetDataSubVolumeAs1DArrayFloats(...
                vMinX, vMinY, vMinZ, vChannel-1, vTime-1, vSizeX, vSizeY, vSizeZ);
          end
          % works on double to allow division and prevent overflows
          eval(sprintf('ch%i = double(vData);', vChannel));
        end

        try
          vData = eval(vAnswer{1});
		  vData(isnan(vData))=0;
  		  vData(isinf(vData))=0;

        catch er
          close(vProgressDisplay);
          msgbox(sprintf(['Error while evaluating the expression.\n\n', ...
            'Possible causes: invalid variable names (ch1, ch2, ...), ', ...
            'invalid operators (use .* instead of *)...\n\n', er.message]));
          return;
        end
        
        try
          if strcmp(vDataSet.GetType,'eTypeUInt8')
              vDataSet.SetDataSubVolumeAs1DArrayBytes(uint8(vData), ...
                vMinX, vMinY, vMinZ, vLastC, vTime-1, vSizeX, vSizeY, vSizeZ);
          elseif strcmp(vDataSet.GetType,'eTypeUInt16')
              vDataSet.SetDataSubVolumeAs1DArrayShorts(uint16(vData), ...
                vMinX, vMinY, vMinZ, vLastC, vTime-1, vSizeX, vSizeY, vSizeZ);
          elseif strcmp(vDataSet.GetType,'eTypeFloat')
              vDataSet.SetDataSubVolumeAs1DArrayFloats(single(vData), ...
                vMinX, vMinY, vMinZ, vLastC, vTime-1, vSizeX, vSizeY, vSizeZ);
          end
        catch er
          close(vProgressDisplay);
          msgbox(sprintf(['The result of the expression is not a valid dataset.\n\n', ...
            'Possible causes: invalid result size.\n\n', er.message]));
          return
        end

        vProgressCount = vProgressCount + 1;
        waitbar(vProgressCount/vProgressTotalCount, vProgressDisplay);

    end
    end
    end
end

vImarisApplication.SetDataSet(vDataSet);
close(vProgressDisplay);

