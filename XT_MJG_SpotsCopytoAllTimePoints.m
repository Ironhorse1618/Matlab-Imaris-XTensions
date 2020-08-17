%
%
%  Spots Copy to all time points
%
%  Copyright Bitplane BPI 2014
%
%
%  Installation:
%
%  - Copy this file into the XTensions folder in the Imaris installation directory
%  - You will find this function in the Image Processing menu
%
%    <CustomTools>
%      <Menu>
%       <Submenu name="Spots Functions">
%        <Item name="Spots Copy to All Time Points" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_SpotsCopytoAllTimePoints(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSpots">
%          <Item name="Spots Copy to All Time Points" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_SpotsCopytoAllTimePoints(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
% 
%
%  Description:
%   
%  This XTension will copy all or selected spots from a single time frame to
%  all of the time points.
% 
%

function XT_MJG_SpotsCopytoAllTimePoints(aImarisApplicationID)

% connect to Imaris interface
if ~isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
  javaaddpath ImarisLib.jar
  vImarisLib = ImarisLib;
  if ischar(aImarisApplicationID)
    aImarisApplicationID = round(str2double(aImarisApplicationID));
  end
  vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
else
  vImarisApplication = aImarisApplicationID;
end


% the user has to create a scene with some spots
vSurpassScene = vImarisApplication.GetSurpassScene;
if isequal(vSurpassScene, [])
  msgbox('Please create some Spots in the Surpass scene!');
  return;
end

% get the spots
vSpots = vImarisApplication.GetFactory.ToSpots(vImarisApplication.GetSurpassSelection);

% search the surfaces if not previously selected
if ~vImarisApplication.GetFactory.IsSpots(vSpots)        
  for vChildIndex = 1:vSurpassScene.GetNumberOfChildren
    vDataItem = vSurpassScene.GetChild(vChildIndex - 1);
    if isequal(vSpots, [])
      if vImarisApplication.GetFactory.IsSpots(vDataItem)
        vSpots = vImarisApplication.GetFactory.ToSpots(vDataItem);
      end
    end
  end
  % did we find the spots?
  if isequal(vSpots, [])
    msgbox('Please create some spots!');
    return;
  end
end

% get all spots
Spots = vImarisApplication.GetFactory.ToSpots(vImarisApplication.GetSurpassSelection);
vSpotsRadius = Spots.GetRadii;
vSpotsPositionXYZ = Spots.GetPositionsXYZ;
vSpotsTime = Spots.GetIndicesT;
vNumberOfSpots = numel(vSpotsTime);

vSpotsName = char(vSpots.GetName);
vSpots.SetVisible(0);


%Generate timeIndex from surfaces
TimeIndexMax=vImarisApplication.GetDataSet.GetSizeT;

vSelectedSpots = vSpots.GetSelectedIds;
vNewSpots = vImarisApplication.GetFactory.CreateSpots;
vSpotsTimeFinal= [];
vSpotsPositionXYZFinal = [];
vSpotsRadiusFinal= [];
%Delete values in column less than 10000000
vSelectedSpots=vSelectedSpots(vSelectedSpots<10000000);
vProgressDisplay= waitbar(0, 'Duplicating spots...','Position',[600 500 285 60]);

if isempty(vSelectedSpots)
        for time=0:TimeIndexMax-1
            vTimeIndex= ones(vNumberOfSpots,1)*time;
            vSpotsTimeFinal =[vSpotsTimeFinal;vTimeIndex];
            vSpotsPositionXYZFinal= [vSpotsPositionXYZFinal;vSpotsPositionXYZ];
            vSpotsRadiusFinal=[vSpotsRadiusFinal;vSpotsRadius];
            waitbar((time+1)/TimeIndexMax,vProgressDisplay);
        end
        vNewSpots.Set(vSpotsPositionXYZFinal, vSpotsTimeFinal, vSpotsRadiusFinal);
else
    vNumberOfSpots = numel(vSelectedSpots);
    vSelectedSpots = int16(vSelectedSpots);
    vSelectedSpots = vSelectedSpots+1;
    vSpotsPositionXYZ=vSpotsPositionXYZ(vSelectedSpots,1:3)
    vSpotsRadius=vSpotsRadius(vSelectedSpots,1);

    for time=0:TimeIndexMax-1
        vTimeIndex= ones(numel(vSpotsRadius),1)*time;
        vSpotsTimeFinal =[vSpotsTimeFinal;vTimeIndex];
        vSpotsPositionXYZFinal= [vSpotsPositionXYZFinal;vSpotsPositionXYZ];
        vSpotsRadiusFinal=[vSpotsRadiusFinal;vSpotsRadius];    
        waitbar((time+1)/TimeIndexMax,vProgressDisplay);
    end 
    vNewSpots.Set(vSpotsPositionXYZFinal, vSpotsTimeFinal, vSpotsRadiusFinal);
end
vSpotsIndices=(0:(TimeIndexMax-1)*vNumberOfSpots-1)';
vTrackEdges=[vSpotsIndices,vSpotsIndices+vNumberOfSpots];
vNewSpots.SetTrackEdges(vTrackEdges);

vNewSpots.SetName([vSpotsName, ' copied to all timepoints']);
vRGBA = vSpots.GetColorRGBA;
vNewSpots.SetColorRGBA(vRGBA);
vSpots.GetParent.AddChild(vNewSpots, -1);
close(vProgressDisplay);




