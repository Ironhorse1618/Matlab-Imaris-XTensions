 %Flatten 3D Spots

%Written by Matthew J. Gastinger, Bitplane Advanced Application Scientist.  
%March 2014.
%
%<CustomTools>
%      <Menu>
%       <Submenu name="Spots Functions">
%        <Item name="Spots Flatten" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_Spots3DFlattened(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSpots">
%          <Item name="Spots Flatten" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_Spots3DFlattened(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description:
%This XTension will flatten all of the spots to a single plane

function XT_MJG_Spots3DFlattened(aImarisApplicationID)
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
    
% get all spots
vSpots = vImarisApplication.GetFactory.ToSpots(vImarisApplication.GetSurpassSelection);
vSpotsRadius = vSpots.GetRadii
vSpotsPositionXYZ = vSpots.GetPositionsXYZ;
vSpotsTime = vSpots.GetIndicesT;
vNumberOfSpots = numel(vSpotsTime);
vSelectedSpots = vSpots.GetSelectedIds
vNewSpots = vImarisApplication.GetFactory.CreateSpots; 

if isempty(vSelectedSpots)
    vSpotsPositionXYZ(:,3)=0;
    vNewSpots.Set(vSpotsPositionXYZ, vSpotsTime, vSpotsRadius);
    vNewSpots.SetTrackEdges(vSpots.GetTrackEdges);
    
else
    vNumberOfSpots = numel(vSelectedSpots);
    vSpotsRadiusSelected=[];
    vSpotsRadiustemp=[];
    vSpotsTimeTemp=[];
    vSpotsTimeSelected=[];
    vSpotsPositionXYZSelected = [];
    vTrackedgesTest = vSpots.GetTrackEdges;
    for s=1:numel(vSelectedSpots)
        vCurrentSpot = vSelectedSpots(s);
        vCurrentSpot = int8(vCurrentSpot);
        vCurrentSpot = vCurrentSpot + 1;
        vSpotsPositionXYZtemp = vSpotsPositionXYZ(vCurrentSpot,:);
        vSpotsTimeTemp = vSpotsTime(vCurrentSpot,:);
        vSpotsPositionXYZSelected = [vSpotsPositionXYZSelected;vSpotsPositionXYZtemp];
        vSpotsTimeSelected = [vSpotsTimeSelected;vSpotsTimeTemp];
        vSpotsRadiusTemp=vSpotsRadius(vCurrentSpot,:);
        vSpotsRadiusSelected=[vSpotsRadiusSelected;vSpotsRadiusTemp];
    end
    vSpotsPositionXYZSelected(:,3)=0;
    vSpotsTime = vSpotsTimeSelected;
    vNewSpots.Set(vSpotsPositionXYZSelected, vSpotsTimeSelected, vSpotsRadius);
end

%Create the a new Spots including tracks if they are there
vRGBA = vSpots.GetColorRGBA;
vNewSpots.SetColorRGBA(vRGBA);
vNewSpots.SetName([char(vSpots.GetName), ' Flattened']);
 
vSpots.SetVisible(0);
vSpots.GetParent.AddChild(vNewSpots, -1);

end