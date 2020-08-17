 %Resize all Spots

%Written by Matthew J. Gastinger, Bitplane Advanced Application Scientist.  
% 2014.
%
%<CustomTools>
%      <Menu>
%       <Submenu name="Spots Functions">
%        <Item name="Resize All Spots" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_SpotsResized(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSpots">
%          <Item name="Resize all Spots" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_SpotsResized(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description:
%This XTension will resize and copy all or selected spots to a new spots
%object, to a size specified by the user.

function XT_MJG_SpotsResized(aImarisApplicationID)
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
vSpotsRadius = vSpots.GetRadii;
vSpotsPositionXYZ = vSpots.GetPositionsXYZ;
vSpotsTime = vSpots.GetIndicesT;
vNumberOfSpots = numel(vSpotsTime);
vSelectedSpots = vSpots.GetSelectedIds;
vNewSpots = vImarisApplication.GetFactory.CreateSpots; 
%Delete values in column less than 10000000
vSelectedSpots=vSelectedSpots(vSelectedSpots<10000000);

%for vId = 0:4
  vSpots.SetLabel(vImarisApplication.GetFactory.CreateObjectLabel(237, 'Group A', 'Label A'));
%end



%Dialog to select new Spot size
vCurrentSpotSize = median(vSpotsRadius);
vCurrentSpotSize = num2str(vCurrentSpotSize);
vQuestion = {sprintf(['Please enter new Spot Radius:'])};
    vAnswer = inputdlg(vQuestion,'Resize Spots',1,{vCurrentSpotSize});
    if isempty(vAnswer), return, end
newRadius = str2double(vAnswer{1}); 

if isempty(vSelectedSpots)
    vSpotsRadius = ones(vNumberOfSpots,1)*newRadius;
    vNewSpots.Set(vSpotsPositionXYZ, vSpotsTime, vSpotsRadius);
    vNewSpots.SetTrackEdges(vSpots.GetTrackEdges);
else
    vNumberOfSpots = numel(vSelectedSpots);
    for s=1:numel(vSelectedSpots)
        vCurrentSpot = vSelectedSpots(s);
        vCurrentSpot = int16(vCurrentSpot);
        vCurrentSpot = vCurrentSpot + 1;
        %adjust selected spot radius
        vSpotsRadius (vCurrentSpot) =  newRadius;
    end
    vNewSpots.Set(vSpotsPositionXYZ, vSpotsTime, vSpotsRadius);
    vNewSpots.SetTrackEdges(vSpots.GetTrackEdges);
end

%Create the a new Spots with new size including tracks if they are there
vRGBA = vSpots.GetColorRGBA;
vNewSpots.SetColorRGBA(vRGBA);
vNewSpots.SetName([char(vSpots.GetName), ' Resized']);
 
vSpots.SetVisible(0);
vSpots.GetParent.AddChild(vNewSpots, -1);

end