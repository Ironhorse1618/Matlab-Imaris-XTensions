%Convert Surface center of homogeneous mass into Spot object

%Written by Matthew J. Gastinger, Bitplane Advanced Application Scientist.  
%March 2014.
%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Center of Mass to Spots" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_SurfacePositionFinal(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Center of Mass to Spots" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_SurfacePositionFinal(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
% 
%Description
%
%This XTension wil collect the XYZ positon of the surface based on the
%center of homogeneous mass, and plot those positions as a Spots object.  It will do
%this in all time points, for all or a selected group of surfaces.

function XT_MJG_SurfacePositionFinal(aImarisApplicationID)
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

% get all surfaces
vSurfaces = vImarisApplication.GetFactory.ToSurfaces(vImarisApplication.GetSurpassSelection);
vNumberOfSurfaces = vSurfaces.GetNumberOfSurfaces;

%%
vSelectedSurfaces = vSurfaces.GetSelectedIndices
if isempty(vSelectedSurfaces)
    %Generate postions from center of mass
    vPositionFinal = vSurfaces.GetCenterOfMass(0)
    for c = 1:vNumberOfSurfaces-1;
        vPositionXYZ = vSurfaces.GetCenterOfMass(c);
        vPositionFinal =[vPositionFinal;vPositionXYZ];
    end;
%Generate timeIndex from surfaces
    vTimeIndexFinal = vSurfaces.GetTimeIndex(0);
    for c = 1:vNumberOfSurfaces-1;
        vTimeIndex = vSurfaces.GetTimeIndex(c);
        vTimeIndexFinal =[vTimeIndexFinal;vTimeIndex];
    end;
%Determine if time laspe and place spots in appropriate timeIndex
    if max(vTimeIndexFinal)>0
        vSpotsTime =vTimeIndexFinal;
    else
        vSpotsTime = zeros(vNumberOfSurfaces,1);
    end;

%%
else
    vPositionFinal=[];
    for s=1:numel(vSelectedSurfaces)
        vPositionXYZ = vSurfaces.GetCenterOfMass(vSelectedSurfaces(s));
        vPositionFinal =[vPositionFinal;vPositionXYZ];
    end
%Generate timeIndex from surfaces
    vTimeIndexFinal = [];
    for c = 1:numel(vSelectedSurfaces);
        vTimeIndex = vSurfaces.GetTimeIndex(vSelectedSurfaces(c));
        vTimeIndexFinal =[vTimeIndexFinal;vTimeIndex];
    end;
%Determine if time laspe and place spots in appropriate timeIndex
    if max(vTimeIndexFinal)>0
        vSpotsTime =vTimeIndexFinal;
    else
        vSpotsTime = zeros(numel(vSelectedSurfaces),1);
    end;

end

%%

%%
%Dialog to select new Spot size
vQuestion = {sprintf(['Please enter new Spot Radius (um):'])};
    vAnswer = inputdlg(vQuestion,'Resize Spots',1,{'1'});
    if isempty(vAnswer), return, end
newRadius = str2double(vAnswer{1});
if isempty(vSelectedSurfaces)
    vSpotsRadius = ones(vNumberOfSurfaces,1)*newRadius;
else 
    vSpotsRadius = ones(numel(vSelectedSurfaces),1)*newRadius;
end
%%
%Create the a new Spots generated from teh center of Mass
vRGBA = vSurfaces.GetColorRGBA;
vNewSpots = vImarisApplication.GetFactory.CreateSpots;
vNewSpots.Set(vPositionFinal, vSpotsTime, vSpotsRadius);
vNewSpots.SetColorRGBA(vRGBA);
vNewSpots.SetName([char(vSurfaces.GetName), ' Center of Mass']);
vNewSpots.SetTrackEdges(vSurfaces.GetTrackEdges)
vSurfaces.SetVisible(0);
vSurfaces.GetParent.AddChild(vNewSpots, -1);

end