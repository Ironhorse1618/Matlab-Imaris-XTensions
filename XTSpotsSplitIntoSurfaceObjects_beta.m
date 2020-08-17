%
%
%  Spots Split Into Surface Objects Function for Imaris 9.3
%
%  Copyright Bitplane 2018 (modified by Igor) - fixed rounding erro for 2D
%  surfaces
%
%
%  Installation:
%
%  - Copy this file into the XTensions folder in the Imaris installation directory.
%  - You will find this function in the Image Processing menu
%
%    <CustomTools>
%      <Menu>
%         <Submenu name="Spots Functions">        
%        <Item name="Split Spots Into Surface Objects_beta" icon="Matlab" tooltip="Spots split into surface objects">
%          <Command>MatlabXT::XTSpotsSplitIntoSurfaceObjects_beta(%i)</Command>
%        </Item>
%         </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSpots">
%          <Item name="Split Into Surface Objects_beta" icon="Matlab" tooltip="Spots split into surface objects">
%            <Command>MatlabXT::XTSpotsSplitIntoSurfaceObjects_beta(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
% 
%
%  Description:
%   
%   For each surfaces object in the same folder of the spots, create a new 
%       spots object containing the spots that lies inside the surfaces.
%
%

function XTSpotsSplitIntoSurfaceObjects_beta(aImarisApplicationID)

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
%%
% get the spots
vSpots = vImarisApplication.GetFactory.ToSpots(vImarisApplication.GetSurpassSelection);
vSurfaces = vImarisApplication.GetFactory.ToSurfaces(vImarisApplication.GetSurpassSelection);

vSurfacesSelected = vImarisApplication.GetFactory.IsSurfaces(vSurfaces);

if vSurfacesSelected
    vScene = vSurfaces.GetParent;
else
    vScene = vImarisApplication.GetSurpassScene;
end
vNumberOfSurfaces = 0;
vSurfacesList{vScene.GetNumberOfChildren} = [];
vNamesList{vScene.GetNumberOfChildren} = [];
for vChildIndex = 1:vScene.GetNumberOfChildren
    vDataItem = vScene.GetChild(vChildIndex - 1);
    if vImarisApplication.GetFactory.IsSurfaces(vDataItem)
        vNumberOfSurfaces = vNumberOfSurfaces+1;
        vSurfacesList{vNumberOfSurfaces} = vImarisApplication.GetFactory.ToSurfaces(vDataItem);
        vNamesList{vNumberOfSurfaces} = char(vDataItem.GetName);
    end
end

if vNumberOfSurfaces==0
    msgbox('Please create at a surfaces object!');
    return;
end

vNamesList = vNamesList(1:vNumberOfSurfaces);
%%
%Choose the surfaces
if vNumberOfSurfaces > 1
    [vPair, vOk] = listdlg('ListString',vNamesList,'SelectionMode','single',...
        'ListSize',[250 150],'Name','Split Spots into Surface Object','InitialValue',1, ...
        'PromptString',{'Please select the Surface'});
    if vOk<1, return, end
    vSurfaces1 = vSurfacesList{vPair(1)}; 
else
    vSurfaces1 = vSurfacesList{1}; 
end
vNumberOfSurfaces1=vSurfaces1.GetNumberOfSurfaces;

%%

if ~vImarisApplication.GetFactory.IsSpots(vSpots)  
  msgbox('Please select some Spots!');
  return;
end

% get the spots coordinates
vSpotsXYZ = vSpots.GetPositionsXYZ;
vSpotsTime = vSpots.GetIndicesT;
vSpotsRadius = vSpots.GetRadiiXYZ;
vSpotsName = char(vSpots.GetName);

vSpots.SetVisible(false);
vTimeInterval = min(vSpotsTime):max(vSpotsTime);
vIndicesSpotsTime = cell(numel(vTimeInterval), 1);
for vTime = vTimeInterval
  vIndicesSpotsTime{vTime - vTimeInterval(1) + 1} = find(vSpotsTime == vTime);
end

% get the parent group
vParentGroup = vSpots.GetParent;
vNewGroup = vImarisApplication.GetFactory.CreateDataContainer;
vNewGroup.SetName([vSpotsName, ' inside Surfaces']);
vParentGroup.AddChild(vNewGroup, -1);

% mask volume
vMin = min(vSpotsXYZ);
vMax = max(vSpotsXYZ);

% add 1% border to be sure to include all spots (avoid edge effects)
vDelta = vMax - vMin;
if ~any(vDelta > 0)
  vDelta = [1, 1, 1];
end
vDim2D = find(vDelta == 0);
vDelta(vDim2D) = mean(vDelta(vDelta > 0));
vMin = vMin - vDelta*0.005;
vMax = vMax + vDelta*0.005;

vMaskSize = 350;
vMaxMaskSize = vMaskSize / max(vMax - vMin);

% spots coordinates on the mask
vSpotsOnMaskXYZ = zeros(size(vSpotsXYZ));
for vDim = 1:3
  vSpotsOnMaskXYZ(:, vDim) = ceil((vSpotsXYZ(:, vDim) - vMin(vDim)) * vMaxMaskSize);
end

% the zeros belongs to the first interval
vSpotsOnMaskXYZ = max(vSpotsOnMaskXYZ, 1); 

vMaskSize = ceil((vMax - vMin) * vMaxMaskSize);

vMaskSize(vDim2D) = 1;
vSpotsOnMaskXYZ(:, vDim2D) = 1;

vProgressDisplay = waitbar(0, 'Splitting spots into surfaces');

% loop through each surface
for vSurfaceIndex = 0:vNumberOfSurfaces1-1
  vAllIndices = vSpotsTime;
  vAllIndicesSize = 0;

  vMask = vSurfaces1.GetSingleMask(vSurfaceIndex, ...
    vMin(1), vMin(2), vMin(3), vMax(1), vMax(2), vMax(3),...
    int32(vMaskSize(1)), int32(vMaskSize(2)), int32(vMaskSize(3)));

  vTimeIndex = vSurfaces1.GetTimeIndex(vSurfaceIndex);
  vMaskImage = vMask.GetDataVolumeAs1DArrayBytes(0, 0);
  vMaskImage = reshape(vMaskImage, vMaskSize);
  
  % search the element of the spot that lies inside the surface
  vIndexSpotsTime = vIndicesSpotsTime{vTimeIndex - vTimeInterval(1) + 1};
  vSpotsCoords = vSpotsOnMaskXYZ(vIndexSpotsTime, :);
  vIndexSpotsInside = vMaskImage( ...
     vSpotsCoords(:, 1) + ...
    (vSpotsCoords(:, 2)-1)*vMaskSize(1) + ...
    (vSpotsCoords(:, 3)-1)*vMaskSize(1)*vMaskSize(2)) == 1;
  vIndexSpotsInside = vIndexSpotsTime(vIndexSpotsInside);

  % copy to complete list
  vSize = numel(vIndexSpotsInside);
  vAllIndices(vAllIndicesSize + (1:vSize)) = vIndexSpotsInside;
  vAllIndicesSize = vAllIndicesSize + vSize;

  vAllIndices = vAllIndices(1:vAllIndicesSize);

  vSpotsInside = vImarisApplication.GetFactory.CreateSpots;
  vSpotsInside.Set(vSpotsXYZ(vAllIndices, :), vSpotsTime(vAllIndices), ...
    zeros(sum(vAllIndices~=0),1));
  vSpotsInside.SetRadiiXYZ(vSpotsRadius(vAllIndices,:));
  vSpotsInside.SetName(sprintf('%s inside %s [%i] t%i', ...
    vSpotsName, char(vSurfaces1.GetName), vSurfaceIndex + 1, vTimeIndex));
  vNumberOfSurfaces = max(vNumberOfSurfaces, 2);
  %vRed = (1-vSurfaceIndex/(vNumberOfSurfaces-1)) * 255;
  %vGreen = (1-vChildIndex/(vNumberOfChildren-1)) * 255;
  %vBlue = (vSurfaceIndex/(vNumberOfSurfaces-1)) * 255;
  vSpotsInside.SetColorRGBA((rand(1, 1)) * 256 * 256 * 256 );

  vNewGroup.AddChild(vSpotsInside, -1);

  waitbar((vChildIndex + vSurfaceIndex/vNumberOfSurfaces1) / ...
    1, vProgressDisplay);
end

close(vProgressDisplay);


