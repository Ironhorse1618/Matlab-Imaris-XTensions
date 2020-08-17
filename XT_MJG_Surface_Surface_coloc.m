%Surface surface Colocalization

%Written by Matthew J. Gastinger, Bitplane Advanced Application Scientist.  
%March 2014.
%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Surface-Surface coloc" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_Surface_Surface_coloc(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Surface-Surface coloc" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_Surface_Surface_coloc(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
% 
%Description
%
%This XTension will mask each of 2 surface scenes.  It will find the voxels
%inside each surface that overlap with each other.  A new channel will be
%made, and a new surface generated from the overlapping regions.

function XT_MJG_Surface_Surface_coloc(aImarisApplicationID)
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



%%
% the user has to create a scene with some surfaces
vSurpassScene = vImarisApplication.GetSurpassScene;
if isequal(vSurpassScene, [])
    msgbox('Please create some Surfaces in the Surpass scene!');
    return;
end

%Create a new folder object for new surfaces
Coloc_surfaces = vImarisApplication.GetFactory;
result = Coloc_surfaces.CreateDataContainer;
result.SetName('Coloc surfaces');

%clone Dataset
vDataSet = vImarisApplication.GetDataSet.Clone;
%%
% get all Surpass surfaces names
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

if vNumberOfSurfaces<2
    msgbox('Please create at least 2 surfaces objects!');
    return;
end

vNamesList = vNamesList(1:vNumberOfSurfaces);
%%
%Choose the surfaces
%Choose how many surfaces to colocalize

vPair = [];
while length(vPair) ~= 2
    [vPair, vOk] = listdlg('ListString',vNamesList,'SelectionMode','multiple',...
        'ListSize',[250 150],'Name','Colocalize surfaces','InitialValue',[1,2], ...
        'PromptString',{'Please select the 2 surfaces to colocalize:'});
    if vOk<1, return, end
    if length(vPair) ~= 2
        vHandle = msgbox(['Please select at least two (2) Surface objects. Use "Control" and left ', ...
            'click to select/unselect an object of the list.']);
        uiwait(vHandle);
    end
end

NumberColocSurfaces=numel(vPair);
vSurfaces1 = vSurfacesList{vPair(1)};
vSurfaces2 = vSurfacesList{vPair(2)};

%%
%Get Image Data parameters
aExtendMaxX = vImarisApplication.GetDataSet.GetExtendMaxX;
aExtendMaxY = vImarisApplication.GetDataSet.GetExtendMaxY;
aExtendMaxZ = vImarisApplication.GetDataSet.GetExtendMaxZ;
aExtendMinX = vImarisApplication.GetDataSet.GetExtendMinX;
aExtendMinY = vImarisApplication.GetDataSet.GetExtendMinY;
aExtendMinZ = vImarisApplication.GetDataSet.GetExtendMinZ;
aSizeX = vImarisApplication.GetDataSet.GetSizeX;
aSizeY = vImarisApplication.GetDataSet.GetSizeY;
aSizeZ = vImarisApplication.GetDataSet.GetSizeZ;
aSizeC = vImarisApplication.GetDataSet.GetSizeC;
aSizeT = vImarisApplication.GetDataSet.GetSizeT;
Xvoxelspacing= (aExtendMaxX-aExtendMinX)/aSizeX;
vSmoothingFactor=Xvoxelspacing*2;

%add additional channel
vDataSet.SetSizeC(aSizeC + 1);
TotalNumberofChannels=aSizeC+1;
vLastChannel=TotalNumberofChannels-1;

%%
%Generate dialog for smoothed or unsmoothed surface
%Have option to chose: 1) no smoothing, 2) default, 3) Custom

vSurfaceSmoothingType = false;
qstring={'Please choose a how to generate Coloc surface '}    
vAnswer = questdlg(qstring, 'Colocalize Surfaces', ...
		'No Smoothing', 'Smoothing', 'Smoothing');
if(isequal(vAnswer, 'Cancel') || isempty(vAnswer))
    return
end
vSurfaceSmoothingType = isequal(vAnswer, 'Smoothing');
if vSurfaceSmoothingType==true
    vSmoothingFactorName = num2str(vSmoothingFactor);
    qstring={'Please set the smoothing factor -- default is double image voxel size'};
    vAnswer2 = inputdlg(qstring,'Smoothing Factor',1,{vSmoothingFactorName});
    if isempty(vAnswer2), return, end
else    
vSmoothingFactor = 0;
end

%%

tic
%Generate surface mask for each surface over time 
for vTimeIndex= 0:aSizeT-1
    vSurfaces1Mask = vSurfaces1.GetMask(aExtendMinX,aExtendMinY,aExtendMinZ,aExtendMaxX,aExtendMaxY,aExtendMaxZ,aSizeX, aSizeY,aSizeZ,vTimeIndex);
    vSurfaces2Mask = vSurfaces2.GetMask(aExtendMinX,aExtendMinY,aExtendMinZ,aExtendMaxX,aExtendMaxY,aExtendMaxZ,aSizeX, aSizeY,aSizeZ,vTimeIndex);
    
    for vIndexZ = 1:aSizeZ
        ch1=vSurfaces1Mask.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1,0,0,aSizeX,aSizeY,1);
        ch2=vSurfaces2Mask.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1,0,0,aSizeX,aSizeY,1);
    
    
    %Determine the Voxels that are colocalized
        Coloc=ch1+ch2;
        Coloc(Coloc<2)=0;
        Coloc(Coloc>1)=1;
        vDataSet.SetDataSubVolumeAs1DArrayBytes(Coloc, ...
            0,0,vIndexZ-1,vLastChannel,vTimeIndex,aSizeX,aSizeY,1);
    end
end    

vDataSet.SetChannelName(vLastChannel,'ColocChannel');
vDataSet.SetChannelRange(vLastChannel,0,1);

vImarisApplication.SetDataSet(vDataSet);

%Run the Surface Creation Wizard on the new channel
ip = vImarisApplication.GetImageProcessing;
Coloc_surfaces1 = ip.DetectSurfaces(vDataSet, [], vLastChannel, vSmoothingFactor, 0, true, 55, '');
Coloc_surfaces1.SetName(sprintf('ColocSurface'));
Coloc_surfaces1.SetColorRGBA((rand(1, 1)) * 256 * 256 * 256 );

%Add new surface to Surpass Scene
vSurfaces1.SetVisible(0);
vSurfaces2.SetVisible(0);

result.AddChild(Coloc_surfaces1, -1);
vImarisApplication.GetSurpassScene.AddChild(result, -1);
toc


