%Surface Surface Contact Area

%Written by Matthew J. Gastinger, Bitplane Advanced Application Scientist.
%Janurary 2018.  Modified to work better in Imaris 9.x
%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Surface-Surface Contact Area" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_SurfaceSurfaceContactArea2(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Surface-Surface Contact Area" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_SurfaceSurfaceContactArea2(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description
%This XTension will find the surface contact area between 2 surfaces.  The
%primary surface is the base, and secondary is one covering the primary.
%
%The result of the XTension will generate a one voxel thick unsmoothed
%surface object above the primary surface representing where the 2 surfaces
%physically overlap.
%
%Two new statistics will be generated.  1)The first will be a total surface
%area of each new surface object.  The measurement will be estimate by
%taking the number of voxels and multiplying by the area a a single (XY
%pixel).  2) The second statistic will be in the "overall" tab, reporting
%the percentage of surface contact area relative to the total surface area
%of the primary surfaces.


function XT_MJG_SurfaceSurfaceContactArea2(aImarisApplicationID)

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

% the user has to create a scene with some surfaces
vSurpassScene = vImarisApplication.GetSurpassScene;
if isequal(vSurpassScene, [])
    msgbox('Please create some Surfaces in the Surpass scene!');
    return;
end

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
vPrimarySurface=1;
vSecondarySurface=1;
vPair = [];
while vPrimarySurface==vSecondarySurface
    %Create Dialog box and allow user to choose the Reference Position
    vPair = [];
    [vPair, vOk] = listdlg('ListString',vNamesList,'SelectionMode','multiple',...
        'ListSize',[250 150],'Name','Surface-Surface Contact Area','InitialValue',[1,1], ...
        'PromptString',{'Please select Primary Surface'});
    vPrimarySurface = vSurfacesList{vPair(1)};
    NumberPrimarySurfaces=vPrimarySurface.GetNumberOfSurfaces;
    %Create Dialog box and allow user to choose the Reference Position
    vPair = [];
    [vPair, vOk] = listdlg('ListString',vNamesList,'SelectionMode','multiple',...
        'ListSize',[250 150],'Name','Surface-Surface Contact Area','InitialValue',[1,1], ...
        'PromptString',{'Please select Secondary coverage surface'});
    vSecondarySurface = vSurfacesList{vPair(1)};
    if vPrimarySurface==vSecondarySurface
        uiwait(msgbox('Please choose 2 different surfaces'));
        
    end
end
str=sprintf('Calculate amount of %s that is in contact with %s', char(vSecondarySurface.GetName),char(vPrimarySurface.GetName));
choice=questdlg(str, 'Surface-Surface contact area analysis',...
    'Analyze','Cancel','Analyze');
%Handle Response
switch choice
    case 'Analyze'
    case 'Cancel'
        return;
end

%%
%Get Image Data parameters
vImarisDataSet = vImarisApplication.GetDataSet.Clone;
vNumberOfChannels = vImarisDataSet.GetSizeC;
vImarisDataSet.SetSizeC(vNumberOfChannels + 3);

vDataMin = [vImarisDataSet.GetExtendMinX, vImarisDataSet.GetExtendMinY, vImarisDataSet.GetExtendMinZ];
vDataMax = [vImarisDataSet.GetExtendMaxX, vImarisDataSet.GetExtendMaxY, vImarisDataSet.GetExtendMaxZ];
vDataSize = [vImarisDataSet.GetSizeX, vImarisDataSet.GetSizeY, vImarisDataSet.GetSizeZ];
aSizeT = vImarisApplication.GetDataSet.GetSizeT;
Xvoxelspacing = (vDataMax(1)-vDataMin(1))/vDataSize(1);
Yvoxelspacing = (vDataMax(2)-vDataMin(2))/vDataSize(2);
Zvoxelspacing = (vDataMax(3)-vDataMin(3))/vDataSize(3);

vSmoothingFactor=Xvoxelspacing*2;
ZLimit=((3*Xvoxelspacing)+100*Xvoxelspacing)/100;%test percent Xvoxelsize

vImarisDataSet = vImarisApplication.GetDataSet.Clone;
vNumberOfChannels = vImarisDataSet.GetSizeC;
vImarisDataSet.SetSizeC(vNumberOfChannels + 3);

%Convert to 32bit
vFloatType = vImarisDataSet.GetType.eTypeFloat;
vImarisDataSet.SetType(vFloatType);

%Convert data to isotropic voxels
if Zvoxelspacing>ZLimit
    endZadj=round(double(Zvoxelspacing/Xvoxelspacing*vDataSize(3)));
    str=sprintf('Warning this volume is does NOT have isotropic voxels.\n\n To properly run this application go to:\n Edit--""Resample3D"" and adjust Z-size to %d steps',endZadj);
    choice2=questdlg(str,'WARNING','Cancel','Cancel');
    %Handle Response
    switch choice2
        case 'Cancel'
            return;
    end
end

vSpotsXYZ = [];
vSpotsTime = [];

vProgressDisplay = waitbar(0, 'Distance Transform: Preparation');

% Create a new channel where the result will be sent
vImarisDataSet.SetChannelName(vNumberOfChannels,['Distance to ', char(vPrimarySurface.GetName)]);
vImarisDataSet.SetChannelColorRGBA(vNumberOfChannels, 255*256*256);
SurfaceIndex=[];
%vImarisApplication.SetVisible(~vImarisApplication.GetVisible);

for vTime = 0:aSizeT-1;
    % Get the mask DataSet
    vMaskDataSetPrimary = vPrimarySurface.GetMask( ...
        vDataMin(1), vDataMin(2), vDataMin(3), ...
        vDataMax(1), vDataMax(2), vDataMax(3), ...
        vDataSize(1), vDataSize(2), vDataSize(3), vTime);
    
    vMaskDataSetSecondary = vSecondarySurface.GetMask( ...
        vDataMin(1), vDataMin(2), vDataMin(3), ...
        vDataMax(1), vDataMax(2), vDataMax(3), ...
        vDataSize(1), vDataSize(2), vDataSize(3), vTime);
    
    for vIndexZ = 1:vDataSize(3)
        %Mask for Primary surface (inside=1)
        vSlice = vMaskDataSetPrimary.GetDataSliceBytes(vIndexZ-1, 0, 0);
        vSlice = vSlice == 1;
        %set set and last slice to all zeros
        if vIndexZ==1 | vIndexZ==vDataSize(3)
            vSlice(vSlice==1)=0;%Replace ones with zeros
        end
        vSlice(:,end)=0;%Replace border values to all zeros
        vSlice(end,:)=0;%Replace border values to all zeros
        vSlice(:,1)=0;%Replace border values to all zeros
        vSlice(1,:)=0;%Replace border values to all zeros

        
       
        vSlice=int8(reshape(vSlice,[numel(vSlice),1]));%covert to 1DArray and int8
        vImarisDataSet.SetDataSubVolumeAs1DArrayBytes(vSlice, ...
            0,0,vIndexZ-1,vNumberOfChannels,vTime,vDataSize(1),vDataSize(2),1);
        
        %Mask for secondary Surface (inside=1)
        vSlice2=vMaskDataSetSecondary.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
        vSlice2 = vSlice2 == 1;% for inside surface only
        vImarisDataSet.SetDataSubVolumeAs1DArrayFloats(vSlice2, ...
            0,0,vIndexZ-1,vNumberOfChannels+1,vTime,vDataSize(1),vDataSize(2),1);
        
        waitbar((vIndexZ+1)/vDataSize(3)/2, vProgressDisplay);
    end
    waitbar((vTime+1)/aSizeT/2, vProgressDisplay);
end

waitbar(0.5, vProgressDisplay, 'Distance Transform: Calculation');
vImarisApplication.GetImageProcessing.DistanceTransformChannel( ...
    vImarisDataSet, vNumberOfChannels, 1, false);
waitbar(1, vProgressDisplay);
close(vProgressDisplay);

% inside secondary surface mask%
%vImarisApplication.SetVisible(~vImarisApplication.GetVisible);
vProgressTotalCount = aSizeT*vDataSize(3);
vProgressDisplay = waitbar(0, 'Surface Contact Area calculation');
vProgressCount = 0;
for vTime2 = 1:aSizeT
    for vSlice3 = 1:vDataSize(3)
        ch1=double(vImarisDataSet.GetDataSubVolumeAs1DArrayFloats(0, 0, vSlice3-1, vNumberOfChannels, vTime2-1, vDataSize(1), vDataSize(2), 1));
        ch2=double(vImarisDataSet.GetDataSubVolumeAs1DArrayFloats(0, 0, vSlice3-1, vNumberOfChannels+1, vTime2-1, vDataSize(1), vDataSize(2), 1));
        %Test slice for surface present
        if ch1>999999999999
            break
        end
        %Multiply distance transform by secondary surface mask
        vData=(ch1.*ch2);
        vImarisDataSet.SetDataSubVolumeAs1DArrayFloats(vData,0,0,vSlice3-1,vNumberOfChannels+2,vTime2-1,vDataSize(1),vDataSize(2),1);
        vProgressCount = vProgressCount + 1;
        waitbar(vProgressCount/vProgressTotalCount, vProgressDisplay)
    end
end

%vImarisApplication.SetVisible(~vImarisApplication.GetVisible);
%vImarisDataSet = vImarisApplication.GetDataSet;
vImarisApplication.SetDataSet(vImarisDataSet);

vPrimarySurface.SetVisible(1);
vSecondarySurface.SetVisible(0);

%Create a new folder object for new surfaces
Newsurfaces = vImarisApplication.GetFactory;
result = Newsurfaces.CreateDataContainer;
result.SetName(sprintf('Surface-Surface contact - %s',char(vPrimarySurface.GetName)));

%Generate Total surface area
ip = vImarisApplication.GetImageProcessing;
UpperThreshold=1.5*Xvoxelspacing;
vAllSurfaceArea=ip.DetectSurfacesWithUpperThreshold(vImarisDataSet,[],vNumberOfChannels,0,0,true,false,0.003,true,false,UpperThreshold,'');
vAllSurfaceArea.SetName(sprintf('Total Surface Area - %s',char(vSecondarySurface.GetName)));
vRGBA=[255,255,255, 0];%for yellow
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); % combine different components (four bytes) into one integer
vAllSurfaceArea.SetColorRGBA(vRGBA);
%Add new surface to Surpass Scene
result.AddChild(vAllSurfaceArea, -1);
vImarisApplication.GetSurpassScene.AddChild(result, -1);

%Get Original Stats from Total surface area
vAllPrimaryStatistics = vAllSurfaceArea.GetStatistics;
vPrimarySurfaceStatNames = cell(vAllPrimaryStatistics.mNames);
vPrimarySurfaceStatValues = vAllPrimaryStatistics.mValues;
vSurfaceIndex1=strmatch('Total Number of Voxels', vPrimarySurfaceStatNames);
vAllSurfacesVoxels = vPrimarySurfaceStatValues(vSurfaceIndex1,:);
vSurfacesArea = vAllSurfacesVoxels*Xvoxelspacing^2;

%Generate contact surface surface
ip = vImarisApplication.GetImageProcessing;
UpperThreshold=1.5*Xvoxelspacing;
vSurfaceContact=ip.DetectSurfacesWithUpperThreshold(vImarisDataSet,[],vNumberOfChannels+2,0,0,true,false,0.003,true,false,UpperThreshold,'');
vSurfaceContact.SetName(sprintf('Surface contact Area - %s',char(vSecondarySurface.GetName)));
vRGBA=[255,255,0, 0];%for yellow
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); % combine different components (four bytes) into one integer
vSurfaceContact.SetColorRGBA(vRGBA);

%Add new surface to Surpass Scene
result.AddChild(vSurfaceContact, -1);
vImarisApplication.GetSurpassScene.AddChild(result, -1);
%Remove/crop the working channels
vImarisApplication.GetDataSet.Crop(0,vDataSize(1),0,vDataSize(2),0,vImarisApplication.GetDataSet.GetSizeZ,0,vNumberOfChannels,0,aSizeT);

%Get surface contact area
vSurfaceContactStatistics = vSurfaceContact.GetStatistics;
vSurfaceContactStatNames = cell(vSurfaceContactStatistics.mNames);
vSurfaceContactStatValues = vSurfaceContactStatistics.mValues;
vSurfaceContactStatIds = vSurfaceContactStatistics.mIds;
vSurfaceIndex1=strmatch('Number of Voxels', vSurfaceContactStatNames);
vSurfaceContactVoxels = vSurfaceContactStatValues(vSurfaceIndex1,:);
[vSurfaceIDssorted indID] = sort(vSurfaceContactStatIds(vSurfaceIndex1,:));
vSurfaceContactArea = vSurfaceContactVoxels*Xvoxelspacing*Yvoxelspacing;

vSurfaceContactAreaFinal=vSurfaceContactArea(indID,:);


if vSurfaceContactStatValues==0
    h=msgbox('No surface contact between selected Surfaces');
else
    
    %Generate Surface contact Area stat
    vInd=1:vSurfaceContact.GetNumberOfSurfaces;
    vIds = vInd - 1;
    vUnits(vInd) = { char(vImarisApplication.GetDataSet.GetUnit) };
    vFactors(vInd) = {'Surface'};
    vNames(vInd) = {' Contact surface area, um^2'};
    vFactorNames = {'Category'};
    vSurfaceContact.AddStatistics(vNames, vSurfaceContactAreaFinal, vUnits, vFactors, vFactorNames, vIds);
    clearvars vInd vNames vFactors vUnits
    %Insert an overall Statistio
    vInd=1;%:vSurfaceContact.GetNumberOfSurfaces;
    vIds = 1;%vInd - 1;
    vUnits(vInd) = { char(vImarisApplication.GetDataSet.GetUnit) };
    vFactors(vInd) = {'Overall'};
    vNames(vInd) = {'% SurfaceArea coverage'};
    vFactorNames = {'Overall'};
    vSurfaceContact.AddStatistics(vNames, sum(vSurfaceContactArea)/vSurfacesArea*100, vUnits, vFactors, vFactorNames, vIds);
end
close(vProgressDisplay);

