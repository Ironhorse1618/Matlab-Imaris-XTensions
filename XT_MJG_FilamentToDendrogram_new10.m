%
%  Filament to Dendrogram - morphological stats
%
%  Copyright Bitplane BPI 2018
%  Matthew Gastinger, Ph.D
%  Advanced Application Scientist
%
%
%  Installation:
%
%  - Copy this file into the XTensions folder in the Imaris installation directory.
%  - You will find this function in the Image Processing menu
%
%    <CustomTools>
%      <Menu>
%       <Submenu name="Filaments Functions">
%        <Item name="MJG_Filament Branch Analysis" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_FilamentToDendrogram_new10(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpFilaments">
%          <Item name="MJG_Filament Branch Analysis" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_FilamentToDendrogram_new10(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description:  SO far this XTension is classify the bifurcation nodes into
%3 categories.  1) Arborization, 2)Continuation, 3) Terminal
%1)Arborization (A) nodes have two bifurcating children. 
%2)Continuation (C) nodes have one bifurcating and one terminating child. 
%3)Termination (T) nodes have two terminating children.
%
%   A map of the dendritic tree will be made from these node
%   classifications, allowing for a advanced level of dendrtie tree
%   morphology

function XT_MJG_FilamentToDendrogram_new10(aImarisApplicationID)

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
tic
% get the filament
vFactory = vImarisApplication.GetFactory;
vFilaments = vFactory.ToFilaments(vImarisApplication.GetSurpassSelection);


% search the filament if not previously selected
vSurpassScene = vImarisApplication.GetSurpassScene;
if ~vFactory.IsFilaments(vFilaments)
    for vChildIndex = 1:vSurpassScene.GetNumberOfChildren
        vDataItem = vSurpassScene.GetChild(vChildIndex - 1);
        if vFactory.IsFilaments(vDataItem)
            vFilaments = vFactory.ToFilaments(vDataItem);
            break;
        end
    end
    % did we find the filament?
    if isequal(vFilaments, [])
        msgbox('Please create some filament!');
        return;
    end
end
%%
%Get Image Data parameters
vDataMin = [vImarisApplication.GetDataSet.GetExtendMinX, vImarisApplication.GetDataSet.GetExtendMinY, vImarisApplication.GetDataSet.GetExtendMinZ];
vDataMax = [vImarisApplication.GetDataSet.GetExtendMaxX, vImarisApplication.GetDataSet.GetExtendMaxY, vImarisApplication.GetDataSet.GetExtendMaxZ];
vDataSize = [vImarisApplication.GetDataSet.GetSizeX, vImarisApplication.GetDataSet.GetSizeY, vImarisApplication.GetDataSet.GetSizeZ];
vNumberOfChannels = vImarisApplication.GetDataSet.GetSizeC;
Xvoxelspacing = (vDataMax(1)-vDataMin(1))/vDataSize(1);
Yvoxelspacing = (vDataMax(2)-vDataMin(2))/vDataSize(2);
Zvoxelspacing = (vDataMax(3)-vDataMin(3))/vDataSize(3);
ZLimit=((3*Xvoxelspacing)+100*Xvoxelspacing)/100;%test percent Xvoxelsize
vThreshold=Xvoxelspacing*2;

% create an array of structures (one element per time point)
vCount = vFilaments.GetNumberOfFilaments;

for FilamentIndex=0:vCount-1
    %%
    % create a structure to hold the data of one filament
    % each element of the structure is a list of XYZR
    vData = CollectFilamentsData(zeros(0, 3), [], zeros(0, 2), [], 0);
    vFieldNames = fieldnames(vData);
    vFieldsCount = numel(vFieldNames);
    vTimes = 1:vCount;
    
    
    for vIndex = 1:vCount
        vTimes(vIndex) = vFilaments.GetTimeIndex(vIndex-1)+1;
    end
    vSizeT = max(vTimes);
    vData = vData(ones(vSizeT, 1));
    %%
    vNumberOfChannels=vImarisApplication.GetDataSet.GetSizeC;
    vRGBA = vFilaments.GetColorRGBA;
    vNewFilament=vImarisApplication.GetFactory.CreateFilaments;
    vTotalShollSpots=0;
    
    FilamentCount=0;
    
    vFilamentsIndexT = vFilaments.GetTimeIndex(FilamentIndex);
    vFilamentsXYZ = vFilaments.GetPositionsXYZ(FilamentIndex);
    vFilamentsRadius = vFilaments.GetRadii(FilamentIndex);
    vFilamentsEdgesSegmentId = vFilaments.GetEdgesSegmentId(FilamentIndex);
    vSegmentIds=unique(vFilamentsEdgesSegmentId);%Idenitfy unique filament segmentIDs
    vFilamentsEdges = vFilaments.GetEdges(FilamentIndex) + 1;
    vTypes = vFilaments.GetTypes(FilamentIndex);
    vBeginningVertex = vFilaments.GetBeginningVertexIndex(FilamentIndex) + 1;
    vIndexT = vFilaments.GetTimeIndex(FilamentIndex);
    
     %Get Original Stats order vNewSpots
%      vAllStatistics = vFilaments.GetStatistics;
%      vFilamentsStatNames = cell(vAllStatistics.mNames);
%      vFilamentsStatValues = vAllStatistics.mValues;
%      vFilamentStatIds = vAllStatistics.mIds;
%      vFilamentsIndex1=strmatch('Dendrite Branch Depth', vFilamentsStatNames);
%      vFilamentsIndex2=strmatch('Dendrite Branch Level', vFilamentsStatNames);
%      [vDendriteIDs vSortIndex]= sort(vFilamentStatIds(vFilamentsIndex1,:));
%      vDendriteBranchDepthpresort = vFilamentsStatValues(vFilamentsIndex1,:);
%      vDendriteBranchLevelpresort = vFilamentsStatValues(vFilamentsIndex2,:);
%      vDendriteBranchDepth=vDendriteBranchDepthpresort(vSortIndex,:);
%      vDendriteBranchLevel=vDendriteBranchLevelpresort(vSortIndex,:);
    
    %Test if the time point has empty filament matrix or filament start
    %point and nothing more
    IsFilamentTest = numel(vFilamentsRadius);
    if IsFilamentTest==1 | isempty(vFilamentsXYZ)
        qIsAFilament='NO';
        continue;
    else
        FilamentCount=FilamentCount+1;
    end
    
    %%
    
    if isempty(vFilamentsXYZ)
        vFilamentsXYZ = zeros(0, 3);
    end
    if isempty(vFilamentsEdges)
        vFilamentsEdges = zeros(0, 2);
    end
    
    vThisData = CollectFilamentsData(vFilamentsXYZ, vFilamentsRadius, ...
        vFilamentsEdges, vTypes, vBeginningVertex);
    
    vTime = vTimes(FilamentIndex + 1);
    for vField = 1:vFieldsCount
        vFieldName = vFieldNames{vField};
        if isempty(vData(vTime).(vFieldName))
            vData(vTime).(vFieldName) = vThisData.(vFieldName);
        elseif ~isempty(vThisData.(vFieldName))
            vData(vTime).(vFieldName) = [vData(vTime).(vFieldName); vThisData.(vFieldName)];
        end
    end
    
    
    %%
    % generate spots
    vParent = vImarisApplication.GetFactory.CreateDataContainer;
    vParent.SetName([char(vFilaments.GetName), '_Points ' num2str(FilamentIndex)]);
    vTracker = vImarisApplication.GetImageProcessing;
    for vField = 1:vFieldsCount
        vFieldName = vFieldNames{vField};
        
        vSizes = 0:vSizeT;
        for vTime = 1:vSizeT
            vSizes(vTime + 1) = vSizes(vTime) + size(vData(vTime).(vFieldName), 1);
        end
        vXYZRT = zeros(vSizes(vSizeT), 5);
        for vTime = 1:vSizeT
            vTimes = vSizes(vTime) + 1: vSizes(vTime + 1);
            vXYZRT(vTimes, 1:4) = vData(vTime).(vFieldName);
            vXYZRT(vTimes, 5) = vTime - 1;
        end
        if vField==1
            vFilamentStartingPoint=vXYZRT(1:3);
        end
        
        if vField==2
            vFilamentBranchPoints=vXYZRT(:,1:3);
        end
        
        if vField==3
            vFilamentTerminalPoints=vXYZRT(:,1:3);
        end
        
        %vSpots = vImarisApplication.GetFactory.CreateSpots;
        %vSpots.Set(vXYZRT(:, 1:3), vXYZRT(:, 5), vXYZRT(:, 4));
        
        %vSpots.SetName([char(vFilaments.GetName), '_', vFieldName]);
        %vColor = hsv2rgb([vField/vFieldsCount, 1, 1]);
        %vSpots.SetColorRGBA(round(vColor * 255) * [1; 256; 256*256]);
        
        % update scene
        %vParent.AddChild(vSpots, -1);
        
    end
    %%
    %Walk thought for each dendrite segment to find
    vNumberOfDendriteBranches = numel(vSegmentIds);%total number dendrite segments
    %Calculate thh new edges for each dendrite segment
    vEdgesNEW = reshape(repmat([0:2:vNumberOfDendriteBranches*2-1,1:2:vNumberOfDendriteBranches*2],2,1),[],2);
    vDendriteIDNEW = repelem((1:vNumberOfDendriteBranches)',2,1);
    
    vDendriteSegmentID=[];
    vDendritesAll=[];
    vDendritesAllRadius=[];
    DendriteCount=1;
    for vBranchIndex=1:(vNumberOfDendriteBranches)
        %Set the ID for dendrite segment
        wSegmentIndex = vSegmentIds(vBranchIndex);
        vDendriteDistance(vBranchIndex)=0;
        %Logical arguement to identify spots in segment
        vSpotsT = vFilamentsEdgesSegmentId == wSegmentIndex;
        vSpotsIndex = find(vSpotsT');
        %Identify position for dendrite segment
        %Test with new method of filtering
        vDendriteEdgesWorking=vFilamentsEdges(vSpotsIndex,:);
        vEdgesUnique=unique(vDendriteEdgesWorking);
        vDendritePositionsWorking=vFilamentsXYZ(vEdgesUnique,1:3);
        vDendriteRadiusWorking=vFilamentsRadius(vEdgesUnique,1);
        vTypesWorking=vTypes(vEdgesUnique,:);
        
        %%
        %%ReOrdering mismatched filament points from filament editing
        %Test working edges to see if reordering is necessary
        %Calculate if sequential and set as true if one or more elements
        %are not sequential it will be false.
        if all(diff(vDendriteEdgesWorking)==1)==true
            vDendritesAll=[vDendritesAll;vDendritePositionsWorking];
            vDendriteSegmentID=[vDendriteSegmentID;repmat(DendriteCount,numel(vDendritePositionsWorking(:,1)),1)];
            %make list of Radius and distance to neuron
            vDendritesAllRadius=[vDendritesAllRadius;vDendriteRadiusWorking];
            for k=2:numel(vDendritePositionsWorking(:,1))
                vX = vDendritePositionsWorking(k-1,1) - vDendritePositionsWorking(k,1);
                vY = vDendritePositionsWorking(k-1,2) - vDendritePositionsWorking(k,2);
                vZ = vDendritePositionsWorking(k-1,3) - vDendritePositionsWorking(k,3);
                vDendriteDistance(vBranchIndex) = vDendriteDistance(vBranchIndex)+...
                    (sqrt(vX.^2 + vY.^2 + vZ.^2));
            end 
            DendriteCount=DendriteCount+1;
            continue
        end;
        
        %find single indices, identifying first and last points of the
        %segment
        vDendriteStartEndIndex=find(accumarray(vDendriteEdgesWorking(:),1)==1);
        %Odd situation when there are not single indices that identify the
        %beginning or end.
        if isempty (vDendriteStartEndIndex)|numel(vDendriteStartEndIndex)==1
            vDendritesAll=[vDendritesAll;vDendritePositionsWorking];
            vDendriteSegmentID=[vDendriteSegmentID;repmat(DendriteCount,numel(vDendritePositionsWorking(:,1)),1)];
            %make list of Radius and distance to neuron
            vDendritesAllRadius=[vDendritesAllRadius;vDendriteRadiusWorking];
            for k=2:numel(vDendritePositionsWorking(:,1))
                vX = vDendritePositionsWorking(k-1,1) - vDendritePositionsWorking(k,1);
                vY = vDendritePositionsWorking(k-1,2) - vDendritePositionsWorking(k,2);
                vZ = vDendritePositionsWorking(k-1,3) - vDendritePositionsWorking(k,3);
                vDendriteDistance(vBranchIndex) = vDendriteDistance(vBranchIndex)+...
                    (sqrt(vX.^2 + vY.^2 + vZ.^2));
            end
            DendriteCount=DendriteCount+1;
            clear SpotsT vDendriteEdgesWorking vEdgesUnique vDendritePositionsWorking vDendriteRadiusWorking
            clear vDendritePositionsNEW vX vY vZ
            continue
        end
        
        for k=1:numel(vDendritePositionsWorking(:,1))
            if k==1 %for first point
                vDendritePositionsNEW(k,:)=vFilamentsXYZ(vDendriteStartEndIndex(1),:);
                vDendriteRadiusWorking(k,:)=vFilamentsRadius(vDendriteStartEndIndex(1),:);
                test=vDendriteEdgesWorking(:,1)==vDendriteStartEndIndex(1);
                if test==false;
                    test=vDendriteEdgesWorking(:,2)==vDendriteStartEndIndex(1);
                end
                %Find next Filament Index
                vDendriteNextFilamentIndex=vDendriteEdgesWorking(test,:);
                vDendriteNextFilamentIndex(vDendriteNextFilamentIndex==vDendriteStartEndIndex(1))=[];
                %remove current edge row
                vDendriteEdgesWorking(test,:)=[];
                test(test,:)=[];               
            else % each additonal point
                vDendritePositionsNEW(k,:)=vFilamentsXYZ(vDendriteNextFilamentIndex,:);
                vDendriteRadiusWorking(k,:)=vFilamentsRadius(vDendriteNextFilamentIndex,:);
                test=vDendriteEdgesWorking(:,1)==vDendriteNextFilamentIndex;
                if test==false;
                    test=vDendriteEdgesWorking(:,2)==vDendriteNextFilamentIndex;
                end
                vDendriteLastFilamentIndex=vDendriteNextFilamentIndex;
                %Find next Filament Index
                vDendriteNextFilamentIndex=vDendriteEdgesWorking(test,:);
                vDendriteNextFilamentIndex(vDendriteNextFilamentIndex==vDendriteLastFilamentIndex)=[];
                %remove current edge row
                vDendriteEdgesWorking(test,:)=[];
                test(test,:)=[];
                %measure length of dendrite
                vX = vDendritePositionsNEW(k-1,1) - vDendritePositionsNEW(k,1);
                vY = vDendritePositionsNEW(k-1,2) - vDendritePositionsNEW(k,2);
                vZ = vDendritePositionsNEW(k-1,3) - vDendritePositionsNEW(k,3);
                vDendriteDistance(vBranchIndex) = vDendriteDistance(vBranchIndex)+...
                    (sqrt(vX.^2 + vY.^2 + vZ.^2));
            end
        end
        %Collate all new reordered positions into new Filament XYZ list
        vDendritesAll=[vDendritesAll;vDendritePositionsNEW];
        vDendritesAllRadius=[vDendritesAllRadius;vDendriteRadiusWorking];
        vDendriteSegmentID=[vDendriteSegmentID;repmat(DendriteCount,numel(vDendritePositionsNEW(:,1)),1)];
        DendriteCount=DendriteCount+1;
        clear SpotsT vDendriteEdgesWorking vEdgesUnique vDendritePositionsWorking vDendriteRadiusWorking
        clear vDendritePositionsNEW vX vY vZ vDendriteRadiusWorking vDendriteDistanceworking
        
    end
%end
%%
%1=(A) nodes have two bifurcating children. 
%2=Continuation (C) nodes have one bifurcating and one terminating child. 
%3=Termination (T) nodes have two terminating children.
vBranchChild=[];
c=1;
test=0;
% test all dendrite branches to set length of the row
for b=1:numel(vFilamentBranchPoints(:,1))
    DendritesPerNode=unique(vDendriteSegmentID(ismember(vDendritesAll,...
        vFilamentBranchPoints(b,:), 'rows')));
    if DendritesPerNode>test
        rowlength=size(DendritesPerNode,1);
    end
end
%%loop each node and encode each node
for vBranchPointIndex=1:numel(vFilamentBranchPoints(:,1))
    %set vBranchChild size with zeros
    vBranchChild=zeros(size(vFilamentBranchPoints,1),rowlength);
    %Identify the dendrite SegmentID that come off this node branch 
    DendritesPerNode=unique(vDendriteSegmentID(ismember(vDendritesAll,...
        vFilamentBranchPoints(vBranchPointIndex,:), 'rows')));
    
    if size(DendritesPerNode,1)==rowlength
        vBranchChild(vBranchPointIndex,:)=DendritesPerNode';
    elseif size(DendritesPerNode,1)==rowlength-1
        vBranchChild(vBranchPointIndex,:)=[DendritesPerNode' 0];
    elseif size(DendritesPerNode,1)==rowlength-2
        vBranchChild(vBranchPointIndex,:)=[DendritesPerNode' 0 0];
    elseif size(DendritesPerNode,1)==rowlength-3
        vBranchChild(vBranchPointIndex,:)=[DendritesPerNode' 0 0 0];
    elseif size(DendritesPerNode,1)==rowlength-4
        vBranchChild(vBranchPointIndex,:)=[DendritesPerNode' 0 0 0 0];
    elseif size(DendritesPerNode,1)==rowlength-5
        vBranchChild(vBranchPointIndex,:)=[DendritesPerNode' 0 0 0 0 0];
    end
    
    %Identify the segmentIDs for dendrites connected to branch point
    workingSegID=vDendriteSegmentID(ismember...
        (vDendriteSegmentID,vBranchChild(vBranchPointIndex,:)),:);
    %Quantify the number Terminal branches comming of the current branches
    vTerminalChild=unique(workingSegID(ismember(vDendritesAll(ismember...
        (vDendriteSegmentID,vBranchChild(vBranchPointIndex,:)),:),vFilamentTerminalPoints,'rows')));
    
    %Set the classification of the Bifurcation Nodes
    if numel(vTerminalChild) == 0
        vNodeClassify(vBranchPointIndex)=1;%A
    elseif numel(vTerminalChild) == 1
        vNodeClassify(vBranchPointIndex)=2;%C
    elseif numel(vTerminalChild) == 2
        vNodeClassify(vBranchPointIndex)=3;%T
    end
    clear workingSegID vTerminalChild DendritesPerNode
end

%Quantify the number of each classifcation of Node per Filament
vNumberAborizationNodes=numel(vNodeClassify(vNodeClassify==1));
vNumberContinuationNodes=numel(vNodeClassify(vNodeClassify==2));
vNumberTerminationNodes=numel(vNodeClassify(vNodeClassify==3));

%Generate Spots object to visualize the different bifurcation nodes.
CreateSpots=true;
%Create a new folder object for new Sholl Spot intersections
vNew_Spots = vImarisApplication.GetFactory;
result2 = vNew_Spots.CreateDataContainer;
result2.SetName(sprintf('Bifurcation Nodes'));
vNewSpots = vImarisApplication.GetFactory.CreateSpots;
vNodeTime= zeros(1,numel(vNodeClassify(vNodeClassify==1)));
vNodeSpotRadius=zeros(1,numel(vNodeClassify(vNodeClassify==1)))+.5;


vNewSpots.Set(vFilamentBranchPoints(vNodeClassify==1,:), vNodeTime,vNodeSpotRadius);
vNewSpots.SetName(sprintf('Arborization Nodes'));
vRGBA=[0 255 0 0];
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); 
vNewSpots.SetColorRGBA(vRGBA);%Set Green Color

%Add new spots to Surpass Scene
result2.AddChild(vNewSpots, -1);

clear vNodeTime vNodeSpotRadius
vNewSpots = vImarisApplication.GetFactory.CreateSpots;
vNodeTime= zeros(1,numel(vNodeClassify(vNodeClassify==2)));
vNodeSpotRadius=zeros(1,numel(vNodeClassify(vNodeClassify==2)))+.5;
vNewSpots.Set(vFilamentBranchPoints(vNodeClassify==2,:), vNodeTime,vNodeSpotRadius);
vNewSpots.SetName(sprintf('Continuation Nodes'));
vRGBA=[255 0 255 0];
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); 
vNewSpots.SetColorRGBA(vRGBA);%Set Blue Color
%Add new spots to Surpass Scene
result2.AddChild(vNewSpots, -1);
clear vNodeTime vNodeSpotRadius
vNewSpots = vImarisApplication.GetFactory.CreateSpots;
vNodeTime= zeros(1,numel(vNodeClassify(vNodeClassify==3)));
vNodeSpotRadius=zeros(1,numel(vNodeClassify(vNodeClassify==3)))+.5;
vNewSpots.Set(vFilamentBranchPoints(vNodeClassify==3,:), vNodeTime,vNodeSpotRadius);
vNewSpots.SetName(sprintf('Termination Nodes'));
vRGBA=[255 0 0 0];
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); 
vNewSpots.SetColorRGBA(vRGBA);%Set Red Color
%Add new spots to Surpass Scene
result2.AddChild(vNewSpots, -1);
%Add to Imaris
vImarisApplication.GetSurpassScene.AddChild(result2, -1);
clear DendritesPerNode vNodeClassify

end








function aData = CollectFilamentsData(aFilamentsXYZ, aFilamentsRadius, ...
    aFilamentsEdges, aFilamentsSpines, aRoot)

[vIsDendriteBranch, vIsDendriteEndPoint] = GetPoints(aFilamentsEdges, ~aFilamentsSpines);
[vIsSpineBranch, vIsSpineEndPoint] = GetPoints(aFilamentsEdges, aFilamentsSpines);

if aRoot > 0
    vIsDendriteBranch(aRoot) = false;
    vIsDendriteEndPoint(aRoot) = false;
else
    aRoot = [];
end

vIsSpineRoot = GetMixedPoints(aFilamentsEdges, ~aFilamentsSpines);
vIsDendriteEndPoint(vIsSpineRoot) = false;
for vIndex = find(vIsSpineEndPoint)
    vCount = numel(find(aFilamentsEdges == vIndex, 2));
    if vCount > 1
        vIsSpineEndPoint(vIndex) = false;
    end
end

aData = struct( ...
    'BeginningPoint', ...
    [aFilamentsXYZ(aRoot, :), aFilamentsRadius(aRoot)], ...
    'DendriteBranch', ...
    [aFilamentsXYZ(vIsDendriteBranch, :), aFilamentsRadius(vIsDendriteBranch)], ...
    'DendriteTerminal', ...
    [aFilamentsXYZ(vIsDendriteEndPoint, :), aFilamentsRadius(vIsDendriteEndPoint)], ...
    'AttachmentPoint', ...
    [aFilamentsXYZ(vIsSpineRoot, :), aFilamentsRadius(vIsSpineRoot)], ...
    'SpineBranch', ...
    [aFilamentsXYZ(vIsSpineBranch, :), aFilamentsRadius(vIsSpineBranch)], ...
    'SpineTerminal', ...
    [aFilamentsXYZ(vIsSpineEndPoint, :), aFilamentsRadius(vIsSpineEndPoint)]);

function [aIsBranch, aIsEndPoint] = GetPoints(aEdges, aIsType)

vPointsSize = numel(aIsType);
vTypeIndices = find(aIsType);
vSize = numel(vTypeIndices);

vIsEdgeLeft = aIsType(aEdges(:, 1));
vIsEdgeRight = aIsType(aEdges(:, 2));
vIsEdge = vIsEdgeLeft & vIsEdgeRight;
vEdges = aEdges(vIsEdge, :);

aIsBranch = false(1, vPointsSize);
aIsEndPoint = false(1, vPointsSize);
for vIndex = 1:vSize
    vTypeIndex = vTypeIndices(vIndex);
    vCount = numel(find(vEdges == vTypeIndex, 3));
    if vCount == 1
        aIsEndPoint(vTypeIndex) = true;
    elseif vCount > 2
        aIsBranch(vTypeIndex) = true;
    end
end

function aMixedPoints = GetMixedPoints(aEdges, aIsType)

vIsEdgeLeft = aIsType(aEdges(:, 1));
vIsEdgeRight = aIsType(aEdges(:, 2));
vEdgesLeft = vIsEdgeLeft & ~vIsEdgeRight;
vEdgesRight = ~vIsEdgeLeft & vIsEdgeRight;
aMixedPoints = [aEdges(vEdgesLeft, 1)', aEdges(vEdgesRight, 2)'];


