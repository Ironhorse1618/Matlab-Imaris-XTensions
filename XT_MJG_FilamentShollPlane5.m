%  Copyright Bitplane BPI 2018
%  Matthew Gastinger, Ph.D
%  Application Support Scientist
%  
%  Modifications:
%  Oran Avivi, Ba.Sc
%  Physics & Electrical Engineer
%  University of Toronto
%  Department of Neuroscience
%
%  Installation:
%
%  - Extract the .m file.   Place into the current active XTensions folder
%    for the current version of Imaris.
%  - This function can be found in the Image Processing menu, as well as in
%       the XTensions (Gear Tab) for the Filament Surpass scene object.
%
%    <CustomTools>
%      <Menu>
%       <Submenu name="Filaments Functions">
%        <Item name="Display Sholl-plane Intersections" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_FilamentShollPlane5(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpFilaments">
%          <Item name="Display Sholl-plane Intersections" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_FilamentShollPlane5(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%
%Description:  This XTension finds and displays the Sholl intersections
%for each Filament that have been traced and contains a StartingPoint.  It
%will be functional for multiple filaments.
%

%For each sholl sphere interval that is defined by the user, the script 
%identifies the dendrite segments that cross that interval at least once.   
%All of the filament points in the segment are measured to the filament
%starting point.  Each intersection is identified and
%classified/reduced to a single position.  Number of Sholl intersections 
%will be displayed as a Spots object for each Sholl pre-defined interval.

%A single Group folder is generated, even if multiple filaments are in the
%Filament Object.  They will be named similarly but differentiated but the
%actual Filament Index in Imaris.  This wil facilitate the export of Sholl
%intersections for mulitple filaments

%This is not an exact replica of the Sholl intersections reported by the
%Filament object (built into Imaris).  The values are very similar, but the rules
%for identifying a single intersection differ slightly, especially
%around where a branch point intersects a Sholl sphere. (See Readme.doc of
%an report of comparison)


%% Set-up
function XT_MJG_FilamentShollPlane5(aImarisApplicationID)

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

% get the filament
vFactory = vImarisApplication.GetFactory;
vFilaments = vFactory.ToFilaments(vImarisApplication.GetSurpassSelection);

%Create a new folder object for new Sholl Spot intersections
vNew_Spots = vImarisApplication.GetFactory;
result2 = vNew_Spots.CreateDataContainer;
result2.SetName(['Sholl Intersections - ',char(vFilaments.GetName)]);


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

%% Preparation
% Step 1: Create an oblique plane in Imaris that is perpendicular to the
% dendrite/astrocyte direction of interest.
% Step 2: place 3 random reference points on the surface of the plane.
% Step 3: Run this code

% Check to see the 3 steps were done (mostly that there are 3 spots.

%Find First MeasurementPoint Object in Surpass scene
%And grab the XYZ positions to define plane
for vChildIndex = 1:vSurpassScene.GetNumberOfChildren
    vDataItem = vSurpassScene.GetChild(vChildIndex - 1);
    if vFactory.IsMeasurementPoints(vDataItem)
        vMeasurementPoints = vFactory.ToMeasurementPoints(vDataItem);
        break;
    end  
end

vMeasurementPointsValues = vMeasurementPoints.Get;
vPositionsXYZ = vMeasurementPointsValues.mPositionsXYZ;
if size(vPositionsXYZ,1)~=3 | isempty(vPositionsXYZ)
    msgbox('Please create 3 measurement points on an Oblique Plane. Ensure the oblique plane is perpendicular to the dendrite/astrocyte growth direction of your dataset');
    return;
end

%% Prompt for choosing plane-plane distances
qQuestion = {'Plane-Plane Distance:', 'Number of Planes:'}; % prompt question
qTitle='Sholl Plane Analysis';
qDefaults = {'10','10'};
vAnswer = inputdlg(qQuestion,qTitle,1,qDefaults);
if isempty(vAnswer), return, end
vShollRadius = str2double(vAnswer{1});
vNumberofSpheres = str2double(vAnswer{2});

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

%% Plane (Oblique Slicer)
% Determine the position of the 3 reference points
vPositionsXYZ1 = vPositionsXYZ(1,:); % Position of first point on plane
vPositionsXYZ2 = vPositionsXYZ(2,:); % Position of second point on plane
vPositionsXYZ3 = vPositionsXYZ(3,:); % Position of third point on plane

% Find the normal vector to the plane
vPlaneVector1 = vPositionsXYZ2 - vPositionsXYZ1; % Find first vector along the plane
vPlaneVector2 = vPositionsXYZ3 - vPositionsXYZ1; % Find second vector along the plane
vNormalPlaneVector = cross(vPlaneVector1, vPlaneVector2); % Find normal vector

% For a plane equation of form a*x + b*y + c*z + d = 0, determine the
% parameters
aa = vNormalPlaneVector(1);
bb = vNormalPlaneVector(2);
cc = vNormalPlaneVector(3);
dd = -(aa*vPositionsXYZ1(1) + bb*vPositionsXYZ1(2) + cc*vPositionsXYZ1(3));

%%
% create a structure to hold the data of one filament
% each element of the structure is a list of XYZR
vData = CollectFilamentsData(zeros(0, 3), [], zeros(0, 2), [], 0);
vFieldNames = fieldnames(vData);
vFieldsCount = numel(vFieldNames);

% create an array of structures (one element per time point)
vCount = vFilaments.GetNumberOfFilaments;
vTimes = 1:vCount;
for vIndex = 1:vCount
    vTimes(vIndex) = vFilaments.GetTimeIndex(vIndex-1)+1;
end
vTimesAll=vTimes;
vSizeT = max(vTimes);
vData = vData(ones(vSizeT, 1));
%%+
vNumberOfChannels=vImarisApplication.GetDataSet.GetSizeC;
vRGBA = vFilaments.GetColorRGBA;
vNewFilament=vImarisApplication.GetFactory.CreateFilaments;
vTotalShollSpots=0;
FilamentCount=0;
vDendriteSegmentID=[];
vDendritesAll=[];
DendriteCount=1;
maxDistanceValuesAll=[];
minDistanceValuesAll=[];
vFilamentBranchPoints=[];
vFilamentTerminalPoints=[];

for FilamentIndex=0:vCount-1
    
    vFilamentsIndexT = vFilaments.GetTimeIndex(FilamentIndex);
    vFilamentsXYZ = vFilaments.GetPositionsXYZ(FilamentIndex);
    vFilamentsRadius = vFilaments.GetRadii(FilamentIndex);
    vFilamentsEdgesSegmentId = vFilaments.GetEdgesSegmentId(FilamentIndex);
    vSegmentIds=unique(vFilamentsEdgesSegmentId);%Idenitfy unique filament segmentIDs
    vFilamentsEdges = vFilaments.GetEdges(FilamentIndex) + 1;
    vTypes = vFilaments.GetTypes(FilamentIndex);
    vBeginningVertex = vFilaments.GetBeginningVertexIndex(FilamentIndex) + 1;
    
    %vFilamentStartingPoint=vFilamentsXYZ(vBeginningVertex,:);
    vIndexT = vFilaments.GetTimeIndex(FilamentIndex);
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
    vTimes=vTimesAll;
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
            %vFilamentStartingPoint=vXYZRT(1:3);
            %vFilamentStartingPoint=vXYZRT(2,1:3)
        end
        
        if vField==2
            vFilamentBranchPointsOne=vXYZRT(:,1:3);
            vFilamentBranchPoints=[vFilamentBranchPoints;vFilamentBranchPointsOne];
        end
        
        if vField==3
            vFilamentTerminalPointsOne=vXYZRT(:,1:3);
            vFilamentTerminalPoints=[vFilamentTerminalPoints;vFilamentTerminalPointsOne];
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
%     vDendriteSegmentID=[];
%     vDendritesAll=[];
%     DendriteCount=1;
    
    for vBranchIndex=1:(vNumberOfDendriteBranches)
        %Set the ID for dendrite segment
        wSegmentIndex = vSegmentIds(vBranchIndex);
        %Logical arguement to identify spots in segment
        vSpotsT = vFilamentsEdgesSegmentId == wSegmentIndex;
        vSpotsIndex = find(vSpotsT');
        %Identify position for dendrite segment
        %Test with new method of filtering+
        vDendriteEdgesWorking=vFilamentsEdges(vSpotsIndex,:);
        vEdgesUnique=unique(vDendriteEdgesWorking);
        vDendritePositionsWorking=vFilamentsXYZ(vEdgesUnique,1:3);
        vDendriteRadiusWorking=vFilamentsRadius(vEdgesUnique,1);
        %vTypesWorking=vTypes(vEdgesUnique,:);
        
        %% Point-Plane Distance Calculations
        %Measure distance of all FilamentPoint position to Filament Startingn point
        %To find out which Denrite segment intersects with each sphere, to
        %hopefully reduce the time of process irrelevant segments
        
        % Define Reference Point Co-ordinates
        x0 = vDendritePositionsWorking(:,1); % X Values
        y0 = vDendritePositionsWorking(:,2); % Y Values
        z0 = vDendritePositionsWorking(:,3); % Z Values
        
        % Determine distance from point to plane
        vShortestDistanceToPlane = abs(aa.*x0 + bb.*y0 + cc.*z0+dd)./sqrt(aa.^2 + bb.^2 + cc.^2); % Equation for shortest distance to plane (along the normal vector)
        maxDistanceValues(DendriteCount,1) = max(vShortestDistanceToPlane);
        minDistanceValues(DendriteCount,1) = min(vShortestDistanceToPlane);
        
%         Previous Code:
%         vX = vFilamentStartingPoint(1) - vDendritePositionsWorking(:,1);
%         vY = vFilamentStartingPoint(2) - vDendritePositionsWorking(:,2);
%         vZ = vFilamentStartingPoint(3) - vDendritePositionsWorking(:,3);
%         maxDistanceValues(vBranchIndex,1) = max(sqrt(vX.^2 + vY.^2 + vZ.^2));
%         minDistanceValues(vBranchIndex,1) = min(sqrt(vX.^2 + vY.^2 + vZ.^2));
%         
        %%
        %%ReOrdering mismatched filament points from filament editing
        %Test working edges to see if reordering is necessary
        %Calculate if sequential and set as true if one or more elements
        %are not sequential it will be false.
        if all(diff(vDendriteEdgesWorking)==1)==true
            vDendritesAll=[vDendritesAll;vDendritePositionsWorking];
            vDendriteSegmentID=[vDendriteSegmentID;repmat(DendriteCount,numel(vDendritePositionsWorking(:,1)),1)];
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
            DendriteCount=DendriteCount+1;
            clear SpotsT vDendriteEdgesWorking vEdgesUnique vDendritePositionsWorking vDendriteRadiusWorking
            clear vDendritePositionsNEW
            continue
        end
        
        for k=1:numel(vDendritePositionsWorking(:,1))
            if k==1 %for first point
                vDendritePositionsNEW(k,:)=vFilamentsXYZ(vDendriteStartEndIndex(1),:);
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
            end
        end
        %Collate all new reordered positions into new Filament XYZ list
        vDendritesAll=[vDendritesAll;vDendritePositionsNEW];
        vDendriteSegmentID=[vDendriteSegmentID;repmat(DendriteCount,numel(vDendritePositionsNEW(:,1)),1)];
        DendriteCount=DendriteCount+1;
        clear SpotsT vDendriteEdgesWorking vEdgesUnique vDendritePositionsWorking vDendriteRadiusWorking
%         clear vDendritePositionsNEW vX vY vZ % This becomes redundant
%         after changing the distance calculations
    end    
    
    clear vXYZRT
    clear vFilamentBranchPointsOne vFilamentTerminalPointsOne vFilamentStartingPoint
end    
    
    %%
%     %Create a new folder object for new Sholl Spot intersections
% vNew_Spots = vImarisApplication.GetFactory;
% result2 = vNew_Spots.CreateDataContainer;
% result2.SetName(['Sholl Intersections - ',char(vFilaments.GetName)]);
    
    
    vShollSpotCount=0;
    vShollNodeCount=0;
    aIntensityLowerThresholdManual=0;
    
    vProgressDisplay = waitbar(0, 'Sholl Analysis');
    for ShollIndex=1:vNumberofSpheres
        %Set lower and upper threshold for the Sholl Mask
        aIntensityLowerThresholdManual=aIntensityLowerThresholdManual+vShollRadius;
        aIntensityUpperThresholdManual=aIntensityLowerThresholdManual+(Xvoxelspacing*2);
        SpotinQuestion=[];
        
        %Process at each branch and find point distance to starting point
        %filter by current sholl plane distance
        %Test if Dendrite fails on the Sholl Plane
        DendritesT=maxDistanceValues>=aIntensityLowerThresholdManual & ...
            minDistanceValues<=aIntensityLowerThresholdManual;
        if DendritesT==false; break;end
        for vBranchIndex=1:(max(DendriteCount)-1)
            if DendritesT(vBranchIndex)==false; continue;end%skip segment that do not intersect
            SpotsT=vDendriteSegmentID==vBranchIndex;
            vDendritePositionsNEW=vDendritesAll(SpotsT,:);
            %Set base for each filament point to false
            vColoc1 = false(numel(vDendritePositionsNEW(:,1)), 1);
            %Measure distance from each point to starting for current
            %segment
            %% Point-Plane Distance Calculations
            %Measure distance of all FilamentPoint position to Filament Startingn point
            %To find out which Denrite segment intersects with each plane, to
            %hopefully reduce the time of process irrelevant segments

            % Define Reference Point Co-ordinates
            x0 = vDendritePositionsNEW(:,1); % new X Values
            y0 = vDendritePositionsNEW(:,2); % new Y Values
            z0 = vDendritePositionsNEW(:,3); % new Z Values

            % Determine distance from point to plane
            vShortestDistanceToPlane2 = abs(aa.*x0 + bb.*y0 + cc.*z0+dd)./sqrt(aa.^2 + bb.^2 + cc.^2); % Equation for shortest distance to plane (along the normal vector)
            vDistanceListValues = vShortestDistanceToPlane2;
            
%             Previous Code:
%             vX = vDendritePositionsNEW(:,1) - vFilamentStartingPoint(1);
%             vY = vDendritePositionsNEW(:,2) - vFilamentStartingPoint(2);
%             vZ = vDendritePositionsNEW(:,3) - vFilamentStartingPoint(3);
%             vDistanceListValues = sqrt(vX.^2 + vY.^2 + vZ.^2);
            
            %Test filament points to find those near current sholl Plane
            ShollIntersectionIndex=find(vDistanceListValues>aIntensityLowerThresholdManual &...
                vDistanceListValues<aIntensityUpperThresholdManual);
            %if no filament point is detected, within these defined limits
            %find the closest point and use that as the lone Sholl
            %intersection
            if isempty(ShollIntersectionIndex)
                [d, ShollIntersectionIndex] = min(abs(vDistanceListValues-aIntensityLowerThresholdManual));   
            end
            %Reset filament point to true if it met above thrshol criteria
            if ~isempty(ShollIntersectionIndex)
                vColoc1(ShollIntersectionIndex)=true;
            else
                continue
            end
            %convert to real numbers
            vDendriteXSphere=double(vColoc1);
            %convert to positionsXYZ
            vDendriteXSpherePositionXYZ=vDendritePositionsNEW(vColoc1,:);
            
            %%
            %Find the Sholl intersections
            if ~isempty(vDendriteXSpherePositionXYZ)
                vDendriteXSphere(end+1)=0;%Need to add one zero to the end of series so that a peak can be identified at the tail end
                [maxtab, mintab] = peakdet(vDendriteXSphere, 0.3);%find dendrite points that overlap
                PeakIndex=maxtab(:,1);%identify index of the peaks overlap
                vDendriteXSphere(end)=[];
                %Test for single gaps in the dendrite peak and fill them. Likely not
                %enough to warrrent a true Sholl intersection cross
                PeakIndexWorking=PeakIndex;
                for q=1:numel(PeakIndexWorking)-1
                    if PeakIndexWorking(q)+1==PeakIndexWorking(q+1)-1
                        vDendriteXSphere(PeakIndexWorking(q)+1)=true;
                        PeakIndex(q+1)=[];
                    end
                end              
                vDendriteXSphere(end+1)=0;
                for i=1:numel(PeakIndex)
                    count=1;
                    b=vDendriteXSphere(PeakIndex(i):1:end);%generated series starting from next peakpoint
                    while b(count)~=0
                        count=count+1;
                    end
                    %Test Peaks if dendrite point intersects with first point
                    %of the dendrite mark it as node for later analysis
                    if PeakIndex(i)==1 & count==2 | PeakIndex(i)==numel(vDendritePositionsNEW(:,1)) & count==2
                        vShollSpotPositionXYZ(vShollSpotCount+1,:)=vDendritePositionsNEW(PeakIndex(i),:);
                        SpotinQuestion(vShollSpotCount+1)=true;%logical count of questionable spots
                        NodeinQuestion(vShollSpotCount+1,:)=vDendritePositionsNEW(PeakIndex(i),:);
                        vShollNodeCount=vShollNodeCount+1;
                        if PeakIndex(i)==1
                            vShollSpotNodeClassify(vShollSpotCount+1,1)=1;%1 sets node a begining of dendrite segment
                        elseif any(ismember(vFilamentTerminalPoints,vDendritePositionsNEW(PeakIndex(i),:),'rows'))==true;
                            %test if member of terminal
                            NodeinQuestion(vShollSpotCount+1,:)=[];%remove
                            SpotinQuestion(vShollSpotCount+1)=[];%remove
                            vShollNodeCount=vShollNodeCount-1;%less one node count
                        else
                            vShollSpotNodeClassify(vShollSpotCount+1,1)=3;%2 set nodel at the end of dendrite segment
                        end
                        vShollSpotCount=vShollSpotCount+1;
                    elseif PeakIndex(i)==1 & count>2 | PeakIndex(i)==numel(vDendritePositionsNEW(:,1)) & count>2
                        vShollSpotPositionXYZ(vShollSpotCount+1,:)=mean(vDendritePositionsNEW(PeakIndex(i):PeakIndex(i)+count-2,:,:));
                        SpotinQuestion(vShollSpotCount+1)=true;%logical count of questionable spots
                        NodeinQuestion(vShollSpotCount+1,:)=vDendritePositionsNEW(PeakIndex(i),:);
                        vShollNodeCount=vShollNodeCount+1;
                        if PeakIndex(i)==1
                            vShollSpotNodeClassify(vShollSpotCount+1,1)=2;%3 sets node a begining of dendrite segment (longer than 2 points)
                        elseif any(ismember(vFilamentTerminalPoints,vDendritePositionsNEW(PeakIndex(i),:),'rows'))==true;
                            %test if member of terminal
                            NodeinQuestion(vShollSpotCount+1,:)=[];%remove
                            SpotinQuestion(vShollSpotCount+1)=[];%remove
                            vShollNodeCount=vShollNodeCount-1;%less one node count
                        else
                            vShollSpotNodeClassify(vShollSpotCount+1,1)=4;%4 set node at the end of dendrite segment (longer than 2 points)
                        end
                        vShollSpotCount=vShollSpotCount+1;
                    elseif count==2
                        vShollSpotPositionXYZ(vShollSpotCount+1,:)=vDendritePositionsNEW(PeakIndex(i),:);
                        vShollSpotCount=vShollSpotCount+1;
                    else
                        vShollSpotPositionXYZ(vShollSpotCount+1,:)=mean(vDendritePositionsNEW(PeakIndex(i):PeakIndex(i)+count-2,:,:));
                        vShollSpotCount=vShollSpotCount+1;
                    end
                end
            end
            clear vShortestDistanceToPlane2 vDistanceListValues
        end
        clear vDendriteXSpherePositionXYZ b count maxtab mintab PeakIndex
        clear vDendritePositionsWorking vDendriteXSphere vDendritePositionsNEW
        clear vDendritePositionIndex q PeakIndexWorking 
        %%
        %Process and remove duplicate spots and reduce spots near node.
        RM=1;
        RemoveSpot=[];
        if any(SpotinQuestion)
            vShollSpotNodeClassify(numel(SpotinQuestion))=0;%add zeros to end to fill to numel(SpotinQuestion)
            
            [q w NodeIndex]=unique(NodeinQuestion,'rows');%identify unique rows in the questionable spot sholl positions
            for z=1:max(NodeIndex)
                SpotsT2=NodeIndex==z;
                [Scenario I]=sort(vShollSpotNodeClassify(SpotsT2));
                if Scenario==0
                    continue
                end
                %Identify certain scenarios of node positions
                switch num2str(Scenario)
                    case ['1';'1';'1']
                        SpotsT3=vShollSpotNodeClassify==1;
                        RemoveSpot(RM,:)=find(SpotsT3, 1,'first');
                        RM=RM+1;
                        RemoveSpot(RM,:)=find(SpotsT3, 1,'last');
                        RM=RM+1;
                    case ['2';'2';'3']%Scenario #2
                        SpotsT3=vShollSpotNodeClassify==3;
                        RemoveSpot(RM,:)=find(SpotsT3);
                        RM=RM+1;
                    case ['1';'1';'4']
                        SpotsT3=vShollSpotNodeClassify==1;
                        RemoveSpot(RM,:)=find(SpotsT3, 1,'first');
                        RM=RM+1;
                        RemoveSpot(RM,:)=find(SpotsT3, 1,'last');
                        RM=RM+1;
                    case ['1';'2';'4']
                        SpotsT3=vShollSpotNodeClassify==4;
                        RemoveSpot(RM,:)=find(SpotsT3);
                        RM=RM+1;
                end
                clear Scenario
            end
            %Remove spots in question idendified in section above
            if ~isempty(RemoveSpot)
                vShollSpotPositionXYZ(RemoveSpot,:)=[];
                vShollSpotCount=vShollSpotCount-numel(RemoveSpot);
            end
        end
        
        vShollSpotRadius=repmat(1,vShollSpotCount,1);
        vShollSpotTime=repmat(vFilamentsIndexT,vShollSpotCount,1);
        
        if vShollSpotCount~=0
            %Create the a new Spots generated from teh center of Mass
            vNewSpots = vImarisApplication.GetFactory.CreateSpots;
            vNewSpots.Set(vShollSpotPositionXYZ, vShollSpotTime',vShollSpotRadius');
            vNewSpots.SetName(sprintf('%s um Sholl Plane (ID:%s)',...
                num2str(aIntensityLowerThresholdManual),num2str(FilamentIndex)));
            vNewSpots.SetColorRGBA((rand(1, 1)) * 256 * 256 * 256 );%Set Random color
            %Add new spots to Surpass Scene
            result2.AddChild(vNewSpots, -1);
            vImarisApplication.GetSurpassScene.AddChild(result2, -1);
        else
            waitbar(vNumberofSpheres/vNumberofSpheres, vProgressDisplay);
            break
        end
        
        %end
        vTotalShollSpots=vTotalShollSpots+vShollSpotCount;
        clear IsNodeSinge NodeinQuestion SpotinQuestion
        waitbar(ShollIndex/vNumberofSpheres, vProgressDisplay);
        vShollSpotCount=0;
        clear vShollSpotPositionXYZ vShollSpotTime vShollSpotRadius    
    end
    close(vProgressDisplay);
%     clear DendritesT maxDistanceValues minDistanceValues vXYZRT
%     clear vFilamentBranchPoints vFilamentTerminalPoints vFilamentStartingPoint
% end









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

function [maxtab, mintab]=peakdet(v, delta, x)
%PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
%        maxima and minima ("peaks") in the vector V.
%        MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%
%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
%        in MAXTAB and MINTAB are replaced with the corresponding
%        X-values.
%
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA.

% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.

maxtab = [];
mintab = [];

v = v(:); % Just in case this wasn't a proper vector

if nargin < 3
    x = (1:length(v))';
else
    x = x(:);
    if length(v)~= length(x)
        error('Input vectors v and x must have same length');
    end
end

if (length(delta(:)))>1
    error('Input argument DELTA must be a scalar');
end

if delta <= 0
    error('Input argument DELTA must be positive');
end

mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;

lookformax = 1;

for i=1:length(v)
    this = v(i);
    if this > mx, mx = this; mxpos = x(i); end
    if this < mn, mn = this; mnpos = x(i); end
    
    if lookformax
        if this < mx-delta
            maxtab = [maxtab ; mxpos mx];
            mn = this; mnpos = x(i);
            lookformax = 0;
        end
    else
        if this > mn+delta
            mintab = [mintab ; mnpos mn];
            mx = this; mxpos = x(i);
            lookformax = 1;
        end
    end
end
