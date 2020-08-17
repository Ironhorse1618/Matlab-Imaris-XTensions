%
%
%
%  Copyright Bitplane BPI 2018
%  Matthew Gastinger, Ph.D
%  Application Support Scientist
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
%        <Item name="Synaptic Bouton Finder" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_FilamentBoutonDetection(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpFilaments">
%          <Item name="Synaptic Bouton Finder" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_FilamentBoutonDetection(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description:




function XT_MJG_FilamentBoutonDetection(aImarisApplicationID)

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
% qQuestion = {'Bouton Detection Sensitivity:'};
% qTitle='Bouton Detection';
% qDefaults = {'0.3'};
% vAnswer = inputdlg(qQuestion,qTitle,1,qDefaults);
% if isempty(vAnswer), return, end
% vPeakDetThreshold = str2double(vAnswer{1});
vPeakDetThreshold = 0.3;


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


%%
% create a structure to hold the data of one filament
% each element of the structure is a list of XYZR
vData = CollectFilamentsData(zeros(0, 3), [], zeros(0, 2), [], 0);
vFieldNames = fieldnames(vData);
vFieldsCount = numel(vFieldNames);

% create an array of structures (one element per time point)
vNumberOfFilaments = vFilaments.GetNumberOfFilaments;
vTimes = 1:vNumberOfFilaments;
for vIndex = 1:vNumberOfFilaments
    vTimes(vIndex) = vFilaments.GetTimeIndex(vIndex-1)+1;
end
vTimesAll=vTimes;
vSizeT = max(vTimes);
vData = vData(ones(vSizeT, 1));
%%
vNumberOfChannels=vImarisApplication.GetDataSet.GetSizeC;
vRGBA = vFilaments.GetColorRGBA;
vNewFilament=vImarisApplication.GetFactory.CreateFilaments;
vNewSpotsBouton = vImarisApplication.GetFactory.CreateSpots;
%Generate Spots-shaped filament to measure intensity and depth code
vNewSpots = vImarisApplication.GetFactory.CreateSpots;
vNewSpotsDendrite = vImarisApplication.GetFactory.CreateSpots;
vNewSpotsSpine = vImarisApplication.GetFactory.CreateSpots;

%Create a new folder object for new Spots
result3 = vFactory.CreateDataContainer;
result3.SetName('New Filament Analysis');
result1 = vFactory.CreateDataContainer;
result1.SetName('Filament displayed as Spots');

TotalBoutons=[];
TotalSegmentIndex=[];
TotalSpotsDendrite=[];
TotalSpotsDendriteTime=[];
TotalSpotsDendriteRadius=[];
FilamentCount=0;
for FilamentIndex=0:vNumberOfFilaments-1
    
    vFilamentsIndexT = vFilaments.GetTimeIndex(FilamentIndex);
    vFilamentsXYZ = vFilaments.GetPositionsXYZ(FilamentIndex);
    vFilamentsRadius = vFilaments.GetRadii(FilamentIndex);
    vFilamentsEdgesSegmentId = vFilaments.GetEdgesSegmentId(FilamentIndex);
    vSegmentIds=unique(vFilamentsEdgesSegmentId);%Idenitfy unique filament segmentIDs
    vFilamentsEdges = vFilaments.GetEdges(FilamentIndex) + 1;
    vTypes = vFilaments.GetTypes(FilamentIndex);
    vBeginningVertex = vFilaments.GetBeginningVertexIndex(FilamentIndex) + 1;
    vFilamentStartingPoint=vFilamentsXYZ(vBeginningVertex,:);
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
    %Get Filament length stats
    if FilamentIndex==0
        vAllStatisticsFilaments = vFilaments.GetStatistics;
        vFilamentStatNames = cell(vAllStatisticsFilaments.mNames);
        vFilamentStatValues = vAllStatisticsFilaments.mValues;
        vFilamentStatIds = vAllStatisticsFilaments.mIds;
        %Find dendrite lengths all filament from scene object
        vDendriteLengthIndex=strmatch('Dendrite Length', vFilamentStatNames);
        vDendriteLength = vFilamentStatValues(vDendriteLengthIndex,:);
        vDendriteLengthIDs = vFilamentStatIds(vDendriteLengthIndex,:);
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
    vDendriteSegmentID=[];
    vDendritesAll=[];
    vDendritesTypesAll=[];
    vDendritesRadiusAll=[];
    vBoutonPositionXYZFinal=[];
    vBoutonRadiusFinal=[];
    vBoutonTimeFinal=[];
    vBoutonNumberFinal=[];
    vBoutonPositionAll=[];
    vBoutonTimeAll=[];
    vBoutonRadiusAll=[];
    wAllSegmentIndex=[];
    DendriteCount=1;
    for vBranchIndex=1:(vNumberOfDendriteBranches)
        %Set the ID for dendrite segment
        wSegmentIndex = vSegmentIds(vBranchIndex);
        %Logical arguement to identify spots in segment
        vSpotsT = vFilamentsEdgesSegmentId == wSegmentIndex;
        vSpotsIndex = find(vSpotsT');
        %Identify position for dendrite segment
        %Test with new method of filtering
        vDendriteEdgesWorking=vFilamentsEdges(vSpotsIndex,:);
        vEdgesUnique=unique(vDendriteEdgesWorking);
        vDendritePositionsNEW=vFilamentsXYZ(vEdgesUnique,1:3);
        vDendriteRadiusNEW=vFilamentsRadius(vEdgesUnique,1);
        vTypesWorking=vTypes(vEdgesUnique,:);
        
        %%
        %find single indices, identifying first and last points of the
        %segment
        vDendriteStartEndIndex=find(accumarray(vDendriteEdgesWorking(:),1)==1);
        
        %%ReOrdering mismatched filament points from filament editing
        %Test working edges to see if reordering is necessary
        %Calculate if sequential and set as true if one or more elements
        %are not sequential it will be false.
        if all(diff(vDendriteEdgesWorking)==1)==true
            vDendritesAll=[vDendritesAll;vDendritePositionsNEW];
            vDendritesRadiusAll=[vDendritesRadiusAll;vDendriteRadiusNEW];
            vDendritesTypesAll=[vDendritesTypesAll;vTypesWorking];            
            vDendriteSegmentID=[vDendriteSegmentID;repmat(DendriteCount,numel(vDendritePositionsNEW(:,1)),1)];
            wAllSegmentIndex=[wAllSegmentIndex;wSegmentIndex];
            DendriteCount=DendriteCount+1;
        elseif isempty (vDendriteStartEndIndex)|numel(vDendriteStartEndIndex)==1
            %Odd situation when there are not single indices that identify the
            %beginning or end.
            vDendritesAll=[vDendritesAll;vDendritePositionsNEW];
            vDendritesRadiusAll=[vDendritesRadiusAll;vDendriteRadiusNEW];
            vDendritesTypesAll=[vDendritesTypesAll;vTypesWorking];
            vDendriteSegmentID=[vDendriteSegmentID;repmat(DendriteCount,numel(vDendritePositionsNEW(:,1)),1)];
            wAllSegmentIndex=[wAllSegmentIndex;wSegmentIndex];
            DendriteCount=DendriteCount+1;
        else
            for k=1:numel(vDendritePositionsNEW(:,1))
                if k==1 %for first point
                    vDendritePositionsNEW(k,:)=vFilamentsXYZ(vDendriteStartEndIndex(1),:);
                    vDendriteRadiusNEW(k,:)=vFilamentsRadius(vDendriteStartEndIndex(1),:);
                    vDendriteTypesNEW(k,:)=vTypes(vDendriteStartEndIndex(1),:);
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
                    vDendriteRadiusNEW(k,:)=vFilamentsRadius(vDendriteNextFilamentIndex,:);
                    vDendriteTypesNEW(k,:)=vTypes(vDendriteNextFilamentIndex,:);
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
            vDendritesRadiusAll=[vDendritesRadiusAll;vDendriteRadiusNEW];
            vDendritesTypesAll=[vDendritesTypesAll;vDendriteTypesNEW];
            vDendriteSegmentID=[vDendriteSegmentID;repmat(DendriteCount,numel(vDendritePositionsNEW(:,1)),1)];
            wAllSegmentIndex=[wAllSegmentIndex;wSegmentIndex];
            DendriteCount=DendriteCount+1;
        end;
%%
        %Quantify Boutons
        %figure; plot(vDendriteRadiusWorking);
        %plot(mintab(:,1), mintab(:,2), 'g*');
        %hold on; plot(maxtab(:,1), maxtab(:,2), 'r*');
        if vTypesWorking==0
            %Identify the peaks in dendrite radius compared to surrounding radii
            [maxtab, mintab] = peakdet(vDendriteRadiusNEW, vPeakDetThreshold);
            if  isempty(maxtab)
                vBoutonNumberFinal=[vBoutonNumberFinal;0];
            else
                vBoutonIndex= maxtab(:,1);
                vBoutonPositionXYZ=vDendritePositionsNEW(vBoutonIndex,:);
                vBoutonRadius=vDendriteRadiusNEW(vBoutonIndex,:);
                vBoutonNumber=size(vBoutonRadius,1);
                vBoutonTime=zeros(size(vBoutonRadius))+vFilamentsIndexT;
                vBoutonPositionXYZFinal=[vBoutonPositionXYZFinal;vBoutonPositionXYZ];
                vBoutonRadiusFinal=[vBoutonRadiusFinal;vBoutonRadius];
                vBoutonTimeFinal=[vBoutonTimeFinal;vBoutonTime];
                vBoutonNumberFinal=[vBoutonNumberFinal;vBoutonNumber];
            end
        end
        clear SpotsT vDendriteEdgesWorking vEdgesUnique vDendritePositionsWorking vDendriteRadiusWorking
        clear vDendritePositionsNEW vDendriteRadiusNEW vX vY vZ vDendritePositionsFinal
        clear vBoutonPositionXYZ vBoutonRadius vBoutonNumber vBoutonTime 
    end
    clear vSegmentIds
    clear vFilamentsXYZ vFilamentsRadius vFilamentsEdgesSegmentId
    %Add Bouton count per dendrite segment
    vNewName2= sprintf(' Dendrite Bouton Number');
    vNewName3= sprintf(' Dendrite Bouton Density');
    
    %Concatonate all boutons from all filaments
    TotalBoutons=[TotalBoutons;vBoutonNumberFinal];
    TotalSegmentIndex=[TotalSegmentIndex;wAllSegmentIndex];
    %Add Spot for each bouton
    vBoutonPositionAll=[vBoutonPositionAll;vBoutonPositionXYZFinal];
    vBoutonTimeAll=[vBoutonTimeAll;vBoutonTimeFinal];
    vBoutonRadiusAll=[vBoutonRadiusAll;vBoutonRadiusFinal];
    %Add the Bouton Spots
    vNewSpotsBouton.Set(vBoutonPositionAll,vBoutonTimeAll,vBoutonRadiusAll);
    vNewSpotsBouton.SetColorRGBA((rand(1, 1)) * 256 * 256 * 256 );
    BoutonSpotName=sprintf(' Boutons on %s ID(%s)', char(vFilaments.GetName),num2str(FilamentIndex))
    vNewSpotsBouton.SetName([BoutonSpotName]);
    result3.AddChild(vNewSpotsBouton, -1);
    clear vNewSpotsBouton
    vNewSpotsBouton = vImarisApplication.GetFactory.CreateSpots;
    if vNumberOfFilaments==FilamentIndex+1
        vInd=1:size(TotalSegmentIndex,1);
        vUnits(vInd) = { char(vImarisApplication.GetDataSet.GetUnit) };
        vFactors2(vInd) = {'Dendrite'};
        vFactorNames = {'Category'};
        TotalBoutons=single(TotalBoutons);
        vNames2(vInd) = {vNewName2};
        vFilaments.AddStatistics(vNames2, TotalBoutons, vUnits, vFactors2, vFactorNames, TotalSegmentIndex);
        vNames2(vInd) = {vNewName3};
        vBoutonDensity=TotalBoutons./vDendriteLength;
        vFilaments.AddStatistics(vNames2, vBoutonDensity, vUnits, vFactors2, vFactorNames, TotalSegmentIndex);
        vBoutonTimeAll=vBoutonTimeAll'; 
    end
end
vImarisApplication.GetSurpassScene.AddChild(result3, -1);




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

function D=bwdistsc(bw,aspect)
% D=BWDISTSC(BW,ASPECT)
% BWDISTSC computes Euclidean distance transform of a binary 3D image BW.
% Distance transform assigns to each pixel in BW a number that is the
% distance from that pixel to the nearest nonzero pixel in BW. BWDISTSC
% can accept a regular 2D image, a 3D array, and a cell array of 2D slices.
% ASPECT is a 3-component vector defining the voxel-aspect-ratio for BW.
% If ASPECT is not given, [1 1 1] isotropic aspect ratio is assumed.
%
% BWDISTSC uses fast optimized scan algorithm and cell-arrays to
% represent internal data, and is less demanding to physical memory as
% well as in many cases up to 10 times faster than MATLAB's native bwdist.
%
% Example:
% bw=zeros(100,100,100);
% bw(40:60,40:60,40:60)=1;
% tic;D=bwdist(bw);toc
% tic;D=bwdistsc(bw);toc
%
% BWDISTSC tries to use MATLAB bwdist from image processing toolbox for 2D
% scans if possible, which is faster, otherwise BWDISTSC will use its own
% algorithm to also perform 2D scans. Own algorithm is also used if x- and
% y-anisotropy scales are not equal; therefore, if your data has only one
% axis that is anisotropic, it is always advantageous to feed it to
% BWDISTSC so that the anisotropic axis is z.
%
%(c) Yuriy Mishchenko HHMI JFRC Chklovskii Lab JUL 2007
% Updated Yuriy Mishchenko (Toros University) SEP 2013

% This implementation uses optimized forward-backward scan version of the
% algorithm of the original bwdistsc (2007), which substantially improves
% its speed and simplifies the code. The improvement is described in the
% part on the selection initial point in the SIVP paper below. The original
% implementation is still used in bwdistsc1, since forward-backward scan
% does not allow limiting computation to a fixed distance value MAXVAL.

% This code is free for use or modifications, just please give credit
% where appropriate. And if you modify code or fix bugs, please drop
% me a message at gmyuriy@hotmail.com.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan algorithms below use following Lema:                     %
% LEMA: let F(X,z) be lower envelope of a family of parabola:   %
% F(X,z)=min_{i} [G_i(X)+(z-k_i)^2];                            %
% and let H_k(X,z)=A(X)+(z-k)^2 be a parabola.                  %
% Then for H_k(X,z)==F(X,z) at each X there exist at most       %
% two solutions k1<k2 such that H_k12(X,z)=F(X,z), and          %
% H_k(X,z)<F(X,z) is restricted to at most k1<k2.               %
% Here X is any-dimensional coordinate.                         %
%                                                               %
% Thus, simply scan away from any z such that H_k(X,z)<F(X,z)   %
% in either direction as long as H_k(X,z)<F(X,z) and update     %
% F(X,z). Note that need to properly choose starting point;     %
% starting point is any z such that H_k(X,z)<F(X,z); z==k is    %
% usually, but not always the starting point!!!                 %
% usually, but not always the starting point!                   %
%                                                               %
% Citation:                                                     %
% Mishchenko Y. (2013) A function for fastcomputation of large  %
% discrete Euclidean distance transforms in three or more       %
% dimensions in Matlab. Signal, Image and Video Processing      %
% DOI: 10.1007/s11760-012-0419-9.                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse inputs
if(nargin<2 || isempty(aspect)) aspect=[1 1 1]; end

% determine geometry of the data
if(iscell(bw)) shape=[size(bw{1}),length(bw)]; else shape=size(bw); end

% correct this for 2D data
if(length(shape)==2) shape=[shape,1]; end
if(length(aspect)==2) aspect=[aspect,1]; end

% allocate internal memory
D=cell(1,shape(3)); for k=1:shape(3) D{k}=zeros(shape(1:2)); end

vProgressDisplay = waitbar(0, 'Distance Transform.....Step1');
%%%%%%%%%%%%% scan along XY %%%%%%%%%%%%%%%%
for k=1:shape(3)
    if(iscell(bw)) bwXY=bw{k}; else bwXY=bw(:,:,k); end
    
    % initialize arrays
    DXY=zeros(shape(1:2));
    D1=zeros(shape(1:2));
    
    % if can, use 2D bwdist from image processing toolbox
    if(exist('bwdist') && aspect(1)==aspect(2))
        D1=aspect(1)^2*bwdist(bwXY).^2;
    else    % if not, use full XY-scan
        %%%%%%%%%%%%%%% X-SCAN %%%%%%%%%%%%%%%
        % reference for nearest "on"-pixel in bw in x direction down
        
        %  scan bottow-up (for all y), copy x-reference from previous row
        %  unless there is "on"-pixel in that point in current row, then
        %  that's the nearest pixel now
        xlower=repmat(Inf,shape(1:2));
        
        xlower(1,find(bwXY(1,:)))=1;    % fill in first row
        for i=2:shape(1)
            xlower(i,:)=xlower(i-1,:);  % copy previous row
            xlower(i,find(bwXY(i,:)))=i;% unless there is pixel
        end
        
        % reference for nearest "on"-pixel in bw in x direction up
        xupper=repmat(Inf,shape(1:2));
        
        xupper(end,find(bwXY(end,:)))=shape(1);
        for i=shape(1)-1:-1:1
            xupper(i,:)=xupper(i+1,:);
            xupper(i,find(bwXY(i,:)))=i;
        end
        
        % build (X,Y) for points for which distance needs to be calculated
        idx=find(~bwXY); [x,y]=ind2sub(shape(1:2),idx);
        
        % update distances as shortest to "on" pixels up/down in the above
        DXY(idx)=aspect(1)^2*min((x-xlower(idx)).^2,(x-xupper(idx)).^2);
        
        %%%%%%%%%%%%%%% Y-SCAN %%%%%%%%%%%%%%%
        % this will be the envelop of parabolas at different y
        D1=repmat(Inf,shape(1:2));
        
        p=shape(2);
        for i=1:shape(2)
            % some auxiliary datasets
            d0=DXY(:,i);
            
            % selecting starting point for x:
            % * if parabolas are incremented in increasing order of y,
            %   then all below-envelop intervals are necessarily right-
            %   open, which means starting point can always be chosen
            %   at the right end of y-axis
            % * if starting point exists it should be below existing
            %   current envelop at the right end of y-axis
            dtmp=d0+aspect(2)^2*(p-i)^2;
            L=D1(:,p)>dtmp;
            idx=find(L);
            D1(idx,p)=dtmp(L);
            
            
            % these will keep track along which X should
            % keep updating distances
            map_lower=L;
            idx_lower=idx;
            
            % scan from starting points down in increments of 1
            for ii=p-1:-1:1
                % new values for D
                dtmp=d0(idx_lower)+aspect(2)^2*(ii-i)^2;
                
                % these pixels are to be updated
                L=D1(idx_lower,ii)>dtmp;
                D1(idx_lower(L),ii)=dtmp(L);
                
                % other pixels are removed from scan
                map_lower(idx_lower)=L;
                idx_lower=idx_lower(L);
                
                if(isempty(idx_lower)) break; end
            end
        end
    end
    D{k}=D1;
    waitbar(k/shape(3), vProgressDisplay);
    
end
close(vProgressDisplay);
vProgressDisplay = waitbar(0, 'Distance Transform.....Step2');

%%%%%%%%%%%%% scan along Z %%%%%%%%%%%%%%%%
D1=cell(size(D));
for k=1:shape(3)
    D1{k}=repmat(Inf,shape(1:2));
    waitbar(k/shape(3), vProgressDisplay);
end
close(vProgressDisplay);
vProgressDisplay = waitbar(0, 'Distance Transform.....Step3');
% start building the envelope
p=shape(3);
for k=1:shape(3)
    % if there are no objects in this slice, nothing to do
    if(isinf(D{k}(1,1)))
        waitbar(k/shape(3), vProgressDisplay);
        continue;
    end
    
    % selecting starting point for (x,y):
    % * if parabolas are incremented in increasing order of k, then all
    %   intersections are necessarily at the right end of the envelop,
    %   and so the starting point can be always chosen as the right end
    %   of the axis
    
    % check which points are valid starting points, & update the envelop
    dtmp=D{k}+aspect(3)^2*(p-k)^2;
    L=D1{p}>dtmp;
    D1{p}(L)=dtmp(L);
    
    % map_lower keeps track of which pixels can be yet updated with the
    % new distance, i.e. all such XY that had been under the envelop for
    % all Deltak up to now, for Deltak<0
    map_lower=L;
    
    % these are maintained to keep fast track of whether map is empty
    idx_lower=find(map_lower);
    
    % scan away from the starting points in increments of -1
    for kk=p-1:-1:1
        % new values for D
        dtmp=D{k}(idx_lower)+aspect(3)^2*(kk-k)^2;
        
        % these pixels are to be updated
        L=D1{kk}(idx_lower)>dtmp;
        map_lower(idx_lower)=L;
        D1{kk}(idx_lower(L))=dtmp(L);
        
        % other pixels are removed from scan
        idx_lower=idx_lower(L);
        
        if(isempty(idx_lower)) break; end
    end
    waitbar(k/shape(3), vProgressDisplay);
end

close(vProgressDisplay);

% prepare the answer
if(iscell(bw))
    D=cell(size(bw));
    for k=1:shape(3) D{k}=sqrt(D1{k}); end
else
    D=zeros(shape);
    for k=1:shape(3) D(:,:,k)=sqrt(D1{k}); end
end

