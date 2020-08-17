%Kiss and run analysis


%Written by Matthew J. Gastinger, Bitplane Advanced Application Scientist.
%July 2018.  Editted to work in Imaris 9.2 version.  Same functionality,
%same statistics and workflow.


%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Kiss and Run Analysis" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_KissandRun7(%i)</Command>
%        </Item>
%       </Submenu>
%       <Submenu name="Spots Functions">
%        <Item name="Kiss and Run Analysis" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_KissandRun7(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Kiss and Run Analysis" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_KissandRun7(%i)</Command>
%          </Item>
%        </SurpassComponent>
%        <SurpassComponent name="bpSpots">
%          <Item name="Kiss and Run Analysis" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_KissandRun7(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%



%Description
%
%This Xtension calculates the contact events between: 1)Two
%different surface objects or 2) Spots object and a surface object.  It
%is calculated by one of 2 methods:


%Method #1: Using the a Distance Transformation outside of the target
%surface object, the closest surface to surface distance is determined of
%tracked object.  Based on user-selected threshold to define a contact event
%the various new statistics will be generated.  This method is typically
%faster and the  added option to chose proximity of the contact event that
%are not necessarily overlapping.


%Method#2:  Using a surface mask for the target and the tracked surfaces,
%the region of surface overlap is determined for each surface object.  In
%order to be defined as a contact event the surface must have at least one
%overlapping voxel in each mask.  Using this method, the process may be
%slower for large number of surfaces and long time lapses.  However, you
%will have the added statistic that will provide the volume of overlap at
%each timepoint for each surface.



function XT_MJG_KissandRun7(aImarisApplicationID)
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


%Launch the GUI
%Collect of the data for the XTension
%Set some variables to global
global vImarisApplicationGlobal f count
global checkSurfaceOverlap checkCreateColocSurface checkDistanceTransform
global checkDistanceTransformChannel editboxDistanceThreshold
global checkNumberofContacts checkNumberofProlongedContacts checkPercentNumberofContacts
global checkTotalContactTime checkPercentContactTime checkLongestContactEvent checkMeanContactEvent
global checkTotalTrackDurationAway qTotalTrackDurationAway
global qSurfaceOverlap qCreateColocSurface qDistanceTransform qCreateDistanceTransformChannel
global qNumberofContacts qNumberofProlongedContacts qPercentNumberofContacts
global qTotalContactTime qPercentContactTime qLongestContactEvent qMeanContactEvent


vImarisApplicationGlobal = vImarisApplication;
count=0;
% Create figure
f = figure('units','pixels','position',[500,500,385,400],...
    'toolbar','none','menu','none',...
    'numbertitle','off','name','Kiss and Run Analysis');
%Option to merge all filamentSpots into a single Spots object
%checkbox Mean Thickness
text = uicontrol('Parent', f,'style','text',...
    'position',[15,350,300,40],...
    'string','Measurement (choose one)',...
    'fontangle','italic','fontweight','bold','fontSize', 15);
checkSurfaceOverlap = uicontrol('Parent', f,'style','checkbox',...
    'position', [65,332,200,27],...
    'string', 'Surface Overlap (slower)','fontweight','bold','fontSize', 11);
text = uicontrol('Parent', f,'style','text',...
    'position',[83,308,200,27],...
    'string','Contact events defined when the two surface objects overlap in 3D',...
    'fontangle','italic','fontangle','italic','horizontalalignment','left');
checkCreateColocSurface = uicontrol('Parent', f,'style','checkbox',...
    'position', [100,283,200,20], 'string', 'Create ColocSurface');

checkDistanceTransform = uicontrol('Parent', f,'style','checkbox',...
    'position', [65,255,250,25],...
    'string', 'Distance Transform (default)','Value', 1,'fontweight','bold','fontSize', 11);
text = uicontrol('Parent', f,'style','text',...
    'position',[83,230,200,27],...
    'string','Contact events defined by closest surface to surface distance measure',...
    'fontangle','italic','fontangle','italic', 'horizontalalignment','left');
editboxDistanceThreshold = uicontrol('Parent',f,'style','edit',...
    'position', [100,205,40,20],...
    'string','0');
text3 = uicontrol('Parent', f,'style','text',...
    'position',[140,192,130,30],...
    'string','Distance Threshold (um)',...
    'fontangle','italic');
checkDistanceTransformChannel = uicontrol('Parent', f,'style','checkbox',...
    'position', [105,180,200,25],...
    'string', 'Create Distance Transform channel');

%Statistics LIST
text = uicontrol('Parent', f,'style','text',...
    'position',[35,155,150,25],...
    'string','Track Statistics',...
    'fontangle','italic','fontweight','bold','fontSize', 15);
checkNumberofContacts = uicontrol('Parent', f,'style','checkbox',...
    'position', [40,130,150,25],...
    'string', 'Number of contacts','Value',1);
checkPercentNumberofContacts = uicontrol('Parent', f,'style','checkbox',...
    'position', [40,110,170,25],...
    'string', 'Percent Number of contacts','Value',1);
checkNumberofProlongedContacts = uicontrol('Parent', f,'style','checkbox',...
    'position', [40,90,170,25],...
    'string', 'Number of Prolonged contacts');
checkTotalContactTime = uicontrol('Parent', f,'style','checkbox',...
    'position', [220,130,150,25],...
    'string', 'Total contact time');
checkTotalTrackDurationAway = uicontrol('Parent', f,'style','checkbox',...
    'position', [220,110,150,25],...
    'string', 'Total non-contact time');
checkLongestContactEvent = uicontrol('Parent', f,'style','checkbox',...
    'position', [220,90,150,25],...
    'string', 'Longest contact event');
checkMeanContactEvent = uicontrol('Parent', f,'style','checkbox',...
    'position', [220,70,150,25],...
    'string', 'Mean contact Event');

% Create ANALYZE and Cancel pushbuttons
analyzebutton = uicontrol('Parent', f,'style','pushbutton','units','pixels',...
    'position',[50,10,130,40],'string','ANALYZE',...
    'callback',@Analyze,...
    'fontweight','bold','fontSize', 17);
cancelbutton = uicontrol('Parent', f,'style','pushbutton','units','pixels',...
    'position',[200,20,70,20],'string','CANCEL',...
    'callback',@Cancel);
CheckAllbutton = uicontrol('Parent', f,'style','pushbutton','units','pixels',...
    'position',[200,156,120,20],'string','CheckAll Stats',...
    'callback',@CheckAll,'fontangle','italic');

function Cancel(varargin)
global f
delete (f);
return

function CheckAll(varargin)
global vImarisApplicationGlobal f count
global checkNumberofContacts checkNumberofProlongedContacts checkPercentNumberofContacts
global checkTotalContactTime checkPercentContactTime checkLongestContactEvent checkMeanContactEvent
global checkTotalTrackDurationAway

%Counter to determine state of the checkall button
if mod(count,2) == 0
    count=count+1;
    N=1;
else
    count=count+1;
    N=0;
end

checkNumberofContacts = uicontrol('Parent', f,'style','checkbox',...
    'position', [40,130,150,25],...
    'string', 'Number of contacts','Value',N);
checkPercentNumberofContacts = uicontrol('Parent', f,'style','checkbox',...
    'position', [40,110,170,25],...
    'string', 'Percent Number of contacts','Value',N);
checkNumberofProlongedContacts = uicontrol('Parent', f,'style','checkbox',...
    'position', [40,90,170,25],...
    'string', 'Number of Prolonged contacts','Value',N);
checkTotalContactTime = uicontrol('Parent', f,'style','checkbox',...
    'position', [220,130,150,25],...
    'string', 'Total contact time','Value',N);
checkTotalTrackDurationAway = uicontrol('Parent', f,'style','checkbox',...
    'position', [220,110,150,25],...
    'string', 'Total non-contact time','Value',N);
checkLongestContactEvent = uicontrol('Parent', f,'style','checkbox',...
    'position', [220,90,150,25],...
    'string', 'Longest contact event','Value',N);
checkMeanContactEvent = uicontrol('Parent', f,'style','checkbox',...
    'position', [220,70,150,25],...
    'string', 'Mean contact Event','Value',N);

function Analyze(varargin)
global vImarisApplicationGlobal f
global checkSurfaceOverlap checkCreateColocSurface checkDistanceTransform
global checkDistanceTransformChannel editboxDistanceThreshold
global checkNumberofContacts checkNumberofProlongedContacts checkPercentNumberofContacts
global checkTotalContactTime checkPercentContactTime checkLongestContactEvent checkMeanContactEvent
global checkTotalTrackDurationAway qTotalTrackDurationAway
global qSurfaceOverlap qCreateColocSurface qDistanceTransform qCreateDistanceTransformChannel
global qNumberofContacts qNumberofProlongedContacts qPercentNumberofContacts
global qTotalContactTime qPercentContactTime qLongestContactEvent qMeanContactEvent

vImarisApplication = vImarisApplicationGlobal;
% connect to Imaris interface
if isa(vImarisApplication, 'Imaris.IApplicationPrxHelper')
    
    %checkSurfaceOverlap
    if get(checkSurfaceOverlap,'Value')==1;
        qSurfaceOverlap='YES';
    else
        qSurfaceOverlap='NO';
    end
    %checkCreateColocSurface
    if get(checkCreateColocSurface,'Value')==1;
        qCreateColocSurface='YES';
    else
        qCreateColocSurface='NO';
    end
    %checkDistanceTransform
    if get(checkDistanceTransform,'Value')==1;
        qDistanceTransform='YES';
    else
        qDistanceTransform='NO';
    end
    %checkDistanceTransform
    if get(checkDistanceTransformChannel,'Value')==1;
        qCreateDistanceTransformChannel='YES';
    else
        qCreateDistanceTransformChannel='NO';
    end
    %
    if isequal(qSurfaceOverlap, 'YES') && isequal(qDistanceTransform, 'YES')
        checkSurfaceOverlap = uicontrol('Parent', f,'style','checkbox',...
            'position', [65,332,200,27],'string', 'Surface Overlap',...
            'fontweight','bold','fontSize', 11,'Value',0);
        checkDistanceTransform = uicontrol('Parent', f,'style','checkbox',...
            'position', [65,255,250,25],'string', 'Distance Transform',...
            'Value', 1,'fontweight','bold','fontSize', 11);
        
        msgbox('Please only chooose one method of analysis')
        
        return
    end
    %
    %Statistics
    %checkNumberofContacts
    if get(checkNumberofContacts,'Value')==1;
        qNumberofContacts='YES';
    else
        qNumberofContacts='NO';
    end
    %checkNumberofProlongedContacts
    if get(checkNumberofProlongedContacts,'Value')==1;
        qNumberofProlongedContacts='YES';
    else
        qNumberofProlongedContacts='NO';
    end
    %checkPercentNumberofContacts
    if get(checkPercentNumberofContacts,'Value')==1;
        qPercentNumberofContacts='YES';
    else
        qPercentNumberofContacts='NO';
    end
    %checkTotalContactTime
    if get(checkTotalContactTime,'Value')==1;
        qTotalContactTime='YES';
    else
        qTotalContactTime='NO';
    end
    %checkTotalTrackDurationAway
    if get(checkTotalTrackDurationAway,'Value')==1;
        qTotalTrackDurationAway='YES';
    else
        qTotalTrackDurationAway='NO';
    end
    %checkLongestContactEvent
    if get(checkLongestContactEvent,'Value')==1;
        qLongestContactEvent='YES';
    else
        qLongestContactEvent='NO';
    end
    %checkMeanContactEvent
    if get(checkMeanContactEvent,'Value')==1;
        qMeanContactEvent='YES';
    else
        qMeanContactEvent='NO';
    end
    
    
    %%
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
    vSpots=vImarisApplication.GetFactory.ToSpots(vImarisApplication.GetSurpassSelection);
    vSpotsSelected = vImarisApplication.GetFactory.IsSpots(vSpots);
    
    if vSpotsSelected==1
        vScene = vSpots.GetParent;
        Selected='NO';
        if isequal(qSurfaceOverlap,'YES')
            msgbox('Surface Overlap method NOT valid for Spots object!');
            checkSurfaceOverlap = uicontrol('Parent', f,'style','checkbox',...
                'position', [65,332,200,27],'Value', 0,...
                'string', 'Surface Overlap (slower)','fontweight','bold','fontSize', 11);
            checkCreateColocSurface = uicontrol('Parent', f,'style','checkbox',...
                'position', [100,283,200,20], 'string', 'Create ColocSurface');
            checkDistanceTransform = uicontrol('Parent', f,'style','checkbox',...
                'position', [65,255,250,25],...
                'string', 'Distance Transform (default)','Value', 1,'fontweight','bold','fontSize', 11);
            return
        end
    elseif vSurfacesSelected
        vScene = vSurfaces.GetParent;
        CurrentSurface=char(vSurfaces.GetName);
        Selected='YES';
    else
        vScene = vImarisApplication.GetSurpassScene;
        Selected='NO';
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
    
    Rnew = vNamesList(~cellfun(@isempty, vNamesList));%Remove Surpass scene object preselected
    
    if vNumberOfSurfaces<2 && vSpotsSelected==0
        msgbox('Please create at least 2 surfaces objects!');
        return;
    end
    
    vNamesList = vNamesList(1:vNumberOfSurfaces);
    %%
    %Choose the surfaces
    
    if isequal (Selected, 'NO') && vSpotsSelected==0
        vTargetSurface=1;
        vSurfaces=1;
        
        while vTargetSurface==vSurfaces
            %Create Dialog box and allow user to choose the Reference Position
            vPair = [];
            [vPair, vOk] = listdlg('ListString',vNamesList,'SelectionMode','multiple',...
                'ListSize',[250 150],'Name','Target Surface','InitialValue',[1,1], ...
                'PromptString',{'Please select Target Surface'});
            vTargetSurface = vSurfacesList{vPair(1)};
            %Create Dialog box and allow user to choose the Reference Position
            vPair = [];
            [vPair, vOk] = listdlg('ListString',vNamesList,'SelectionMode','multiple',...
                'ListSize',[250 150],'Name','Primary Tracked Surface Selection','InitialValue',[2,2], ...
                'PromptString',{'Please select Primary tracked surface'});
            vSurfaces = vSurfacesList{vPair(1)};
            if vTargetSurface==vSurfaces
                uiwait(msgbox('Please choose 2 different surfaces'));
            end
        end
    elseif isequal (Selected, 'YES') && vSpotsSelected==0
        PrimarySceneIndex = strmatch(CurrentSurface,Rnew);%Identify index of selected surface
        if numel(PrimarySceneIndex)==1
            vNamesList{1,PrimarySceneIndex}='DO NOT USE';%Remove tracked surface from list
        end
        vPair = [];
        check1='NO';
        while isequal(check1,'NO')%Test is surface surface was selected for both
            [vPair, vOk] = listdlg('ListString',vNamesList,'SelectionMode','multiple',...
                'ListSize',[250 150],'Name','Target Surface','InitialValue',[1,1], ...
                'PromptString',{'Please select Target Surface'});
            check1=char(vNamesList (vPair(1)));
            if isequal (check1,'DO NOT USE')
                uiwait(msgbox('Please choose Target Surface again','modal'));
                check1='NO';
                vPair = [];
            end
        end
        vTargetSurface = vSurfacesList{vPair(1)};
        
    else
        vPair = [];
        [vPair, vOk] = listdlg('ListString',vNamesList,'SelectionMode','multiple',...
            'ListSize',[250 150],'Name','Target Surface','InitialValue',[1,1], ...
            'PromptString',{'Please select Target Surface'});
        vTargetSurface = vSurfacesList{vPair(1)};
    end
    
    %%
    %Get Image Data parameters
    vDataMin = [vImarisApplication.GetDataSet.GetExtendMinX, vImarisApplication.GetDataSet.GetExtendMinY, vImarisApplication.GetDataSet.GetExtendMinZ];
    vDataMax = [vImarisApplication.GetDataSet.GetExtendMaxX, vImarisApplication.GetDataSet.GetExtendMaxY, vImarisApplication.GetDataSet.GetExtendMaxZ];
    vDataSize = [vImarisApplication.GetDataSet.GetSizeX, vImarisApplication.GetDataSet.GetSizeY, vImarisApplication.GetDataSet.GetSizeZ];
    aSizeC = vImarisApplication.GetDataSet.GetSizeC;
    aSizeT = vImarisApplication.GetDataSet.GetSizeT;
    
    Xvoxelspacing= (vDataMax(1)-vDataMin(1))/vDataSize(1);
    Zvoxelspacing = (vDataMax(3)-vDataMin(3))/vDataSize(3);
    vSmoothingFactor=Xvoxelspacing*2;
    
    
    if isequal (qCreateColocSurface,'YES') && vSpotsSelected==0;
        %add additional channel
        TotalNumberofChannels=aSizeC+1;
        vLastChannel=TotalNumberofChannels-1;
        %clone Dataset
        vDataSet = vImarisApplication.GetDataSet.Clone;
        vDataSet.SetSizeC(aSizeC + 1);
        
        %Generate dialog for smoothed or unsmoothed surface
        %Have option to chose: 1) no smoothing, 2) default, 3) Custom
        %vSurfaceSmoothingType = false;
        %qstring={'Please choose a how to generate Coloc surface '}
        %vAnswer = questdlg(qstring, 'Colocalize Surfaces', ...
        %        'No Smoothing', 'Smoothing', 'Smoothing');
        %if(isequal(vAnswer, 'Cancel') || isempty(vAnswer))
        %    return
        %end
        %vSurfaceSmoothingType = isequal(vAnswer, 'Smoothing');
        %if vSurfaceSmoothingType==true
        %    vSmoothingFactorName = num2str(vSmoothingFactor);
        %    qstring={'Please set the smoothing factor -- default is double image voxel size'};
        %    vAnswer2 = inputdlg(qstring,'Smoothing Factor',1,{vSmoothingFactorName});
        %    if isempty(vAnswer2), return, end
        %else
        vSmoothingFactor = 0;
        %end
    end
    %%
    %Identify ID for each track
    if vSpotsSelected==1
        vAllIds = vSpots.GetIds;
        edges = vSpots.GetTrackEdges + 1;
        TrackIds=vSpots.GetTrackIds;
    else
        vAllIds = vSurfaces.GetIds;
        edges = vSurfaces.GetTrackEdges + 1;
        TrackIds=vSurfaces.GetTrackIds;
    end
    
    if ~isempty(edges)
        vTrackID=zeros(size(vAllIds,1),1);
        uniqueTrackIDs=unique(TrackIds);
        % loop each trackID fin
        for nextTrackID=1:size(uniqueTrackIDs)
            ObjectIndex=unique(edges(find(TrackIds==uniqueTrackIDs(nextTrackID)),:));
            vtrackID(ObjectIndex)=uniqueTrackIDs(nextTrackID);
        end
        
        
        
        %old code broken in 9.1
        %     edges_forspots = 1:size(vAllIds);%size(vTracks.GetPositionsXYZ, 1);
        %     edges_forspots( : ) = size(edges, 1) + 1; % initialize array to fictive edge
        %     edges_forspots(edges(:, 1)) = 1:size(edges, 1);
        %     edges_forspots(edges(:, 2)) = 1:size(edges, 1);
        %     trackid_foredges = [TrackIds; 0]; % add fictive track id
        %     trackid_forspots = trackid_foredges(edges_forspots);
        %convert TrackID into a single integer
        %OriginalTrackID = double(trackid_forspots);
        vtrackID = vtrackID-1000000000;
        vtrackIDmax = size(uniqueTrackIDs,1);
    else
        vtrackIDmax=0;
        vtrackID=zeros(numel(vAllIds),1);
    end
    %Set the size of the new variable for overlap volume
    vAllTimeIndices=[];
    vNumberofContactsTotal=[];
    
    
    %%
    
    if isequal (qDistanceTransform, 'YES')
        
        vImarisDataSet = vImarisApplication.GetDataSet.Clone;
        vNumberOfChannels = vImarisDataSet.GetSizeC;
        vImarisDataSet.SetSizeC(vImarisDataSet.GetSizeC + 1);%Add new channel
        %Convert to 32bit
        vImarisDataSet.SetType(vImarisDataSet.GetType.eTypeFloat);
        % Create a new channel where the result will be sent
        vImarisDataSet.SetChannelName(vNumberOfChannels,['Distance to ', char(vTargetSurface.GetName)]);
        vImarisDataSet.SetChannelColorRGBA(vNumberOfChannels, 255*256*256);
        %Get Distance threshold from user input
        DistanceThreshold=str2double(get(editboxDistanceThreshold,'String'));%
        
        %Identify the vAllTimeIndices points that contain tracked surfaces
        if vSpotsSelected==1
            vAllTimeIndices=vSpots.GetIndicesT;
        else
            for GetTimeIndex = 0:size(vAllIds,1)-1;
                vAllTimeIndices=[vAllTimeIndices;vSurfaces.GetTimeIndex(GetTimeIndex)];
            end
        end
        vValidTimePoints=unique (vAllTimeIndices);
        
        vProgressDisplay = waitbar(0, 'Distance Transform: Preparation');
        %Calculate distance transform for each vAllTimeIndices point
        for vTime = 0:size(vValidTimePoints,1)-1
            CurrentTimePoint=vValidTimePoints(vTime+1);
            vMaskDataSetTarget = vTargetSurface.GetMask( ...
                vDataMin(1), vDataMin(2), vDataMin(3), ...
                vDataMax(1), vDataMax(2), vDataMax(3), ...
                vDataSize(1), vDataSize(2), vDataSize(3), CurrentTimePoint);
            for vIndexZ = 1:vDataSize(3)
                vSlice=vMaskDataSetTarget.GetDataSubVolumeAs1DArrayBytes(...
                    0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                vSlice = vSlice == 1;%Outside Mask
                
                %Remove border voxels and set to zero
                %vSlice(1:5,:)=0;%set first row to zero
                %vSlice((end-(vDataSize(1)-1)):end,:)=0;%Set last row to zero
                %vSlice(1:vDataSize(1):end)=0;%set left left to zero
                %vSlice(vDataSize(1):vDataSize(1):end)=0;%Set right side to zero
                
                vImarisDataSet.SetDataSubVolumeAs1DArrayFloats(vSlice, ...
                    0,0,vIndexZ-1,vNumberOfChannels,CurrentTimePoint,vDataSize(1),vDataSize(2),1);
            end
            waitbar((vTime+1+1)/(size(vValidTimePoints,1)), vProgressDisplay);
        end
        waitbar(0.5, vProgressDisplay, 'Distance Transform: Calculation');
        vImarisApplication.GetImageProcessing.DistanceTransformChannel( ...
            vImarisDataSet, vNumberOfChannels, 1, false);
        waitbar(1, vProgressDisplay);
        close (vProgressDisplay);
        vImarisApplication.SetDataSet(vImarisDataSet);
        
        %Get Original Stats order vNewSpots
        if vSpotsSelected==1
            vAllStatistics = vSpots.GetStatistics;
            vSpotStatNames = cell(vAllStatistics.mNames);
            vSpotStatValues = vAllStatistics.mValues;
            vSpotsStatIds = vAllStatistics.mIds;
            vSpotIndex=strmatch('Intensity Min', vSpotStatNames);
            vAllMinDist = vSpotStatValues(vSpotIndex,:);
            vObjectIDsTemp=vSpotsStatIds(vSpotIndex,:)
        else
            vAllStatistics = vSurfaces.GetStatistics;
            vSurfaceStatNames = cell(vAllStatistics.mNames);
            vSurfaceStatValues = vAllStatistics.mValues;
            vSurfaceStatAllIds = vAllStatistics.mIds;
            vSurfaceIndex=strmatch('Intensity Min', vSurfaceStatNames);
            vAllMinDist = vSurfaceStatValues(vSurfaceIndex,:);
            vObjectIDsTemp=vSurfaceStatAllIds(vSurfaceIndex,:);
        end
        %Question to Distance transform channel
        
        if isequal (qCreateDistanceTransformChannel, 'NO')
            vImarisDataSet.SetSizeC(vImarisDataSet.GetSizeC - 1);%Add new channel
        end
        vImarisApplication.SetDataSet(vImarisDataSet);
        
        %Separate out the new Distance Transform channel value
        vMinDist=[];
        vMinDistObjID=[];
        Ids=unique(vObjectIDsTemp);
        for c=1:size(Ids,1)
            x=(vAllMinDist(vObjectIDsTemp==Ids(c)));
            vMinDist=[vMinDist;x(end)];
            vMinDistObjID=[vMinDistObjID;Ids(c)];
        end
        
        %Calculate the Contact intervals for each timepoint
        if aSizeT>1
            if vSpotsSelected==1
                t(1,:) = datetime(cell(vSpots.GetTimePoint(0)),...
                    'Format','yyyy-MM-dd HH:mm:ss.SSS');
                tInt(1,:)=datevec(between(t(1),t(1)));
                for w=1:aSizeT-1
                    t(w+1,:) = datetime(cell(vSpots.GetTimePoint(w)),...
                        'Format','yyyy-MM-dd HH:mm:ss.SSS');
                    tInt(w+1,:)=datevec(between(t(w),t(w+1)));
                    vContactIntervalsec=(tInt(w,6)+tInt(w+1,6))/2;
                    vContactIntervalmin=(tInt(w,5)+tInt(w+1,5))/2*60;
                    vContactIntervalhour=(tInt(w,4)+tInt(w+1,4))/2*60*60;
                    vContactInterval(w,:)=vContactIntervalsec+vContactIntervalmin+vContactIntervalhour;
                    
                end
            else
                t(1,:) = datetime(cell(vSurfaces.GetTimePoint(0)),...
                    'Format','yyyy-MM-dd HH:mm:ss.SSS');
                tInt(1,:)=datevec(between(t(1),t(1)));
                for w=1:aSizeT-1
                    t(w+1,:) = datetime(cell(vSurfaces.GetTimePoint(w)),...
                        'Format','yyyy-MM-dd HH:mm:ss.SSS');
                    tInt(w+1,:)=datevec(between(t(w),t(w+1)));
                    vContactIntervalsec=(tInt(w,6)+tInt(w+1,6))/2;
                    vContactIntervalmin=(tInt(w,5)+tInt(w+1,5))/2*60;
                    vContactIntervalhour=(tInt(w,4)+tInt(w+1,4))/2*60*60;
                    vContactInterval(w,:)=vContactIntervalsec+vContactIntervalmin+vContactIntervalhour;
                end
            end
            vContactIntervalsec=(tInt((w+1),6))/2;
            vContactIntervalmin=(tInt((w+1),5))/2*60;
            vContactIntervalhour=(tInt((w+1),4))/2*60*60;
            vContactInterval(w+1,:)=vContactIntervalsec+vContactIntervalmin+vContactIntervalhour;
            
            %Calculate Kiss and run events. loop all objects in each track
            for trackloop = 1:vtrackIDmax;
                
                %Generate TrackID sequence for new Statistics
                TrackStatID(trackloop)=uniqueTrackIDs(trackloop);
                vSurfacesT = vtrackID == uniqueTrackIDs(trackloop)-1000000000;
                vSurfacesIndex = find(vSurfacesT);%Generate SurfaceIndex based on the logical matrix
                vWorkingMinDist=vMinDist(vSurfacesIndex,:);%Find DistMin for each object in track
                NumberofContacts=sum(vWorkingMinDist<=DistanceThreshold);%Count the number of contacts that meet threshold
                vNumberofContactsTotal=[vNumberofContactsTotal;NumberofContacts];
                vPercentContactTotal(trackloop)=NumberofContacts/numel(vSurfacesIndex)*100;
                
                %Find the indices for contact and non-contact events
                ContactNO=find(vWorkingMinDist>DistanceThreshold);
                ContactYES=find(vWorkingMinDist<=DistanceThreshold);
                %Calculate the number of events, find which are sequential
                SeqGroups = ([0 cumsum(diff(ContactYES')~=1)])'; % find consecutive sequential contact
                NumberofSeqGroups=max(SeqGroups)+1;%Total number
                
                NumberProlongedKissEvents=0;
                TrackDurationAllKissEvents=0;
                if ~isempty (ContactYES)
                    for q=0:NumberofSeqGroups-1
                        vSeqIndex=find(SeqGroups==q);
                        if  numel(vSeqIndex)>1
                            NumberProlongedKissEvents=NumberProlongedKissEvents+1;
                            TrackDurationAllKissEvents(q+1,:)=sum(vContactInterval(vAllTimeIndices(vSurfacesIndex(ContactYES(vSeqIndex)))+1));
                        else
                            TrackDurationAllKissEvents(q+1,:)=sum(vContactInterval(vAllTimeIndices(vSurfacesIndex(ContactYES(vSeqIndex)))+1));
                        end
                    end
                    %Number of prolonged kiss events greater than 2 timepoints
                    %Extract the Longest kiss event
                    LongestKissEventTime(trackloop)=max(TrackDurationAllKissEvents);
                    MeanKissEventTime(trackloop)=mean(TrackDurationAllKissEvents);
                    TotalTrackDurationAllKissEvents(trackloop)=sum(TrackDurationAllKissEvents);
                    TotalNumberofProlongedKissEvents(trackloop)=NumberProlongedKissEvents;
                    TotalTrackDurationAway(trackloop)=sum(vContactInterval(vAllTimeIndices(vSurfacesIndex([ContactNO;ContactYES]))+1))-sum(TrackDurationAllKissEvents);
                else
                    LongestKissEventTime(trackloop)=0;
                    MeanKissEventTime(trackloop)=0;
                    TotalNumberofProlongedKissEvents(trackloop)=0;
                    TotalTrackDurationAllKissEvents(trackloop)=0;
                    TotalTrackDurationAway(trackloop)=sum(vContactInterval(vAllTimeIndices(vSurfacesIndex(ContactNO))+1));
                end
            end
        else
            %when there is only one timepoint
            TrackStatID=0;
            vSurfacesT = vtrackID == 0;
            vSurfacesIndex = find(vSurfacesT);%Generate SurfaceIndex based on the logical matrix
            vNumberofContactsTotal=sum(vMinDist<=DistanceThreshold);%Count the number of contacts that meet threshold
            vPercentContactTotal=vNumberofContactsTotal/numel(vSurfacesIndex)*100;
            TotalNumberofProlongedKissEvents=0;
        end
        
        if vSpotsSelected==1
            %Add new statistic to Spots
            vInd=1:size(vAllIds);
            vIds=vAllIds;
            vUnits(vInd) = {'um'};
            vFactors(vInd) = {'Spot'};
            vFactors(2, vInd) = num2cell(vAllTimeIndices+1);
            vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
            vFactorNames = {'Category','Time'};
            vNames(vInd) = {sprintf(' Distance to %s',char(vTargetSurface.GetName))};
            vSpots.AddStatistics(vNames, vMinDist, vUnits, vFactors, vFactorNames, vIds);
            
            %Overall Statistics
            for i=1:aSizeT
                OverallNumberOfContactsperTimpoint(i)=sum(vMinDist(vAllTimeIndices+1==i)<=DistanceThreshold);
                OverallTotalNumberOfSurfacesperTimpoint(i)=sum(vMinDist(vAllTimeIndices+1==i)>=0);
            end
            PercentageContactsperTimpoint=round(OverallNumberOfContactsperTimpoint./OverallTotalNumberOfSurfacesperTimpoint*100,2);
            PercentageContactsperTimpoint(isnan(PercentageContactsperTimpoint))=0;%set NAN to zero
            
            clear vInd vIds vUnits vFactors vNames
            vInd=1:aSizeT;
            vIds(vInd)=0;
            vUnits(vInd) = {'%'};%{ char(vImarisApplication.GetDataSet.GetUnit) };
            Indices=1:aSizeT;
            vFactors(vInd) = {'Overall'};
            vFactors(2, vInd) = num2cell(Indices);
            vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
            vFactorNames = {'Overall','Time'};
            vNames(vInd) = {sprintf(' Percent Surface Contacts per Timepoint with %s',char(vTargetSurface.GetName))};
            vSpots.AddStatistics(vNames, PercentageContactsperTimpoint', vUnits, vFactors, vFactorNames, vIds);
            vUnits(vInd) = {''};
            vNames(vInd) = {sprintf(' Number of Contacts per Timepoint with %s',char(vTargetSurface.GetName))};
            vSpots.AddStatistics(vNames, OverallNumberOfContactsperTimpoint', vUnits, vFactors, vFactorNames, vIds);
            
            %Set Track Statistics
            
            if aSizeT~=1
    
                vIndT=1:vtrackIDmax;%Total number of tracks
                vUnitsT(vIndT) = {'sec'};
                vFactorsT(vIndT) = {'Track'};
                
                if isequal (qNumberofContacts,'YES')
                    vNamesT(vIndT) = {sprintf(' Number of contacts with %s',char(vTargetSurface.GetName))};
                    vFactorNamesT = {'Category'};
                    vSpots.AddStatistics(vNamesT, vNumberofContactsTotal, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qNumberofProlongedContacts,'YES') && aSizeT>1
                    vNamesT(vIndT) = {sprintf(' Number prolonged contact events with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {''};
                    vSpots.AddStatistics(vNamesT, TotalNumberofProlongedKissEvents, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qPercentNumberofContacts,'YES')
                    vNamesT(vIndT) = {sprintf(' Percent Surface contact with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {'%'};
                    vSpots.AddStatistics(vNamesT, vPercentContactTotal, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qTotalContactTime,'YES') && aSizeT>1
                    vNamesT(vIndT) = {sprintf(' Total Time in contact with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {'sec'};
                    vSpots.AddStatistics(vNamesT, TotalTrackDurationAllKissEvents, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qTotalTrackDurationAway,'YES') && aSizeT>1
                    vNamesT(vIndT) = {sprintf(' Total Time without contact with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {'sec'};
                    vSpots.AddStatistics(vNamesT, TotalTrackDurationAway, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qLongestContactEvent,'YES') && aSizeT>1
                    vNamesT(vIndT) = {sprintf(' Longest contact event with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {'sec'};
                    vSpots.AddStatistics(vNamesT, LongestKissEventTime, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qMeanContactEvent,'YES') && aSizeT>1
                    vNamesT(vIndT) = {sprintf(' Mean Length contact event with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {'sec'};
                    vSpots.AddStatistics(vNamesT, MeanKissEventTime, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                
            else
            end
            %Append original surface name in Surpass Scene
            vSpots.SetName(sprintf('Analyzed Distance threshold - %s',char(vSpots.GetName)));
            vImarisApplication.GetSurpassScene.AddChild(vSpots, -1);
            
        else
            %Add new statistic to Surfaces
            vInd=1:size(vAllIds);
            vIds=vAllIds;
            vUnits(vInd) = {'um'};
            vFactors(vInd) = {'Surface'};
            vFactors(2, vInd) = num2cell(vAllTimeIndices+1);
            vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
            vFactorNames = {'Category','Time'};
            vNames(vInd) = {sprintf(' Distance to %s',char(vTargetSurface.GetName))};
            vSurfaces.AddStatistics(vNames, vMinDist, vUnits, vFactors, vFactorNames, vIds);
            
            %Overall Statistics
            for i=1:aSizeT
                OverallNumberOfContactsperTimpoint(i)=sum(vMinDist(vAllTimeIndices+1==i)<=DistanceThreshold);
                OverallTotalNumberOfSurfacesperTimpoint(i)=sum(vMinDist(vAllTimeIndices+1==i)>=0);
            end
            PercentageContactsperTimpoint=round(OverallNumberOfContactsperTimpoint./OverallTotalNumberOfSurfacesperTimpoint*100,2);
            PercentageContactsperTimpoint(isnan(PercentageContactsperTimpoint))=0;
            
            clear vInd vIds vUnits vFactors vNames
            vInd=1:aSizeT;
            vIds(vInd)=0;
            vUnits(vInd) = {'%'};%{ char(vImarisApplication.GetDataSet.GetUnit) };
            Indices=1:aSizeT;
            vFactors(vInd) = {'Overall'};
            vFactors(2, vInd) = num2cell(Indices);
            vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
            vFactorNames = {'Overall','Time'};
            vNames(vInd) = {sprintf(' Percent Surface Contacts per Timepoint with %s',char(vTargetSurface.GetName))};
            vSurfaces.AddStatistics(vNames, PercentageContactsperTimpoint', vUnits, vFactors, vFactorNames, vIds);
            vUnits(vInd) = {''};
            vNames(vInd) = {sprintf(' Number of Contacts per Timepoint with %s',char(vTargetSurface.GetName))};
            vSurfaces.AddStatistics(vNames, OverallNumberOfContactsperTimpoint', vUnits, vFactors, vFactorNames, vIds);
            
            %Set Track Statistics
            
            if aSizeT~=1
                
                vIndT=1:vtrackIDmax;%Total number of tracks
                vUnitsT(vIndT) = {'sec'};
                vFactorsT(vIndT) = {'Track'};
      
                if isequal (qNumberofContacts,'YES')
                    vNamesT(vIndT) = {sprintf(' Number of contacts with %s',char(vTargetSurface.GetName))};
                    vFactorNamesT = {'Category'};
                    vSurfaces.AddStatistics(vNamesT, vNumberofContactsTotal, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qNumberofProlongedContacts,'YES') && aSizeT>1
                    vNamesT(vIndT) = {sprintf(' Number prolonged contact events with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {''};
                    vSurfaces.AddStatistics(vNamesT, TotalNumberofProlongedKissEvents, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qPercentNumberofContacts,'YES')
                    vNamesT(vIndT) = {sprintf(' Percent Surface contact with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {'%'};
                    vSurfaces.AddStatistics(vNamesT, vPercentContactTotal, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qTotalContactTime,'YES') && aSizeT>1
                    vNamesT(vIndT) = {sprintf(' Total time in contact with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {'sec'};
                    vSurfaces.AddStatistics(vNamesT, TotalTrackDurationAllKissEvents, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qTotalTrackDurationAway,'YES') && aSizeT>1
                    vNamesT(vIndT) = {sprintf(' Total Time without Contact with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {'sec'};
                    vSurfaces.AddStatistics(vNamesT, TotalTrackDurationAway, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qLongestContactEvent,'YES') && aSizeT>1
                    vNamesT(vIndT) = {sprintf(' Longest contact event with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {'sec'};
                    vSurfaces.AddStatistics(vNamesT, LongestKissEventTime, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
                if isequal (qMeanContactEvent,'YES') && aSizeT>1
                    vNamesT(vIndT) = {sprintf(' Mean Length contact event with %s',char(vTargetSurface.GetName))};
                    vUnitsT(vIndT) = {'sec'};
                    vSurfaces.AddStatistics(vNamesT, MeanKissEventTime, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
                end
            else
            end
            %Rename new surface to Surpass Scene
            vSurfaces.SetName(sprintf('Analyzed Distance threshold - %s',char(vSurfaces.GetName)));
            vImarisApplication.GetSurpassScene.AddChild(vSurfaces, -1);
        end   
    end
     
    if isequal (qSurfaceOverlap, 'YES');
        
        WorkingSumColoc=0;
        WorkingSumPrimary=0;
        vProgressDisplay = waitbar(0, 'Surface Overlap Calculation...');
        %vImarisApplication.SetVisible(~vImarisApplication.GetVisible);
        vOverlapVolume=zeros(size(vAllIds));
        vSurfaceMaskVolume=zeros(size(vAllIds));
        StatIdsAll=[];
        TotalColoc=0;
        
        %Calculate vAllTimeIndices intervals
        t(1,:) = datetime(cell(vSurfaces.GetTimePoint(0)),...
            'Format','yyyy-MM-dd HH:mm:ss.SSS');
        tInt(1,:)=datevec(between(t(1),t(1)));
        for w=1:aSizeT-1
            t(w+1,:) = datetime(cell(vSurfaces.GetTimePoint(w)),...
                'Format','yyyy-MM-dd HH:mm:ss.SSS');
            tInt(w+1,:)=datevec(between(t(w),t(w+1)));
            vContactInterval(w,:)=(tInt(w,6)+tInt(w+1,6))/2;
        end
        vContactInterval(w+1,:)=(tInt((w+1),6))/2;
        
        %%
        
        %Calculate vAllTimeIndices Indices
        for i=0:size(vAllIds,1)-1
            vAllTimeIndices(i+1,:)=vSurfaces.GetTimeIndex(i)+1;
        end
        c=1;
        for TimeIndex=1:aSizeT
            vSurfacesIndex = vAllIds(find(vAllTimeIndices==TimeIndex));
            
            vMaskTarget = vTargetSurface.GetMask( ...
                vDataMin(1), vDataMin(2), vDataMin(3), ...
                vDataMax(1), vDataMax(2), vDataMax(3), ...
                vDataSize(1), vDataSize(2), vDataSize(3),TimeIndex-1);
            
            for q=0:size(vSurfacesIndex,1)-1
                
                vMaskPrimary = vSurfaces.GetSingleMask(c-1, ...
                    vDataMin(1), vDataMin(2), vDataMin(3), ...
                    vDataMax(1), vDataMax(2), vDataMax(3), ...
                    vDataSize(1), vDataSize(2), vDataSize(3));
                %Cycle thru stack and find overlap for current surface
                for vIndexZ = 1:vDataSize(3)
                    ch1=vMaskPrimary.GetDataSubVolumeAs1DArrayBytes( ...
                        0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                    ch2=vMaskTarget.GetDataSubVolumeAs1DArrayBytes( ...
                        0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                    %Determine the Voxels that are colocalized
                    Coloc=ch1+ch2;
                    Coloc(Coloc<2)=0;
                    Coloc(Coloc>1)=1;
                    WorkingSumColoc=WorkingSumColoc+sum(Coloc);
                    WorkingSumPrimary=WorkingSumPrimary+sum(ch1);
                end
                vOverlapVolume(c)=WorkingSumColoc*Xvoxelspacing^2*Zvoxelspacing;
                vSurfaceMaskVolume(c)=WorkingSumPrimary*Xvoxelspacing^2*Zvoxelspacing;
                WorkingSumColoc=0;
                WorkingSumPrimary=0;
                
                waitbar(c/size(vAllIds,1), vProgressDisplay);
                c=c+1;
            end
        end
        %%
        if aSizeT>1
            for trackloop = 1:vtrackIDmax;
                
                %Generate TrackID sequence for new Statistics
                TrackStatID(trackloop)=uniqueTrackIDs(trackloop);
                
                vSurfacesT = vtrackID == uniqueTrackIDs(trackloop)-1000000000;
                vSurfacesIndex = find(vSurfacesT);%Generate SurfaceIndex based on the logical matrix
                vWorkingOverlapVolume=vOverlapVolume(vSurfacesIndex,:);%Find DIstMin for each object in track
                NumberofContacts=sum(vWorkingOverlapVolume>0);%Count the number of contacts that meet threshold
                vNumberofContactsTotal=[vNumberofContactsTotal;NumberofContacts];
                vPercentContactTotal(trackloop)=NumberofContacts/numel(vSurfacesIndex)*100;
                
                ContactNO=find(vWorkingOverlapVolume==0);
                ContactYES=find(vWorkingOverlapVolume>0);
                SeqGroups = ([0 cumsum(diff(ContactYES')~=1)])'; % find consecutive seq groups
                NumberofSeqGroups=max(SeqGroups)+1;
                
                NumberProlongedKissEvents=0;
                TrackDurationAllKissEvents=0;
                if ~isempty (ContactYES)
                    for q=0:NumberofSeqGroups-1
                        vSeqIndex=find(SeqGroups==q);
                        if  numel(vSeqIndex)>1
                            NumberProlongedKissEvents=NumberProlongedKissEvents+1;
                            TrackDurationAllKissEvents(q+1,:)=sum(vContactInterval(vAllTimeIndices(vSurfacesIndex(ContactYES(vSeqIndex)))));
                        else
                            TrackDurationAllKissEvents(q+1,:)=sum(vContactInterval(vAllTimeIndices(vSurfacesIndex(ContactYES(vSeqIndex)))));
                        end
                    end
                    %Number of prolonged kiss events greater than 2 timepoints
                    %Convert calendar duration to number and sum all vAllTimeIndices values
                    %Extract the Longest kiss eventTotalTrackDurationofKissEvents(trackloop+1)=sum(TrackDurationAllKissEvents);
                    LongestKissEventTime(trackloop)=max(TrackDurationAllKissEvents);
                    MeanKissEventTime(trackloop)=mean(TrackDurationAllKissEvents);
                    TotalTrackDurationAllKissEvents(trackloop)=sum(TrackDurationAllKissEvents);
                    TotalNumberofProlongedKissEvents(trackloop)=NumberProlongedKissEvents;
                    TotalTrackDurationAway(trackloop)=sum(vContactInterval(vAllTimeIndices(vSurfacesIndex([ContactNO;ContactYES]))))-sum(TrackDurationAllKissEvents);
                else
                    LongestKissEventTime(trackloop)=0;
                    MeanKissEventTime(trackloop)=0;
                    TotalNumberofProlongedKissEvents(trackloop)=0;
                    TotalTrackDurationAllKissEvents(trackloop)=0;
                    TotalTrackDurationAway(trackloop)=sum(vContactInterval(vAllTimeIndices(vSurfacesIndex(ContactNO))));
                end
            end
        else
            
            TrackStatID=0;
            %when there is only one timepoint
            vSurfacesT = vtrackID == 0;
            vSurfacesIndex = find(vSurfacesT);%Generate SurfaceIndex based on the logical matrix
            vNumberofContactsTotal=sum(vOverlapVolume>0);%Count the number of contacts that meet threshold
            vPercentContactTotal=vNumberofContactsTotal/numel(vSurfacesIndex)*100;
        end
        
        %%
        close(vProgressDisplay);
        
        if isequal (qCreateColocSurface,'YES');
            
            vProgressDisplay = waitbar(0, 'Creating overlap Surface reconstruction...');
            waitbar(0.5, vProgressDisplay);
            vDataSet.SetChannelName(vLastChannel,'ColocChannel');
            vDataSet.SetChannelRange(vLastChannel,0,1);
            
            %Generate surface mask for each surface over Time
            for vTimeIndex= 0:aSizeT-1
                vMaskPrimary = vSurfaces.GetMask( ...
                    vDataMin(1), vDataMin(2), vDataMin(3), ...
                    vDataMax(1), vDataMax(2), vDataMax(3), ...
                    vDataSize(1), vDataSize(2), vDataSize(3),vTimeIndex);
                vMaskTarget = vTargetSurface.GetMask( ...
                    vDataMin(1), vDataMin(2), vDataMin(3), ...
                    vDataMax(1), vDataMax(2), vDataMax(3), ...
                    vDataSize(1), vDataSize(2), vDataSize(3),vTimeIndex);
                for vIndexZ = 1:vDataSize(3)
                    ch1=vMaskPrimary.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                    ch2=vMaskTarget.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                    %Determine the Voxels that are colocalized
                    Coloc=ch1+ch2;
                    Coloc(Coloc<2)=0;
                    Coloc(Coloc>1)=1;
                    vDataSet.SetDataSubVolumeAs1DArrayBytes(Coloc, ...
                        0,0,vIndexZ-1,vLastChannel,vTimeIndex,vDataSize(1),vDataSize(2),1);
                end
            end
            
            
            vImarisApplication.SetDataSet(vDataSet);
            %Run the Surface Creation Wizard on the new channel
            ip = vImarisApplication.GetImageProcessing;
            Coloc_surfaces1 = ip.DetectSurfaces(vDataSet, [], vLastChannel, vSmoothingFactor, 0, true, 55, '');
            Coloc_surfaces1.SetName(sprintf('ColocSurface'));
            Coloc_surfaces1.SetColorRGBA((rand(1, 1)) * 256 * 256 * 256 );
            
            %Add new surface to Surpass Scene
            vTargetSurface.SetVisible(0);
            vSurfaces.SetVisible(0);
            vImarisApplication.GetSurpassScene.AddChild(Coloc_surfaces1, -1);
            waitbar(1, vProgressDisplay);
            close(vProgressDisplay);
        end
     
        %Add new statistic to Surfaces
        vInd=1:size(vAllIds);
        vIds=vAllIds;
        vUnits(vInd) = {'um^3'};
        vFactors(vInd) = {'Surface'};
        vFactors(2, vInd) = num2cell(vAllTimeIndices);
        vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
        vFactorNames = {'Category','Time'};
        vNames(vInd) = {sprintf(' Surface Overlap Volume with %s',char(vTargetSurface.GetName))};
        vSurfaces.AddStatistics(vNames, vOverlapVolume, vUnits, vFactors, vFactorNames, vIds);
        vNames(vInd) = {' Surface Mask Volume'};
        vSurfaces.AddStatistics(vNames, vSurfaceMaskVolume, vUnits, vFactors, vFactorNames, vIds);
      
        %Overall Statistics
        for i=1:aSizeT
            OverallNumberOfContactsperTimpoint(i)=sum(vOverlapVolume(vAllTimeIndices==i)>0);
            OverallTotalNumberOfSurfacesperTimpoint(i)=sum(vOverlapVolume(vAllTimeIndices==i)>=0);
        end
        PercentageContactsperTimpoint=round(OverallNumberOfContactsperTimpoint./OverallTotalNumberOfSurfacesperTimpoint*100,2);
        PercentageContactsperTimpoint(isnan(PercentageContactsperTimpoint))=0;
        
        %Overall statistics adding
        clear vInd vIds vUnits vFactors vNames
        vInd=1:aSizeT;
        vIds(vInd)=0;
        vUnits(vInd) = {'%'};
        Indices=1:aSizeT;
        vFactors(vInd) = {'Overall'};
        vFactors(2, vInd) = num2cell(Indices);
        vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
        vFactorNames = {'Overall','Time'};
        vNames(vInd) = {sprintf(' Percent Surface Contacts per Timepoint with %s',char(vTargetSurface.GetName))};
        vSurfaces.AddStatistics(vNames, PercentageContactsperTimpoint', vUnits, vFactors, vFactorNames, vIds);
        vUnits(vInd) = {''};
        vNames(vInd) = {sprintf(' Number of Contacts per Timepoint with %s',char(vTargetSurface.GetName))};
        vSurfaces.AddStatistics(vNames, OverallNumberOfContactsperTimpoint', vUnits, vFactors, vFactorNames, vIds);
 
        %Set Track Statistics
        if aSizeT~=1
            vIndT=1:vtrackIDmax;%Total number of tracks
            vUnitsT(vIndT) = {''};
            vFactorsT(vIndT) = {'Track'};
            
            if isequal (qNumberofContacts,'YES')
                vNamesT(vIndT) = {sprintf(' Number of contacts with %s',char(vTargetSurface.GetName))};
                vFactorNamesT = {'Category'};
                vSurfaces.AddStatistics(vNamesT, vNumberofContactsTotal, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
            end
            if isequal (qNumberofProlongedContacts,'YES') && aSizeT>1
                vNamesT(vIndT) = {sprintf(' Number prolonged contact events with %s',char(vTargetSurface.GetName))};
                vUnitsT(vIndT) = {''};
                vSurfaces.AddStatistics(vNamesT, TotalNumberofProlongedKissEvents, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
            end
            if isequal (qPercentNumberofContacts,'YES')
                vNamesT(vIndT) = {sprintf(' Percent Surface contact with %s',char(vTargetSurface.GetName))};
                vUnitsT(vIndT) = {'%'};
                vSurfaces.AddStatistics(vNamesT, vPercentContactTotal, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
            end
            if isequal (qTotalContactTime,'YES') && aSizeT>1
                vNamesT(vIndT) = {sprintf(' Total Time in contact with %s',char(vTargetSurface.GetName))};
                vUnitsT(vIndT) = {'sec'};
                vSurfaces.AddStatistics(vNamesT, TotalTrackDurationAllKissEvents, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
            end
            if isequal (qTotalTrackDurationAway,'YES') && aSizeT>1
                vNamesT(vIndT) = {sprintf(' Total Time without contact with %s',char(vTargetSurface.GetName))};
                vUnitsT(vIndT) = {'sec'};
                vSurfaces.AddStatistics(vNamesT, TotalTrackDurationAway, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
            end
            if isequal (qLongestContactEvent,'YES') && aSizeT>1
                vNamesT(vIndT) = {sprintf(' Longest contact event with %s',char(vTargetSurface.GetName))};
                vUnitsT(vIndT) = {'sec'};
                vSurfaces.AddStatistics(vNamesT, LongestKissEventTime, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
            end
            if isequal (qMeanContactEvent,'YES') && aSizeT>1
                vNamesT(vIndT) = {sprintf(' Mean length contact event with %s',char(vTargetSurface.GetName))};
                vUnitsT(vIndT) = {'sec'};
                vSurfaces.AddStatistics(vNamesT, MeanKissEventTime, vUnitsT, vFactorsT, vFactorNamesT, TrackStatID);
            end
            
        else
        end
        
        %Rename new surface to Surpass Scene
        vSurfaces.SetName(sprintf('Analyzed Surface Overlap - %s',char(vSurfaces.GetName)));
        vImarisApplication.GetSurpassScene.AddChild(vSurfaces, -1);
        %vImarisApplication.SetVisible(~vImarisApplication.GetVisible);
        
    end
end
delete (f);

