% <CustomTools>
%   <Menu>
%     <Item name="AddIDStatisticsMJG" icon="Matlab" tooltip="Adds a new Statistics value ID to the selected spots or surface component">
%       <Command>MatlabXT::XT_MJG_AddIDStatistics2(%i)</Command>
%     </Item>
%   </Menu>
% </CustomTools>

function XT_MJG_AddIDStatistics2(aImarisApplicationID)
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


%if isequal(vImarisApplication, [])
%  display('Could not connect to Imaris')
%  return
%end
%%
vSelection = vImarisApplication.GetSurpassSelection;
vSpots = vImarisApplication.GetFactory.ToSpots(vSelection);
vSurfaces = vImarisApplication.GetFactory.ToSurfaces(vSelection);
vFilaments = vImarisApplication.GetFactory.ToFilaments(vSelection);
vCells = vImarisApplication.GetFactory.ToCells(vSelection);

if isequal(vSpots, []) && isequal(vSurfaces, []) && isequal(vFilaments, []) && isequal(vCells, []);
  display('No spots or surfaces selected')
  return
end
%%
%Is a spot there
if ~isequal(vSpots, [])
  vCategory = 'Spot';
  vEdges = vSpots.GetTrackEdges + 1; % indices start from 1 here (matlab)
  vAllIds = vSpots.GetIds;
  if isempty(vEdges)%Check to see if there are tracked items;
    vItemIndicesT = vSpots.GetIndicesT;
    vNumberOfItems = size(vItemIndicesT,1);
    vInd = 1:vNumberOfItems;
    vIds = vInd - 1;
    vUnits(1:vNumberOfItems) = {''};%create an array of cells containing empty string
    vFactors(vInd) = {'Spot'};
    vFactors(2, vInd) = {'1'};
    vFactorNames = {'Category', 'Time'};
    vNames(vInd) = {' IDSpots'};
    vSpots.AddStatistics(vNames, vInd', vUnits, vFactors, vFactorNames, vAllIds);
    return
  else
    vedges_forspots = 1:size(vSpots.GetPositionsXYZ, 1);
    vedges_forspots( : ) = size(vEdges, 1) + 1; % initialize array to fictive edge
    vedges_forspots(vEdges(:, 1)) = 1:size(vEdges, 1);
    vedges_forspots(vEdges(:, 2)) = 1:size(vEdges, 1);
    trackid_foredges = [vSpots.GetTrackIds; 0]; % add fictive track id
    trackid_forspots = trackid_foredges(vedges_forspots);
    tempx = double(trackid_forspots);
    vtrackID = tempx-1000000000;
	trackIDmax = max(vtrackID+1);
    vNumberOfItems = trackIDmax;
%Set the Track ID for Spot statistics    
    vSpotID=(vtrackID+1)'; 
    vItemIndicesT = vSpots.GetIndicesT;
    vNumberOfItems = size(vItemIndicesT,1);
    vInd = 1:vNumberOfItems;
    vIDs = vInd - 1;
    vUnits(1:vNumberOfItems) = {''};%create an array of cells containing empty string
    vFactorNames = {'Category', 'Time'};
    vFactors(vInd) = {'Spot'};
    vFactors(2, vInd) = num2cell((vItemIndicesT+1)');
    vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
    vNames(vInd) = {' IDSpots'};
    vSpots.AddStatistics(vNames, vSpotID, vUnits, vFactors, vFactorNames, vAllIds);
    clear vUnits
    clear vFactors
    clear vNames
    clear vFactorNames
    clear vIDs
%Add ID as a new statistics       
    vValues = (1:trackIDmax);
    vInd = 1:trackIDmax;
    vFactorNames = {'Category'};
    vNames(vInd) = {' IDSpotsTrack'};
    vUnits(1:trackIDmax) = {''};%create an array of cells containing empty string
    vIDs = vInd - 1;
    vFactors(vInd) = {'Track'};
    vSpots.AddStatistics(vNames, vValues, vUnits, vFactors, vFactorNames, vIDs);
    return
  end
end
%%
%Is there a surface
if ~isequal(vSurfaces, [])
   vCategory = 'Surface';
   vEdges = vSurfaces.GetTrackEdges + 1; % indices start from 1 here (matlab)
   vNumberOfItems = vSurfaces.GetNumberOfSurfaces;
   vIndicesT = vSurfaces.GetIds;

   if isempty(vEdges)%Check to see if there are tracked items;
    vInd = 1:size(vIndicesT,1);   
    vIds = vInd - 1;
    vUnits(1:size(vIndicesT,1)) = {''};%create an array of cells containing empty string
    vFactors(vInd) = {'Surface'};
    vFactors(2, vInd) = {'1'};
    vFactorNames = {'Category', 'Time'};
    vNames(vInd) = {' IDSurfaces'};
    vSurfaces.AddStatistics(vNames, vInd', vUnits, vFactors, vFactorNames, vIndicesT);
    return
   else
   vTrackID=zeros(size(vIndicesT,1),1);
   uniqueTrackIDs=unique(vSurfaces.GetTrackIds);
   for nextTrackID=1:size(uniqueTrackIDs)    
       ObjectIndex=unique(vEdges(find(vSurfaces.GetTrackIds==uniqueTrackIDs(nextTrackID)),:));
       vtrackID(ObjectIndex)=uniqueTrackIDs(nextTrackID);  
   end
           
%     vedges_forspots = 1:vNumberOfItems;
%     vedges_forspots( : ) = size(vEdges, 1) + 1; % initialize array to fictive edge
%     vedges_forspots(vEdges(:, 1)) = 1:size(vEdges, 1);
%     vedges_forspots(vEdges(:, 2)) = 1:size(vEdges, 1);
%     trackid_foredges = [; 0]; % add fictive track id
%     trackid_forspots = trackid_foredges(vedges_forspots);
%     tempx = double(trackid_forspots);
%     vtrackID = tempx-1000000000;
    vtrackID = vtrackID-1000000000;
	trackIDmax = max(vtrackID+1);
    vNumberOfTracks = trackIDmax;
%Set the Track ID for Surface statistics    
    vSurfaceID=(vtrackID)';
    for i=0:vNumberOfItems-1
        vAllTimeIndices(i+1)=vSurfaces.GetTimeIndex(i);
    end  
    vIndicesT = vSurfaces.GetIds;
    vNumberOfItems = size(vIndicesT,1);
    %vInd = 1:vNumberOfItems;
    vInd = 1:size(vIndicesT,1);   
    vIDs = vInd - 1;
    vUnits(1:vNumberOfItems) = {''};%create an array of cells containing empty string
    vFactorNames = {'Category', 'Time'};
    vFactors(vInd) = {'Surface'};
    vFactors(2, vInd) = num2cell(sort(vAllTimeIndices));
    vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
    vNames(vInd) = {' IDSurfaces'};
    vSurfaces.AddStatistics(vNames, vInd', vUnits, vFactors, vFactorNames, vIndicesT);
    clear vUnits
    clear vFactors
    clear vNames
    clear vFactorNames
    clear vIDs
%Add ID as a new statistic    
    vValues = (1:size(uniqueTrackIDs,1));
    vInd = 1:size(uniqueTrackIDs,1);
    vFactorNames = {'Category'};
    vNames(vInd) = {' IDSurfacesTrack'};
    vUnits(1:size(uniqueTrackIDs,1)) = {''};%create an array of cells containing empty string
    vIDs = uniqueTrackIDs - 1000000000;
    vFactors(vInd) = {'Track'};
    vSurfaces.AddStatistics(vNames, vValues, vUnits, vFactors, vFactorNames, vIDs);
    return
   end
end
%%
%Is there a filament and add ID as new statistic
if ~isequal(vFilaments, [])
    vCategory = 'Filament';
    vInd = 1:vFilaments.GetNumberOfFilaments;
    vIds = vInd - 1 + 100000000;
    vUnits(vInd) = {''};
    vFactors(vInd) = {'Filament'};
    vFactors(2, vInd) = {'1'};
    vFactorNames = {'Category', 'Time'};
    vNames(vInd) = {' IDFilaments'};
    vFilaments.RemoveStatistics('IDFilaments');
    vFilaments.AddStatistics(vNames, vInd', vUnits, vFactors, vFactorNames, vIds);
    return
end
%%
%is there an Imaris Cell
if ~isequal(vCells, [])
    vCategory = 'Cell';
    vEdges = vCells.GetTrackEdges + 1; % indices start from 1 here (matlab)
    vNumberOfItems = vCells.GetNumberOfCells;
  if isempty(vEdges)%Check to see if there are track data
    vInd = 1:vNumberOfItems;
    vIds = vInd - 1;
    vUnits(1:vNumberOfItems) = {''};%create an array of cells containing empty string
    vFactors(vInd) = {'Surface'};
    vFactors(2, vInd) = {'1'};
    vFactorNames = {'Category', 'Time'};
    vNames(vInd) = {' IDCells'};
    vCells.AddStatistics(vNames, vInd', vUnits, vFactors, vFactorNames, vIds); 
  else
    vedges_forspots = 1:vNumberOfItems;
    vedges_forspots( : ) = size(vEdges, 1) + 1; % initialize array to fictive edge
    vedges_forspots(vEdges(:, 1)) = 1:size(vEdges, 1);
    vedges_forspots(vEdges(:, 2)) = 1:size(vEdges, 1);
    trackid_foredges = [vCells.GetTrackIds; 0]; % add fictive track id
    trackid_forspots = trackid_foredges(vedges_forspots);
    tempx = double(trackid_forspots);
    vtrackID = tempx-1000000000;
	trackIDmax = max(vtrackID+1);
    vNumberOfItems = trackIDmax;
%Set the Track ID for Spot statistics    
    vCellID=(vtrackID+1)'; 
    vItemIndicesT = vCells.GetIndicesT;
    vNumberOfItems = size(vItemIndicesT,1);
    vInd = 1:vNumberOfItems;
    vIDs = vInd - 1;
    vUnits(1:vNumberOfItems) = {''};%create an array of cells containing empty string
    vFactorNames = {'Category', 'Time'};
    vFactors(vInd) = {'Cell'};
    vFactors(2, vInd) = num2cell((vItemIndicesT+1)');
    vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
    vNames(vInd) = {' IDCells'};
    vSpots.AddStatistics(vNames, vCellID, vUnits, vFactors, vFactorNames, vIDs);
    clear vUnits
    clear vFactors
    clear vNames
    clear vFactorNames
    clear vIDs
%Add ID as a new statistic    
    vValues = (1:trackIDmax);
    vInd = 1:trackIDmax;
    vFactorNames = {'Category'};
    vNames(vInd) = {' IDCellsTrack'};
    vUnits(1:trackIDmax) = {''};%create an array of cells containing empty string
    vIDs = vInd - 1;
    vFactors(vInd) = {'Track'};
    vCells.AddStatistics(vNames, vValues, vUnits, vFactors, vFactorNames, vIDs);
  end 
%Calculate Nucleus track IDs.    
    vEdges = vCells.GetTrackEdges + 1; % indices start from 1 here (matlab)
    vNumberOfNuclei = sum(vCells.GetNumberOfNuclei);
  if isempty(vEdges)%Check to see if there are track data
    vInd = 1:vNumberOfNuclei;
    vIds = vInd - 1;
    vUnits(1:vNumberOfNuclei) = {''};%create an array of cells containing empty string
    vFactors(vInd) = {'Surfaces'};
    vFactors(2, vInd) = {'1'};
    vFactorNames = {'Category', 'Time'};
    vNames(vInd) = {' IDNuclei'};
    vCells.AddStatistics(vNames, vInd', vUnits, vFactors, vFactorNames, vIds); 
  else
    vedges_forspots = 1:vNumberOfItems;
    vedges_forspots( : ) = size(vEdges, 1) + 1; % initialize array to fictive edge
    vedges_forspots(vEdges(:, 1)) = 1:size(vEdges, 1);
    vedges_forspots(vEdges(:, 2)) = 1:size(vEdges, 1);
    trackid_foredges = [vCells.GetTrackIds; 0]; % add fictive track id
    trackid_forspots = trackid_foredges(vedges_forspots);
    tempx = double(trackid_forspots);
    vtrackID = tempx-1000000000;
	trackIDmax = max(vtrackID+1);
    vNumberOfItems = trackIDmax;
%Add ID as a new statistic    
    vValues = (1:trackIDmax);
    vInd = 1:trackIDmax;
    vFactorNames = {'Category'};
    vNames(vInd) = {' IDNuclei'};
    vUnits(1:trackIDmax) = {''};%create an array of cells containing empty string
    vIDs = vInd - 1;
    vFactors(vInd) = {'Track'};
    vCells.AddStatistics(vNames, vValues, vUnits, vFactors, vFactorNames, vIDs);    
    
  end

end