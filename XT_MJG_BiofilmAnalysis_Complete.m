%Biofilm Analysis

%Written by Matthew J. Gastinger, PhD
%Bitplane Advanced Application Scientist
%2016 May 20
%
%<CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="BiofilmAnalysis" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_BiofilmAnalysis_Complete(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="BiofilmAnalysis" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_BiofilmAnalysis_Complete(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description
%This XTension will measure the thickness and Roughness of isosurfaces,
%overall values and/or individual measurements using 3 methods:
%1. Surface Mask
%2. BioVolume (no gaps, to substratum)
%3. Surface to Surface edge (subvoxel)
%Three types of measurments will be calculated: Mean thickness, Max
%thickness, and a Roughness coefficient (variance).  For detailed description
%of these value, please see PDF attached to the zip file.


%Additional overall measurements include:
%--color map of mean thickness values
%--Biomass (um3/um2)
%--Biovolume
%--Continuity Ratio (porosity)






function XT_MJG_BiofilmAnalysis_Complete(aImarisApplicationID)


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


%Launch the GUI
%Collect of the data for the XTension
%Set some variables to global
global vImarisApplicationGlobal qHeatMap editboxReductionFactor checkAutoDetectSubstratumStart
global f qOverall qFlip qSurfaceMask qCreateSpots checkHeatMap checkCreateAllSpots qRoughness
global qSelectedSurfaces  qHighResThickness qThicknessDistribution qCreateAllSpots
global checkOverall checkSurfaceMask checkThicknessDistribution checkCreateSpots
global checkHighResolution checkSelectedSurfaces qHighResolution qBioMass checkBiofilmHeatMap
global checkBioMass checkBioVolume qBioVolume qContinuityRatio checkContinuityRatio checkRoughness
global editboxSubstratumLayer checkIsBiofilm qIsBiofilm qAutoDetectSubstratumStart qBiofilmHeatMap
vImarisApplicationGlobal = vImarisApplication;

% Create figure
f = figure('units','pixels','position',[500,500,385,400],...
    'toolbar','none','menu','none',...
    'numbertitle','off','name','Biofilm Analysis');

%Option to merge all filamentSpots into a single Spots object
%checkbox Mean Thickness
checkOverall = uicontrol('Parent', f,'style','checkbox',...
    'position', [25,350,200,25],...
    'string', 'Overall Calculations');
checkSelectedSurfaces = uicontrol('Parent', f,'style','checkbox',...
    'position', [25,330,200,25],...
    'string', 'Individual Surface measurements');
%Text for check box1
text2 = uicontrol('Parent', f,'style','text',...
    'position',[35,300,280,25],...
    'string','Choose the Method of Measurement:',...
    'fontangle','italic','fontweight','bold','fontSize', 11);

checkSurfaceMask = uicontrol('Parent', f,'style','checkbox',...
    'position', [65,280,200,25],...
    'string', 'Surface Mask');
checkThicknessDistribution = uicontrol('Parent', f,'style','checkbox',...
    'position', [65,260,230,25],...
    'string', 'BioVolume (no gaps, to substratum)');
editboxSubstratumLayer = uicontrol('Parent',f,'style','edit',...
    'position', [95,240,20,20],...
    'string','1');
text3 = uicontrol('Parent', f,'style','text',...
    'position',[118,227,120,30],...
    'string','Substratum start slice',...
    'fontangle','italic');
checkAutoDetectSubstratumStart = uicontrol('Parent', f,'style','checkbox',...
    'position', [240,238,100,25],...
    'string', 'AutoDetect');
checkHighResolution = uicontrol('Parent', f,'style','checkbox',...
    'position', [65,215,280,25],...
    'string', 'Surface to surface edge (subvoxel)');

editboxReductionFactor = uicontrol('Parent',f,'style','edit',...
    'position', [95,195,30,20],...
    'string','5');
text2 = uicontrol('Parent', f,'style','text',...
    'position',[125,182,100,30],...
    'string','Reducing Factor',...
    'fontangle','italic');
checkCreateSpots = uicontrol('Parent', f,'style','checkbox',...
    'position', [225,192,170,25],...
    'string', 'show vertices');
checkCreateAllSpots = uicontrol('Parent', f,'style','checkbox',...
    'position', [315,192,200,25],...
    'string', 'ALL');
%checkNoReduction = uicontrol('Parent', f,'style','checkbox',...
%    'position', [225,192,200,25],...
%    'string', 'show spot vertices');
checkIsBiofilm = uicontrol('Parent', f,'style','checkbox',...
    'position', [95,168,250,25],...
    'string', 'AutoReduction for large surface(biofilm)');
text2 = uicontrol('Parent', f,'style','text',...
    'position',[20,139,230,30],...
    'string','Other Overall Measurements/Features',...
    'fontangle','italic','fontweight','bold');
checkHeatMap = uicontrol('Parent', f,'style','checkbox',...
    'position', [65,130,250,25],...
    'string', 'Surface mask thickness heatmap');
checkBiofilmHeatMap = uicontrol('Parent', f,'style','checkbox',...
    'position', [250,130,250,25],...
    'string', 'BioVolume heatmap');
checkBioVolume = uicontrol('Parent', f,'style','checkbox',...
    'position', [65,110,230,25],...
    'string', 'BioVolume');
checkBioMass = uicontrol('Parent', f,'style','checkbox',...
    'position', [65,90,200,25],...
    'string', 'Biomass(um3/um2)');
checkRoughness = uicontrol('Parent', f,'style','checkbox',...
    'position', [250,90,200,25],...
    'string', 'Roughness (StdDev)');
checkContinuityRatio = uicontrol('Parent', f,'style','checkbox',...
    'position', [65,70,230,25],...
    'string', 'Continuity Ratio');
text1 = uicontrol('Parent', f,'style','text',...
    'position',[75,40,240,28],...
    'string','Measurements: Mean thickness, Max thickness, Roughness Coefficent (Variance)',...
    'fontangle','italic');

% Create ANALYZE and Cancel pushbuttons
analyzebutton = uicontrol('Parent', f,'style','pushbutton','units','pixels',...
    'position',[50,10,100,30],'string','ANALYZE',...
    'callback',@Analyze,...
    'fontweight','bold','fontSize', 14);
cancelbutton = uicontrol('Parent', f,'style','pushbutton','units','pixels',...
    'position',[170,15,70,20],'string','CANCEL',...
    'callback',@Cancel);
%CheckAllbutton = uicontrol('Parent', f,'style','pushbutton','units','pixels',...
%                'position',[170,15,70,20],'string','Check All Analyze',...
%                'callback',@CheckAll);

%End the Script after clicking CANCEL button
function Cancel(varargin)
global f
delete (f);
return



% Analyze function - runs the XT script
function Analyze(varargin)
% GET data from the dialog menus
global vImarisApplicationGlobal checkBioMass checkBioVolume checkContinuityRatio
global f qOverall qFlip checkHeatMap qHeatMap qBioMass qContinuityRatio checkAutoDetectSubstratumStart
global qSelectedSurfaces qSurfaceMask qHighResolution qBioVolume qAutoDetectSubstratumStart qRoughness
global qThicknessDistribution qCreateSpots editboxReductionFactor qIsBiofilm checkCreateAllSpots checkRoughness
global checkOverall checkSurfaceMask checkThicknessDistribution checkIsBiofilm qCreateAllSpots qBiofilmHeatMap
global checkHighResolution checkSelectedSurfaces checkCreateSpots editboxSubstratumLayer checkBiofilmHeatMap
vImarisApplication = vImarisApplicationGlobal;
% connect to Imaris interface
if isa(vImarisApplication, 'Imaris.IApplicationPrxHelper')
    qOverall='NO';
    qFlip='NO';
    vSpotReductionFactor=str2double(get(editboxReductionFactor,'String'));
    if vSpotReductionFactor==0
        vSpotReductionFactor=1;
    end
    %%
    %checkOverallThickness
    if get(checkOverall,'Value')==1;
        qOverall='YES';
    else
        qOverall='NO';
    end
    %checkSelectedSurfaces
    if get(checkSelectedSurfaces,'Value')==1;
        qSelectedSurfaces='YES';
    else
        qSelectedSurfaces='NO';
    end
    %checkqThickness Distribution
    if get(checkThicknessDistribution,'Value')==1;
        qThicknessDistribution='YES';
    else
        qThicknessDistribution='NO';
    end
    %checkqAutoDetectSubstratum start
    if get(checkAutoDetectSubstratumStart,'Value')==1;
        qAutoDetectSubstratumStart='YES';
    else
        qAutoDetectSubstratumStart='NO';
    end
    %checkSurfaceMask
    if get(checkSurfaceMask,'Value')==1;
        qSurfaceMask='YES';
    else
        qSurfaceMask='NO';
    end
    %checkHigh Resolution
    if get(checkHighResolution,'Value')==1;
        qHighResolution='YES';
    else
        qHighResolution='NO';
    end
    %checkCreateSPots
    if get(checkCreateSpots,'Value')==1;
        qCreateSpots='YES';
    else
        qCreateSpots='NO';
    end
    %checkCreateAllSPots
    if get(checkCreateAllSpots,'Value')==1;
        qCreateAllSpots='YES';
    else
        qCreateAllSpots='NO';
    end
    %checkHeatMap
    if get(checkBiofilmHeatMap,'Value')==1;
        qBiofilmHeatMap='YES';
    else
        qBiofilmHeatMap='NO';
    end
    %checkHeatMap
    if get(checkHeatMap,'Value')==1;
        qHeatMap='YES';
    else
        qHeatMap='NO';
    end
    %checkBioMass
    if get(checkBioMass,'Value')==1;
        qBioMass='YES';
    else
        qBioMass='NO';
    end
    %checkRoughness
    if get(checkRoughness,'Value')==1;
        qRoughness='YES';
    else
        qRoughness='NO';
    end
    %checkBioVolume
    if get(checkBioVolume,'Value')==1;
        qBioVolume='YES';
    else
        qBioVolume='NO';
    end
    %checkContinuityRatio
    if get(checkContinuityRatio,'Value')==1;
        qContinuityRatio='YES';
    else
        qContinuityRatio='NO';
    end
    vSubtstratumStartSlice=str2double(get(editboxSubstratumLayer,'String'));
    %checkIsBiofilm
    if get(checkIsBiofilm,'Value')==1;
        qIsBiofilm='YES';
        vSpotReductionFactor=50;
    else
        qIsBiofilm='NO';
    end
    %clearvars checkbox1 checkbox2 checkbox3 checkbox4 checkbox5
    
    
    %%
    % the user has to create a scene with some surfaces
    vSurpassScene = vImarisApplication.GetSurpassScene;
    if isequal(vSurpassScene, [])
        msgbox('Please create some Surfaces in the Surpass scene!');
        return;
    end
    %type of measurement
    if isequal (qOverall,'NO') && isequal (qSelectedSurfaces,'NO')
        msgbox('Please choose method of Analysis!');
        return;
    end
    %type of measurement2
    if isequal (qSurfaceMask,'NO') && isequal (qThicknessDistribution,'NO') && isequal (qHighResolution,'NO')
        msgbox('Please choose method of Measurement!');
        return;
    end
    %%
    % get all Surpass surfaces names
    vSurfaces = vImarisApplication.GetFactory.ToSurfaces(vImarisApplication.GetSurpassSelection);
    vSurfacesSelected = vImarisApplication.GetFactory.IsSurfaces(vSurfaces);
    vNumberOfSurfaces=vSurfaces.GetNumberOfSurfaces;
    
    if vSurfacesSelected
        vScene = vSurfaces.GetParent;
    else
        msgbox('Please create a Surface object!');
        return;
    end
    tic
    vImarisDataSet = vImarisApplication.GetDataSet.Clone;
    
    %%
    %Get Image Data parameters
    vDataMin = [vImarisDataSet.GetExtendMinX, vImarisDataSet.GetExtendMinY, vImarisDataSet.GetExtendMinZ];
    vDataMax = [vImarisDataSet.GetExtendMaxX, vImarisDataSet.GetExtendMaxY, vImarisDataSet.GetExtendMaxZ];
    vDataSize = [vImarisDataSet.GetSizeX, vImarisDataSet.GetSizeY, vImarisDataSet.GetSizeZ];
    
    aSizeC = vImarisApplication.GetDataSet.GetSizeC;%Dataset number of channels
    aSizeT = vImarisApplication.GetDataSet.GetSizeT;%Dataset number of time points
    Xvoxelspacing = (vDataMax(1)-vDataMin(1))/vDataSize(1);
    Yvoxelspacing = (vDataMax(2)-vDataMin(2))/vDataSize(2);
    Zvoxelspacing = (vDataMax(3)-vDataMin(3))/vDataSize(3);
    vSmoothingFactor=Xvoxelspacing*2;
    
    if isequal (qHeatMap,'YES') | isequal (qBiofilmHeatMap,'YES')
        %Convert dataset to 32bit float
        vImarisDataSet = vImarisApplication.GetDataSet.Clone;
        vFloatType = vImarisDataSet.GetType.eTypeFloat;
        vImarisDataSet.SetType(vFloatType);
        % Create a new channel where the result will be sent
        if isequal (qHeatMap,'YES') && isequal (qBiofilmHeatMap,'YES')
            vImarisDataSet.SetSizeC(aSizeC + 2);
            vImarisDataSet.SetChannelName(aSizeC,['Surface thickness of ', char(vSurfaces.GetName)]);
            vImarisDataSet.SetChannelName(aSizeC+1,['Biofilm thickness of ', char(vSurfaces.GetName)]);
        elseif isequal (qBiofilmHeatMap,'YES')
            vImarisDataSet.SetSizeC(aSizeC + 1);
            vImarisDataSet.SetChannelName(aSizeC,['Biofilm thickness of ', char(vSurfaces.GetName)]);
        else
            vImarisDataSet.SetSizeC(aSizeC + 1);
            vImarisDataSet.SetChannelName(aSizeC,['Surface thickness of ', char(vSurfaces.GetName)]);
        end
        
    end
    
    %Test for selected surfaces
    vSurfacesSelected = vImarisApplication.GetFactory.ToSurfaces(vImarisApplication.GetSurpassSelection);
    vSurfacesSelectedIds = vSurfacesSelected.GetSelectedIds;
    vSurfacesSelectedIndices = vSurfacesSelected.GetSelectedIndices;
    vNumberOfSelectedSurfaces = size(vSurfacesSelectedIndices,1);
    vSurfaceIndicesAll=vSurfaces.GetIds;
    vSurfacesSelectedIds=vSurfacesSelectedIds(vSurfacesSelectedIds<90000000);
    vSurfaceIndicesAll=vSurfaces.GetIds;
    
    vSumSlice=zeros(vDataSize(1)*vDataSize(2),1);
    vSumSliceSubstratum=zeros(vDataSize(1)*vDataSize(2),1);
    vSumSliceSurface=zeros(vDataSize(1)*vDataSize(2),1);
    SubstratumTest=zeros(1,vDataSize(3));
    vSubstratumArea=zeros(1,vDataSize(3));
    vAllPopulationStdDev=[];
    vAllSampleStdDev=[];
    vAllVariance=[];
    vAllPopulationStdDevTD=[];
    vAllSampleStdDevTD=[];
    vAllVarianceTD=[];
    vSliceAllTD=[];
    count=0;
    if  isequal(qOverall, 'YES') && isequal(qThicknessDistribution, 'YES') | isequal(qSurfaceMask, 'YES')
        vProgressDisplay2 = waitbar(0, 'Biofilm Overall calculations...',...
            'position',[620,500,270,50]);
        if  isequal(qBioVolume, 'YES') | isequal(qBioMass, 'YES')
            %Get Original Stats from surfaces
            vAllSurfaceStatistics = vSurfaces.GetStatistics;
            vSurfacesStatNames = cell(vAllSurfaceStatistics.mNames);
            vSurfacesStatValues = vAllSurfaceStatistics.mValues;
            vSurfaceVolumeIndex=strmatch('Volume', vSurfacesStatNames);
            vSurfaceTimeIndex=strmatch('Time Index', vSurfacesStatNames);
            vSurfacesVolume = vSurfacesStatValues(vSurfaceVolumeIndex,:);
            vSurfacesTimeIndex = vSurfacesStatValues(vSurfaceTimeIndex,:);
        end
        
        for vTime = 0:aSizeT-1;
            % Get the mask DataSet
            vMaskDataSetPrimary = vSurfaces.GetMask( ...
                vDataMin(1), vDataMin(2), vDataMin(3), ...
                vDataMax(1), vDataMax(2), vDataMax(3), ...
                vDataSize(1), vDataSize(2), vDataSize(3), vTime);
            %Concatonate mean Biovolume per time point
            if  isequal(qBioVolume, 'YES') | isequal(qBioMass, 'YES')
                vOverallBioVolume(vTime+1) = sum(vSurfacesVolume(vSurfacesTimeIndex==vTime+1));
            end
            %Loop through volume, slice by slice, and sums the mass
            for vIndexZ = 1:vDataSize(3)
                if isequal (qFlip,'NO')
                    vSlice = vMaskDataSetPrimary.GetDataSubVolumeAs1DArrayBytes(0,0,1+vDataSize(3)-vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                else
                    vSlice = vMaskDataSetPrimary.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                end
                vSlice = vSlice == 1;
                vSumSlice=vSumSlice+double(vSlice);%generate summed projection of surface mask
                %Calculate the area of the surface contacting the substratum
                vSubstratumArea(1+vDataSize(3)-vIndexZ)=sum(vSlice)*Xvoxelspacing*Yvoxelspacing;
                %AutoDetect Substratum StartSlice
                SubstratumTest(1+vDataSize(3)-vIndexZ)=nnz(vSlice);
                
                %Thickness distribution
                CurrentSliceNumber=1+vDataSize(3)-vIndexZ;
                vSliceTD=double(vSlice);
                vSliceTD(vSliceTD==1)=CurrentSliceNumber;%-vSubstratumStartSlice;%Replace ones with zeros
                vSliceAllTD=[vSliceAllTD  vSliceTD];
                
                count=count+1;
                waitbar(count/((aSizeT-1)*2*vDataSize(3)), vProgressDisplay2);
            end
            
            %AutoDetect Substratum Start slice
            if isequal (qAutoDetectSubstratumStart, 'YES')
                %Find the min slice index where the substratum starts, when the slice
                %has nonzero elements
                vSubtstratumStartSlice=min(find(SubstratumTest));%-2;
                if vSubtstratumStartSlice<0
                    vSubtstratumStartSlice=1;
                end
                vSliceAllTD=vSliceAllTD-vSubtstratumStartSlice+1;
                vSliceAllTD(vSliceAllTD<0)=0;
            else%No autodetect but user selects the substratum start slice
                vSliceAllTD=vSliceAllTD-vSubtstratumStartSlice+1;
                vSliceAllTD(vSliceAllTD<0)=0;%set all non zero values from substratum correction to zero
            end
            
            %Calculate the Thickness distribution in microns (corrected for
            %substratum startslice
            vThicknessDistribution=max(vSliceAllTD,[],2)*Zvoxelspacing;
            
            %Concatonate mean Biovolume per time point
            if  isequal(qBioVolume, 'YES') | isequal(qBioMass, 'YES')
                vOverallBioVolume(vTime+1) = sum(vSurfacesVolume(vSurfacesTimeIndex==vTime+1));
                %Total volume of the biomass, this includes all voxels in surface
                %mask.  The substratum start slice is not a factor
                vTotalBioVolume(vTime+1)=sum(vSumSlice)*Zvoxelspacing*Xvoxelspacing*Yvoxelspacing;%Total volume of Biomass
                vTotalSubstratumArea(vTime+1)=vSubstratumArea(vSubtstratumStartSlice);
                %Calculate substratum area by a average first 5 slices from the
                %substratum start slice
                vMeanSubstratumArea(vTime+1)=mean(vSubstratumArea(vSubtstratumStartSlice:vSubtstratumStartSlice+5));
                
            end
            vSumSlice=vSumSlice*Zvoxelspacing;%Thickness at each max projected voxel
            vmaxThickness(vTime+1)=max(vSumSlice);%max thickness at each time point
            %loop through volume, slice by slice, set 1's to thickness measure
            for vIndexZ = 1:vDataSize(3)
                vSlice = vMaskDataSetPrimary.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                vSlice = vSlice == 1;
                %Define the volume from the substratum start Slice
                if vIndexZ>vSubtstratumStartSlice && isequal(qBioVolume, 'YES') | isequal(qBioMass, 'YES')
                    vSumSliceSubstratum=vSumSliceSubstratum+double(vSlice);%generate summed projection of surface mask from substratum
                end
                
                
                vSliceThickness=double(vSlice).*vSumSlice;%set thickness transform channel
                if isequal (qHeatMap,'YES') && isequal (qBiofilmHeatMap,'YES')
                    vImarisDataSet.SetDataSubVolumeAs1DArrayBytes(vSumSlice, ...
                        0,0,vIndexZ-1,aSizeC,vTime,vDataSize(1),vDataSize(2),1);
                    vImarisDataSet.SetDataSubVolumeAs1DArrayBytes(vThicknessDistribution, ...
                        0,0,vIndexZ-1,aSizeC+1,vTime,vDataSize(1),vDataSize(2),1);
                elseif isequal (qBiofilmHeatMap,'YES')
                    vImarisDataSet.SetDataSubVolumeAs1DArrayBytes(vThicknessDistribution, ...
                        0,0,vIndexZ-1,aSizeC,vTime,vDataSize(1),vDataSize(2),1);
                elseif isequal (qHeatMap, 'YES')
                    vImarisDataSet.SetDataSubVolumeAs1DArrayBytes(vSumSlice, ...
                        0,0,vIndexZ-1,aSizeC,vTime,vDataSize(1),vDataSize(2),1);
                end
                
                
                count=count+1;
                waitbar(count/((aSizeT)*2*vDataSize(3)), vProgressDisplay2);
            end
            if  isequal(qBioVolume, 'YES') | isequal(qBioMass, 'YES')
                vTotalBioVolumeSubstratum(vTime+1)=sum(vSumSliceSubstratum)*Zvoxelspacing*Xvoxelspacing*Yvoxelspacing;%Total volume of Biovolume from substratum
            end
            
            %Calculate mean thickness of surface mask and Variance per time point
            vSumSliceFinal = vSumSlice(vSumSlice~=0);
            if isempty (vSumSliceFinal)
                vAllMeanThickness(vTime+1)=0;
                vAllMaxThickness(vTime+1)=0;
                vAllPopulationStdDev=[vAllPopulationStdDev;0];
                vAllSampleStdDev=[vAllSampleStdDev;0];
                vAllVariance=[vAllVariance;0];
                vAllMeanThicknessDistribution(vTime+1)=0;
                vAllMaxThicknessDistribution(vTime+1)=0;
                vAllPopulationStdDevTD=[vAllPopulationStdDevTD;0];
                vAllSampleStdDevTD=[vAllSampleStdDevTD;0];
                vAllVarianceTD=[vAllVarianceTD;0];
                continue
            end
            vAllMeanThickness(vTime+1)=mean(vSumSliceFinal);
            vAllMaxThickness(vTime+1)=max(vSumSliceFinal);
            vMeanDifferenceSqrd=(vAllMeanThickness(vTime+1)-vSumSliceFinal).^2;
            %Calculate Population Standard Deviation
            SampleStdDev=sqrt(sum(vMeanDifferenceSqrd)/(numel(vSumSliceFinal)-1));
            %Calcualte Sample Standard Deviation
            PopulationStdDev=sqrt(sum(vMeanDifferenceSqrd)/numel(vSumSliceFinal));
            %Calculate Variance
            %vVariance=sum(vMeanDifferenceSqrd)/numel(vSumSliceFinal);
            vVariance=sum(vMeanDifferenceSqrd/mean(vSumSliceFinal))/numel(vSumSliceFinal);
            
            %Compile all data from STD and Variance
            vAllPopulationStdDev=[vAllPopulationStdDev;PopulationStdDev];
            vAllSampleStdDev=[vAllSampleStdDev;SampleStdDev];
            vAllVariance=[vAllVariance;vVariance];
            
            %Calculate Thicknes Distribution and Variance per time point
            vThicknessDistributionFinal=vThicknessDistribution(vThicknessDistribution~=0);
            vAllMeanThicknessDistribution(vTime+1)=mean(vThicknessDistributionFinal);
            vAllMaxThicknessDistribution(vTime+1)=max(vThicknessDistributionFinal);
            vMeanDifferenceSqrdTD=(vAllMeanThicknessDistribution(vTime+1)-vThicknessDistributionFinal).^2;
            
            %Calculate Sample Standard Deviation
            SampleStdDevTD=sqrt(sum(vMeanDifferenceSqrdTD)/(numel(vThicknessDistributionFinal)-1));
            %Calcualte Population Standard Deviation
            PopulationStdDevTD=sqrt(sum(vMeanDifferenceSqrdTD)/numel(vThicknessDistributionFinal));
            %Calculate Variance
            vVarianceTD=sum(vMeanDifferenceSqrdTD/mean(vThicknessDistributionFinal))/numel(vThicknessDistributionFinal);
            %Compile all data from STD and Variance
            vAllPopulationStdDevTD=[vAllPopulationStdDevTD;PopulationStdDevTD];
            vAllSampleStdDevTD=[vAllSampleStdDevTD;SampleStdDevTD];
            vAllVarianceTD=[vAllVarianceTD;vVarianceTD];
            
            clear vSliceThickness vMeanDifference vSumSliceFinal vSumSlice
            clear vThicknessDistributionFinal vThicknessDistribution vSliceAllTD
            vSumSlice=zeros(vDataSize(1)*vDataSize(2),1);
            vThicknessDistribution=zeros(vDataSize(1)*vDataSize(2),1);
            vSliceAllTD=[];
        end
        
        close(vProgressDisplay2);
        toc
        %Set display range for the Thickness Transform channel
        if isequal (qHeatMap,'YES') && isequal (qBiofilmHeatMap,'YES')
            vMaxThicknessinTime=max(vmaxThickness);
            vImarisDataSet.SetChannelRange(aSizeC,0,vMaxThicknessinTime);
            vMaxThicknessinTime=max(vAllMaxThicknessDistribution);
            vImarisDataSet.SetChannelRange(aSizeC+1,0,vMaxThicknessinTime);
        elseif isequal (qBiofilmHeatMap,'YES')
            vMaxThicknessinTime=max(vAllMaxThicknessDistribution);
            vImarisDataSet.SetChannelRange(aSizeC,0,vMaxThicknessinTime);
        elseif isequal (qHeatMap, 'YES')
            vMaxThicknessinTime=max(vmaxThickness);
            vImarisDataSet.SetChannelRange(aSizeC,0,vMaxThicknessinTime);
        end
        
        
        
        %Get Imaris version
        aVersion = char(vImarisApplication.GetVersion);
        aImarisFolderEnd = strfind(aVersion, ' [');
        if numel(aImarisFolderEnd) ~= 1
            msgbox('Invalid Imaris version')
            return
        end
        aImarisFolder = strtrim(aVersion(1:aImarisFolderEnd));
        aDelimiters = strfind(aImarisFolder, '-');
        if numel(aDelimiters) == 2
            aImarisFolder(aDelimiters(2)) = [];
            aImarisFolder(aDelimiters(1)) = ' ';
        end
        %Find/remove x64 --- fix for Imaris 9.8 and later
        aImarisNumberVersion = str2num(aImarisFolder(end-4:end-2));
        if aImarisNumberVersion <= 9.8
           aImarisFolder = erase(aImarisFolder,' x64');
        end
        
        
        if isequal(qHeatMap, 'YES') | isequal(qBiofilmHeatMap, 'YES') && isequal (qSurfaceMask, 'YES') | isequal (qThicknessDistribution, 'YES')%will only work for surface and observed area
            
            %Set Colot LUT
            if ismac
                % Manually generate the Spectrum LUT
                a=[5;5;5;4];
                a=repmat(a, 2,1);
                a=[a;5;5;4];
                a=[a;a;5;5;5;4;5];
                
                b=[-1024;-1280];
                b=repmat(b,9,1);
                b=[b;-1024];
                b=[b;b;b];
                
                c=[262144;327680];
                c=repmat(c,9,1);
                c=[c;262144];
                c=[c;c;c];
                
                d=[-4;-5];
                d=repmat(d,9,1);
                d=[d;-4];
                d=[d;d;d];
                
                e=[1024;1280];
                e=repmat(e,9,1);
                e=[e;1024];
                e=[e;e;e];
                SpectrumColorCreate=[a;b;c;d;e];
                colorTable=16711808;
                temp=16711808;
                for i =1:254
                    colorTable(i+1)=temp-SpectrumColorCreate(i);
                    temp=colorTable(i+1);
                end
                colorTable=colorTable';
                colorTable=[colorTable;255];
            else
                file=sprintf('C:\\Program Files\\Bitplane\\%s\\colorTables\\Spectrum.pal',aImarisFolder);
                RGB = load(file);
                % convert for Imaris
                colorTable = RGB(:, 1) + RGB(:, 2) * 256 + RGB(:, 3) * 256 * 256;
            end
            
            % set first color in LUT to black for the background
            colorTable(1)=0;
            
            if isequal (qHeatMap,'YES') && isequal (qBiofilmHeatMap,'YES')
                vImarisDataSet.SetChannelColorTable(aSizeC, colorTable, 0);
                vImarisDataSet.SetChannelColorTable(aSizeC+1, colorTable, 0);
                
            else
                vImarisDataSet.SetChannelColorTable(aSizeC, colorTable, 0);
                
            end
            vImarisApplication.SetDataSet(vImarisDataSet);
        end
    end
    
    
    %%
    vSurfaceThickness = vImarisApplication.GetFactory.ToSurfaces(vImarisApplication.GetSurpassSelection);
    %Create a new folder object for new Surface thickness
    %Newsurfaces = vImarisApplication.GetFactory;
    %result = Newsurfaces.CreateDataContainer;
    %result.SetName(sprintf('Biofilm analysis',char(vSurfaces.GetName)));
    %vSurfaceThickness=vSurfaces;%duplicate surfaces?
    %vSurfaceThickness.SetName(sprintf('Biofilm analysis - %s',char(vSurfaces.GetName)));
    %Add new surface to Surpass Scene
    %result.AddChild(vSurfaceThickness, -1);
    vImarisApplication.GetSurpassScene.AddChild(vSurfaceThickness, -1);
    
    %%
    if isequal (qCreateAllSpots, 'YES')& isequal (qHighResolution, 'NO')
        TotalAllvVertices=[];
        TotalAllTime=[];
        TotalAllRadii=[];
        vSpotsAll = vImarisApplication.GetFactory.CreateSpots;
        % create new group
        vNormalsGroup = vImarisApplication.GetFactory.CreateDataContainer;
        vNormalsGroup.SetName(['Vertices surface ', char(vSurfaces.GetName)]);
        % add group
        vSurpassScene.AddChild(vNormalsGroup, -1);
        
        for vSurfaceIndex2 = 0:vNumberOfSurfaces-1
            % There is one normal per vertex
            vNormalsAll = vSurfaces.GetNormals(vSurfaceIndex2);
            vVerticesAll = vSurfaces.GetVertices(vSurfaceIndex2);
            vNormalsAllTest = vNormalsAll(:,3)>0;
            vTestAll=[vNormalsAllTest vNormalsAllTest vNormalsAllTest];
            %vVerticesAll(~vTestAll)=0;
            %vVerticesAll(all(vVerticesAll==0,2),:)=[];
            
            if size(vVerticesAll(:,1),1)<50000
                idx = round((size(vVerticesAll(:,1),1)-0).*rand(round(size(vVerticesAll(:,1),1)*0.9),1) + 0);
                %rIndex = round((size(vVerticesAll(:,1),1)-0).*rand(100000,1) + 0);
            elseif size(vVerticesAll(:,1),1)<100000
                idx = round((size(vVerticesAll(:,1),1)-0).*rand(round(size(vVerticesAll(:,1),1)*0.85),1) + 0);
            elseif size(vVerticesAll(:,1),1)<300000
                idx = round((size(vVerticesAll(:,1),1)-0).*rand(round(size(vVerticesAll(:,1),1)*0.75),1) + 0);
            elseif size(vVerticesAll(:,1),1)>=300000
                idx = round((size(vVerticesAll(:,1),1)-0).*rand(round(size(vVerticesAll(:,1),1)*0.5),1) + 0);
            end
            idx(idx==0)=[];%Remove any zeros
            vVerticesAllFinal=vVerticesAll(idx,:);
            AllTime=ones(size(vVerticesAllFinal(:,1),1),1)*vSurfaces.GetTimeIndex(vSurfaceIndex2);
            AllRadii=zeros(size(vVerticesAllFinal(:,1),1),1)+0.01;
            
            TotalAllvVertices=[TotalAllvVertices;vVerticesAllFinal];
            TotalAllTime=[TotalAllTime;AllTime];
            TotalAllRadii=[TotalAllRadii;AllRadii];
            
            
            if vSurfaceIndex2 == vNumberOfSurfaces-1
                
                vSpotsAll.Set(TotalAllvVertices, TotalAllTime, TotalAllRadii);
                %vSpotsAll.SetColorRGBA(vRed + vGreen*256 + vBlue*256*256);
                vSpotsAll.SetColorRGBA(150 + 0*256 + 255*256*256);
                vSpotsAll.SetName(['All vertices ',num2str(vSurfaceIndex2+1)]);
                vNormalsGroup.AddChild(vSpotsAll, -1);
            end
            %Idx=[];
            clear idx vVerticesAllFinal
        end
        vNormalsGroup.SetVisible(0);
    end
    
    
    if isequal(qHighResolution, 'YES')
        if isequal (qCreateSpots, 'YES') | isequal (qCreateAllSpots, 'YES')
            % create spots factory
            vSpots1 = vImarisApplication.GetFactory.CreateSpots;
            vSpots2 = vImarisApplication.GetFactory.CreateSpots;
            vSpotsAll = vImarisApplication.GetFactory.CreateSpots;
            
            % create new group
            vNormalsGroup = vImarisApplication.GetFactory.CreateDataContainer;
            vNormalsGroup.SetName(['Vertices surface ', char(vSurfaces.GetName)]);
            % add group
            vSurpassScene.AddChild(vNormalsGroup, -1);
        end
        vAllSampleStdDev2=[];
        vAllPopulationStdDev2=[];
        vAllVariance2=[];
        TotalSpotsTop=[];
        TotalSpotsBottom=[];
        TotalTimeTop=[];
        TotalTimeBottom=[];
        TotalRadiiTop=[];
        TotalRadiiBottom=[];
        TotalAllvVertices=[];
        TotalAllTime=[];
        TotalAllRadii=[];
        
        
        count=1;
        vProgressDisplay = waitbar(0, 'Calculating Thickness (subvoxel)...');
        
        if ~isempty(vSurfacesSelectedIndices) & isequal (qSelectedSurfaces, 'YES')
            vNumberOfSurfaces=vNumberOfSelectedSurfaces;
        else
            
        end
        
        for vSurfaceIndex2 = 0:vNumberOfSurfaces-1
            if isempty(vSurfacesSelectedIndices) %& isequal (qOverall, 'YES')
                vSurfaceIndex=vSurfaceIndex2;
            else
                vSurfaceIndex=vSurfacesSelectedIndices(vSurfaceIndex2+1);
            end
            % get the vertices and normals
            vVertices1 = vSurfaces.GetVertices(vSurfaceIndex);
            vVertices2 = vSurfaces.GetVertices(vSurfaceIndex);
            if isequal (qCreateAllSpots, 'YES')
                vVerticesAll=vVertices1;
                vNormalsAll = vSurfaces.GetNormals(vSurfaceIndex2);
                vNormalsAllTest = vNormalsAll(:,3)>0;
                vTestAll=[vNormalsAll vNormalsAll vNormalsAll];
                vVerticesAll(~vTestAll)=0;
                vVerticesAll(all(vVerticesAll==0,2),:)=[];
            end
            
            % There is one normal per vertex
            vNormals1 = vSurfaces.GetNormals(vSurfaceIndex);
            vNormals2 = vSurfaces.GetNormals(vSurfaceIndex);
            
            %Loop for top and bottom vertices
            %Filter the normal Vectors
            vNormalsTopTest = vNormals1(:,3)<-0.75;
            vNormalsBottomTest = vNormals2(:,3)>0.55;
            
            vTestTop=[vNormalsTopTest vNormalsTopTest vNormalsTopTest];
            vTestBottom=[vNormalsBottomTest vNormalsBottomTest vNormalsBottomTest];
            vNormals1(~vTestTop)=0;
            vNormals1(all(vNormals1==0,2),:)=[];
            vNormalsStats1=vNormals1;
            vNormals2(~vTestBottom)=0;
            vNormals2(all(vNormals2==0,2),:)=[];
            vNormalsStats2=vNormals2;
            
            vVertices1(~vTestTop)=0;
            vVertices1(all(vVertices1==0,2),:)=[];
            vVertices1Final=vVertices1;
            vNumberOfVertices1 = size(vVertices1, 1);
            
            vVertices2(~vTestBottom)=0;
            vVertices2(all(vVertices2==0,2),:)=[];
            %vVertices2Final=vVertices2;
            vNumberOfVertices2 = size(vVertices2, 1);
            
            if vNumberOfVertices1<3 |vNumberOfVertices2<3| isempty(vNumberOfVertices1)|isempty(vNumberOfVertices2)
                vMeanThickness2(vSurfaceIndex+1)=0;
                vTimeIndices(vSurfaceIndex+1)=vSurfaces.GetTimeIndex(vSurfaceIndex);%Get time Index for each surface
                continue
            end
            %Filtering of spots to speed up analysis
            P=vVertices2;
            %idx=(randperm(n,nr))';
            %idx=randi([1 n],nr,1);
            
            if vNumberOfVertices1<1000
                out=vVertices1;
            else if isequal(checkIsBiofilm, 'YES')
                    nr=5000;
                    %idx=(randperm(n,nr))';
                    idx = round((size(vVertices1(:,1),1)-0).*rand(nr,1) + 0);
                    idx(idx==0)=[];%Remove zeros
                    out=vVertices1(idx(1:nr),:);%generate new random sample of all 1/3 total vertices
                else
                    n=length(vVertices1);
                    nr=round(n/vSpotReductionFactor);%Reduce the total number by 1/3
                    %idx=(randperm(n,nr))';
                    idx = round((size(vVertices1(:,1),1)-0).*rand(nr,1) + 0);
                    idx(idx==0)=[];%Remove zeros
                    out=vVertices1(idx(1:nr),:);%generate new random sample of all 1/3 total vertices
                end
            end
            vNumberOfVertices1 = size(out, 1);%number of vertices to sample
            
            % Test the number of vertices and bring up dialog to continue
            if vNumberOfVertices1>20000
                str=sprintf('Due to a large number of vertices in the surface\n\nContinuing with this analysis will take a long time \n\nPlease consider readjusting the Reducing factor\n\nOr USE the option for large surfaces (biofilm) ');
                choice=questdlg(str, 'WARNING !!!', 'Continue','Cancel','Reset Analysis','Cancel');
                
                %Handle Response
                switch choice
                    case 'Analyze'
                    case 'Cancel'
                        delete(f);
                        close(vProgressDisplay);
                        return;
                    case 'Reset Analysis'
                        close(vProgressDisplay);
                        return;
                end
            end
            
            
            %Create Test Spots to verify normal vertice choice
            if isequal (qCreateAllSpots, 'YES')
                if size(vVerticesAll(:,1),1)<100000
                    vVerticesAllFinal=vVerticesAll;
                else
                    rIndex = round((size(vVerticesAll(:,1),1)-0).*rand(100000,1) + 0);
                    vVerticesAllFinal=vVerticesAll(rIndex,:);
                end
                
                AllTime=ones(size(vVerticesAllFinal(:,1),1),1)*vSurfaces.GetTimeIndex(vSurfaceIndex);
                AllRadii=zeros(size(vVerticesAllFinal(:,1),1),1)+0.01;
                TotalAllvVertices=[TotalAllvVertices;vVerticesAllFinal];
                TotalAllTime=[TotalAllTime;AllTime];
                TotalAllRadii=[TotalAllRadii;AllRadii];
                
                
                if vSurfaceIndex2 == vNumberOfSurfaces-1
                    vSpotsAll.Set(TotalAllvVertices, TotalAllTime, TotalAllRadii);
                    %vSpotsAll.SetColorRGBA(vRed + vGreen*256 + vBlue*256*256);
                    vSpotsAll.SetColorRGBA(150 + 0*256 + 255*256*256);
                    vSpotsAll.SetName(['All vertices ',num2str(vSurfaceIndex+1)]);
                    vNormalsGroup.AddChild(vSpotsAll, -1);
                end
                
            end
            
            if isequal (qCreateSpots, 'YES')
                %Spots parameters
                vNbPointsPerNormal = 1;
                vPointRadius = 0.01;
                %Create spots
                vTimeTop=zeros(vNumberOfVertices1,1)+vSurfaces.GetTimeIndex(vSurfaceIndex);
                vSpotCountTop=size(vNumberOfVertices1,1);
                vRadiiTop=zeros(vNumberOfVertices1,1);
                
                %vSpotCountTop=size(vNormals1,1);
                vSpotCountBottom=size(vNormals2,1);
                %vTimeTop=zeros(vSpotCountTop,1);
                vTimeBottom=zeros(vSpotCountBottom,1)+vSurfaces.GetTimeIndex(vSurfaceIndex);
                %vRadiiTop=zeros(vSpotCountTop,1);
                vRadiiBottom=zeros(vSpotCountBottom,1);
                %vRadiiTop = vRadiiTop + 0.01;
                vRadiiBottom = vRadiiBottom + 0.01;
                
                vRed = 0;
                vGreen = 0;
                vBlue = 255;
                %vSpots1.Set(vVertices1Final, vTimeTop, vRadiiTop);
                qSpotsMerge='YES';
                if isequal (qSpotsMerge,'YES');
                    TotalSpotsTop=[TotalSpotsTop;out];
                    TotalSpotsBottom=[TotalSpotsBottom;vVertices2];
                    TotalTimeTop=[TotalTimeTop;vTimeTop];
                    TotalTimeBottom=[TotalTimeBottom;vTimeBottom];
                    TotalRadiiTop=[TotalRadiiTop;vRadiiTop];
                    TotalRadiiBottom=[TotalRadiiBottom;vRadiiBottom];
                    if vSurfaceIndex2 == vNumberOfSurfaces-1
                        vSpots1.Set(TotalSpotsTop, TotalTimeTop, TotalRadiiTop);
                        vSpots1.SetColorRGBA(vRed + vGreen*256 + vBlue*256*256);
                        vSpots1.SetName(['Vertices at Substratum ',num2str(vSurfaceIndex+1)]);
                        vNormalsGroup.AddChild(vSpots1, -1);
                        
                        vSpots2.Set(TotalSpotsBottom, TotalTimeBottom, TotalRadiiBottom);
                        vSpots2.SetColorRGBA(255 + vGreen*256 + vBlue*256*256);
                        vSpots2.SetName(['Vertices on upper edge ',num2str(vSurfaceIndex+1)]);
                        vNormalsGroup.AddChild(vSpots2, -1);
                    end
                else
                    vSpots1.Set(out, vTimeTop, vRadiiTop);
                    vSpots1.SetColorRGBA(vRed + vGreen*256 + vBlue*256*256);
                    vSpots1.SetName(['Spots on top ',num2str(vSurfaceIndex+1)]);
                    vNormalsGroup.AddChild(vSpots1, -1);
                    
                    vSpots2.Set(vVertices2, vTimeBottom, vRadiiBottom);
                    vSpots2.SetColorRGBA(255 + vGreen*256 + vBlue*256*256);
                    vSpots2.SetName(['Spots on bottom ',num2str(vSurfaceIndex+1)]);
                    vNormalsGroup.AddChild(vSpots2, -1);
                end
            end
            vertexcount=0;
            vertexcount2=1;
            vMinDistance=zeros(1,vNumberOfVertices1);
            
            %This will loop through each of the vertices1 on the bottom, filtered and
            %renamed 'out'.  the total vVertices1 have been reduced by 1/3
            for vIndex = 1:vNumberOfVertices1
                vertexcount=vertexcount+1;
                if vertexcount==1000
                    vertexcount2=vertexcount2+1;
                    vertexcount=0;
                    
                end
                %Q1=vVertices1(vIndex,:);
                %Q2=vVertices1(vIndex,:);
                %Q2(:,3)=1000;
                Q1=out(vIndex,:);
                Q2=out(vIndex,:);
                Q2(:,3)=1000;
                %Calculate norm(vVertices2-Q1)^2
                test=[vVertices2(:,1)-Q1(:,1) vVertices2(:,2)-Q1(:,2) vVertices2(:,3)-Q1(:,3)];
                N = (sqrt(sum(abs(test).^2,2))).^2;
                
                %Calculate (norm(Q2-Q1)^2)
                NORMQ2_Q1=norm(Q2-Q1);
                NORMQ2_Q1=(repmat(NORMQ2_Q1,size(vVertices2, 1),1)).^2;
                
                %Calculate dot(Q2-Q1,vVertices2-Q1)^2)
                DIFF=Q2-Q1;
                DIFF=repmat(DIFF,size(vVertices2, 1),1);
                DOT=sum(conj(DIFF).*test,2).^2;
                
                %d2 = sqrt(norm(Q2-Q1)^2*norm(vVertices21-Q1)^2-dot(Q2-Q1,vVertices21-Q1)^2)/norm(Q2-Q1);
                d1 = sqrt(NORMQ2_Q1.*N-DOT)./sqrt(NORMQ2_Q1);
                
                %Find minimum distance to the orthogonal line
                [M,I] = min(d1);
                %Measure distance from indexed spot
                
                %vMinDistance(vIndex)=sqrt((vVertices1(vIndex,1)-vVertices2(I,1))^2+(vVertices1(vIndex,2)-vVertices2(I,2))^2+(vVertices1(vIndex,3)-vVertices2(I,3))^2);
                %vMinDists(vIndex)=sqrt((vVertices1(1,1)-vVertices2(I,1))^2+(vVertices1(1,2)-vVertices2(I,2))^2+(vVertices1(1,3)-vVertices2(I,3))^2);
                vMinDistance(vIndex)=sqrt((out(vIndex,1)-vVertices2(I,1))^2+(out(vIndex,2)-vVertices2(I,2))^2+(out (vIndex,3)-vVertices2(I,3))^2);
                
                clear M I
                if vNumberOfVertices1<4000
                    waitbar (vSurfaceIndex/vNumberOfSurfaces);
                else
                    waitbar (vertexcount2/(vNumberOfVertices1/1000));
                end
                
            end
            %
            vMeanThickness2(count)=mean(vMinDistance);
            vMaxThickness2(count)=max(vMinDistance);
            vMeanDifferenceSqrd2=(vMeanThickness2(count)-vMinDistance).^2;
            %Calculate Population Standard Deviation
            vPopulationStdDev2=sqrt(sum(vMeanDifferenceSqrd2)/size(vMinDistance,2));
            %Calcualte Sample Standard Deviation
            vSampleStdDev2=sqrt(sum(vMeanDifferenceSqrd2)/(size(vMinDistance,2)-1));
            %Calculate Variance
            vVariance2=sqrt(sum(vMeanDifferenceSqrd2/mean(vMinDistance))/size(vMinDistance,2));
            
            
            
            vAllPopulationStdDev2(count)=vPopulationStdDev2;%all surfaces in this time point
            vAllSampleStdDev2(count)=vSampleStdDev2;%all surfaces in this time point
            vAllVariance2(count)=vVariance2;%all surfaces in this time point
            
            vTimeIndices(count)=vSurfaces.GetTimeIndex(vSurfaceIndex);%Get time Index for each surface
            waitbar (count/vNumberOfSurfaces);
            count=count+1;
            clear vTimeTop vTimeBottom vRadiiTop vRadiiBottom
            clear idx
            
        end
        %Loop time to generate time lapse for mean thickness and variance
        if isequal (qOverall, 'YES')
            for timeindex=1:aSizeT
                vAllMeanThicknessSurfaces2(timeindex) = mean(vMeanThickness2(vTimeIndices==timeindex-1));
                vAllMeanVarianceSurfaces2(timeindex) = mean(vAllVariance2(vTimeIndices==timeindex-1));
                vAllMaxThicknessSurfaces2(timeindex) = mean(vMaxThickness2(vTimeIndices==timeindex-1));
            end
            %Convert NAN to zeros
            vAllMeanThicknessSurfaces2(isnan(vAllMeanThicknessSurfaces2))=0;
            vAllMeanVarianceSurfaces2(isnan(vAllMeanVarianceSurfaces2))=0;
        end
        if isequal(qCreateSpots,'YES') | isequal(qCreateAllSpots,'YES')
            vNormalsGroup.SetVisible(0);
        end
        close(vProgressDisplay);
    end
    
    %%
    %Get individual mean thicknesses for each surface
    if  isequal(qSelectedSurfaces,'YES') && isequal(qSurfaceMask,'YES') | isequal(qThicknessDistribution,'YES')
        vProgressDisplay = waitbar(0, 'Biofilm: Thickness of Selected Surfaces',...
            'position',[620,400,270,50]);
        
        count2=0;
        %Test length of processing
        if (vNumberOfSelectedSurfaces*vDataSize(3))>500 | (size(vSurfaceIndicesAll(1))*vDataSize(3))>500 |(vNumberOfSurfaces*vDataSize(3))>500
            vImarisApplication.SetVisible(~vImarisApplication.GetVisible);
        end
        %%
        if isempty(vSurfacesSelectedIndices)
            %Processing all surfaces
            if isequal (qSurfaceMask, 'YES') | isequal (qThicknessDistribution, 'YES')
                for W=0:vNumberOfSurfaces-1
                    % Get the mask DataSet
                    vMaskDataSetSurface = vSurfaces.GetSingleMask(W,...
                        vDataMin(1), vDataMin(2), vDataMin(3), ...
                        vDataMax(1), vDataMax(2), vDataMax(3), ...
                        vDataSize(1), vDataSize(2), vDataSize(3));
                    %Loop through volume, slice by slice, and sums the mask
                    for vIndexZ = 1:vDataSize(3)
                        if isequal (qFlip,'NO')
                            vSlice = vMaskDataSetSurface.GetDataSubVolumeAs1DArrayBytes(0,0,1+vDataSize(3)-vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                        else
                            vSlice = vMaskDataSetSurface.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                        end
                        vSlice = vSlice == 1;
                        
                        vSumSliceSurface=vSumSliceSurface+double(vSlice);%generate summed projection of surface mask
                        
                        %Thickness distribution
                        if isequal (qThicknessDistribution, 'YES')
                            CurrentSliceNumber=1+vDataSize(3)-vIndexZ;
                            vSliceTD=double(vSlice);
                            if max(vSliceTD)>0,
                                vSliceTD(vSliceTD==1)=CurrentSliceNumber-vSubtstratumStartSlice;%Replace ones with zeros
                                vSliceAllTD=[vSliceAllTD  vSliceTD];
                            end
                        end
                        count2=count2+1;
                        waitbar(count2/(vNumberOfSurfaces*vDataSize(3)), vProgressDisplay);
                    end
                    vSliceAllTD(vSliceAllTD<0)=0;
                    
                    %Calculate Thickness Distribution and Variance per time point
                    if isequal (qThicknessDistribution, 'YES')
                        vThicknessDistributionSurfaces=max(vSliceAllTD,[],2)*Zvoxelspacing;%Calculate the Thickness distribution in microns
                        vThicknessDistributionFinal=vThicknessDistributionSurfaces(vThicknessDistributionSurfaces~=0);
                        vMeanThicknessDistribution(W+1)=mean(vThicknessDistributionFinal);
                        vMaxThicknessTD(W+1)=max(vThicknessDistributionFinal);%max thickness at each time point
                        vMeanDifferenceSqrdTD=(vMeanThicknessDistribution(W+1)-vThicknessDistributionFinal).^2;
                        
                        %Population Standard Deviation
                        vPopulationStdDevTD(W+1)=sqrt(sum(vMeanDifferenceSqrdTD)/numel(vThicknessDistributionFinal));
                        %Sample Standard Deviation
                        vSampleStdDevTD(W+1)=sqrt(sum(vMeanDifferenceSqrdTD)/(numel(vThicknessDistributionFinal)-1));
                        %Variance
                        vVarianceTD(W+1)=sum(vMeanDifferenceSqrdTD/mean(vThicknessDistributionFinal))/numel(vThicknessDistributionFinal);
                        
                    end
                    vSumSliceSurface=vSumSliceSurface*Zvoxelspacing;%Thickness at each max projected voxel
                    vMaxSurfaceThickness(W+1)=max(vSumSliceSurface);%max thickness at each time point
                    vSumSliceSurfaceFinal = vSumSliceSurface(vSumSliceSurface~=0);%remove zeros for mean calculation
                    vMeanSurfaceThickness(W+1)=mean(vSumSliceSurfaceFinal);%calculate mean thickness
                    vMeanDifferenceSqrdSurface=(vMeanSurfaceThickness(W+1)-vSumSliceSurfaceFinal).^2;
                    
                    %Population Standard Deviation
                    vPopulationStdDevSurface(W+1)=sqrt(sum(vMeanDifferenceSqrdSurface)/numel(vSumSliceSurfaceFinal));
                    %Sample Standard Deviation
                    vSampleStdDevSurface(W+1)=sqrt(sum(vMeanDifferenceSqrdSurface)/(numel(vSumSliceSurfaceFinal)-1));
                    %Variance
                    vVarianceSurface(W+1)=sum(vMeanDifferenceSqrdSurface/mean(vSumSliceSurfaceFinal))/numel(vSumSliceSurfaceFinal);
                    
                    vSumSliceSurface=zeros(vDataSize(1)*vDataSize(2),1);%rezero the Summedslice projection
                    vTimeIndices(W+1)=vSurfaces.GetTimeIndex(W);%Get time Index for each surface
                    clear vSumSliceFinal
                end
            end
            %info for new stats all surfaces
            vInd=1:vNumberOfSurfaces;
            vIds = vSurfaceIndicesAll;
            vUnits(vInd) = { char(vImarisApplication.GetDataSet.GetUnit) };
            vFactors(vInd) = {'Surface'};
            vFactors(2, vInd) = num2cell(vTimeIndices+1);
            vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
            close(vProgressDisplay);
            
            %%
            
            
        else
            %%    %Processing all selected surfaces
            for W=0:vNumberOfSelectedSurfaces-1
                vSurfaceIndex=vSurfacesSelectedIndices(W+1);
                % Get the mask DataSet
                vMaskDataSetSurface = vSurfaces.GetSingleMask(vSurfaceIndex,...
                    vDataMin(1), vDataMin(2), vDataMin(3), ...
                    vDataMax(1), vDataMax(2), vDataMax(3), ...
                    vDataSize(1), vDataSize(2), vDataSize(3));
                %Loop through volume, slice by slice, and sums the mask
                for vIndexZ = 1:vDataSize(3)
                    if isequal (qFlip,'NO')
                        vSlice = vMaskDataSetSurface.GetDataSubVolumeAs1DArrayBytes(0,0,1+vDataSize(3)-vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                    else
                        vSlice = vMaskDataSetSurface.GetDataSubVolumeAs1DArrayBytes(0,0,vIndexZ-1,0,0,vDataSize(1),vDataSize(2),1);
                    end
                    vSlice = vSlice == 1;
                    vSumSliceSurface=vSumSliceSurface+double(vSlice);%generate summed projection of surface mask
                    
                    %Thickness distribution
                    CurrentSliceNumber=1+vDataSize(3)-vIndexZ;
                    vSliceTD=double(vSlice);
                    if max(vSliceTD)>0,
                        vSliceTD(vSliceTD==1)=CurrentSliceNumber-vSubtstratumStartSlice;%Replace ones with zeros
                        vSliceAllTD=[vSliceAllTD  vSliceTD];
                    end
                    
                    count2=count2+1;
                    waitbar(count2/(vNumberOfSelectedSurfaces*vDataSize(3)), vProgressDisplay);
                end
                vSliceAllTD(vSliceAllTD<0)=0;
                %Calculate Thickness Distribution and Variance per time point
                vThicknessDistributionSurfaces=max(vSliceAllTD,[],2)*Zvoxelspacing;%Calculate the Thickness distribution in microns
                vThicknessDistributionFinal=vThicknessDistributionSurfaces(vThicknessDistributionSurfaces~=0);
                vMeanThicknessDistribution(W+1)=mean(vThicknessDistributionFinal);
                vMaxThicknessTD(W+1)=max(vThicknessDistributionFinal);%max thickness at each time point
                vMeanDifferenceSqrdTD=(vMeanThicknessDistribution(W+1)-vThicknessDistributionFinal).^2;
                %Population Standard Deviation
                vPopulationStdDevTD(W+1)=sqrt(sum(vMeanDifferenceSqrdTD)/numel(vThicknessDistributionFinal));
                %Sample Standard Deviation
                vSampleStdDevTD(W+1)=sqrt(sum(vMeanDifferenceSqrdTD)/(numel(vThicknessDistributionFinal)-1));
                %Variance
                vVarianceTD(W+1)=sum(vMeanDifferenceSqrdTD/mean(vThicknessDistributionFinal))/numel(vThicknessDistributionFinal);
                
                vSumSliceSurface=vSumSliceSurface*Zvoxelspacing;%Thickness at each max projected voxel
                vMaxSurfaceThickness(W+1)=max(vSumSliceSurface);%max thickness at each time point
                vSumSliceSurfaceFinal = vSumSliceSurface(vSumSliceSurface~=0);%remove zeros for mean calculation
                vMeanSurfaceThickness(W+1)=mean(vSumSliceSurfaceFinal);%calculate mean thickness
                vMeanDifferenceSqrdSurface=(vMeanSurfaceThickness(W+1)-vSumSliceSurfaceFinal).^2;
                
                %Population Standard Deviation
                vPopulationStdDevSurface(W+1)=sqrt(sum(vMeanDifferenceSqrdSurface)/numel(vSumSliceSurfaceFinal));
                %Sample Standard Deviation
                vSampleStdDevSurface(W+1)=sqrt(sum(vMeanDifferenceSqrdSurface)/(numel(vSumSliceSurfaceFinal)-1));
                %Variance
                vVarianceSurface(W+1)=sum(vMeanDifferenceSqrdSurface/mean(vSumSliceSurfaceFinal))/numel(vSumSliceSurfaceFinal);
                
                %vMeanThicknessOverall(vTime+1)=mean(vMeanSurfaceThickness);
                
                vSumSliceSurface=zeros(vDataSize(1)*vDataSize(2),1);%rezero the Summedslice projection
                vTimeIndices(W+1)=vSurfaces.GetTimeIndex(vSurfaceIndex);%Get time Index for each surface
                clear vSumSliceFinal
            end
            
            close(vProgressDisplay);
        end
        if (vNumberOfSelectedSurfaces*vDataSize(3))>500 | (size(vSurfaceIndicesAll(1))*vDataSize(3))>500 |(vNumberOfSurfaces*vDataSize(3))>500
            vImarisApplication.SetVisible(~vImarisApplication.GetVisible);
        end
        
        
    end
    if  isequal(qSelectedSurfaces,'YES')
        %Info for new stats from selected surfaces
        if isempty (vSurfacesSelectedIds)
            vInd=1:vSurfaces.GetNumberOfSurfaces;
            vIds=vSurfaces.GetIds;
        else
            vInd=1:vNumberOfSelectedSurfaces;
            vIds = vSurfacesSelectedIds;
        end
        
        
        vUnits(vInd) = { char(vImarisApplication.GetDataSet.GetUnit) };
        vFactors(vInd) = {'Surface'};
        vFactors(2, vInd) = num2cell(vTimeIndices+1);
        vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
        vFactorNames = {'Category','Time'};
        if isequal (qThicknessDistribution, 'YES')
            vNames(vInd) = {' Mean thickness (um)(no gaps to substratum)'};
            vSurfaces.AddStatistics(vNames, vMeanThicknessDistribution, vUnits, vFactors, vFactorNames, vIds);
            vNames(vInd) = {' Max thickness (um)(no gaps to substratum)'};
            vSurfaces.AddStatistics(vNames, vMaxThicknessTD, vUnits, vFactors, vFactorNames, vIds);
            vNames(vInd) = {' StdDev Population (no gaps to substratum)'};
            vSurfaces.AddStatistics(vNames, vPopulationStdDevTD, vUnits, vFactors, vFactorNames, vIds);
        end
        if  isequal(qSurfaceMask, 'YES')
            vNames(vInd) = {' Mean Thickness surfacemask'};
            vSurfaces.AddStatistics(vNames, vMeanSurfaceThickness, vUnits, vFactors, vFactorNames, vIds);
            vNames(vInd) = {' Population Standard Deviation (surface mask)'};
            vSurfaces.AddStatistics(vNames, vPopulationStdDevSurface, vUnits, vFactors, vFactorNames, vIds);
            vNames(vInd) = {' Max Thickness surfacemask'};
            vSurfaces.AddStatistics(vNames, vMaxSurfaceThickness, vUnits, vFactors, vFactorNames, vIds);
        end
        if  isequal(qHighResolution, 'YES')
            vNames(vInd) = {' Mean Thickness (subvoxel)'};
            vSurfaces.AddStatistics(vNames, vMeanThickness2, vUnits, vFactors, vFactorNames, vIds);
            vNames(vInd) = {' StdDev Population (subvoxel)'};
            vSurfaces.AddStatistics(vNames, vAllPopulationStdDev2, vUnits, vFactors, vFactorNames, vIds);
            vNames(vInd) = {' Max Thickness (subvoxel)'};
            vSurfaces.AddStatistics(vNames, vMaxThickness2, vUnits, vFactors, vFactorNames, vIds);
        end
        
        clear vNames;
    end
    %%
    
    
    for aIndexC = 0:aSizeC-1
        vImarisApplication.SetChannelVisibility(aIndexC,0);
    end
    %Setup values for overall statistic
    if isequal(qOverall,'YES')
        %vSurfaces.SetName(sprintf('Biofilm Analysis - %s',char(vSurfaces.GetName)));
        clear vInd vIds vUnits vFactors
        vInd=1:aSizeT;
        %vIds = vInd-1;
        vIds(vInd)=-1;
        vUnits(vInd) = {''};%{ char(vImarisApplication.GetDataSet.GetUnit) };
        Indices=1:aSizeT
        ;
        vFactors(vInd) = {'Overall'};
        vFactors(2, vInd) = num2cell(Indices);
        vFactors(2, vInd) = cellfun(@num2str, vFactors(2, vInd), 'UniformOutput', false);
        vFactorNames = {'Overall','Time'};
        
        
        if isequal(qBioMass, 'YES')&& isequal(qSurfaceMask, 'YES') | isequal(qThicknessDistribution, 'YES')
            %Overall Statistic for BioMass
            %vOverallBioMass = vTotalBioVolumeSubstratum/vSubstratumArea(vSubtstratumStartSlice);
            vOverallBioMass = vTotalBioVolumeSubstratum./vMeanSubstratumArea;
            vNames(vInd) = {' Biomass um3/um2'};
            vSurfaces.AddStatistics(vNames, vOverallBioMass, vUnits, vFactors, vFactorNames, vIds);
        end
        if isequal(qBioVolume, 'YES') && isequal(qSurfaceMask, 'YES') | isequal(qThicknessDistribution, 'YES')
            %Overall Statistic for BioVolume
            vNames(vInd) = {' BioVolume ImarisSurface (um3)'};
            vSurfaces.AddStatistics(vNames, vOverallBioVolume, vUnits, vFactors, vFactorNames, vIds);
            vNames(vInd) = {' BioVolume (surface mask)(um3)'};
            vSurfaces.AddStatistics(vNames, vTotalBioVolume, vUnits, vFactors, vFactorNames, vIds);
            vNames(vInd) = {' BioVolume (surface mask) from Substratum(um3)'};
            vSurfaces.AddStatistics(vNames, vTotalBioVolumeSubstratum, vUnits, vFactors, vFactorNames, vIds);
        end
        if isequal(qContinuityRatio, 'YES')  && isequal(qSurfaceMask, 'YES') && isequal(qThicknessDistribution, 'YES')
            %Overall Statistic for Surface Continuity
            vNames(vInd) = {' Biofilm Continuity Coefficient (1=no gaps in surface)'};
            vSurfaces.AddStatistics(vNames, vAllMeanThickness'./vAllMeanThicknessDistribution', vUnits, vFactors, vFactorNames, vIds);
        end
        if isequal(qSurfaceMask, 'YES')
            if isequal (qRoughness,'YES')
                %Overall Statistic for Population Standard Deviation
                vNames(vInd) = {' StdDev Population '};
                vSurfaces.AddStatistics(vNames, vAllPopulationStdDev, vUnits, vFactors, vFactorNames, vIds);
                %Overall Statistic for Sample Standard Deviation
                %vNames(vInd) = {' StdDev Sample '};
                %vSurfaces.AddStatistics(vNames, vAllSampleStdDev, vUnits, vFactors, vFactorNames, vIds);
                %Overall Statistic for Variance
                vNames(vInd) = {' Roughness Coefficient(surface mask)'};
                vSurfaces.AddStatistics(vNames, vAllVariance, vUnits, vFactors, vFactorNames, vIds);
            end
            %Stat for Overall mean thickness
            vNames(vInd) = {' Mean surface thickness bioVolume(um)'};
            vSurfaces.AddStatistics(vNames, vAllMeanThickness', vUnits, vFactors, vFactorNames, vIds);
            %Stat for Overall max thickness
            vNames(vInd) = {' Max surface thickness bioVolume(um)'};
            vSurfaces.AddStatistics(vNames, vAllMaxThickness', vUnits, vFactors, vFactorNames, vIds);
        end
        if isequal(qThicknessDistribution, 'YES')
            %Stat for Overall mean thickness Distribution
            vNames(vInd) = {' Mean thickness (um)(no gaps to substratum)'};
            vSurfaces.AddStatistics(vNames, vAllMeanThicknessDistribution, vUnits, vFactors, vFactorNames, vIds);
            %Stat for Overall max thickness Distribution
            vNames(vInd) = {' Max thickness  (um)(no gaps to substratum)'};
            vSurfaces.AddStatistics(vNames, vAllMaxThicknessDistribution, vUnits, vFactors, vFactorNames, vIds);
            if isequal (qRoughness,'YES')
                %Overall Statistic for Population Standard Deviation of Thickness Distribution
                vNames(vInd) = {' StdDev Population(no gaps to substratum)'};
                vSurfaces.AddStatistics(vNames, vAllPopulationStdDevTD, vUnits, vFactors, vFactorNames, vIds);
                %Overall Statistic for Sample Standard Deviation of Thickness Distribution
                %vNames(vInd) = {' StdDev Sample(no gaps to substratum)'};
                %vSurfaces.AddStatistics(vNames, vAllSampleStdDevTD, vUnits, vFactors, vFactorNames, vIds);
                %Stat for Overall Variance of Thickness Distribution
                vNames(vInd) = {' Roughness Coefficient(no gaps to substratum)'};
                vSurfaces.AddStatistics(vNames, vAllVarianceTD, vUnits, vFactors, vFactorNames, vIds);
            end
        end
        if isequal(qHighResolution, 'YES')
            %Stat for Overall mean thickness high Resolution
            vNames(vInd) = {' Mean thickness (subvoxel)(um)'};
            vSurfaces.AddStatistics(vNames, vAllMeanThicknessSurfaces2, vUnits, vFactors, vFactorNames, vIds);
            %Stat for Overall Max thickness of high Resolution
            vNames(vInd) = {' Max thickness (subvoxel)'};
            vSurfaces.AddStatistics(vNames, vAllMaxThicknessSurfaces2, vUnits, vFactors, vFactorNames, vIds);
            if isequal (qRoughness,'YES')
                %Overall Statistic for Population Standard Deviation
                vNames(vInd) = {' StdDev Population (subvoxel)'};
                vSurfaces.AddStatistics(vNames, vAllPopulationStdDev2, vUnits, vFactors, vFactorNames, vIds);
                %Overall Statistic for Sample Standard Deviation
                vNames(vInd) = {' StdDev Sample (subvoxel)'};
                vSurfaces.AddStatistics(vNames, vAllSampleStdDev2, vUnits, vFactors, vFactorNames, vIds);
                %Stat for Overall Variance of high Resolution
                vNames(vInd) = {' Roughness Coefficient(subvoxel)'};
                vSurfaces.AddStatistics(vNames, vAllMeanVarianceSurfaces2, vUnits, vFactors, vFactorNames, vIds);
            end
        end
        
    end
    
    %%
    
    %Rename new surface to Surpass Scene
    if isequal (qAutoDetectSubstratumStart, 'YES')
        vSurfaces.SetName(sprintf('Analyzed - %s subtratum startslice#%d',char(vSurfaces.GetName),vSubtstratumStartSlice));
        vImarisApplication.GetSurpassScene.AddChild(vSurfaces, -1);
    else
        vSurfaces.SetName(sprintf('Analyzed - %s subtratum startslice#%d',char(vSurfaces.GetName),vSubtstratumStartSlice));
        vImarisApplication.GetSurpassScene.AddChild(vSurfaces, -1);
    end
    
end
delete(f);





