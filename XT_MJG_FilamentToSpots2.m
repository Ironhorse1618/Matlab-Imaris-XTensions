%
%
%  filament to Spots
%
%  Copyright Bitplane BPI 2015
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
%        <Item name="Filaments To Spots" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_FilamentToSpots2(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpFilaments">
%          <Item name="Filaments To Spots" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_FilamentToSpots2(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
% 
%
%  Description:
%       This XTension will do several things:
%       1)Place spots at each point along the filament, allowing for
%           visualizing and measuring the instensity along the dendrites



function XT_MJG_FilamentToSpots2(aImarisApplicationID)

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

%Create a new folder object for new Spots
result = vFactory.CreateDataContainer;
result.SetName('Filament displayed as Spots');
vNewSpots1 = vImarisApplication.GetFactory.CreateSpots;
vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;

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
%
vDendritePositionsWorking=[];
vDendriteRadiusWorking=[];
%vRGBA = vFilaments.GetColorRGBA
vCount = vFilaments.GetNumberOfFilaments;
for FilamentIndex=0:vCount-1
    vFilamentsIndexT = vFilaments.GetTimeIndex(FilamentIndex);
    vFilamentsXYZ = vFilaments.GetPositionsXYZ(FilamentIndex);
    vFilamentsRadius = vFilaments.GetRadii(FilamentIndex);
    vTypes = vFilaments.GetTypes(FilamentIndex);
    vFilamentEdges = vFilaments.GetEdges(FilamentIndex);
    vNewSpots1 = vImarisApplication.GetFactory.CreateSpots;
    vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;

%Logical arguement to identify spots in segment
    for c = 0:1 %loop twice for each filamant, 0=dendrite 1=spine, and generate a 
        %separate spots object for dendrites and spines.
        vSpotsT= vTypes == c;    
        vSpotsIndex = find(vSpotsT');
        vDendritePositionsWorking=vFilamentsXYZ(vSpotsIndex,:);
        vDendriteRadiusWorking=vFilamentsRadius(vSpotsIndex,:);
        vDendritevTypesWorking=vTypes(vSpotsIndex,:);
        vNumberOfFilamentPoints= numel(vDendriteRadiusWorking);%including starting point
        vTimeIndex(1:vNumberOfFilamentPoints)=vFilamentsIndexT;%for filament spot creation
        if c==0 %Do first look for dendrites
            vNewSpots1.Set(vDendritePositionsWorking, vTimeIndex, vDendriteRadiusWorking);
            vNewSpots1.SetName([char(vFilaments.GetName), '_dendrites_ID:' num2str(FilamentIndex+1)]);
            vNewSpots1.SetColorRGBA((rand(1, 1)) * 256 * 256 * 256 );
            result.AddChild(vNewSpots1, -1);
        elseif ~isempty(vDendritevTypesWorking)%test second loop if spines exist, if not not make spine spots object
            vNewSpots2.Set(vDendritePositionsWorking, vTimeIndex, vDendriteRadiusWorking);
            vNewSpots2.SetName([char(vFilaments.GetName), '_spines_ID:' num2str(FilamentIndex+1)]);
            vNewSpots2.SetColorRGBA((rand(1, 1)) * 256 * 256 * 256 );
            result.AddChild(vNewSpots2, -1);
        end
        
        vImarisApplication.GetSurpassScene.AddChild(result, -1);

        clearvars vTimeIndex
    end
    
end

 