%
%
%  Visualizes Surface of surface via spot at each vertex
%
%  Copyright Bitplane BPI 2014
%
%
%  Installation:
%
%  - Copy this file into the XTensions folder in the Imaris installation directory
%  - You will find this function in the Image Processing menu
%
%    <CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Visualize Surface Spots" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_VisualizeSurfaceSpots(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Visualize Surface Spots" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_VisualizeSurfaceSpots(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
% 
%
%  Description:
%   
%   Places a spot at each vertex on the surface of the surface.  This
%   XTension only identifies those vertices with a normal z-vector < -0.25.
%   Thus only identifying the top surface
% 
%

function XT_MJG_VisualizeSurfaceSpots(aImarisApplicationID)

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
%make Imaris invisible
%vImarisApplication.SetVisible(~vImarisApplication.GetVisible);
vProgressDisplay = waitbar(0, 'Creating normals');

% the user has to create a scene with some surfaces
vSurpassScene = vImarisApplication.GetSurpassScene;
if isequal(vSurpassScene, [])
  msgbox('Please create some Surfaces in the Surpass scene!');
  return;
end

% get the surfaces
vSurfaces = vImarisApplication.GetFactory.ToSurfaces(vImarisApplication.GetSurpassSelection);
if ~vImarisApplication.GetFactory.IsSurfaces(vSurfaces)  
  msgbox('Please select some Surfaces!');
  return;
end

% create new group
vNormalsGroup = vImarisApplication.GetFactory.CreateDataContainer;
vNormalsGroup.SetName(['Spots on surface of ', char(vSurfaces.GetName)]);

% add group
vSurpassScene.AddChild(vNormalsGroup, -1);

%Spots parameters
vNbPointsPerNormal = 1;
vPointRadius = 0.01;
vRed = 0;
vGreen = 0;
vBlue = 255;


vNumberOfSurfaces = vSurfaces.GetNumberOfSurfaces;



%%
for vSurfaceIndex = 0:vNumberOfSurfaces-1
	% create spot
    vSpots = vImarisApplication.GetFactory.CreateSpots;

    % get the vertices and normals
	vVertices = vSurfaces.GetVertices(vSurfaceIndex);
	% There is one normal per vertex
	vNormals = vSurfaces.GetNormals(vSurfaceIndex);
   
   
%%
    %Filter the normal Vectors
    vNormalsTopTest = vNormals;%(:,3)<-0.25;
    
    %vTest=[vNormalsTopTest vNormalsTopTest vNormalsTopTest];
    %vNormals(~vTest)=0;
    %vNormals(all(vNormals==0,2),:)=[];
    vNormalsStats=vNormals;
    
    %vVertices(~vTest)=0;
    vVertices(all(vVertices==0,2),:)=[];
    vVerticesFinal=vVertices;
    vNumberOfVertices = size(vVertices, 1);

    vSpotCountTop=size(vNormals,1);
    vTimeTop=zeros(vSpotCountTop,1);
    vRadiiTop=zeros(vSpotCountTop,1);
    vRadiiTop = vRadiiTop + 0.01;
    vNbPointsPerNormal=1;
    
	% Normalise the normal
	% Norme = sqrt(x^2 + y^2 + z^2)
	vNormalNorme = vNormals.^2;
	% Double transpose of the matrix to sum each row and no each column
	vNormalNorme = (sum(vNormalNorme'))';
	vNormalNorme = sqrt(vNormalNorme);
	% Divide each row of the matrix by the vector element at the same row
	vNormals = bsxfun(@rdivide, vNormals, vNormalNorme);
	
	% Duplicate each row on five rows
	% vNumberOfVertices is used because the number of vertices = the number of normals
	vNormals=vNormals(ceil((1:vNbPointsPerNormal*vSpotCountTop)/vNbPointsPerNormal), :);
	
	% Creation of a points sequence for one normal
	vPointsSequence = (0:vNbPointsPerNormal-1)';
	vPointsSequence = vPointsSequence * vPointRadius * 2;
	% Repeat this sequence for each normal
	vPointsSequence = repmat(vPointsSequence, vSpotCountTop, 1);
	vPointsSequencePos = bsxfun(@times, vNormals, vPointsSequence);
	
	% Duplicate each row on five rows
	vVertices=vVertices(ceil((1:vNbPointsPerNormal*vSpotCountTop)/vNbPointsPerNormal), :);
    vNormals = vPointsSequencePos + vVertices;
    
    vSpots.Set(vNormals, vTimeTop, vRadiiTop);
    vSpots.SetColorRGBA(vRed + vGreen*256 + vBlue*256*256);
    vSpots.SetName(['Spots on surface ',num2str(vSurfaceIndex+1)]);
    
    %Create a new NormalVector Statistic
    vPointsCount = size(vNormalsStats,1);
    vNormalsStats=vNormalsStats(:);
    vStatsCount = numel(vNormalsStats);
    
    vNames = cell(vStatsCount, 1);
    vNames(1:vPointsCount) = {'Normal X'};
    vNames(vPointsCount+1:vPointsCount*2) = {'Normal Y'};
    vNames(vPointsCount*2+1:vStatsCount) = {'Normal Z'};
    vIds = [1:vPointsCount, 1:vPointsCount, 1:vPointsCount]' - 1;
    vUnits = cell(vStatsCount, 1);
    vUnits(:) = {'um'};
    vFactors = cell(3, vStatsCount);
    vFactors(1, :) = {'Spot'};
    vFactors(2, :) = {'1'};
    vFactors(3, :) = {'Normal Vector'};
    vFactorNames = {'Category','Time','Collection'};
    
    vSpots.AddStatistics(vNames, vNormalsStats, vUnits, vFactors, vFactorNames, vIds);
    
    
	vNormalsGroup.AddChild(vSpots, -1);
    

    
    waitbar(vSurfaceIndex/(vNumberOfSurfaces-1));
end
close(vProgressDisplay);
%make Imaris visible
%vImarisApplication.SetVisible(~vImarisApplication.GetVisible);



