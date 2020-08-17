%
%
%  Surface to Convex Hull for Imaris 8.3.0
%
%  Copyright Bitplane BPI2016
%  Matthew J. Gastinger, PhD
%
%  Installation:
%
%  - Copy this file into the XTensions folder in the Imaris installation directory.
%  - You will find this function in the Image Processing menu
%
%    <CustomTools>
%      <Menu>
%       <Submenu name="Surfaces Functions">
%        <Item name="Surface to Convex Hull" icon="Matlab" tooltip="Create a Surface which contains the convex hull of the Spots points.">
%          <Command>MatlabXT::XT_MJG_SurfaceToConvexHull(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="Surface Convex Hull">
%            <Command>MatlabXT::XT_MJG_SurfaceToConvexHull(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
% 
%
%  Description:
%
%   Create a Surface which contains the convex hull of the surface points.
%
%

function XT_MJG_SurfaceToConvexHull(aImarisApplicationID)

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

% get the surfaces
vSurfaces = vImarisApplication.GetFactory.ToSurfaces(vImarisApplication.GetSurpassSelection);
if ~vImarisApplication.GetFactory.IsSurfaces(vSurfaces)  
  msgbox('Please select some Surfaces!');
  return;
end

% create new group
%vNormalsGroup = vImarisApplication.GetFactory.CreateDataContainer;
%vNormalsGroup.SetName(['Spots on surface of ', char(vSurfaces.GetName)]);

% add group
%vSurpassScene.AddChild(vNormalsGroup, -1);

%Spots parameters
vNbPointsPerNormal = 1;
vPointRadius = 0.01;
%vRed = 0;
%vGreen = 0;
%vBlue = 255;

vFactory = vImarisApplication.GetFactory;
vNumberOfSurfaces = vSurfaces.GetNumberOfSurfaces;
result3 = vFactory.CreateDataContainer;
result3.SetName('Convex Hulls');


%%
for vSurfaceIndex = 0:vNumberOfSurfaces-1
	% create spot
    %vSpots = vImarisApplication.GetFactory.CreateSpots;
    % get the vertices and normals
	vVertices = vSurfaces.GetVertices(vSurfaceIndex);
	% There is one normal per vertex
	vNormals = vSurfaces.GetNormals(vSurfaceIndex);
   
   
%%
    %Filter the normal Vectors
    vNormalsTopTest = vNormals;
    
    vTest=[vNormalsTopTest vNormalsTopTest vNormalsTopTest];
    vNormals(~vTest)=0;
    vNormals(all(vNormals==0,2),:)=[];
    vNormalsStats=vNormals;
    
    vVertices(~vTest)=0;
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
    vNormals = vPointsSequencePos + vVertices; % Vertice Positions
    
    nr=5000;
    %idx=(randperm(n,nr))';
    idx = round((size(vNormals(:,1),1)-0).*rand(nr,1) + 0);
    idx(idx==0)=[];%Remove zeros
    vNormalsReduced=vNormals(idx,:);%Reduced set of vertices    
    %vNormalsReduced=vNormals(idx(1:nr),:);

vCountDegeneratedData = 0;
% create result object
vSurfaceHull = vFactory.CreateSurfaces;
vIndexT = vSurfaces.GetTimeIndex(vSurfaceIndex);
 
  try
    % yes, matlab does most of the work!
    %vConvexHull = convhull(double(vSpotsXYZ(:,1)',double(vSpotsXYZ(:,2)')));
    vConvexHull = convhulln(double(vNormalsReduced));
  catch er
      er.message; % suppress warning
      vCountDegeneratedData = vCountDegeneratedData + 1;
      %continue
  end
  % select the necessary points
  vNumberOfPoints = size(vNormalsReduced, 1);
  vPoints = false(vNumberOfPoints, 1);
  vPoints(vConvexHull(:)) = true;
  vPoints = find(vPoints);
  vVertices = vNormalsReduced(vPoints, :);

  % remap vertex indices to our selection
  % and reorder triangle vertices (clockwise to counter)
  vPointsMap = zeros(vNumberOfPoints, 1);
  vPointsMap(vPoints) = 1:numel(vPoints);
  vTriangles = vPointsMap(vConvexHull(:, [1, 3, 2]));

  % calculate normals (do not normalize them, imaris will do it)
  % follow rays from center to vertices
  
  vNbTriangles = size(vTriangles,1);
  vTrianglesNormals = zeros(size(vTriangles));
  vNbVertices = size(vVertices,1);
  vNormals2 = zeros(size(vVertices));
  
  % Vectors containing the first, second and third vertices of the triangles
  vTriangleVertices1 = vVertices(vTriangles(:,1),:);
  vTriangleVertices2 = vVertices(vTriangles(:,2),:);
  vTriangleVertices3 = vVertices(vTriangles(:,3),:);
  
  % Calculate the cross product for each triangle --> give the normals per triangle
  vTrianglesNormals = cross(vTriangleVertices2-vTriangleVertices1,vTriangleVertices3-vTriangleVertices1);
  
  
  % Pair triangle number / vertice number
  vNbTrianglesElements = 1:numel(vTriangles);  
  vTrianglesElementsIndices = mod(vNbTrianglesElements-1,vNbTriangles)+1;
  
  % Map representing in which triangle each vertice appears ("1" if it appears and "0" otherwise)
  % The third dimension is used to store a normal (on X, Y and Z axis)
  vMappingTrianglesVertices = zeros(vNbVertices,vNbTriangles, 3);
  vMappingTrianglesVertices(sub2ind(size(vMappingTrianglesVertices),vTriangles(vNbTrianglesElements),vTrianglesElementsIndices)) = 1;
  
  % Copy the same map for the Y and the Z
  vMappingTrianglesVertices(:,:,2) = vMappingTrianglesVertices(:,:,1);
  vMappingTrianglesVertices(:,:,3) = vMappingTrianglesVertices(:,:,1);
  
  % Set the triangle normal for each triangle in which the vertice appears (on X, Y and Z axis)
  vMappingTrianglesVertices(:,:,1) = bsxfun(@times, vMappingTrianglesVertices(:,:,1), vTrianglesNormals(:,1)');
  vMappingTrianglesVertices(:,:,2) = bsxfun(@times, vMappingTrianglesVertices(:,:,2), vTrianglesNormals(:,2)');
  vMappingTrianglesVertices(:,:,3) = bsxfun(@times, vMappingTrianglesVertices(:,:,3), vTrianglesNormals(:,3)');
  
  % Sum all the triangle normals to have a mean normal
  vNormals3D = sum(vMappingTrianglesVertices,2);

  % Reshape into a 2D matrice to match with Imaris Interface
  vNormals2(:,1) = vNormals3D(:,:,1);
  vNormals2(:,2) = vNormals3D(:,:,2);
  vNormals2(:,3) = vNormals3D(:,:,3);

  vTriangles = vTriangles - 1;
  vSurfaceHull.AddSurface(vVertices, vTriangles, vNormals2, vIndexT);
  vSurfaceHull.SetName(['Convex Hull of ', char(vSurfaces.GetName)]);
  vSurfaceHull.SetColorRGBA(vSurfaces.GetColorRGBA);
  result3.AddChild(vSurfaceHull, -1);
  vImarisApplication.GetSurpassScene.AddChild(result3, -1);
end



if vCountDegeneratedData > 0
  msgbox(['Could not create convex hull for some of the filaments. ', ...
    'Filament points must not lie on a plane to build a valid convex hull.'])
end
