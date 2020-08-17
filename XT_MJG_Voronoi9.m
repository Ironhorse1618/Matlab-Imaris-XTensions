%2D Voronoi diagram

%Written by Matthew J. Gastinger, Bitplane Advanced Application Scientist.
%Jan 2019.
%
%<CustomTools>
%      <Menu>
%       <Submenu name="Spots Functions">
%        <Item name="2D Voronoi Diagram" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_Voronoi9(%i)</Command>
%        </Item>
%       </Submenu>
%       <Submenu name="Surface Functions">
%        <Item name="2D Voronoi Diagram" icon="Matlab">
%          <Command>MatlabXT::XT_MJG_Voronoi9(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSpots">
%          <Item name="2D Voronoi Diagram" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_Voronoi9(%i)</Command>
%          </Item>
%        </SurpassComponent>
%        <SurpassComponent name="bpSurfaces">
%          <Item name="2D Voronoi Diagram" icon="Matlab">
%            <Command>MatlabXT::XT_MJG_Voronoi9(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%Description: 
%This XTension takes the 2D spot positions and generates a
%Voronoi diagram in Matlab.  It will generate a new channel in Imaris with
%the Voronoi polygons lines.
%
%For each time point, a random distribution of spots (same number as the
%original) generates a random Voronoi diagram
%
%Using Surfaces it generates a 2D surface of the polygons.  You must be in
%the "Surface view" to see them properly


function XT_MJG_Voronoi9(aImarisApplicationID)

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


%%
%Get Image Data parameters
vDataMin = [vImarisApplication.GetDataSet.GetExtendMinX, vImarisApplication.GetDataSet.GetExtendMinY, vImarisApplication.GetDataSet.GetExtendMinZ];
vDataMax = [vImarisApplication.GetDataSet.GetExtendMaxX, vImarisApplication.GetDataSet.GetExtendMaxY, vImarisApplication.GetDataSet.GetExtendMaxZ];
vDataSize = [vImarisApplication.GetDataSet.GetSizeX, vImarisApplication.GetDataSet.GetSizeY, vImarisApplication.GetDataSet.GetSizeZ];
aSizeC = vImarisApplication.GetDataSet.GetSizeC;
aSizeT = vImarisApplication.GetDataSet.GetSizeT;

Xvoxelspacing= (vDataMax(1)-vDataMin(1))/vDataSize(1);
Yvoxelspacing= (vDataMax(2)-vDataMin(2))/vDataSize(2);
Zvoxelspacing = (vDataMax(3)-vDataMin(3))/vDataSize(3);
vSmoothingFactor=Xvoxelspacing*2;

%Spot information
if isempty (vSurfaces)
    vSpotsRadius = vSpots.GetRadii;
    vObjectPositionXYZ = vSpots.GetPositionsXYZ;
    vObjectTime = vSpots.GetIndicesT;
    vTotalNumberOfObjects = numel(vSpotsRadius);
else
    vTotalNumberOfObjects=numel(vSurfaces.GetIds);
    for SurfaceIndex=0:vTotalNumberOfObjects-1
        vObjectPositionXYZ(SurfaceIndex+1,:)=vSurfaces.GetCenterOfMass(SurfaceIndex);
        vObjectTime(SurfaceIndex+1,:)=vSurfaces.GetTimeIndex(SurfaceIndex);
    end
end
%Add 2 new channels to Imaris
vImarisDataSet = vImarisApplication.GetDataSet.Clone;
vNumberOfChannels = vImarisDataSet.GetSizeC;
vImarisDataSet.SetSizeC(vNumberOfChannels + 2);
vIndexZ=1;

for vTime=0:max(vObjectTime)
    %Find spots per each time point
    vValid = find(vObjectTime == vTime);
    vObjectsWorking=vObjectPositionXYZ(vValid,:);
    vNumberOfSpots=size(vObjectsWorking,1);

    %Set positions in Pixels relative to image extensions
    x=round((double(vObjectsWorking(:,1)-vDataMin(1)))/Xvoxelspacing,0);
    %Invert Yposition to match image in Imaris
    y=round((vDataMin(2)+vDataMax(2)-double(vObjectsWorking(:,2))-vDataMin(2))/Yvoxelspacing,0);
    
    %Generate Random spots in field
    xRandom=(randi([0 vDataSize(1)-1],1,vNumberOfSpots))';
    yRandom=(randi([0 vDataSize(2)-1],1,vNumberOfSpots))';
    
    %Set image dimensions for ROI Voronoi mask
    img=[vDataSize(2) vDataSize(1)];
    %Calculate mask of data
    mask = voronoi2mask(x',y',img,vDataSize);
    %mask = label2rgb(mask, 'jet', 'c', 'shuffle');
    maskRand = voronoi2mask(xRandom',yRandom',img,vDataSize);
    
    %set borders to 255
    mask(mask>0)=255;
    maskRand(maskRand>0)=255;
    
%% show results
%     subplot(2,1,1)
%     imagesc(img);colormap('gray');
%     hold on;
%     h = voronoi(x,y);
%     set(h(:),'Color',[0 1 0]);
%     axis image;
%     title('original image with voronoi diagram');
%     
%     subplot(2,1,2)
%     imshow(mask);
%     hold on;
%     h = voronoi(x,y);
%     set(h(:),'Color',[0 0 0]);
%     axis image;
%     title('output image with voronoi diagram');
%%    
    %convert mask to grayscale
    %maskGreyscale=rgb2gray(mask);
    %maskGreyscale=mask;
    %convert to 1D array
    %vPolygons=(reshape(mask(end:-1:1,:).', 1, []))';
    %vPolygonsRandom=(reshape(maskRand(end:-1:1,:).', 1, []))';
    
    %invert channels
    vPolygons=double(~logical((reshape(mask(end:-1:1,:).', 1, []))'));
    vPolygons(vPolygons==1)=255;    
    vPolygonsRandom=double(~logical((reshape(maskRand(end:-1:1,:).', 1, []))'));
    vPolygonsRandom(vPolygons==1)=255;    

    
    vImarisDataSet.SetDataSubVolumeAs1DArrayShorts(vPolygons,0,0,vIndexZ-1,...
        vNumberOfChannels,vTime,vDataSize(1),vDataSize(2),1);
    vImarisDataSet.SetDataSubVolumeAs1DArrayShorts(vPolygonsRandom,0,0,vIndexZ-1,...
        vNumberOfChannels+1,vTime,vDataSize(1),vDataSize(2),1);
    clear vPolygons vPolygonsRandom vSpotsWorking
end
% Create a new channel where the result will be sent
vImarisDataSet.SetChannelName(vNumberOfChannels,['Voronoi polygons']);
vImarisDataSet.SetChannelColorRGBA(vNumberOfChannels, 255*256*256);
vImarisDataSet.SetChannelName(vNumberOfChannels+1,['Random Voronoi polygons']);
vRGBA=[0 255 0 0];%green
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]);
vImarisDataSet.SetChannelColorRGBA(vNumberOfChannels+1, vRGBA);
vImarisApplication.SetDataSet(vImarisDataSet);


%Create a new folder object for new surfaces
Voronoi_Surfaces = vImarisApplication.GetFactory;
result = Voronoi_Surfaces.CreateDataContainer;
result.SetName('Voronoi Surfaces');

%Run the Surface Creation Wizard on the Voronoi channel
ip = vImarisApplication.GetImageProcessing;
Voronoi_Surfaces = ip.DetectSurfaces(vImarisDataSet, [], vNumberOfChannels, 0, 0, true, 55, '');
Voronoi_Surfaces.SetName(sprintf('Voronoi Surface'));
Voronoi_Surfaces.SetColorRGBA((rand(1, 1)) * 256 * 256 * 256 );

%Run the Surface Creation Wizard on the random channel
ip = vImarisApplication.GetImageProcessing;
VoronoiRandom_Surfaces = ip.DetectSurfaces(vImarisDataSet, [], vNumberOfChannels+1, 0, 0, true, 55, '');
VoronoiRandom_Surfaces.SetName(sprintf('Random Voronoi Surface'));
VoronoiRandom_Surfaces.SetColorRGBA((rand(1, 1)) * 256 * 256 * 256 );

%Add new surface to Surpass Scene
VoronoiRandom_Surfaces.SetVisible(0);
Voronoi_Surfaces.SetVisible(1);

result.AddChild(Voronoi_Surfaces, -1);
result.AddChild(VoronoiRandom_Surfaces, -1);
vImarisApplication.GetSurpassScene.AddChild(result, -1);

vImarisApplication.SetChannelVisibility(vNumberOfChannels+1,0);
vImarisApplication.SetChannelVisibility(vNumberOfChannels ,0);


function mask = voronoi2mask(x,y,szImg,vDataSize)
%
% voronoi2mask Convert Voronoi cells to region mask
%
%   mask = voronoi2mask(x,y,szImg) computes a mask of the Voronoi cells
%   given points,'x' and 'y' and a 2d-image size 'szImg', which
%   the points are extracted from. The voronoi diagram is created
%   using Matlab's voronoi function.
%
%     Example
%     -------
%
%         % get image from steve's blog at mathworks.com
%         if ~exist('nuclei.png','file')
%            img = imread('http://blogs.mathworks.com/images/steve/60/nuclei.png');
%            imwrite(img,'nuclei.png');
%         else
%            img = imread('/nuclei.png');
%         end
%         img = double(img);
%
%         % crop image
%         img = img(1:300,1:350);
%
%         % "blur" image with imopen
%         se = strel('disk', 15);
%         imgo = imopen(img, se);
%
%         % find regional max
%         imgPros = imregionalmax(imgo,4);
%
%         % get centroids of regional max
%         objects = regionprops(imgPros,{'Centroid', 'BoundingBox','Image'});
%
%         % save centroids to array
%         centroids = nan([numel(objects),2]);
%         for i = 1:numel(objects)
%             centroids(i,:) = objects(i).Centroid;
%         end
%
%         % based on the centroids, create the voronoi diagram
%         % and transform the Voronoi cells to an image.
%         mask = voronoi2mask(centroids(:,1),centroids(:,2),size(img));
%         mask = label2rgb(mask, 'jet', 'c', 'shuffle');
%
%         % show results
%         subplot(2,1,1)
%         imagesc(img);colormap('gray');
%         hold on;
%         h = voronoi(centroids(:,1),centroids(:,2));
%         set(h(:),'Color',[0 1 0]);
%         axis image;
%         title('original image with voronoi diagram');
%
%         subplot(2,1,2)
%         imshow(mask);
%         hold on;
%         h = voronoi(centroids(:,1),centroids(:,2));
%         set(h(:),'Color',[0 0 0]);
%         axis image;
%         title('output image with voronoi diagram');
%
%     See also poly2mask, roipoly.
%
%
% $Created: 1.0 $ $Date: 2013/08/11 20:00$ $Author: Pangyu Teng $
%

if nargin < 3
    display('requires 3 inputs. (voronoi2mask.m)');
    return;
end

% format x, y to be column wise
if size(x,1) < size(x,2)
    x = x';
end

if size(y,1) < size(y,2)
    y = y';
end

% create voronoi diagram and get its finite vertices
[vx, vy] = voronoi(x,y);

% create a mask to draw the vertices
border = logical(false(szImg));
mask = zeros(szImg);

% draw vertices on mask
for i = 1:size(vx,2)
    
    % create line function between 2 points
    f = makelinefun(vy(1,i),vx(1,i),vy(2,i),vx(2,i),2);
    
    % get distance between 2 points
    dist = round(1.5*sqrt(diff(vx(:,i)).^2+diff(vy(:,i)).^2));
    
    % create 'dist' points on the line
    [vxLine, vyLine] = f(dist);
    
    % round the line
    vxLine = round(vxLine);
    vyLine = round(vyLine);
    
    % contrain line to be within the image
    validInd = vxLine >= 1 & vxLine <= szImg(1) & vyLine >= 1 & vyLine <= szImg(2);
    vxLine = vxLine(validInd);
    vyLine = vyLine(validInd);
    
    % draw the line to an image
    newInd = sub2ind(szImg,vxLine,vyLine);
    border(newInd) = true;
end

% round xs and yx
x = round(x);
y = round(y);
x(x==0)=1;
y(y==0)=1;
% number each region based on the index of the centroids
% (xs and ys are "flipped" ...)
mask=zeros(vDataSize(2),vDataSize(1));
for i = 1:numel(x)
    if i>1
        bwOld=bw;
    end
    bw = imfill(border,sub2ind(szImg,y(i),x(i)));
    %mask(bw(:)==1) = i;
    
    %test and show overlapping positions
    if i>1
        mask=mask+bw+bwOld;
        mask(mask==1)=0;%Set all 1's to zero
        %Coloc(Coloc<2)=0;
        %Coloc(Coloc>1)=1;
    end
end








