function h = displayGrid(background,dispFieldX,dispFieldY,gridSpacing,figNr, figPosition)
% Function for displaying a grid that has been displaced by the dispacement
% field. 
%
% Inputs:   - background:  Can be a rgb image, a black and white image or a
%                          vector that defines the background color. For 
%                          example, the vector [1 0 0] defines a red 
%                          background. If kept empty, the background is set
%                          to gray. 
%
%           - dispFieldX:  Displacement field contaning the x component.
%    
%           - dispFieldY:  Displacement field contaning the y component.
% 
%           - gridSpacing: The distance between the grid lines. Can be
%                          defined as vector defining the distance in y and
%                          x direction or as an integer resutling in the
%                          same distance in x and y direction.
%
%           - figNr:       The figure number. If kept empty, [], or not
%                          given as input argument, a new figure with the
%                          lowest unused numbering is created. 
%
%           - figPosition: The position and the size of the figure. Can be 
%                          left empty. It can for example be set to 
%                          [100, 100, 512, 512] to get a figure that is
%                          rectangular.
%
% Output:   - h:           Handle to the figure.
%



%% Settings

color     = 'y';
lineWidth = 1;
backgroundColor = [0.8 0.8 0.8];

%% Test the inarguments


if nargin < 4
    gridSpacing   = [10,10];
elseif isempty(gridSpacing)
    gridSpacing   = [10,10];
elseif size(gridSpacing,2) == 1
    gridSpacing = [gridSpacing, gridSpacing];
end

if isempty(background) 
    background = backgroundColor;
end

if nargin < 5
    h = figure;
elseif isempty(figNr) || figNr < 1 
    h = figure;
else    
    h = figure(round(figNr));
end


if nargin == 6
    set(h,'Position',figPosition);%[100, 100, 512, 512]);
end


%% Plot the background

if(ndims(background)>1 && size(background,2) > 4)
    imagesc(background);
    colormap gray; axis image;
elseif(ndims(background)>1 && size(background,2) == 3)
    backgroundColor = background; 
    set(gca,'Color',backgroundColor);
else
    error(['The third input argument can be an image, a 1x3'...
          ' vector defining the background color or be left' ...
          ' empty where the color will be set to gray']);
end

colormap gray; %axis image; axis off
%set(gca,'position',[0 0 1 1],'units','normalized')

%% Plot the grid
[sy,sx] = size(dispFieldY);
% 
[x,y] = meshgrid(1:sx,1:sy);

rows = 1:gridSpacing(2):sx;
cols = 1:gridSpacing(1):sy;

dispFieldY = y + dispFieldY;
dispFieldX = x + dispFieldX;

hold on
p1 = plot(dispFieldX(:,rows),dispFieldY(:,rows));
p2 = plot(dispFieldX(cols,:)',dispFieldY(cols,:)');
hold off
%axis([1 sx 1 sy])

ph = [p1;p2];
set(ph,'color',color,'linewidth',lineWidth);
set(ph,'linewidth',lineWidth);

end