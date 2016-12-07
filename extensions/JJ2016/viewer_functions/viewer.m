function viewer(varargin)
% VIEWER is a GUI for visualizing the individual image slices in a
% medical data set. It contains sliders for selecting the image slice and
% adjusting the image contrast. The image data can also be visualized
% together with a binary mask, which can be used to evaluate the results
% of segmenentation algorithms.
%
% Input:    - Medical data volume, i.e. a volume of CT slices.
%           - A binary mask showing the segmentation result. (Optional)
%           - Initial window for the contrast. If the initial contrast  
%             window is left empty or set to [] the max and the min volume 
%             of the medical data volume are set as initial window values.   
%           - 'xy', 'xz' or 'yz' plane. Default is set to 'xy'. 



% global val, global w1, global w2, global alpha;
val = 1; alpha = 0;
width = 1000;

minSize = 500;
minValue = min(varargin{1,1}(:));
maxValue = max(varargin{1,1}(:));

% If only the medical data is used, create a binary volume containing zeros
if nargin == 1
    imageData = varargin{1,1};
    binaryMap = zeros(size(imageData));
    useBinaryVol = 0;
elseif isempty(varargin{1,2})
    imageData = varargin{1,1};
    binaryMap = zeros(size(imageData));
    useBinaryVol = 0;
else
    imageData = varargin{1,1};
    binaryMap = double(varargin{1,2});
    useBinaryVol = 1;
end

% Setting the direction of the slices
if numel(varargin) == 4
    switch varargin{1,4}
        case 'xy'
            imageData = imageData;
            binaryMap = double(binaryMap);
            binaryMap(binaryMap<0) = 0;
        case 'xz'
            imageData = permute(imageData,[3 2 1]);
            binaryMap = double(permute(binaryMap,[3 2 1])); % 2 3 1
            binaryMap(binaryMap<0) = 0;
        case 'yz'
            imageData = permute(imageData,[3 1 2]);
            binaryMap = double(permute(binaryMap,[3 1 2]));
            binaryMap(binaryMap<0) = 0;
        otherwise
            disp('4th inargument needs to be: xy, xz or yz')
    end
end

% Get the size of the image.
[rows, cols,number] = size(imageData);

if(rows > minSize)
    sizeY = rows;
else
    sizeY = minSize;
end

if(cols > minSize)
    sizeX = cols;
else
    sizeX = minSize;
end

% Set the initial windowing
if numel(varargin) >= 3
    if ~isempty(varargin{1,3})
        w1 = varargin{1,3}(1);
        w2 = varargin{1,3}(2);
    else
        w1 = max(double(minValue),1); 
        w2 = double(maxValue)-1;
    end
else
    w1 = max(double(minValue),1); 
    w2 = double(maxValue)-1;
end


%% Create sliders and Windows
f = figure('Visible','on','Position',[width,500,sizeX + 150,sizeY+150]);
% position, storlek

sld1 = uicontrol('Style', 'slider',...
    'Min',1,'Max',number,'Value',1,...
    'Position', [20 40 120 20],...    % position, storlek
    'SliderStep',[1/number 0.2],'Callback', @sliderCont_Callback);

wind1 = uicontrol('Style', 'slider',...
     'Min',1,'Max',maxValue,'Value',w1,...
     'Position', [((sizeY-100+50)/2) 40 120 20],...    % position, storlek
     'Callback', @sliderthr1_Callback); 
 
wind2 = uicontrol('Style', 'slider',...
     'Min',1,'Max',maxValue,'Value',w2,...
     'Position', [((sizeY-100+50)/2) 20 120 20],...    % position, storlek
     'Callback', @sliderthr2_Callback); 

if useBinaryVol 
    blend = uicontrol('Style', 'slider',...
        'Min',0,'Max',1,'Value',0,...
        'Position', [(sizeY-40) 40 120 20],...    % position, storlek
        'Callback', @sliderblend_Callback); 
end

txt = uicontrol('Style','text',...
    'Position',[60 15 40 20],...
    'String',num2str(val));


haxes1 = axes('Units','pixels','Position',[50,100,sizeX+25,sizeY+25]); % [50,100,size_x+100,size_y]);

%% Choose the color for the vizualization of the segmentation
o = ones( rows, cols);
z = zeros(rows, cols);
segmentColor = cat(3, 0.8*o, 0.8*o, z);  

%% Defining the calback functions

    function sliderCont_Callback(source,eventdata)
        val = round(get(source,'Value'));
        set(txt,'String',num2str(val))
        image_uppdate(val,haxes1,w1,w2,alpha)
    end
    function sliderthr1_Callback(source,eventdata)
        w1 = double(round(get(source,'Value')));
        %set(txt,'String',num2str(val))
        image_uppdate(val,haxes1,w1,w2,alpha)
    end
    function sliderthr2_Callback(source,eventdata)
        w2 = double(round(get(source,'Value')));
        %set(txt,'String',num2str(val))
        image_uppdate(val,haxes1,w1,w2,alpha)
    end
    function sliderblend_Callback(source,eventdata)
        alpha = get(source,'Value');
        %set(txt,'String',num2str(val))
        image_uppdate(val,haxes1,w1,w2,alpha)
    end


    function image_uppdate(val,ha,w1,w2,alpha)
%        warning off;
        im1 = double(imageData(:,:,val));
%        im2 = double(varargin{1,2}(:,:,val));
%        show_im = (1-alpha)*im1 + w2*alpha*im2; 
        imagesc(im1,[w1,w2]);
        colormap gray;
        axis image
        
        % Manual placement of colorbar
        warning('off','MATLAB:colorbar:InvalidRulerLocation')
        hc = colorbar('location','eastoutside','position',[0.9 0.15 0.05 0.81]);
        set(hc,'xaxisloc','top');
        
        
        if useBinaryVol
            hold on; 
            h1 = imagesc(segmentColor);
            hold off
            set(h1, 'AlphaData', min(alpha*binaryMap(:,:,val),1))
        end
        
%        warning on;
     end

 
end