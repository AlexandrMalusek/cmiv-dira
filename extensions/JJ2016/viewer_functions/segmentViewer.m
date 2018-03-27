function segmentViewer(varargin)
% SEGMENTVIEWER is a GUI for visualizing the segmentation results of the 
% JJ2016 algorithm. 
%
% Input:    - Medical data volume, i.e. a volume of CT slices.
%
%           - Binary volume containing the bone tissue. 
%           - Binary volume containing the adipose tissue. 
%           - Binary volume containing the rectume tissue.
%           - Binary volume containing the prostate tissue.
%           - Binary volume containing the air inside the body. 
%
%           - Initial window for the contrast. If the initial contrast  
%             window is left empty or set to [] the max and the min volume 
%             of the medical data volume are set as initial window values.   
%
%           - 'xy', 'xz' or 'yz' plane. Default is set to 'xy'. 
%
% 
% EXAMPLE: segmentViewer(imageData, boneMask, adiposeMask, rectumMask, ...
%                        prostateMask, airMask);
%
%


% global val, global w1, global w2, global alpha;
val = 1; alpha = 0;
width = 1000;

minSize = 500;
minValue = min(varargin{1,1}(:));
maxValue = max(varargin{1,1}(:));


% bone, adipose, rectum, prostate, air)

% Tissues
imageData = varargin{1,1};
boneMap = double(varargin{1,2});
adiposeMap = double(varargin{1,3});
rectumMap = double(varargin{1,4});
prostateMap = double(varargin{1,5});
airMap = double(varargin{1,6});



% Setting the direction of the slices
if numel(varargin) == 8
    switch varargin{1,8}
        case 'xy'

        case 'xz'
            imageData = permute(imageData,[3 2 1]);
            boneMap = double(permute(boneMap,[3 2 1]));
            adiposeMap = double(permute(adiposeMap,[3 2 1]));
            rectumMap = double(permute(rectumMap,[3 2 1]));
            prostateMap = double(permute(prostateMap,[3 2 1]));
            airMap = double(permute(airMap,[3 2 1]));
        case 'yz'
            imageData = permute(imageData,[3 1 2]);
            boneMap = double(permute(boneMap,[3 1 2]));
            adiposeMap = double(permute(adiposeMap,[3 1 2]));
            rectumMap = double(permute(rectumMap,[3 1 2]));
            prostateMap = double(permute(prostateMap,[3 1 2]));
            airMap = double(permute(airMap,[3 1 2]));
        otherwise
            disp('4th inargument needs to be: xy, xz or yz')
    end
end

boneMap(boneMap<0) = 0;
adiposeMap(adiposeMap<0) = 0;
rectumMap(rectumMap<0) = 0;
prostateMap(prostateMap<0) = 0;
airMap(airMap<0) = 0;
            
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
if numel(varargin) >= 7
    if ~isempty(varargin{1,3})
        w1 = varargin{1,7}(1);
        w2 = varargin{1,7}(2);
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

blend = uicontrol('Style', 'slider',...
     'Min',0,'Max',1,'Value',0,...
     'Position', [(sizeY-40) 40 120 20],...    % position, storlek
     'Callback', @sliderblend_Callback); 


txt = uicontrol('Style','text',...
    'Position',[60 15 40 20],...
    'String',num2str(val));


haxes1 = axes('Units','pixels','Position',[50,100,sizeX+25,sizeY+25]); % [50,100,size_x+100,size_y]);

%% Choose the color for the vizualization of the segmentation
o = ones( rows, cols);
z = zeros(rows, cols);

% colors
red = cat(3,o, z, z);
green = cat(3,  z,0.5*o, z);
blue = cat(3, z, z, o); 
yellow = cat(3, 0.8*o, 0.8*o, z);   
margenta = cat(3, o, z, o);
cyan = cat(3,  z,o, o);





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
        im1 = double(imageData(:,:,val));
%        im2 = double(varargin{1,2}(:,:,val));
%        show_im = (1-alpha)*im1 + w2*alpha*im2; 
        imagesc(im1,[w1,w2]);
        colormap gray;
        axis image
        
        % Manual placement of colorbar
%         colorbar('location','manual','position',[0.85 0.15 0.15 0.8], ...
%                  'plotboxaspectratiomode','manual', ...
%                  'plotboxaspectratio',[0.3 5 1]);
        warning('off','MATLAB:colorbar:InvalidRulerLocation')
        hc = colorbar('location','eastoutside','position',[0.9 0.15 0.05 0.81]);
        set(hc,'xaxisloc','top');
        
             
        % Segmentation results     
        hold on; h1 = imagesc(yellow); hold off
        set(h1, 'AlphaData', min(alpha*boneMap(:,:,val),1))
             
        hold on; h1 = imagesc(green); hold off
        set(h1, 'AlphaData', min(alpha*adiposeMap(:,:,val),1))
             
        hold on; h1 = imagesc(red); hold off
        set(h1, 'AlphaData', min(alpha*rectumMap(:,:,val),1))
             
        hold on; h1 = imagesc(blue); hold off
        set(h1, 'AlphaData', min(alpha*prostateMap(:,:,val),1))
             
        hold on; h1 = imagesc(margenta); hold off
        set(h1, 'AlphaData', min(alpha*airMap(:,:,val),1))
             
             
    end


end