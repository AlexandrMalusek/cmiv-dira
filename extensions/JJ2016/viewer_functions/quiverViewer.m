function quiverViewer(varargin)
% GUI for visualizing the displacement field as arrows. The GUI allows the
% user to go throught the slices in the medical data, change the contrast
% and change the arrow size and distance using sliders.
%
%
% Input:    - Volume containing the medical data. 
%
%           - Displacment field in y direction.
%           - Displacment field in x direction.
%           - Displacment field in z direction.
%
%           - Direction of slices. Options are: 'xy', 'xz' and 'yz'.
%             At default it is set to 'xy'. (Optional) 
%
%           - Initial window for the contrast of the image slices. It needs 
%             to be an interval [lowBoundary, highBoundary]. (Optional)
%
%           - Color of the arrows using matlabs standard colors. The 
%             default setting is yellow ('y'). (Optional)
%
%           - Initial scale of the arrows, If set to zero, the autoscaling
%             of the quiver function is used. (Optional) 
%
% Example: 
%   quiverViewer(ctVolume, displacementFieldX, displacementFieldY, ...
%             displacementFieldZ, 'xy', [10 1250], 'r', 0)
%
%


quiverScale = 0;
quiverDist  = 1;

if(nargin == 8)
    if(isempty(varargin{1,8}))
        quiverScale =  10;
    elseif varargin{1,8} == 0
        quiverScale = 0;
    else
    
        quiverScale =  varargin{1,8};
    end
end


if(nargin >= 7)
    if(isempty(varargin{1,7}))
        yellowColors = 1;
    else
        yellowColors = 0;
        color     = varargin{1,7};
    end
else
    yellowColors = 1;
end

% Setting the direction of the slices
if numel(varargin) >= 5
    if ~isempty(varargin{1,5})
        switch varargin{1,5}
            case 'xz'
                ctVol = permute(varargin{1,1},[3 2 1]);
                dispFieldY = permute(varargin{1,3},[3 2 1]);
                dispFieldX = permute(varargin{1,4},[3 2 1]);
            case 'yz'
                ctVol = permute(varargin{1,1},[3 1 2]);
                dispFieldY = permute(varargin{1,2},[3 1 2]);
                dispFieldX = permute(varargin{1,4},[3 1 2]);
            otherwise
                ctVol = varargin{1,1};
                dispFieldY = varargin{1,2};
                dispFieldX = varargin{1,3}; 
        end
    else
        ctVol = varargin{1,1};
        dispFieldY = varargin{1,2};
        dispFieldX = varargin{1,3};
    end
else
    ctVol = varargin{1,1};
    dispFieldY = varargin{1,2};
    dispFieldX = varargin{1,3};
end


%%
% global val, global w1, global w2;
val = 1; 
minSize = 500;
minVal = min(ctVol(:));
maxVal = max(ctVol(:));

width = 1000;
[rows, cols,number] = size(ctVol);

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
if numel(varargin) >= 6
    if ~isempty(varargin{1,6})
        w1 = varargin{1,6}(1);
        w2 = varargin{1,6}(2);
    else
        w1 = max(double(minVal),1); 
        w2 = double(maxVal)-1;
    end
else
    w1 = max(double(minVal),1); 
    w2 = double(maxVal)-1;
end

[sy,sx,sz] = size(ctVol);


%% Create sliders and Windows
f = figure('Visible','on','Position',[width,500,sizeX + 150,sizeY+150]);
% position, storlek


% Slider text
txtSlice = uicontrol('Style','text',...
    'Position',[45 60 80 20],...
    'String','Image slice');

txtContrast = uicontrol('Style','text',...
    'Position',[250 60 100 20],...
    'String','Image contrast');

txtArrow = uicontrol('Style','text',...
    'Position',[483 60 130 20],...
    'String','Arrow size & distance');





% Sliders
sld1 = uicontrol('Style', 'slider',...
    'Min',1,'Max',number,'Value',1,...
    'Position', [20 40 120 20],...    % position, storlek
    'SliderStep',[1/number 0.2],'Callback', @sliderCont_Callback);

wind1 = uicontrol('Style', 'slider',...
     'Min',1,'Max',maxVal,'Value',w1,...
     'Position', [((sizeY-100+50)/2) 40 120 20],...    % position, storlek
     'Callback', @sliderthr1_Callback); 
 
wind2 = uicontrol('Style', 'slider',...
     'Min',1,'Max',maxVal,'Value',w2,...
     'Position', [((sizeY-100+50)/2) 20 120 20],...    % position, storlek
     'Callback', @sliderthr2_Callback); 
 
arrowSize = uicontrol('Style', 'slider',...
      'Min',0,'Max',5,'Value',0,...
      'Position', [(sizeY-40) 40 120 20],...    % position, storlek
      'Callback', @sliderQuiverSize_Callback); 
  
quivDist = uicontrol('Style', 'slider',...
      'Min',1,'Max',25,'Value',10,...
      'Position', [(sizeY-40) 20 120 20],...    % position, storlek
      'Callback', @sliderQuiverDist_Callback); 

% Slice number  
txt = uicontrol('Style','text',...
    'Position',[60 15 40 20],...
    'String',num2str(val));




haxes1 = axes('Units','pixels','Position',[50,100,sizeX+25,sizeY+25]);



%% Defining the calback functions

    function sliderCont_Callback(source,eventdata)
        val = round(get(source,'Value'));
        set(txt,'String',num2str(val))
        imageUppdate(val,varargin,haxes1,w1,w2,quiverScale,quiverDist)
    end
    function sliderthr1_Callback(source,eventdata)
        w1 = double(round(get(source,'Value')));
        %set(txt,'String',num2str(val))
        imageUppdate(val,varargin,haxes1,w1,w2,quiverScale,quiverDist)
    end
    function sliderthr2_Callback(source,eventdata)
        w2 = double(round(get(source,'Value')));
        %set(txt,'String',num2str(val))
        imageUppdate(val,varargin,haxes1,w1,w2,quiverScale,quiverDist)
    end
    function sliderQuiverSize_Callback(source,eventdata) % DO I NEED ROUND?
        quiverScale = get(source,'Value');
        %set(txt,'String',num2str(val))
        imageUppdate(val,varargin,haxes1,w1,w2,quiverScale,quiverDist)
    end

    function sliderQuiverDist_Callback(source,eventdata)
        quiverDist = round(get(source,'Value'));
        %set(txt,'String',num2str(val))
        imageUppdate(val,varargin,haxes1,w1,w2,quiverScale,quiverDist)
    end


    function imageUppdate(val,varargin,ha,w1,w2,quiverScale,quiverDist)   
        im1 = double(ctVol(:,:,val));
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
        
        
        xGrid = 1:quiverDist:cols;
        yGrid = 1:quiverDist:rows;

        axis([1 sx 1 sy])
        hold on
        if(yellowColors)
            quiver(xGrid, yGrid, dispFieldX(yGrid,xGrid,val), ...
                   dispFieldY(yGrid,xGrid,val), quiverScale,'y');
        else
            quiver(xGrid, yGrid, dispFieldX(yGrid,xGrid,val), ...
                   dispFieldY(yGrid,xGrid,val),quiverScale,color);
        end
        hold off

     end

 
end