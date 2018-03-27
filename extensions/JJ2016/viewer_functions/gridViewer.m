function gridViewer(varargin)
%
% Input:    - CT Volume
%           - Displacment field in y direction.
%           - Displacment field in x direction.
%           - Displacment field in z direction.
%           - Direction of slices (Optional)
%           - Initial window (Optional)
%           - Grid color (Optional). When left empty it is set to multiple colors 
%           - Distance for gridlines (Optional). 

spacing =  10;%
%color     = 'y';
lineWidth = 1;

if(nargin == 8)
    if(isempty(varargin{1,8}))
        spacing =  10;
    elseif varargin{1,8} == 0
        spacing = 10;
    else
    
        spacing =  varargin{1,8};
    end
end


if(nargin >= 7)
    if(isempty(varargin{1,7}))
        multipleColors = 1;
    else
        multipleColors = 0;
        color     = varargin{1,7};
    end
else
    multipleColors = 1;
end

% Setting the direction of the slices
if numel(varargin) >= 5
    if ~isempty(varargin{1,5})
        switch varargin{1,5}
            case 'xz'
                ctVol = permute(varargin{1,1},[3 2 1]);
                yDisp = permute(varargin{1,3},[3 2 1]);
                xDisp = permute(varargin{1,4},[3 2 1]);
            case 'yz'
                ctVol = permute(varargin{1,1},[3 1 2]);
                yDisp = permute(varargin{1,2},[3 1 2]);
                xDisp = permute(varargin{1,4},[3 1 2]);
            otherwise
                ctVol = varargin{1,1};
                yDisp = varargin{1,2};
                xDisp = varargin{1,3}; 
        end
    else
        ctVol = varargin{1,1};
        yDisp = varargin{1,2};
        xDisp = varargin{1,3};
    end
else
    ctVol = varargin{1,1};
    yDisp = varargin{1,2};
    xDisp = varargin{1,3};
end


%%
% global val, global w1, global w2, global alpha;
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
[gridRows,gridCols, ...
dispFieldY,dispFieldX] = makeGrid(ctVol,spacing);

%% Create sliders and Windows
f = figure('Visible','on','Position',[width,500,sizeX + 150,sizeY+150]);
% position, storlek

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
 
blend = uicontrol('Style', 'slider',...
      'Min',2,'Max',30,'Value',10,...
      'Position', [(sizeY-40) 40 120 20],...    % position, storlek
      'Callback', @sliderLineDist_Callback); 


txt = uicontrol('Style','text',...
    'Position',[60 15 40 20],...
    'String',num2str(val));


haxes1 = axes('Units','pixels','Position',[50,100,sizeX+25,sizeY+25]);



%% Defining the calback functions

    function sliderCont_Callback(source,eventdata)
        val = round(get(source,'Value'));
        set(txt,'String',num2str(val))
        imageUppdate(val,varargin,haxes1,w1,w2,...
                      gridRows,gridCols,dispFieldY,dispFieldX)
    end
    function sliderthr1_Callback(source,eventdata)
        w1 = double(round(get(source,'Value')));
        %set(txt,'String',num2str(val))
        imageUppdate(val,varargin,haxes1,w1,w2,...
                      gridRows,gridCols,dispFieldY,dispFieldX)
    end
    function sliderthr2_Callback(source,eventdata)
        w2 = double(round(get(source,'Value')));
        %set(txt,'String',num2str(val))
        imageUppdate(val,varargin,haxes1,w1,w2,...
                      gridRows,gridCols,dispFieldY,dispFieldX)
    end
    function sliderLineDist_Callback(source,eventdata)
        spacing = round(get(source,'Value'));
        [gridRows,gridCols, ...
         dispFieldY,dispFieldX] = makeGrid(ctVol,spacing);
        %set(txt,'String',num2str(val))
        imageUppdate(val,varargin,haxes1,w1,w2,...
                      gridRows,gridCols,dispFieldY,dispFieldX)
    end


    function [gridRows,gridCols, ...
              dispFieldY,dispFieldX] = makeGrid(ctVol,spacing)
        %% Make the grid
        [sy,sx,sz] = size(ctVol);
        %
        [x,y] = meshgrid(1:sx,1:sy,1:sz);
        
        gridRows = 1:spacing:sx;
        gridCols = 1:spacing:sy;
        
        dispFieldY = y + yDisp;
        dispFieldX = x + xDisp;
        
    end



    function imageUppdate(val,varargin,ha,w1,w2,...
                          gridRows,gridCols,dispFieldY,dispFieldX)   
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
        
       
        hold on
        p1 = plot(dispFieldX(:,gridRows,val),dispFieldY(:,gridRows,val));
        p2 = plot(dispFieldX(gridCols,:,val)',dispFieldY(gridCols,:,val)');
        hold off
        axis([1 sx 1 sy])
        
        ph = [p1;p2];
        if(multipleColors)
           set(ph,'linewidth',lineWidth);
        else
           set(ph,'color',color,'linewidth',lineWidth);set(ph,'linewidth',lineWidth);
        end
%         hold on; 
%         h1 = imagesc(segment_col);
%         hold off
%         set(h1, 'AlphaData', alpha*varargin{1,2}(:,:,val))
     end

 
end