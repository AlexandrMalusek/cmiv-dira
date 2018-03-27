function SegmentationVisualizer(varargin)
% Function for visualizing the segmentation. The first input is the ct
% volume, the second is the first segmentation result(red) and the third is
% the second segmentation result (green). In areas where the segmentation
% overlap the segmentation is visulaized using the color yellow.


%% Get sliders right, etc.
global val, global w1, global w2, global alpha;
val = 1; alpha = 0.6;
min_size = 500;
min_val = min(varargin{1,1}(:));
max_val = max(varargin{1,1}(:));

w1 = double(min_val); 
w2 = double(max_val);

width = 1000;
[rows, cols,number] = size(varargin{1,1});

if(rows > min_size)
    size_y = rows;
else
    size_y = min_size;
end

if(cols > min_size)
    size_x = cols;
else
    size_x = min_size;
end

% Set the initial windowing
if numel(varargin) == 4
    w1 = varargin{1,4}(1);
    w2 = varargin{1,4}(2);
else
    w1 = max(double(min_val),1); 
    w2 = double(max_val)-1;
end

%% Create a volume with the overlap and one for the differeces, ... 
segm1 = varargin{1,2};
segm2 = varargin{1,3};

% 
overlap   = segm1 & segm2;
seg1_diff = segm1 - overlap;
seg2_diff = segm2 - overlap;

% Colors
o = ones( rows, cols); 
z = zeros(rows, cols);

red    = cat(3,     o,     z, z);
green  = cat(3,     z, 0.5*o, z);
blue   = cat(3,     z,     z, o);
yellow = cat(3, 0.8*o, 0.8*o, z);  
%%
f = figure('Visible','on','Position',[width,500,size_x + 150,size_y+150]);
% position, storlek

sld1 = uicontrol('Style', 'slider',...
    'Min',1,'Max',number,'Value',1,...
    'Position', [20 40 120 20],...    % position, storlek
    'SliderStep',[1/number 0.2],'Callback', @sliderCont_Callback);

wind1 = uicontrol('Style', 'slider',...
     'Min',1,'Max',max_val,'Value',w1,...
     'Position', [((size_y-100+50)/2) 40 120 20],...    % position, storlek
     'Callback', @sliderthr1_Callback); 
 
wind2 = uicontrol('Style', 'slider',...
     'Min',1,'Max',max_val,'Value',w2,...
     'Position', [((size_y-100+50)/2) 20 120 20],...    % position, storlek
     'Callback', @sliderthr2_Callback); 
 
blend = uicontrol('Style', 'slider',...
     'Min',0,'Max',1,'Value',alpha,...
     'Position', [(size_y-40) 40 120 20],...    % position, storlek
     'Callback', @sliderblend_Callback); 


txt = uicontrol('Style','text',...
    'Position',[60 15 40 20],...
    'String',num2str(val));


haxes1 = axes('Units','pixels','Position',[50,100,size_x+100,size_y]);


    function sliderCont_Callback(source,eventdata)
        val = round(get(source,'Value'));
        set(txt,'String',num2str(val))
        image_uppdate(val,varargin,haxes1,w1,w2,alpha)
    end
    function sliderthr1_Callback(source,eventdata)
        w1 = double(round(get(source,'Value')));
        %set(txt,'String',num2str(val))
        image_uppdate(val,varargin,haxes1,w1,w2,alpha)
    end
    function sliderthr2_Callback(source,eventdata)
        w2 = double(round(get(source,'Value')));
        %set(txt,'String',num2str(val))
        image_uppdate(val,varargin,haxes1,w1,w2,alpha)
    end
    function sliderblend_Callback(source,eventdata)
        alpha = get(source,'Value');
        %set(txt,'String',num2str(val))
        image_uppdate(val,varargin,haxes1,w1,w2,alpha)
    end


    function image_uppdate(val,varargin,ha,w1,w2,alpha)
        im1 = double(varargin{1,1}(:,:,val));
        %im2 = double(varargin{1,2}(:,:,val));
        %show_im = (1-alpha)*im1 + w2*alpha*im2; 
        imagesc(im1,[w1,w2]);%,'Parent',ha);
        colorbar; colormap gray;
        axis image
 
        hold on; h1 = imagesc(yellow);hold off
        set(h1, 'AlphaData', alpha*overlap(:,:,val))
        hold on; h2 = imagesc(red);hold off
        set(h2, 'AlphaData', alpha*seg1_diff(:,:,val))
        hold on; h3 = imagesc(blue);hold off
        set(h3, 'AlphaData', alpha*seg2_diff(:,:,val))
        %hold on
        %h = imagesc(true_color); axis image; axis off
        %hold off

        %set(h, 'AlphaData', influence_map)
        
     end

 
end