function Viewer(varargin)
% This function contains a Gui for testing different thresholds. The Gui 
% provides three slides, one for the selection of the image slice in a 
% volume containing medical data, and two for setting an upper and a lower
% threshold for threshold segmentation. 
%
% Input: Medical volume containing image slices. 


global val thresh1 thresh2;
val = 1;

width = 1000;
[rows, cols,number] = size(varargin{1,1});
max_intensity = max(max(max(varargin{1,1}))) ;%
min_intensity = min(min(min(varargin{1,1}))) ; % plus one to not get out of bounds

f = figure('Visible','on','Position',[width,500,2*cols + 150,rows+150]); 
% position, storlek

 sld1 = uicontrol('Style', 'slider',...
     'Min',1,'Max',number,'Value',1,...
     'Position', [((cols-120+50)/2) 40 120 20],...    % position, storlek
     'Callback', @sliderCont_Callback);
 
 thr1 = uicontrol('Style', 'slider',...
     'Min',min_intensity,'Max',max_intensity,'Value',min_intensity,...
     'Position', [((width + cols-120+200)/2) 45 120 20],...    % position, storlek
     'Callback', @sliderthr1_Callback); 
 
 thr2 = uicontrol('Style', 'slider',...
     'Min',min_intensity,'Max',max_intensity,'Value',min_intensity,...
     'Position', [((width + cols-120+200)/2) 15 120 20],...    % position, storlek
     'Callback', @sliderthr2_Callback); 
 
 
 txt = uicontrol('Style','text',...
     'Position',[((rows-40+50)/2) 15 40 20],...
     'String',num2str(val));
 thresh_txt1 = uicontrol('Style','text',...
                'Position',[((width + cols-120+200)/2+200) 45 40 20],...
                'String',num2str(thresh1)); 
 thresh_txt2 = uicontrol('Style','text',...
                'Position',[((width + cols-120+200)/2+200) 15 40 20],...
                'String',num2str(thresh2)); 
 
 
 haxes1 = axes('Units','pixels','Position',[50,100,rows,cols]);
 haxes2 = axes('Units','pixels','Position',[cols + 100,100,rows,cols]); 

 image_uppdate(1,varargin,haxes1);
 
%     function image_Callback(source,eventdata, varargin)

    function image_uppdate(val,varargin,ha)
        imagesc(varargin{1,1}(:,:,val),'Parent',ha);
        colormap gray;
     end
    
    function sliderCont_Callback(source,eventdata)
        val = round(get(source,'Value'));
        set(txt,'String',num2str(val))
        image_uppdate(val,varargin,haxes1)
        thresh_uppdate(val,thresh1, thresh2,varargin,haxes2)
    end


    function sliderthr1_Callback(source,eventdata)
        thresh1 = round(get(source,'Value'));
        if(thresh1 > thresh2)
           set(thr2,'Value',thresh1);
           thresh2 = thresh1+1;
        end
        set(thresh_txt1,'String',num2str(thresh1))
        set(thresh_txt2,'String',num2str(thresh2))
        thresh_uppdate(val, thresh1, thresh2,varargin,haxes2)
    end

    function sliderthr2_Callback(source,eventdata)
        thresh2 = round(get(source,'Value'));
        if(thresh2 < thresh1)
           set(thr1,'Value',thresh2);
           thresh1 = thresh2-1;
        end
        set(thresh_txt2,'String',num2str(thresh2))
        set(thresh_txt1,'String',num2str(thresh1))
        thresh_uppdate(val, thresh1, thresh2,varargin,haxes2)
    end

    function thresh_uppdate(val,thresh1, thresh2,varargin,ha)
        im = varargin{1,1}(:,:,val);
        %[~, imageThreshold]=histc(im,[thresh1 thresh2]);
        imageThreshold = thresh1 <= im & im <= thresh2;
        imageThreshold(imageThreshold~=1)=0;
        imagesc(imageThreshold,'Parent',ha);
        colormap gray;
    end
end