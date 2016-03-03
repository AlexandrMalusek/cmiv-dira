function SingleViewer(varargin)
% Function for testing different thresholds on a gray scale image. The
% function is creating a figure with two image windows. In the first window
% is the original image shown and in the second a thresholded version of the
% image. The thresholds are set by adjusting two sliders, one for the upper
% and one for the lower threshold.
%
% Input: Gray scale image
%
% Written by Julius Jeuthe, 2015 



global thresh1 thresh2;


width = 1000;
[rows, cols] = size(varargin{1,1});
max_intensity = max(max(max(varargin{1,1}))) ;%
min_intensity = min(min(min(varargin{1,1}))) ; % plus one to not get out of bounds

f = figure('Visible','on','Position',[width,500,2*cols + 150,rows+150]); 
% position, storlek

 
 thr1 = uicontrol('Style', 'slider',...
     'Min',min_intensity,'Max',max_intensity,'Value',min_intensity,...
     'Position', [((width + cols-120+200)/2) 45 120 20],...    % position, storlek
     'Callback', @sliderthr1_Callback); 
 
 thr2 = uicontrol('Style', 'slider',...
     'Min',min_intensity,'Max',max_intensity,'Value',min_intensity,...
     'Position', [((width + cols-120+200)/2) 15 120 20],...    % position, storlek
     'Callback', @sliderthr2_Callback); 
 
 

 thresh_txt1 = uicontrol('Style','text',...
                'Position',[((width + cols-120+200)/2+200) 45 40 20],...
                'String',num2str(thresh1)); 
 thresh_txt2 = uicontrol('Style','text',...
                'Position',[((width + cols-120+200)/2+200) 15 40 20],...
                'String',num2str(thresh2)); 
 
 
 haxes1 = axes('Units','pixels','Position',[50,100,rows,cols]);
 haxes2 = axes('Units','pixels','Position',[cols + 100,100,rows,cols]); 
 
 
%     function image_Callback(source,eventdata, varargin)

    %function image_uppdate(varargin,ha)
    imagesc(varargin{1,1},'Parent',haxes1);
    colormap gray;
    % end
    
    function sliderthr1_Callback(source,eventdata)
        thresh1 = round(get(source,'Value'));
        if(thresh1 > thresh2)
           set(thr2,'Value',thresh1);
           thresh2 = thresh1+1;
        end
        set(thresh_txt1,'String',num2str(thresh1))
        set(thresh_txt2,'String',num2str(thresh2))
        thresh_uppdate( thresh1, thresh2,varargin,haxes2)
    end

    function sliderthr2_Callback(source,eventdata)
        thresh2 = round(get(source,'Value'));
        if(thresh2 < thresh1)
           set(thr1,'Value',thresh2);
           thresh1 = thresh2-1;
        end
        set(thresh_txt2,'String',num2str(thresh2))
        set(thresh_txt1,'String',num2str(thresh1))
        thresh_uppdate( thresh1, thresh2,varargin,haxes2)
    end

    function thresh_uppdate(thresh1, thresh2,varargin,ha)
        im = varargin{1,1};
        %[~, imageThreshold]=histc(im,[thresh1 thresh2]);
        imageThreshold = thresh1 <= im & im <= thresh2;
        imageThreshold(imageThreshold~=1)=0;
        imagesc(imageThreshold,'Parent',ha);
        colormap gray;
    end
end