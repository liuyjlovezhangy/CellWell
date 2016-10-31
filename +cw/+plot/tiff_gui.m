function fighdl = tiff_gui( im, title_string )

    if ~iscell(im)
        im_cell = {im};
    else
        im_cell = im;
    end

    fighdl = figure('Units','normalized','Position',[0.1,0.1,0.7,0.7]);
    
    frame_scrollbar=uicontrol('style','slider','units','normalized','position',[.1 0.05 0.8 .02],'callback',@framescroll_Callback);   
    channel_scrollbar=uicontrol('style','slider','units','normalized','position',[.9 0.05 0.02 .2],'callback',@channelscroll_Callback);   
    
    handles.im_cell = im_cell;
    handles.title_string = title_string;
    
    handles.frame_scrollbar = frame_scrollbar;
    handles.channel_scrollbar = channel_scrollbar;
    
    guidata(fighdl,handles);
    
    update_scrollbar(im_cell{1},frame_scrollbar);
    
    set(channel_scrollbar, 'min', 1);
    set(channel_scrollbar, 'max', numel(im_cell));
    set(channel_scrollbar, 'Value', 1); % Somewhere between max and min.
    set(channel_scrollbar, 'SliderStep',1/numel(im_cell) * [1 1]);
    
    set(fighdl,'CurrentAxes',axes('Position',[0.1,0.1,0.8,0.8]))

    %%% MAKE FRAME SCROLLBAR UPDATE PLOT WHILE SLIDING
    % because MATLAB guis are a MESS

    hhSlider = handle(frame_scrollbar);
    hProp = findprop(hhSlider,'Value');
    try    % R2014b and newer
       % hProp is a matlab.graphics.internal.GraphicsMetaProperty object
       addlistener(frame_scrollbar,hProp,'PostSet',@framescroll_Callback);
    catch  % R2014a and older
        try
            % hProp is a schema.prop object
            hListener = handle.listener(hhSlider,hProp,'PropertyPostSet',@framescroll_Callback);
            setappdata(frame_scrollbar,'sliderListener',hListener);  % this is important - read above
        catch e
            try    % R2013b and older
                addlistener(frame_scrollbar,'ActionEvent',@framescroll_Callback);
            catch  % R2014a and newer
                addlistener(frame_scrollbar,'ContinuousValueChange',@framescroll_Callback);
            end
        end
    end

    update_gui(im_cell,1,1,title_string);
end

function update_scrollbar(im, frame_scrollbar)
    set(frame_scrollbar, 'min', 1);
    set(frame_scrollbar, 'max', size(im,3));
    set(frame_scrollbar, 'Value', 1); % Somewhere between max and min.
    set(frame_scrollbar, 'SliderStep',1/size(im,3) * [1 10]);
end

function update_gui(im_cell,channel,frame_idx,title_string)
    im = im_cell{channel};
    
    if size(im,3) == 1
        im = mat2gray(im);
        
%         colormap gray
    end
    
    imagesc(squeeze(im(:,:,frame_idx,:)))
    
    axis image
    set(gca,'Ydir','Reverse')
    axis off
    
    colormap gray
    
    
    
    title([title_string ' Channel: ' num2str(channel) ' Frame: ' num2str(frame_idx) ' of ' num2str(size(im,3)) ' [Close fig to move on]'])    
    
    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');
end

function framescroll_Callback(varargin)
    handles = guidata(gcf);

    channel = round(get(handles.channel_scrollbar,'value'));
    frame_idx = round(get(handles.frame_scrollbar,'value'));

    update_gui(handles.im_cell,channel,frame_idx,handles.title_string);
end

function channelscroll_Callback(varargin)
    handles = guidata(gcf);

    channel = round(get(handles.channel_scrollbar,'value'));
    frame_idx = round(get(handles.frame_scrollbar,'value'));

    update_gui(handles.im_cell,channel,frame_idx,handles.title_string);
end
