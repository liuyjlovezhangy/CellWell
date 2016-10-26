function fighdl = tiff_gui( im, title_string )

    fighdl = figure('Units','normalized','Position',[0.1,0.1,0.7,0.7]);
    
    hscrollbar=uicontrol('style','slider','units','normalized','position',[.1 0.05 0.8 .02],'callback',@hscroll_Callback);   
    
    handles.im = im;
    handles.title_string = title_string;
    handles.hscrollbar = hscrollbar;
    
    guidata(fighdl,handles);
    
    set(hscrollbar, 'min', 1);
    set(hscrollbar, 'max', size(im,3));
    set(hscrollbar, 'Value', 1); % Somewhere between max and min.
    set(hscrollbar, 'SliderStep',1/size(im,3) * [1 10]);
    
    set(fighdl,'CurrentAxes',axes('Position',[0.1,0.1,0.8,0.8]))

    % because MATLAB guis are a MESS

    hhSlider = handle(hscrollbar);
    hProp = findprop(hhSlider,'Value');
    try    % R2014b and newer
       % hProp is a matlab.graphics.internal.GraphicsMetaProperty object
       addlistener(hscrollbar,hProp,'PostSet',@hscroll_Callback);
    catch  % R2014a and older
        try
            % hProp is a schema.prop object
            hListener = handle.listener(hhSlider,hProp,'PropertyPostSet',@hscroll_Callback);
            setappdata(hscrollbar,'sliderListener',hListener);  % this is important - read above
        catch e
            try    % R2013b and older
                addlistener(hscrollbar,'ActionEvent',@hscroll_Callback);
            catch  % R2014a and newer
                addlistener(hscrollbar,'ContinuousValueChange',@hscroll_Callback);
            end
        end
    end

    update_gui(im,1,title_string);
end

function update_gui(im,frame_idx,title_string)
    imagesc(im(:,:,frame_idx))
    
    axis image
    set(gca,'Ydir','Reverse')
    axis off
    colormap gray
    
    title([title_string ' Frame: ' num2str(frame_idx) ' of ' num2str(size(im,3)) ' [Close fig to move on]'])    
    
    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');
end

function hscroll_Callback(varargin)

    handles = guidata(gcf);
    im = handles.im;
    title_string = handles.title_string;

    frame_idx = round(get(handles.hscrollbar,'value'));

    update_gui(im,frame_idx,title_string);
end

