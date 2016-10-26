function answer = confirm_results( im, titlestr )
    fighdl = cw.plot.tiff_gui(im,titlestr);
    waitfor(fighdl);
    
    answer = questdlg('Are you satisfied with the results?');
end

