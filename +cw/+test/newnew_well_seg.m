function [ output_args ] = newnew_well_seg( input_args )
%NEWNEW_WELL_SEG Summary of this function goes here
%   Detailed explanation goes here

    im = imread('C5-Microwells_4color-04(10)_s10_BRIGHT-1.tif');
    
    im = mat2gray(im);
    
    im = wiener2(im);
    
    d=5;
    order=2;
    [r, c]=size(im);
    
    im_res = homofil(im,d,r,c,order);
    
%     im_res = mat2gray(im_res);
    
    figure(13048)
    clf

        subtightplot(1,2,1)
            hold all

            imagesc(im)

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
        subtightplot(1,2,2)
            hold all

            imagesc(im_res)

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
end

