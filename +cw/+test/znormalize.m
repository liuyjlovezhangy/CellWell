function [ output_args ] = test_znormalize( input_args )
%TEST_ZNORMALIZE Summary of this function goes here
%   Detailed explanation goes here

    im = zloadim('movies/Timelapse01.ome.tiff');
    
    im1 = znormalize(im);
    im2 = znormalize(im,4);
    
    
end

