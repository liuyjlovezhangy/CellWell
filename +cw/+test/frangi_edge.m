function [ output_args ] = Untitled( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    im = imread('em_watershed_test.tif');
    im = mat2gray(im);

    %   Options : Struct with input options,
%       .FrangiScaleRange : The range of sigmas used, default [1 8]
%       .FrangiScaleRatio : Step size between sigmas, default 2
%       .FrangiBetaOne : Frangi correction constant, default 0.5
%       .FrangiBetaTwo : Frangi correction constant, default 15
%       .BlackWhite : Detect black ridges (default) set to true, for
%                       white ridges set to false.
%       .verbose : Show debug information, default true
    
    options.FrangiScaleRange = [1 8];
    options.FrangiScaleRatio = 0.12;
    options.FrangiBetaOne = 1;
    options.FrangiBetaTwo = 40;
    options.BlackWhite = 0;

    [outIm,whatScale,Direction] = FrangiFilter2D(im, options);
    
    figure(1033);clf;imagesc(im);colormap gray;axis image
    figure(1034);clf;imagesc(outIm);colormap gray;axis image
end

