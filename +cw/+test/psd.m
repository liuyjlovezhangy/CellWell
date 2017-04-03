function [ output_args ] = Untitled( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     im = imread('em_watershed_test.tif');
    im = imread('C2-Timelapse01.ome.tiff kept stack.tif');
%     im = imread('C3-Timelapse01.ome.tiff kept stack.tif');
    
    im = mat2gray(im);
    
%     im = wiener2(im);

%     back_sz = 10;
%     N = 2*back_sz + 1;
%     bx = zeros(1,N)  + 1/N;
%     by = bx';
%     
%     box1 = conv2(im,bx,'same');
%     box_final = conv2(box1,by,'same');
% 
%     im = im - box_final;
    
    [pstruct, mask, imgLM, imgLoG] = pointSourceDetection(im, 3, 'Alpha', 0.0001);
    
    imgLoG = mat2gray(imgLoG);
    
    pstruct.A
    
    im(imgLM>0) = 0;    
    imgLoG(imgLM>0) = 0;   
    imgLM(imgLM>0) = 1;
    
    comp_im = cat(3,imgLM,im,zeros(size(im)));
    comp_im_log = cat(3,imgLM,imgLoG,zeros(size(im)));
    
    figure(1033);clf;hold all;imagesc(im);colormap gray;axis image;plot(pstruct.x,pstruct.y,'mo','MarkerSize',20)
    figure(1034);clf;hold all;imagesc(mask);colormap gray;axis image;plot(pstruct.x,pstruct.y,'mo','MarkerSize',20)
    figure(1035);clf;hold all;imagesc(comp_im);colormap gray;axis image;plot(pstruct.x,pstruct.y,'mo','MarkerSize',20)
    figure(1036);clf;hold all;imagesc(comp_im_log);colormap gray;axis image;plot(pstruct.x,pstruct.y,'mo','MarkerSize',20)
end

