function [ output_args ] = well_seg_3( input_args )
%WELL_SEG_3 Summary of this function goes here
%   Detailed explanation goes here

% I = imread('hard_brightfield_rot.tif');
I = zloadim('C4-Farid_test-02(1)-2_BF.tif');
I = I(:,:,220);

It1 = imread('well_templates/wt1.tif');
It2 = imread('well_templates/wt2.tif');
It3 = imread('well_templates/wt4.tif');

I = mat2gray(I);

It1 = mat2gray(It1);
It2 = mat2gray(It2);
It3 = mat2gray(It3);

[~,I_NCC1]=template_matching(It1,I);
[~,I_NCC2]=template_matching(It2,I);
[~,I_NCC3]=template_matching(It3,I);

%%%%%%%

processparam = [1 3];
pad = 3*processparam(2);

im_pad = padarray(I_NCC1,[pad pad],NaN);
im_bpass1 = bpass(im_pad,processparam(1),processparam(2));
im_bpass1 = mat2gray(im_bpass1(pad+1:end-pad,pad+1:end-pad));

im_pad = padarray(I_NCC2,[pad pad],NaN);
im_bpass2 = bpass(im_pad,processparam(1),processparam(2));
im_bpass2 = mat2gray(im_bpass2(pad+1:end-pad,pad+1:end-pad));

im_pad = padarray(I_NCC3,[pad pad],NaN);
im_bpass3 = bpass(im_pad,processparam(1),processparam(2));
im_bpass3 = mat2gray(im_bpass3(pad+1:end-pad,pad+1:end-pad));

im_bpass_final = mat2gray(im_bpass1 + im_bpass2 + im_bpass3);

im_peaks = imclose(im_bpass_final,strel('disk',30));

centers=FastPeakFind(im_peaks,0);            
centers = reshape(centers,2,[])';
idcs = sub2ind(size(im_peaks), centers(:,2), centers(:,1));
values = im_peaks(idcs);
[~,top_idcs] = sort(values,1,'descend');
top_idcs = top_idcs(1:36);
selected_centers = centers(top_idcs,:);

%%%%%%%%%%%%%

% [asdf,ncc_bp1_ncc]=template_matching(nccbp1,im_bpass_final);
% [~,ncc_bp2_ncc]=template_matching(nccbp2,im_bpass_final);
% 
% figure;imagesc(imcomplement(asdf));colormap(gray)
% 
% return



%%%%%%%%%%%%%

% imwrite(im_bpass_final,'ncc_bandpass.tif');

% return
% 
% mask = im_bpass_final > 0.01;
% 
% ibpc = imclose(mask,strel('line',5,0));
% ibpc = imclose(ibpc,strel('line',10,90));
% 
% ibh = imopen(ibpc,strel('line',12,0));
% ibv = imopen(ibpc,strel('line',8,90));
% figure;imagesc(ibh)
% figure;imagesc(ibv)
% figure;imagesc(ibh & ibv)
% return
% im_bpass_final_opened = imopen(,vse);

%%%%%%

figure(13042)
    clf
            hold all

            imagesc(I)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            

figure(13043)
    clf
    
        subtightplot(2,3,1)
            hold all

            imagesc(I_NCC1)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
        subtightplot(2,3,2)
            hold all

            imagesc(I_NCC2)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
        subtightplot(2,3,3)
            hold all

            imagesc(I_NCC3)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
        subtightplot(2,3,4)
            hold all

            imagesc(im_bpass1)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
        subtightplot(2,3,5)
            hold all

            imagesc(im_bpass2)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
        subtightplot(2,3,6)
            hold all

            imagesc(im_bpass3)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
   figure(134265);clf;
    hold all
        subtightplot(1,4,1)
            imagesc(im_bpass_final)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
        subtightplot(1,4,2)
            hold all
            imagesc(im_peaks)
            plot(centers(:,1),centers(:,2),'mo','MarkerSize',20,'LineWidth',4)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
        subtightplot(1,4,3)
            hold all
            imagesc(im_peaks)
            plot(selected_centers(:,1),selected_centers(:,2),'mo','MarkerSize',20,'LineWidth',4)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
        subtightplot(1,4,4)
            hold all
            imagesc(I)
            plot(selected_centers(:,1),selected_centers(:,2),'mo','MarkerSize',20,'LineWidth',4)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
%     figure(13423);clf;hist(im_bpass1(im_bpass1(:)~=0),1000)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

I = imread('hard_brightfield.tif');
I_rot = imrotate(I,2);

imwrite(I_rot,'hard_brightfield_rot.tif')

figure(13042)
    clf
    
        subtightplot(1,2,1)
            hold all

            imagesc(I)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title('pre rot')
            
        subtightplot(1,2,2)
            hold all

            imagesc(I_rot)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
end

