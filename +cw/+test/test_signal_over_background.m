function [ output_args ] = test_signal_over_background( input_args )
%TEST_SIGNAL_OVER_BACKGROUND Summary of this function goes here
%   Detailed explanation goes here

    im_sig1 = imread('signal1.tif');
    im_sig2 = imread('signal2.tif');
    im_sig3 = imread('signal3.tif');
    im_sig4 = imread('signal4.tif');
    
    im_bknd1 = imread('background1.tif');
    im_bknd2 = imread('background2.tif');
    
    min_int = double(min([im_sig1(:)', im_sig2(:)', im_sig3(:)', im_sig4(:)', im_bknd1(:)', im_bknd2(:)']))
    max_int = double(max([im_sig1(:)', im_sig2(:)', im_sig3(:)', im_sig4(:)', im_bknd1(:)', im_bknd2(:)']))
    
    im_sig1 = mat2gray(im_sig1,[min_int,max_int]);
    im_sig2 = mat2gray(im_sig2,[min_int,max_int]);
    im_sig3 = mat2gray(im_sig3,[min_int,max_int]);
    im_sig4 = mat2gray(im_sig4,[min_int,max_int]);
    
    im_bknd1 = mat2gray(im_bknd1,[min_int,max_int]);
    im_bknd2 = mat2gray(im_bknd2,[min_int,max_int]);
    
    
    [sig1_hist,sig1_hist_loc] = histogram(im_sig1(:));
    [sig2_hist,sig2_hist_loc] = histogram(im_sig2(:));
    [sig3_hist,sig3_hist_loc] = histogram(im_sig3(:));
    [sig4_hist,sig4_hist_loc] = histogram(im_sig4(:));
    
    [bknd1_hist,bknd1_hist_loc] = histogram(im_bknd1(:));
    [bknd2_hist,bknd2_hist_loc] = histogram(im_bknd2(:));
    
%     p = kruskalwallis([im_sig1(:)',im_bknd1(:)'],[ones(1,numel(im_sig1)),2*ones(1,numel(im_bknd1))])
%     p = kruskalwallis([im_bknd1(:)',im_bknd2(:)'],[ones(1,numel(im_bknd1)),2*ones(1,numel(im_bknd2))])
    
    p1 = ranksum(im_sig1(:),im_bknd1(:))
    p2 = ranksum(im_sig2(:),im_bknd1(:))
    p3 = ranksum(im_sig3(:),im_bknd1(:))
    p4 = ranksum(im_sig4(:),im_bknd1(:))
    
    pb = ranksum(im_bknd1(:),im_bknd2(:))
    
    figure(21342)
    clf
    
        subtightplot(2,6,1)
            hold all
            
            imagesc(im_sig1)
            
            axis image
            set(gca,'Ydir','Reverse')

            set(gca,'Color','white')
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            
        subtightplot(2,6,7)
            hold all
            
            plot(sig1_hist_loc,sig1_hist,'LineWidth',4)
            plot(bknd1_hist_loc,bknd1_hist,'LineWidth',4)
            
            %xlim([0 1])
            
        %%%
            
        subtightplot(2,6,2)
            hold all
            
            imagesc(im_sig2)
            
            axis image
            set(gca,'Ydir','Reverse')

            set(gca,'Color','white')
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            
        subtightplot(2,6,8)
            hold all
            
            plot(sig2_hist_loc,sig2_hist,'LineWidth',4)
            plot(bknd1_hist_loc,bknd1_hist,'LineWidth',4)
            
            %xlim([0 1])
            
        %%%
            
        subtightplot(2,6,3)
            hold all
            
            imagesc(im_sig3)
            
            axis image
            set(gca,'Ydir','Reverse')

            set(gca,'Color','white')
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            
        subtightplot(2,6,9)
            hold all
            
            plot(sig3_hist_loc,sig3_hist,'LineWidth',4)
            plot(bknd1_hist_loc,bknd1_hist,'LineWidth',4)
            
            %xlim([0 1])
            
        %%%
            
        subtightplot(2,6,4)
            hold all
            
            imagesc(im_sig4)
            
            axis image
            set(gca,'Ydir','Reverse')

            set(gca,'Color','white')
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            
        subtightplot(2,6,10)
            hold all
            
            plot(sig4_hist_loc,sig4_hist,'LineWidth',4)
            plot(bknd1_hist_loc,bknd1_hist,'LineWidth',4)
            
            %xlim([0 1])
            
        %%%
            
        subtightplot(2,6,5)
            hold all
            
            imagesc(im_bknd1)
            
            axis image
            set(gca,'Ydir','Reverse')

            set(gca,'Color','white')
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            
        subtightplot(2,6,11)
            hold all
            
            plot(bknd1_hist_loc,bknd1_hist,'LineWidth',4)
            plot(bknd1_hist_loc,bknd1_hist,'LineWidth',4)
            
            %xlim([0 1])
            
        %%%
            
        subtightplot(2,6,6)
            hold all
            
            imagesc(im_bknd2)
            
            axis image
            set(gca,'Ydir','Reverse')

            set(gca,'Color','white')
            set(gca,'XTick',[])
            set(gca,'YTick',[])
        
        subtightplot(2,6,12)
            hold all
            
            plot(bknd2_hist_loc,bknd2_hist,'LineWidth',4)
            plot(bknd1_hist_loc,bknd1_hist,'LineWidth',4)
            
            %xlim([0 1])
            
        %%%
            
    colormap gray
            
    set(findall(gcf,'type','text'),'fontSize',16,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',16,'fontWeight','bold','LineWidth',5)
    set(gcf, 'color', 'white');
end

