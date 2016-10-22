function test_detect_noise()

    opts.thresh = 1.3;
    opts.do_wiener = 0;
    opts.do_blur = 1;
    opts.blur_sigma = 20;
    opts.blur_hsize = [20 20];

    noise_im = zloadim('movies/single_molecule_background.tif');
    signal_im = zloadim('movies/single_molecule_signal.tif');
    
    % correctly normalize
    maxval = max([noise_im(:)',signal_im(:)']);
    minval = min([noise_im(:)',signal_im(:)']);
    limits = double([minval maxval]);
    
    noise_im = mat2gray(noise_im,limits);
    signal_im = mat2gray(signal_im,limits);
    
    [noise_im_is_noise,noise_im_filtered_list,noise_im_filtered] = detect_noise( noise_im, opts );
    [signal_im_is_noise,signal_im_filtered_list,signal_im_filtered] = detect_noise( signal_im, opts );
    
    max_ent_noise = max(noise_im_filtered_list(:));
    max_ent_signal = max(signal_im_filtered_list(:));
    max_ent = max([max_ent_noise max_ent_signal]);
    
    delete('noise_demo.avi')
    
    writerObj = VideoWriter('noise_demo.avi');
    writerObj.FrameRate = 20;
    writerObj.Quality = 100;
    open(writerObj);

    figure(1342)
    clf
        subplot(2,1,1)
            hold all

            plot(-1,-1,'LineStyle','-','Color','r','LineWidth',6)
            plot(-1,-1,'LineStyle','-','Color','g','LineWidth',6)

            xlim([0 1])
            ylim([0 1])

            axis off

            legend('Pure background','Signal present')
    
        subplot(2,1,2)
            hold all

            plot(-1,-1,'LineStyle','-','Color','r','LineWidth',6)
            plot(-1,-1,'LineStyle','-','Color','g','LineWidth',6)
            plot(-1,-1,'LineStyle','-','Color','b','LineWidth',6)
            plot(-1,-1,'LineStyle','-','Color','m','LineWidth',6)
            xlim([0 1])
            ylim([0 1])

            axis off

            legend('Threshold failure','Threshold success','Signal level','Background level')
    
	set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');
        
    %%%
    
    background_means = [];
    signal_means = [];
    background_stdevs = [];
    signal_stdevs = [];
    
    figure(130892)
    
    for frame_idx = 1:size(noise_im,3)

        cur_noise_filtered = noise_im_filtered(:,:,frame_idx);
        cur_signal_filtered = signal_im_filtered(:,:,frame_idx);

        clf
        
        % Original images
        
            subtightplot(3,2,1,0.05)
                hold all

                imagesc(noise_im(:,:,frame_idx))
                
                colormap gray
                axis image
                set(gca,'Ydir','Reverse')
                axis off
                
                plot_noise_box(noise_im_is_noise(frame_idx))
                
                title('Background image')
                
                freezeColors

            subtightplot(3,2,2,0.05)
                hold all

                imagesc(signal_im(:,:,frame_idx))
                
                colormap gray
                axis image
                set(gca,'Ydir','Reverse')
                axis off
                
                plot_noise_box(signal_im_is_noise(frame_idx))
                
                title('Signal image')
                
                freezeColors

        % Entropy images
              
            subtightplot(3,2,3,0.05)
                hold all

                imagesc(cur_noise_filtered)

                dist = cur_noise_filtered(:);
                if all(dist==0)
%                     colormap([0 0 0])
                else
                    colormap gray
                end
                
                axis image
                set(gca,'Ydir','Reverse')
                axis off
                
                plot_noise_box(noise_im_is_noise(frame_idx))
                
                title('Entropy-filtered background image')
                
                freezeColors

            subtightplot(3,2,4,0.05)
                hold all

                imagesc(cur_signal_filtered)
                
                dist = cur_signal_filtered(:);
                if all(dist==0)
                    
%                     colormap([0 0 0])
                else
                    colormap gray
                end
                
                axis image
                set(gca,'Ydir','Reverse')
                axis off
                
                plot_noise_box(signal_im_is_noise(frame_idx))

                title('Entropy-filtered signal image')
                
                freezeColors
                
        % Entropy histogram

            subtightplot(3,1,3,0.05)
                hold all
                
                [noise_hist,xnoise] = hist(noise_im_filtered_list(frame_idx,:),50);
                [signal_hist,xsignal] = hist(signal_im_filtered_list(frame_idx,:),50);

                noise_hist = noise_hist ./ max(noise_hist);
                signal_hist = signal_hist ./ max(signal_hist);
                
                mean_noise = mean(noise_im_filtered_list(frame_idx,:));
                mean_signal = mean(signal_im_filtered_list(frame_idx,:));
                
                std_noise = std(noise_im_filtered_list(frame_idx,:));
                std_signal = std(signal_im_filtered_list(frame_idx,:));
                
                background_means = [background_means, mean_noise];
                signal_means = [signal_means, mean_signal];
                
                background_stdevs = [background_stdevs, std_noise];
                signal_stdevs = [signal_stdevs, std_signal];
                
                plot(xnoise,noise_hist,'Color','m','LineWidth',3)
                plot(xsignal,signal_hist,'Color','b','LineWidth',3)
                
                patch_min = min([mean_noise,mean_signal]);
                patch_max = max([mean_noise,mean_signal]);
                
                patch_x = [patch_min,patch_max,patch_max,patch_min];
                patch_y = [0,0,1,1];
                
                patch(patch_x,patch_y,[0.7 0.7 0.7],'EdgeColor','none');
                
                if mean_noise > opts.thresh || mean_signal < opts.thresh
                    line_color = 'r';
                else
                    line_color = 'g';
                end
                
                line([opts.thresh opts.thresh],ylim,'LineStyle','-','Color',line_color,'LineWidth',6)
                line([mean_noise mean_noise],ylim,'LineStyle','--','Color','m','LineWidth',4)
                line([mean_signal mean_signal],ylim,'LineStyle','--','Color','b','LineWidth',4)
                
                xlim([0 max_ent])
                ylim([0 1]);
                
                xlabel('Entropy [bits]')
                ylabel('Normalized counts')
                title('Histogram of entropy')
                
                pos = get(gca,'OuterPosition');
                set(gca,'OuterPosition',[pos(1) + 0.01, pos(2) + 0.01, pos(3), pos(4)])
                
                box on
                grid on
                
                set(gca,'layer','top')
                
        suptitle(['Signal vs. background discrimination test. Frame: ' num2str(frame_idx) ' of ' num2str(size(noise_im,3))])
        
        set(findall(gcf,'type','text'),'fontSize',16,'fontWeight','bold')
        set(findall(gcf,'type','axes'),'fontSize',16,'fontWeight','bold','LineWidth',3)
        set(gcf, 'color', 'white');
        
        drawnow
        
        F=getframe(gcf);

        writeVideo(writerObj,F);
    end
    
    close(writerObj);
    
    figure(13542)
    clf
        subplot(1,2,1)
            hold all
            
            best_mean_thresh_dist = mean([background_means;signal_means],1);
            best_mean_thresh = mean(best_mean_thresh_dist);
            [best_thresh_hist,x0] = hist(best_mean_thresh_dist,35);

            best_thresh_hist = best_thresh_hist ./ max(best_thresh_hist);
            
            plot(x0,best_thresh_hist,'Color','k','LineWidth',3)
%             line([max(background_means) max(background_means)],ylim,'LineStyle','--','Color','r','LineWidth',3)
%             line([min(signal_means) min(signal_means)],ylim,'LineStyle','--','Color','r','LineWidth',3)
            line([best_mean_thresh,best_mean_thresh],ylim,'LineStyle','--','Color','g','LineWidth',3)
            
            xlabel('Threshold level')
            ylabel('Normalized counts')
            title('Optimal mean threshold')
            
            box on
            grid on
            
            pos = get(gca,'OuterPosition');
            set(gca,'OuterPosition',[pos(1), pos(2), pos(3), 0.9*pos(4)])
            
        subplot(1,2,2)
            hold all

            best_mean_thresh_dist = mean([background_stdevs;signal_stdevs],1);
            best_mean_thresh = mean(best_mean_thresh_dist);
            [best_thresh_hist,x0] = hist(best_mean_thresh_dist,35);
            
            best_thresh_hist = best_thresh_hist ./ max(best_thresh_hist);

            plot(x0,best_thresh_hist,'Color','k','LineWidth',3)
            line([best_mean_thresh,best_mean_thresh],ylim,'LineStyle','--','Color','g','LineWidth',3)
            
            xlabel('Threshold level')
            ylabel('Normalized counts')
            title('Optimal stdev threshold')
            
            box on
            grid on
            
            pos = get(gca,'OuterPosition');
            set(gca,'OuterPosition',[pos(1), pos(2), pos(3), 0.9*pos(4)])
        
    suptitle('Optimal threshold levels')
        
    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');
end
