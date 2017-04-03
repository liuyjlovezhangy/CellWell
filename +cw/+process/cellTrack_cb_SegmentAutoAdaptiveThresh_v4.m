% Author: Simon Gordonov 
% Date: 10/16
% Automated cell body segmentation using adaptive thresholding

%Version 2 addting alternative ways to get Totsu

function [processedFrames,repeat] = ...
    cellTrack_cb_SegmentAutoAdaptiveThresh_v4(processedFrames,IcbNoTophat,params,currFrameNum,frameNumsToSeg,repeat)

    %Define the nuclei centers
    maskNucShrink = logical(full(processedFrames.nucMask{currFrameNum}));
    
    %Scale image
    IcbNoTophat = mat2gray(IcbNoTophat);
    
    %Perform background subtraction using tophat:
    if params.cbSeg.tophatDiskR ~= 0
        Icb = imtophat(IcbNoTophat,strel('disk',params.cbSeg.tophatDiskR));
    else
        Icb = IcbNoTophat;
    end
       
    %%%%%%%%%%%% PARAMETER SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ismember(currFrameNum,frameNumsToSeg) || ~processedFrames.goodFrame(currFrameNum-1) %if start of pre-treat, treat, or washout, scane a larger range of params, or previous frame bad
        
        %Find elbow of cumulative histogram of intensity values
        [h,binCent] = hist(Icb(:),1000);
        hCS = cumsum(h);
        elbowIntHist = binCent(elbowOfCurve(hCS));
%         cbMask = im2bw(Icb,elbowIntHist);
%         cbMask = imreconstruct(maskNucShrink,cbMask);
%         figure(1)
%         imshowpair(Icb,bwperim(cbMask)|maskNucShrink)
%         
        %Find initial Otsu threshold and set range around it for segmentation
%         TotsuActual = graythresh(Icb);
        
        %Set starting value of Otsu threshold
        TotsuMin = 0.001;
        TotsuMax = 0.8;
        %Set Otsu threshold increment
%         TotsuIncr = 0.005;

        %Set number of Otsu thresholds to sweep
        TotsuNumVals = 32;
        
        a = linspace(0,1,TotsuNumVals/2);
        y = exp(-5*a);
    	y = (y-min(y))/range(y);
        
        y1 = y.*(abs(elbowIntHist-TotsuMin));
        y1 = elbowIntHist-y1;
        y2 = y.*(abs(TotsuMax-elbowIntHist));
        y2 = elbowIntHist+y2;
        
        TotsuVec = [y1(1:end-1),fliplr(y2)];
        
%                 plot(TotsuVec,ones(1,numel(TotsuVec)),'.')

%                 plot(y1,ones(1,numel(y)),'.')

        
        %Set vector of Totsu
%         TotsuVec = linspace(Totsu,TotsuActual+(0.1*TotsuActual),TotsuNumVals);%Totsu:TotsuIncr:(Totsu+(TotsuIncr*(TotsuNumVals-1)));
        
        TotsuNumVals = numel(TotsuVec);

        %EDGES
        
        %Find initial Canny edge threshold and set range around it for segmentation
        [~,TedgeActual] = edge(IcbNoTophat,'canny');
        %Set starting values of Canny edge detection thresholds
        TedgeMin = [TedgeActual(1)-0.05*TedgeActual(1),TedgeActual(2)-0.05*TedgeActual(2)];
        TedgeMax = [TedgeActual(1)+0.2*TedgeActual(1),TedgeActual(2)+0.2*TedgeActual(2)];
%         Tedge = [0.01,0.03];
        
        %Set edge threshold increment
%         TedgeIncr = 0.01;

        %Set number of edge thresholds to sweep
        TedgeNumVals = 5;
        
        %Set vector of Tedge
        TedgeVecLow = linspace(TedgeMin(1),TedgeMax(1),TedgeNumVals);%Totsu:TotsuIncr:(Totsu+(TotsuIncr*(TotsuNumVals-1)));
        TedgeVecHigh = linspace(TedgeMin(2),TedgeMax(2),TedgeNumVals);%Totsu:TotsuIncr:(Totsu+(TotsuIncr*(TotsuNumVals-1)));
        TedgeVec = [TedgeVecLow',TedgeVecHigh'];
       
%         TedgeVecLow = Tedge(1):TedgeIncr:Tedge(1)+(TedgeIncr*TedgeNumVals);
%         TedgeVecHigh = Tedge(2):TedgeIncr:Tedge(2)+(TedgeIncr*TedgeNumVals);
%         TedgeVec = [TedgeVecLow',TedgeVecHigh'];

        
    else
        
        %Set starting value of Otsu threshold
%         Totsu = 0.001;
% 
%         %Set Otsu threshold increment
%         TotsuIncr = 0.005;

        %Set number of Otsu thresholds to sweep
        TotsuNumVals = 9;
        
%         if processedFrames.TotsuCB(currFrameNum-1) < Totsu
            d = processedFrames.TotsuCB(currFrameNum-1)/TotsuNumVals;
            TotsuVec = processedFrames.TotsuCB(currFrameNum-1)-(d*floor(TotsuNumVals/2)):d:processedFrames.TotsuCB(currFrameNum-1)+(d*TotsuNumVals/2);%Totsu:TotsuIncr:(Totsu+(TotsuIncr*(TotsuNumVals-1)));
%         else
%             TotsuVec = (processedFrames.TotsuCB(currFrameNum-1) - (TotsuIncr*floor(TotsuNumVals/2))):TotsuIncr:(processedFrames.TotsuCB(currFrameNum-1) + (TotsuIncr*floor(TotsuNumVals/2)));
            TotsuVec(TotsuVec > 1) = [];
            TotsuVec(TotsuVec < 0) = [];
%         end
        
        TotsuNumVals = numel(TotsuVec);

%         %Set starting values of Canny edge detection thresholds
%         TedgeLow = max([0.01,(processedFrames.TedgeCB{currFrameNum-1}(1)-(0.01))]);
%         TedgeHigh = max([0.03,(processedFrames.TedgeCB{currFrameNum-1}(2)-(0.01))]);
%         Tedge = [TedgeLow,TedgeHigh];
%         
%         %Set edge threshold increment
%         TedgeIncr = 0.01;
% 
%         %Set number of edge thresholds to sweep
%         TedgeNumVals = 5;

       %Set number of Edge thresholds to sweep
%         TedgeNumVals = 5;
%         
%         dLow = processedFrames.TedgeCB{currFrameNum-1}(1)/TedgeNumVals;
%         TedgeVecLow = processedFrames.TedgeCB{currFrameNum-1}(1)-(dLow*floor(TedgeNumVals/2)):dLow:processedFrames.TedgeCB{currFrameNum-1}(1)+(dLow*floor(TedgeNumVals/2));
%         dHigh = processedFrames.TedgeCB{currFrameNum-1}(2)/TedgeNumVals;
%         TedgeVecHigh= processedFrames.TedgeCB{currFrameNum-1}(2)-(dHigh*floor(TedgeNumVals/2)):dHigh:processedFrames.TedgeCB{currFrameNum-1}(2)+(dHigh*floor(TedgeNumVals/2));
%         
%         TedgeVec = [TedgeVecLow',TedgeVecHigh'];
%         
%         TedgeNumVals = numel(TedgeVecHigh);
        
%         %Set starting values of Canny edge detection thresholds
%         Tedge = [0.01,0.03];
% 
%         %Set edge threshold increment
%         TedgeIncr = 0.01;
% 
%         %Set number of edge thresholds to sweep
%         TedgeNumVals = 5;
        
%         %Set starting values of Canny edge detection thresholds
        [~,TedgeActual] = edge(IcbNoTophat,'canny');
        TedgeMin = [TedgeActual(1)-0.1*TedgeActual(1),TedgeActual(2)-0.1*TedgeActual(2)];
        TedgeMax = [TedgeActual(1)+0.2*TedgeActual(1),TedgeActual(2)+0.2*TedgeActual(2)];

        %Set number of edge thresholds to sweep
        TedgeNumVals = 5;
        
        %Set vector of Tedge
        TedgeVecLow = linspace(TedgeMin(1),TedgeMax(1),TedgeNumVals);%Totsu:TotsuIncr:(Totsu+(TotsuIncr*(TotsuNumVals-1)));
        TedgeVecHigh = linspace(TedgeMin(2),TedgeMax(2),TedgeNumVals);%Totsu:TotsuIncr:(Totsu+(TotsuIncr*(TotsuNumVals-1)));
        TedgeVec = [TedgeVecLow',TedgeVecHigh'];
%                 
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Allocate cell arrays to store masks of Otsu only, Edge only, combo
    %masks, area of foreground, intnsities of foreground
    otsuMasks = cell(TotsuNumVals,1);
    edgeMasks = cell(TedgeNumVals,1);
    comboMasks = cell(TotsuNumVals,TedgeNumVals);
    areas = zeros(TotsuNumVals,TedgeNumVals);
    intens = zeros(TotsuNumVals,TedgeNumVals);
    bkgMeanIntens = zeros(TotsuNumVals,TedgeNumVals);
    diffAreaPrevFrame = zeros(TotsuNumVals,TedgeNumVals);
    
    %Perform Otsu segmentations over the threshold ranges specified
    for i = 1:TotsuNumVals
        otsuMasks{i} = im2bw(Icb,TotsuVec(i)); 
    end
    
    %Perform Edge detections over the threshold ranges specified
    for i = 1:TedgeNumVals
        edgeMasks{i} = edge(IcbNoTophat,'canny',TedgeVec(i,:));%Tedge + ((i-1)*TedgeIncr));
    end
    
    %Perform final segmentation using the Otsu and Edge masks
    
    for i = 1:TotsuNumVals
        for j = 1:TedgeNumVals
            
%             %Combine the Otsu foreground with the Canny edges that touch the the
%             %foreground:
%             cbMask = otsuMasks{i} | edgeMasks{j};
% 
%             %Perform reconstruction of cell body channel objects based on cell nuclei points
%             %used as markers. Use actin channel with filled holes to ensure the
%             %Euler number is due to nuclei present and not cells creating
%             %"closed rings":
%             cbMask = imreconstruct(maskNucShrink,cbMask);
%             
%             %Perform additional tweaks as needed.
%             comboMasks{i,j} = imclose(cbMask,strel('disk',2));
% 
%             %Compute area of foreground pixels
%             areas(i,j) = sum(comboMasks{i,j}(:));
%             
%             %Compute intensities of foreground pixels
%             intens(i,j) = sum(Icb(comboMasks{i,j}));
%             
%             %Mean intensity of the background
%             bkgMeanIntens(i,j) = mean(Icb(~comboMasks{i,j}));
            
            % OLD
            %Combine the Otsu foreground with the Canny edges that touch the the
            %foreground:
            cbMask = otsuMasks{i} | imreconstruct(otsuMasks{i},edgeMasks{j});

            %Perform additional tweaks as needed.
            cbMask = imclose(cbMask,strel('disk',2));

            %Create a filled holes cell body channel used for identifying
            %single nontouching cells later:
            cbMask = imfill(cbMask,'holes');
            
            %Perform additional tweaks as needed.
%             cbMask = imerode(cbMask,strel('disk',1));

            %Remove objects that are smaller than minObjSize # of pixels:
%             cbMaskFillHoles = bwareaopen(cbMaskFillHoles,params.cbSeg.minObjSize);

            %Perform reconstruction of cell body channel objects based on cell nuclei points
            %used as markers. Use actin channel with filled holes to ensure the
            %Euler number is due to nuclei present and not cells creating
            %"closed rings":
            comboMasks{i,j} = imreconstruct(maskNucShrink,cbMask);
            
            %Compute area of foreground pixels
            areas(i,j) = sum(comboMasks{i,j}(:));
            
            %Compute intensities of foreground pixels
            intens(i,j) = sum(Icb(comboMasks{i,j}));
            
            %Mean intensity of the background
            bkgMeanIntens(i,j) = mean(Icb(~comboMasks{i,j}));
            
            %Compute deviation of area from previous frame
            if currFrameNum > 1
%                 diffAreaPrevFrame(i,j) = abs(1-(processedFrames.Area{currFrameNum-1}/areas(i,j)));
                diffAreaPrevFrame(i,j) = max([processedFrames.Area{currFrameNum-1}/areas(i,j),areas(i,j)/processedFrames.Area{currFrameNum-1}]);
            end

        end
    end
    
    bkgMeanIntens(isnan(bkgMeanIntens)) = 0;
    
    %Compute difference in areas and change intensities of foreground
    %pixels
%     dAreasOtsu = diff(areas,1,1);
%     dAreasOtsu = [dAreasOtsu;zeros(1,TedgeNumVals)];
%     dAreasEdge = diff(areas,1,2);
%     dAreasEdge = [dAreasEdge,zeros(TotsuNumVals,1)];
% 
%     dIntensOtsu = diff(intens,1,1);
%     dIntensOtsu = [dIntensOtsu;zeros(1,TedgeNumVals)];
%     dIntensEdge = diff(intens,1,2);
%     dIntensEdge = [dIntensEdge,zeros(TotsuNumVals,1)];

    %Find elbow point of areas with varying intensity thresholds
    elbowAreas = zeros(1,TedgeNumVals);
    if ismember(currFrameNum,frameNumsToSeg) || ~processedFrames.goodFrame(currFrameNum-1) %if current frame has close previous frame threshold setting, don't use curve elbow
        for i = 1:TedgeNumVals
%             secDer = (areas(3:end,i)-areas(1:end-2,i))./(2*areas(2:end-1,i));
%             secDer(secDer == -Inf) = NaN;
%             secDer(secDer == Inf) = NaN;
%             [~,tempMinLoc] = min(secDer); 
%             tempMinLoc = min(tempMinLoc);
            elbowAreas(i) = elbowOfCurve(areas(:,1)');%tempMinLoc + 1;
        end
    end
    
    %Area changes Otsu
    dAreasOtsu = Inf(TotsuNumVals-1,TedgeNumVals);
    for i = 1:TotsuNumVals-1
        for j = 1:TedgeNumVals
            if i >= elbowAreas(j)
                dAreasOtsu(i,j) = areas(i,j)/areas(i+1,j);
            end
        end
    end
    dAreasOtsu = [dAreasOtsu;Inf(1,TedgeNumVals)];

    %Intensity changes Otsu
    dIntensOtsu = Inf(TotsuNumVals-1,TedgeNumVals);
    for i = 2:TotsuNumVals-1
        for j = 1:TedgeNumVals
            if i >= elbowAreas(j)
                if dAreasOtsu(i,j) == 1 %if no area change
                    dIntensOtsu(i,j) = 0;
                else
%                     dIntensOtsu(i,j) = bkgMeanIntens(i,j)/((intens(i,j)-intens(i+1,j))/(areas(i,j)-areas(i+1,j)));
                    dIntensOtsu(i,j) = ((intens(i,j)-intens(i-1,j))/(areas(i,j)-areas(i-1,j)))/bkgMeanIntens(i-1,j);
                end
            end
        end
    end
    dIntensOtsu = [dIntensOtsu;Inf(1,TedgeNumVals)];

    %Area change Edge
    dAreasEdge = Inf(TotsuNumVals,TedgeNumVals-1);
    for i = 1:TotsuNumVals
        for j = 1:(TedgeNumVals-1)
            if i >= elbowAreas(j)
                dAreasEdge(i,j) = areas(i,j)/areas(i,j+1);
            end
        end
    end
    dAreasEdge = [dAreasEdge,Inf(TotsuNumVals,1)];
    
    %Intensity change Edge
%     dIntensEdge = Inf(TotsuNumVals,TedgeNumVals-1);
%     for i = 1:TotsuNumVals
%         for j = 2:TedgeNumVals
%             if i >= elbowAreas(j)
%                 if dAreasEdge(i,j) == 1 %if no intensity change
%                     dIntensEdge(i,j) = 0;
%                 else
% %                     dIntensEdge(i,j) = bkgMeanIntens(i,j)/((intens(i,j)-intens(i,j+1))/(areas(i,j)-areas(i,j+1)));
%                       dIntensEdge(i,j) = ((intens(i,j)-intens(i,j-1))/(areas(i,j)-areas(i,j-1)))/bkgMeanIntens(i,j-1);
%                 end
%             end
%         end
%     end
%     dIntensEdge = [dIntensEdge,Inf(TotsuNumVals,1)];
    
%     InucMaskStack = false(processedFrames.frameSize(1),processedFrames.frameSize(2),TotsuNumVals);
%     IcbMaskStack = false(processedFrames.frameSize(1),processedFrames.frameSize(2),TotsuNumVals);
%     IcbStack = zeros(processedFrames.frameSize(1),processedFrames.frameSize(2),TotsuNumVals);
%     for i = 1:TotsuNumVals
%         InucMaskStack(:,:,i) = imdilate(maskNucShrink,ones(3,3));
%         IcbMaskStack(:,:,i) = comboMasks{i,2};
%         IcbStack(:,:,i) = Icb;
%     end
%     imseriesmaskshow(IcbStack,{IcbMaskStack,InucMaskStack},'displayRange',[0 1],'maskAlphas',[0.2,0.2],'maskColors',[0 1 0;1 0 0])
% % % % % % %     
%     figure;
%     cmap = linspecer(TedgeNumVals);
%     for i = 1:TedgeNumVals
%         plot(areas(:,i),'o-','color',cmap(i,:));
%         hold on;
%         plot(elbowAreas(i),areas(elbowAreas(i),i),'b.','markersize',20);
%     end
%     
%     figure;
%     cmap = linspecer(TedgeNumVals);
%     for i = 1:TedgeNumVals
%         plot(dIntensOtsu(:,i),'color',cmap(i,:));
%         hold on;
%     end
    
    
    %Compute the score for the combination of Otsu and Edge thresholds
    if currFrameNum == 1 || ~processedFrames.goodFrame(currFrameNum-1)
        
%         figure;
%         cmap = linspecer(TedgeNumVals);
%         for i = 1:TedgeNumVals
%             plot(areas(:,i),'o-','color',cmap(i,:));
%             hold on;
%             plot(elbowAreas(i),areas(elbowAreas(i),i),'b.','markersize',20);
%         end

        score = dAreasOtsu + dAreasEdge;% + dIntensOtsu + dIntensEdge;
        score(isnan(score)) = Inf;
        repeat = 0;
        
    else
        
        score = dAreasOtsu + dAreasEdge + diffAreaPrevFrame; %+ dIntensOtsu;% + dIntensEdge;
        score(diffAreaPrevFrame>1.5) = Inf;
        score(isnan(score)) = Inf;
        
        if repeat == 1 && min(score(:)) == Inf
            score = dAreasOtsu + dIntensOtsu + dAreasEdge; %+ dIntensEdge;
            score(isnan(score)) = Inf;
        end
        repeat = 0;
        if min(score(:)) == Inf
%             disp(['ALL INFINITY IN SCORE MATRIX FOR FRAME #',num2str(currFrameNum)])
            repeat = 1;
            
%             if ~exist(fullfile(params.parentFolderForAnalysis,'ElbowPlotsProblemFields'),'dir')
%                 mkdir(fullfile(params.parentFolderForAnalysis,'ElbowPlotsProblemFields'));
%             end
%             
%             figure;
%             set(gcf,'visible','off')
%             cmap = linspecer(TedgeNumVals);
%             for i = 1:TedgeNumVals
%                 plot(areas(:,i),'o-','color',cmap(i,:));
%                 hold on;
%                 plot(elbowAreas(i),areas(elbowAreas(i),i),'b.','markersize',20);
%             end
%             saveas(gcf,fullfile(params.parentFolderForAnalysis,'ElbowPlotsProblemFields',['ElbowPlot_',num2str(fieldNum),'_',num2str(currFrameNum),'.tif']));
            
            return;
        end
    end
    
        %Find best combo of thresholds
    if currFrameNum == 1 || ~processedFrames.goodFrame(currFrameNum-1) || ismember(currFrameNum,frameNumsToSeg)
        [r,c] = find(imregionalmin(score));
        locMinr = find(r==min(r));
%         locMinc = c(locMinr);
        rBest = r(locMinr(1));
        cBest = c(locMinr(1));
    else
        tempMin = imregionalmin(score);
        score(~tempMin) = Inf;
        [rBest,cBest] = find(score==min(score(:)));
        rBest = rBest(end);
        cBest = cBest(end);
    end
    
    cbMaskNucRec = comboMasks{rBest,cBest};
    
    %Segment image using final set of selected thresholds
%     TotsuBest = Totsu + ((rBest-1)*TotsuIncr);
%     TedgeBest = Tedge + ((cBest-1)*TedgeIncr);
%     cbMask = run_cb_segmentation_Final_AdaptiveThresh(Icb,maskNucShrink,TotsuBest,TedgeBest,params);

    %Remove all spur pixels
    maskOld = cbMaskNucRec;
    maskNew = bwmorph(cbMaskNucRec,'spur');
    while 1
        if all(maskOld(:) == maskNew(:))
            break
        else
            maskOld = maskNew;
            maskNew = bwmorph(maskOld,'spur');
        end
    end
    cbMaskNucRec = maskNew;
    
    %Fill holes
%     cbMaskNucRecFillHoles = imfill(cbMaskNucRec,'holes');
    cbMaskNucRecFillHoles = cbMaskNucRec;
    
    %Remove objects that are smaller than minObjSize # of pixels:
    cbMaskNucRecFillHoles = bwareaopen(cbMaskNucRecFillHoles,params.cbSeg.minObjSize);
    
    %Use nuclei markers as holes in the cell body channel objects, and
    %compute euler number per cell body channel object. If euler number is
    %not zero (i.e. cells are touching), then remove them. This step
    %retains only those cells that do not touch each other.
    cbMaskRecWithNucHoles = cbMaskNucRecFillHoles.*imcomplement(maskNucShrink);
    
    CCcb = bwconncomp(cbMaskRecWithNucHoles);

    %Extract Euler number of each cell body channel object with nuclei
    %holes:
    CCcbEulerNum = regionprops(CCcb,'EulerNumber'); %structure
    CCcbEulerNum = [CCcbEulerNum.EulerNumber]'; %convert to array

    %Loop through the cell body objects and set pixels of cell body to
    %2 if cell DOES NOT touch others, and leave to 1 if it touches
    %others.
    for idx = 1:numel(CCcbEulerNum)
        if CCcbEulerNum(idx) == 0 
            cbMaskRecWithNucHoles(CCcb.PixelIdxList{idx}) = 2;
        end
    end
    
    %Fill holes in objects and output final cell body
    %mask (logical):
    cbMask = imfill(cbMaskRecWithNucHoles,'holes');
%     cbMask = cbMaskRecWithNucHoles;
    
    %Close holes smaller than nuclei size
%     maskBkg = ~cbMask;
%     Itemp = bwareaopen(maskBkg,params.nucSeg.minObjSize);
%     cbMask = cbMask | xor(maskBkg,Itemp);
    
    %Set the objects in cbMask that touch the border to be 3's (recall 1's
    %in cbMask mean objects touch each other, while 2's mean those cells
    %don't touch each other).
    temp = cbMask & ~imclearborder(cbMask);
    cbMask(temp==1) = 3; 
    
%     cbMask = cbMask .* cbMaskNucRec;
        
    %Set data for output
    processedFrames.cbMask{currFrameNum} = cbMask; 
    processedFrames.nucMask{currFrameNum} = double(maskNucShrink) .* processedFrames.cbMask{currFrameNum};
    processedFrames.TotsuCB(currFrameNum) = TotsuVec(rBest);
    processedFrames.TedgeCB{currFrameNum} = TedgeVec(cBest,:);%Tedge + ((cBest-1)*TedgeIncr);
    processedFrames.Area{currFrameNum} = sum(logical(cbMask(:)));
    
%     figure(1);
%     imshowpair(Icb,bwperim(cbMask));
%     title(['Frame number ',num2str(currFrameNum)]);
%     drawnow()

%     

end  









