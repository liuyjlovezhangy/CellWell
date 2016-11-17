function [level, level_low, level_high, xvals, yvals] = adaptiveThresholdFinder(Istack3d, scale, tracking_params, levels_prev)

smooth_vols = 5;

%%%% Parameters %%%%
threshold_density = tracking_params.threshold_density; %density of threshold values to try (on log scale) when generating mean/median volume curves
% start = round(10*prctile(log(Istack3d(:)),90))/10; %assumes that objects take up less than 10% of total pixels
% if isnan(start) || isinf(start)
%     start = min(log(Istack3d(Istack3d>0)));
% end
start = min(log(Istack3d(Istack3d>0)));

levels = start:threshold_density:-0.1;

vols_mean = zeros(1,numel(levels));

for level_idx = 1:numel(levels)
    Ithresh = threshold3D(Istack3d,10^levels(level_idx));
    L=bwlabeln(Ithresh);
    STATS = regionprops(L,Ithresh,'Area');
    vols_mean(level_idx) = mean([STATS.Area]);
end
    


%Find optimal threshold level (first peak of mean volume curve after initial neg slope)
% sm = movingmean(vols_mean',3,[],[],1);
% sm2deriv = movingmean(diff(diff(sm)),3,[],[],1);

sm = smooth(vols_mean,smooth_vols);
sm2deriv = smooth(diff(diff(sm)),smooth_vols);

max2deriv = max(sm2deriv(1:end-2));
if max2deriv > 0, minidx2deriv = find(sm2deriv>max2deriv/3,1); %find beginning of first significant peak in second derivative (i.e. first end of a negative slope)
else minidx2deriv = find(sm2deriv>(max2deriv-sm2deriv(1))/3,1); end
pkidx2deriv = find(diff(sm2deriv(minidx2deriv:end))<0,1) + minidx2deriv - 1; %first significant peak in second derivative (not necessarily the max peak)

if pkidx2deriv
    minidx = find(diff(sm(pkidx2deriv+1:end))>0,1) + pkidx2deriv; %beginning of next positive slope in mean volume curve
else
    minidx = [];
end

if minidx
    [maxvol, maxidx] = max(sm(minidx:end)); %maximum value of mean volume after the neg slope
    maxidx = maxidx(end) + minidx - 1;
    minvol = min(sm(minidx:maxidx)); %minimum dip in mean volume between the beginning of positive slope and the peak
    lowidx = find(sm(minidx:end)>minvol+(maxvol-minvol)/4,1) + minidx - 1; %first point above one third of the max value
    if strcmp(tracking_params.peak_stringency,'high')
        pkidx = length(sm)-find(diff(sm(end:-1:lowidx))<0,1)+1; %last peak
    else
        pkidx = find(diff(sm(lowidx:end))<0,1) + lowidx - 1; %first peak above one third of the max value
    end
    if pkidx
%         disp('    peak')
        level_low = 10^(start + threshold_density*(lowidx-1)) * scale;
        level = 10^(start + threshold_density*(pkidx-1)) * scale;
        level_high = 2*level - level_low;
    else
%         disp('    no peak')
        if ~isempty(levels_prev)
            level_low = levels_prev.level_low;
            level = levels_prev.level;
            level_high = levels_prev.level_high;
        else
            level_low = 10^(start + threshold_density*(lowidx-1)) * scale;
            level = 2*level_low;
            level_high = 2*level - level_low;
        end
    end
else
    [shoulder, shoulderidx] = min(sm2deriv(pkidx2deriv+1:end));
    if shoulder <= 0 
%         disp('    shoulder peak')
        shoulderidx = shoulderidx + pkidx2deriv;
        level_low = 10^(start + threshold_density*(pkidx2deriv-1)) * scale;
        level = 10^(start + threshold_density*(shoulderidx-1)) * scale;
        level_high = 2*level - level_low;
    else
%         disp('    no peak')
        if ~isempty(levels_prev)
            level_low = levels_prev.level_low;
            level = levels_prev.level;
            level_high = levels_prev.level_high;
        else
            if pkidx2deriv
                level_low = 10^(start + threshold_density*(pkidx2deriv-1)) * scale;
            else
                level_low = 10^(start + threshold_density*(minidx2deriv-1)) * scale;
            end
            level = 2*level_low;
            level_high = 2*level - level_low;
        end
    end
end

if ~isempty(levels_prev) && (level/levels_prev.level>3 || levels_prev.level/level>3)   %%% should update this to find the nearest peak to the previous level
    level_low = levels_prev.level_low;
    level = levels_prev.level;
    level_high = levels_prev.level_high;
end

if ~isempty(levels_prev) && strcmp(tracking_params.threshold_smoothing,'on')
    level_low = (level_low + levels_prev.level_low)/2;
    level = (level + levels_prev.level)/2;
    level_high = (level_high + levels_prev.level_high)/2;
end

% disp(level)

xvals = 10.^[start:threshold_density:-0.1];
% yvals = movingmean(vols_mean',3,[],[],1);
yvals = smooth(vols_mean,smooth_vols);

%hold on

% figure(17631)
% clf
% plot(xvals,yvals,'b',[level level] / scale, ylim)
% 
% hold on, line([level_low level_low],[0 50]), line([level level],[0 50]);
    
end