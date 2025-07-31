%% Load mat file 
%load('\\10.165.57.13\Sutter_backup\Hammad\Ephys\LeverTask\Data_for_Figures\M1GridBaseline\WavePropCombined.mat')

%% Make parameters variable 
parameters.Fs = 1000;
parameters.ts = 1/parameters.Fs;
parameters.windowBeforePull = 1.5; % in seconds
parameters.windowAfterPull = 1.5; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds 
parameters.windowAfterCue = 1.5; % in seconds 
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 1.5; % in seconds 

%% Generate evaluation points from wavePresent

% wavesHit
for i = 1:size(wavesHit,2)
    wavesHit(i).evaluationPoints = find(diff(wavesHit(i).wavePresent)==1)+1;
end
% wavesMiss
for i = 1:size(wavesMiss,2)
    wavesMiss(i).evaluationPoints = find(diff(wavesMiss(i).wavePresent)==1)+1;
end
% wavesHitReward
for i = 1:size(wavesHitReward,2)
    wavesHitReward(i).evaluationPoints = find(diff(wavesHitReward(i).wavePresent)==1)+1;
end
% wavesFA
for i = 1:size(wavesFA,2)
    wavesFA(i).evaluationPoints = find(diff(wavesFA(i).wavePresent)==1)+1;
end
% wavesMIHit
for i = 1:size(wavesMIHit,2)
    wavesMIHit(i).evaluationPoints = find(diff(wavesMIHit(i).wavePresent)==1)+1;
end
% wavesFA
for i = 1:size(wavesFA,2)
    wavesMIFA(i).evaluationPoints = find(diff(wavesMIFA(i).wavePresent)==1)+1;
end

%% Getting file start and end points
varNames = {'Hit','Miss', 'HitReward', 'FA', 'MIHit', 'MIFA'};
for k = 1:numel(varNames)
    var = eval(['waves' varNames{k}]); % get the struct, e.g. wavesMiss
    fID = [find(diff(horzcat(var.fileID)) >= 1) size(var,2)];
    fID(2,:) = [1 fID(1:end-1)+1];
    fID = flip(fID, 1);
    fileIDArray.(varNames{k}) = fID;
end

%% Wave speed z-scored to baseline
avgRT = mean(vertcat(wavesHit(1:fileIDArray.Hit(2,end)).RT));
nPoints = 60; 
interval = (parameters.Fs*(parameters.windowAfterCue+parameters.windowBeforeCue))/nPoints;
t = ((interval:interval:interval*nPoints)/parameters.Fs)-1.5;
for i=1:size(fileIDArray.Hit,2)
    [waveSpeedHit(i,:),waveSpeedzHit(i,:)] = getWaveSpeedNorm(wavesHit(fileIDArray.Hit(1,i):fileIDArray.Hit(2,i)),parameters,nPoints);
    [waveSpeedMiss(i,:),waveSpeedzMiss(i,:)] = getWaveSpeedNorm(wavesMiss(fileIDArray.Miss(1,i):fileIDArray.Miss(2,i)),parameters,nPoints);
    [waveSpeedMIHit(i,:),waveSpeedzMIHit(i,:)] = getWaveSpeedNorm(wavesMIHit(fileIDArray.MIHit(1,i):fileIDArray.MIHit(2,i)),parameters,nPoints);
    [waveSpeedMIFA(i,:),waveSpeedzMIFA(i,:)] = getWaveSpeedNorm(wavesMIFA(fileIDArray.MIFA(1,i):fileIDArray.MIFA(2,i)),parameters,nPoints);
end

waveSpeedzHit      = arrayfun(@(col) horzcat(waveSpeedzHit{:,col}),      1:size(waveSpeedzHit,2),      'UniformOutput', false);
waveSpeedzMiss     = arrayfun(@(col) horzcat(waveSpeedzMiss{:,col}),     1:size(waveSpeedzMiss,2),     'UniformOutput', false);
% waveSpeedzHitReward= arrayfun(@(col) horzcat(waveSpeedzHitReward{:,col}),1:size(waveSpeedzHitReward,2),'UniformOutput', false);
% waveSpeedzFA       = arrayfun(@(col) horzcat(waveSpeedzFA{:,col}),       1:size(waveSpeedzFA,2),       'UniformOutput', false);
waveSpeedzMIFA     = arrayfun(@(col) horzcat(waveSpeedzMIFA{:,col}),     1:size(waveSpeedzMIFA,2),     'UniformOutput', false);
waveSpeedzMIHit    = arrayfun(@(col) horzcat(waveSpeedzMIHit{:,col}),    1:size(waveSpeedzMIHit,2),    'UniformOutput', false);

smoothWin = 3; % window size

figure();
subplot(2,1,1)
y = cellfun(@mean, waveSpeedzHit);
err = cellfun(@(x) std(x) / sqrt(numel(x)), waveSpeedzHit);
h1 = plot(t,smoothdata(y,'movmean', smoothWin),'Color', [0.8500 0.3250 0.0980],'LineWidth',2); hold on;
plot(t,smoothdata(y-err,'movmean', smoothWin),'Color', [0.8500 0.3250 0.0980],'LineWidth',1); hold on;
plot(t,smoothdata(y+err,'movmean', smoothWin),'Color', [0.8500 0.3250 0.0980],'LineWidth',1); hold on;
y = cellfun(@mean, waveSpeedzMiss);
err = cellfun(@(x) std(x) / sqrt(numel(x)), waveSpeedzMiss);
h2 = plot(t,smoothdata(y,'movmean', smoothWin),'Color', [0 0 0],'LineWidth',2); hold on;
plot(t,smoothdata(y-err,'movmean', smoothWin),'Color', [0.7 0.7 0.7],'LineWidth',1); hold on;
plot(t,smoothdata(y+err,'movmean', smoothWin),'Color', [0.7 0.7 0.7],'LineWidth',1); hold on;
xline(0,'--r','Cue');xlabel('Time (ms)'); ylabel('Change in wave speed');
xline(avgRT,'--r','Avg RT');
legend([h1 h2],'Hits','Misses','Location','best'); %ylim([5 15]);
title('Wave Speed - M1')
xlim([-0.5 1.5]);box off;set(gca,'TickDir','out','fontsize',14');ylim([-0.25 0.15]);

subplot(2,1,2)
y = cellfun(@mean, waveSpeedzMIHit);
err = cellfun(@(x) std(x) / sqrt(numel(x)), waveSpeedzMIHit);
h1 = plot(t,smoothdata(y,'movmean', smoothWin),'Color', [0.8500 0.3250 0.0980],'LineWidth',2); hold on;
plot(t,smoothdata(y-err,'movmean', smoothWin),'Color', [0.8500 0.3250 0.0980],'LineWidth',1); hold on;
plot(t,smoothdata(y+err,'movmean', smoothWin),'Color', [0.8500 0.3250 0.0980],'LineWidth',1); hold on;
y = cellfun(@mean, waveSpeedzMIFA);
err = cellfun(@(x) std(x) / sqrt(numel(x)), waveSpeedzMIFA);
h2 = plot(t,smoothdata(y,'movmean', smoothWin),'Color', [0 0 0],'LineWidth',2); hold on;
plot(t,smoothdata(y-err,'movmean', smoothWin),'Color', [0.7 0.7 0.7],'LineWidth',1); hold on;
plot(t,smoothdata(y+err,'movmean', smoothWin),'Color', [0.7 0.7 0.7],'LineWidth',1); hold on;
xline(0,'--r','MI');xlabel('Time (ms)'); ylabel('Change in wave speed');
legend([h1 h2],'MIHits','MIFA','Location','best'); %ylim([5 15]);
title('Wave Speed - M1')
xlim([-0.5 1.5]);box off;set(gca,'TickDir','out','fontsize',14');ylim([-0.25 0.15]);
beautifyFigure;

%% Plotting PGD 

trialTime = -1.5:1/parameters.Fs:1.5;
avgRT = mean(vertcat(wavesHit(1:fileIDArray.Hit(2,end-2)).RT));
PGDHits = vertcat(wavesHit.PGD);
PGDMiss = vertcat(wavesMiss.PGD);
PGDMIHit = vertcat(wavesMIHit.PGD);
PGDMIFA = vertcat(wavesMIFA.PGD);

for i=1:6%size(fileIDArray.Hit,2)
    mu = mean(PGDHits(fileIDArray.Hit(1,i):fileIDArray.Hit(2,i),1:1500),'all');
    std1 = std(PGDHits(fileIDArray.Hit(1,i):fileIDArray.Hit(2,i),1:1500),0,[1 2]);
    PGD.Hit{i} = (PGDHits(fileIDArray.Hit(1,i):fileIDArray.Hit(2,i),:)-mu)/std1;
    meanPGD.Hit(i,:) = mean((PGDHits(fileIDArray.Hit(1,i):fileIDArray.Hit(2,i),:)-mu)/std1,1);
    
    mu = mean(PGDMiss(fileIDArray.Miss(1,i):fileIDArray.Miss(2,i),1:1500),'all');
    std1 = std(PGDMiss(fileIDArray.Miss(1,i):fileIDArray.Miss(2,i),1:1500),0,[1 2]);
    PGD.Miss{i} = (PGDMiss(fileIDArray.Miss(1,i):fileIDArray.Miss(2,i),:)-mu)/std1;
    meanPGD.Miss(i,:) = mean((PGDMiss(fileIDArray.Miss(1,i):fileIDArray.Miss(2,i),:)-mu)/std1,1);

    mu = mean(PGDMIHit(fileIDArray.MIHit(1,i):fileIDArray.MIHit(2,i),1:800),'all');
    std1 = std(PGDMIHit(fileIDArray.MIHit(1,i):fileIDArray.MIHit(2,i),1:800),0,[1 2]);
    PGD.MIHit{i} = (PGDMIHit(fileIDArray.MIHit(1,i):fileIDArray.MIHit(2,i),:)-mu)/std1;
    meanPGD.MIHit(i,:) = mean((PGDMIHit(fileIDArray.MIHit(1,i):fileIDArray.MIHit(2,i),:)-mu)/std1,1);

    mu = mean(PGDMIFA(fileIDArray.MIFA(1,i):fileIDArray.MIFA(2,i),1:800),'all');
    std1 = std(PGDMIFA(fileIDArray.MIFA(1,i):fileIDArray.MIFA(2,i),1:800),0,[1 2]);
    PGD.MIFA{i} = (PGDMIFA(fileIDArray.MIFA(1,i):fileIDArray.MIFA(2,i),:)-mu)/std1;
    meanPGD.MIFA(i,:) = mean((PGDMIFA(fileIDArray.MIFA(1,i):fileIDArray.MIFA(2,i),:)-mu)/std1,1);
end

smoothWin = 50; % window size
figure();
subplot(2,1,1);
shadedErrorBar(trialTime, cell2mat(PGD.Hit'), ...
    {@(x) smoothdata(mean(x), 'movmean', smoothWin), ...
     @(x) smoothdata(std(x)/sqrt(size(x,1)), 'movmean', smoothWin)}, ...
     'lineprops', '-b');
hold on;
shadedErrorBar(trialTime, cell2mat(PGD.Miss'), ...
    {@(x) smoothdata(mean(x), 'movmean', smoothWin), ...
     @(x) smoothdata(std(x)/sqrt(size(x,1)), 'movmean', smoothWin)}, ...
     'lineprops', '-k');
title('Baseline-normalized PGD');ylabel('Change in PGD wrt baseline'); xlabel('Time from Cue(s)');
xline(0,'--r','Cue','LabelVerticalAlignment','top');xline(avgRT,'--r','Avg RT','LabelVerticalAlignment','top');box off;  legend('Hits','Misses');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');ylim([-0.2 0.3]);

subplot(2,1,2);
shadedErrorBar(trialTime, cell2mat(PGD.MIHit'), ...
    {@(x) smoothdata(mean(x), 'movmean', smoothWin), ...
     @(x) smoothdata(std(x)/sqrt(size(x,1)), 'movmean', smoothWin)}, ...
     'lineprops', '-b');
hold on;
shadedErrorBar(trialTime, cell2mat(PGD.MIFA'), ...
    {@(x) smoothdata(mean(x), 'movmean', smoothWin), ...
     @(x) smoothdata(std(x)/sqrt(size(x,1)), 'movmean', smoothWin)}, ...
     'lineprops', '-k');
title('Baseline-normalized PGD');ylabel('Change in PGD wrt baseline'); xlabel('Time from MI(s)');
xline(0,'--r','MI','LabelVerticalAlignment','top');box off;  legend('MI Hit','MI FA');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');ylim([-0.2 0.3]);

%% Wave direction

avgRT = mean(vertcat(wavesHit(1:fileIDArray.Hit(2,end)).RT));
nPoints = 15; 
interval = (parameters.Fs*(parameters.windowAfterCue+parameters.windowBeforeCue))/nPoints;
t = ((interval:interval:interval*nPoints)/parameters.Fs)-1.5;

clear waveDirectionality waveDirection 
for i=1:size(fileIDArray.Hit,2)
    [waveDirectionality.Hit(i,:),waveDirection.Hit(i,:)] = getNetWaveDirectionality(wavesHit(fileIDArray.Hit(1,i):fileIDArray.Hit(2,i)),nPoints,parameters);
end

figure();
plot(t,mean(waveDirectionality.Hit(:,:),1));


figure();
subplot(2,1,1);
plot(t,waveDirectionality.Hit(:,:));
subplot(2,1,2);
plot(t,waveDirection.Hit(:,:));

%% Plotting waves

for i=1:size(wavesHit,2)
    [wavesHit(i).wave_id, wavesHit(i).found_flag] = findWaveByCriteria(wavesHit(i), parameters.Fs, 10, 1, 0.9, 40);
end

animateWaves(66,IntanBehaviour.cueHitTrace,Waves.wavesHit,0,17)

animateWaves(66,IntanBehaviour.cueHitTrace,Waves.wavesHit,0,1)


animateWaves(60,IntanBehaviour.cueHitTrace,Waves.wavesHit,0,1)
animateWaves(60,IntanBehaviour.cueHitTrace,Waves.wavesHit,0,8)


animateWavesPhase(66,IntanBehaviour.cueHitTrace,Waves.wavesHit,0,1)

animateWavesAmplitude(66,IntanBehaviour.cueHitTrace,Waves.wavesHit,0,1)
animateWavesAmplitude(60,IntanBehaviour.cueHitTrace,Waves.wavesHit,0,8)
animateWavesAmplitude(66,IntanBehaviour.cueHitTrace,Waves.wavesHit,0,17)

trial = 60;
waveno = 8;
LFPamplitude = reshape(IntanBehaviour.cueHitTrace(trial).xf(:,:,Waves.wavesHit(trial).waveTime(waveno,1)-30:Waves.wavesHit(trial).waveTime(waveno,2)+30),32,[]);
% LFPamplitude = reshape(IntanBehaviour.cueHitTrace(trial).xf(:,:,:),32,[]);
figure();
plot(LFPamplitude');

figure();
stack_plot(LFPamplitude,0,1,1000);



animateWavesAmplitude_PNG(66,IntanBehaviour.cueHitTrace,Waves.wavesHit,0,17)
%% LFP with GP

%% look at file
LFP_with_GP

%% Phase gradient map and rho

trial = 66;
waveno = 17;
timepoint = round((Waves.wavesHit(trial).waveTime(waveno,1)+Waves.wavesHit(trial).waveTime(waveno,2))/2);
p = IntanBehaviour.cueHitTrace(trial).xgp; 
[pm,pd,dx,dy] = phase_gradient_complex_multiplication(p, parameters.xspacing, parameters.yspacing );
sourcepoint = find_source_points(timepoint,X,Y,dx, dy );

figure();
X = parameters.X;Y = parameters.Y;
imagesc(angle(p(:,:,timepoint)));hold on; axis off;
% map = colorcet( 'C2' ); colormap( map ); %caxis([-pi pi]);
colormap(viridis);
cb = colorbar;
% h2=quiver(X,Y,dx(:,:,timepoint),dy(:,:,timepoint));
h2=quiver(X,Y,cos(pd(:,:,timepoint)),sin(pd(:,:,timepoint)));
h2.Color = 'k';
h2.LineWidth = 1.5;
h2.ShowArrowHead = 'on';
h2.MaxHeadSize = 2;
h2.AutoScale = 1;
h2.AutoScaleFactor = 0.4;
plot(sourcepoint(1), sourcepoint(2), '.', 'markersize', 35, 'color', [.7 .7 .7]);


figure();
imagesc(angle(p(:,:,timepoint))-angle(p(sourcepoint(1), sourcepoint(2),timepoint)));hold on; axis off;
colormap(viridis);
cb = colorbar;

del_phase = reshape(angle(p(:,:,timepoint))-angle(p(sourcepoint(1), sourcepoint(2),timepoint)),32,[]);
[cc,pv,D] = phase_correlation_distance( angle(p(:,:,timepoint)), sourcepoint, parameters.xspacing, parameters.yspacing);

figure();
scatter(D,del_phase);hold on;

coeffs = polyfit(D, del_phase, 1);    % [slope, intercept]
y_fit = polyval(coeffs, D);

% Plot linear fit
plot(D, y_fit, 'r-', 'LineWidth', 2);