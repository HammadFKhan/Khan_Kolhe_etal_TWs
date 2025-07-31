%% Cooled Spiking data

files = dir(fullfile('D:\M1Cooling\SpikesGSP\','*.mat'));
M1DynamicsCooled = struct();
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    Spikes = makeSpikeGPFA(Spikes);
    Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
    for n = 1:IntanBehaviour.nCueHit%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
        Spikes.GPFA.HitMiss.dat(n).trialId = n;
    end
    Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
    for n = length(IntanBehaviour.MIHitTrace)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
        Spikes.GPFA.MIHitFA.dat(n).trialId = n;
    end
    %%%
    addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
    addpath(genpath('mat_results'));
    if exist('mat_results','dir'),rmdir('mat_results','s'),end
    [Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
    [Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
    [Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
    [Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
    [Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
    [Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
    close all
    [M1DynamicsCooled(fileNum).neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
    M1DynamicsCooled(fileNum).IntanBehaviour = IntanBehaviour;
    M1DynamicsCooled(fileNum).filename = files(fileNum).name;
end
%% Calculate Speed and trajectory dynamics for cooled and uncooled conditions for HIT trials
% We mainly only care about the hit conditions in this analysis since that
% is the metric we want to track over this experiment. That's not to say we
% should check the other conditions. We do; but it is a supplemental
% finding.

% Plot speed over time, highlighting different states
if ~exist('M1DynamicsCooled','var')
    load('D:\M1Cooling\M1DynamicsCooled.mat');
end

dynamics = M1DynamicsCooled;
dimension = 1;
speedTotBaseline = [];
speedTotCooled = [];
tempCutoff = -9; % Cuttoff of temperature cooling
rtbaseline = [];
rtcooled = [];
for n = 1:length(dynamics)
    if ~exist('dynamics(n).IntanBehaviour.hitTemp')
        dynamics(n).IntanBehaviour = grabTemp(dynamics(n).IntanBehaviour);
    end
    temperatureId = (dynamics(n).IntanBehaviour.hitTemp<tempCutoff);
    speed_data = dynamics(n).neuralDynamics.hit.speed;
    rtbaseline{n} = dynamics(n).neuralDynamics.rawreactionTime(~temperatureId);
    rtcooled{n} = dynamics(n).neuralDynamics.rawreactionTime(temperatureId);
    
    speedTotBaseline{n} = squeeze(speed_data.speed(dimension,2:end,~temperatureId));
    speedTotCooled{n} = squeeze(speed_data.speed(dimension,2:end,temperatureId));
end

speedTotBaseline = horzcat(speedTotBaseline{:});
speedTotCooled = speedTotCooled(~cellfun(@isempty, speedTotCooled));
speedTotCooled = horzcat(speedTotCooled{:});
rtbaseline = horzcat(rtbaseline{:});
rtcooled = horzcat(rtcooled{:});


colors = [166/255 14/255 90/255;40/255 153/255 196/255];
time = -1499:20:1500;
figure;
plot(time(2:end),mean(speedTotBaseline,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTotBaseline,2)+std(speedTotBaseline,[],2)/sqrt(size(speedTotBaseline,2)),'color',colors(1,:),'linewidth',2)
plot(time(2:end),mean(speedTotBaseline,2)-std(speedTotBaseline,[],2)/sqrt(size(speedTotBaseline,2)),'color',colors(1,:),'linewidth',2)
hold on;

plot(time(2:end),mean(speedTotCooled,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),mean(speedTotCooled,2)+std(speedTotCooled,[],2)/sqrt(size(speedTotCooled,2)),'color',colors(2,:),'linewidth',2)
plot(time(2:end),mean(speedTotCooled,2)-std(speedTotCooled,[],2)/sqrt(size(speedTotCooled,2)),'color',colors(2,:),'linewidth',2)
hold on;
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.005 0.1])
xline(0, '--r', 'Cue');
xline(mean(rtbaseline)*1000, '--g', 'MI');
xline(nanmean(rtcooled(rtcooled>.400))*1000, '--b', 'MI');
xlabel('Time (s)');
ylabel('Average Speed');
%% Do bootstrap
% Shuffle baseline and cooled and see how different the effect is to
% baseline
%   speedTotBaseline: [n_baseline × time] matrix
%   speedTotCooled: [n_cooled × time] matrix

% 1. Combine data
tot = [speedTotBaseline, speedTotCooled]';
n_baseline = size(speedTotBaseline, 2);
n_cooled = size(speedTotCooled, 2);

% 2. Run test
n_permutations = 1000;
[p_values, obs_diff, perm_diffs] = permutation_test(tot, n_baseline, n_cooled, n_permutations);

% 3. Interpret results
significant_cue_mov = mean(p_values(75:90));
disp(['Significant val ', num2str((significant_cue_mov))]);

%% Miss
dimension = 1;
speedTotBaseline = [];
speedTotCooled = [];
tempCutoff = -12; % Cuttoff of temperature cooling
for n = 1:length(dynamics)
    temperatureId = (dynamics(n).IntanBehaviour.missTemp<tempCutoff);
    speed_data = dynamics(n).neuralDynamics.miss.speed;
    
    speedTotBaseline{n} = squeeze(speed_data.speed(dimension,2:end,~temperatureId));
    speedTotCooled{n} = squeeze(speed_data.speed(dimension,2:end,temperatureId));
end
speedTotBaseline = horzcat(speedTotBaseline{:});
speedTotCooled = speedTotCooled(~cellfun(@isempty, speedTotCooled));
speedTotCooled = horzcat(speedTotCooled{:});

colors = [166/255 14/255 90/255;40/255 153/255 196/255];
time = -1499:20:1500;
figure;
plot(time(2:end),nanmean(speedTotBaseline,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),nanmean(speedTotBaseline,2)+nanstd(speedTotBaseline,[],2)/sqrt(size(speedTotBaseline,2)),'color',colors(1,:),'linewidth',2)
plot(time(2:end),nanmean(speedTotBaseline,2)-nanstd(speedTotBaseline,[],2)/sqrt(size(speedTotBaseline,2)),'color',colors(1,:),'linewidth',2)
hold on;

plot(time(2:end),nanmean(speedTotCooled,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),nanmean(speedTotCooled,2)+nanstd(speedTotCooled,[],2)/sqrt(size(speedTotCooled,2)),'color',colors(2,:),'linewidth',2)
plot(time(2:end),nanmean(speedTotCooled,2)-nanstd(speedTotCooled,[],2)/sqrt(size(speedTotCooled,2)),'color',colors(2,:),'linewidth',2)
hold on;
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.005 0.055])
xlabel('Time (s)');
ylabel('Average Speed');
%% FA
dimension = 1;
speedTotBaseline = [];
speedTotCooled = [];
tempCutoff = -13; % Cuttoff of temperature cooling
for n = 1:length(dynamics)
    temperatureId = (dynamics(n).IntanBehaviour.FATemp<tempCutoff);
    speed_data = dynamics(n).neuralDynamics.MIFA.speed;
    speedTotBaseline{n} = squeeze(speed_data.speed(dimension,2:end,~temperatureId));
    speedTotCooled{n} = squeeze(speed_data.speed(dimension,2:end,temperatureId));
end
speedTotBaseline = horzcat(speedTotBaseline{:});
speedTotCooled = speedTotCooled(~cellfun(@isempty, speedTotCooled));
speedTotCooled = horzcat(speedTotCooled{:});

colors = [166/255 14/255 90/255;40/255 153/255 196/255];
time = -1499:20:1500;
figure;
plot(time(2:end),nanmean(speedTotBaseline,2),'color',colors(1,:),'linewidth',2),hold on
plot(time(2:end),nanmean(speedTotBaseline,2)+nanstd(speedTotBaseline,[],2)/sqrt(size(speedTotBaseline,2)),'color',colors(1,:),'linewidth',2)
plot(time(2:end),nanmean(speedTotBaseline,2)-nanstd(speedTotBaseline,[],2)/sqrt(size(speedTotBaseline,2)),'color',colors(1,:),'linewidth',2)
hold on;

plot(time(2:end),nanmean(speedTotCooled,2),'color',colors(2,:),'linewidth',2),hold on
plot(time(2:end),nanmean(speedTotCooled,2)+nanstd(speedTotCooled,[],2)/sqrt(size(speedTotCooled,2)),'color',colors(2,:),'linewidth',2)
plot(time(2:end),nanmean(speedTotCooled,2)-nanstd(speedTotCooled,[],2)/sqrt(size(speedTotCooled,2)),'color',colors(2,:),'linewidth',2)
hold on;
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([-500 1500]),ylim([0.005 0.055])
xlabel('Time (s)');
ylabel('Average Speed');
%% Check speed cooling vs baseline over dimensions
totalSpeedDimensionBaseline = [];
totalSpeedDimensionCooled = [];

tempCutoff = -12; % Cuttoff of temperature cooling
speed_data = [];
for dimension = 1:6
    speedHitCueMovement = [];
    speedMissCueMovement = [];
    speedMIhitCueMovement = [];
    speedFACueMovement = [];
    for n = 1:length(dynamics)
        temperatureId = (dynamics(n).IntanBehaviour.hitTemp<tempCutoff);
        speed_data = dynamics(n).neuralDynamics.hit.speed;
        speedHitCueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,~temperatureId), [2]));
        
        speed_data = dynamics(n).neuralDynamics.MIhit.speed;
        speedMIhitCueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,~temperatureId), [2]));
        
        temperatureId = (dynamics(n).IntanBehaviour.missTemp<tempCutoff);
        
        speed_data = dynamics(n).neuralDynamics.miss.speed;
        speedMissCueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,~temperatureId), [2]));
        
        temperatureId = (dynamics(n).IntanBehaviour.FATemp<tempCutoff);
        
        speed_data = dynamics(n).neuralDynamics.MIFA.speed;
        speedFACueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,~temperatureId), [2]));
        
    end
    speedHitCueMovement = vertcat(speedHitCueMovement{:});
    speedMissCueMovement = vertcat(speedMissCueMovement{:});
    speedMIhitCueMovement = vertcat(speedMIhitCueMovement{:});
    speedFACueMovement = vertcat(speedFACueMovement{:});
    
    temp = nan(max([length(speedHitCueMovement),length(speedMissCueMovement),length(speedMIhitCueMovement),length(speedFACueMovement)]),4);
    temp(1:length(speedHitCueMovement),1) = speedHitCueMovement;
    temp(1:length(speedMissCueMovement),2) = speedMissCueMovement;
    temp(1:length(speedMIhitCueMovement),3) = speedMIhitCueMovement;
    temp(1:length(speedFACueMovement),4) = speedFACueMovement;
    
    totalSpeedDimensionBaseline{dimension} = temp;
    
    speedHitCueMovement = [];
    speedMissCueMovement = [];
    speedMIhitCueMovement = [];
    speedFACueMovement = [];
    for n = 1:length(dynamics)
        temperatureId = (dynamics(n).IntanBehaviour.hitTemp<tempCutoff);
        speed_data = dynamics(n).neuralDynamics.hit.speed;
        speedHitCueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,temperatureId), [2]));
        
        speed_data = dynamics(n).neuralDynamics.MIhit.speed;
        speedMIhitCueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,temperatureId), [2]));
        
        temperatureId = (dynamics(n).IntanBehaviour.missTemp<tempCutoff);
        
        speed_data = dynamics(n).neuralDynamics.miss.speed;
        speedMissCueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,temperatureId), [2]));
        
        temperatureId = (dynamics(n).IntanBehaviour.FATemp<tempCutoff);
        
        speed_data = dynamics(n).neuralDynamics.MIFA.speed;
        speedFACueMovement{n} = squeeze(mean(speed_data.preMovement(dimension,:,temperatureId), [2]));
        
    end
    speedHitCueMovement = vertcat(speedHitCueMovement{:});
    speedMissCueMovement = vertcat(speedMissCueMovement{:});
    speedMIhitCueMovement = vertcat(speedMIhitCueMovement{:});
    speedFACueMovement = vertcat(speedFACueMovement{:});
    
    temp = nan(max([length(speedHitCueMovement),length(speedMissCueMovement),length(speedMIhitCueMovement),length(speedFACueMovement)]),4);
    temp(1:length(speedHitCueMovement),1) = speedHitCueMovement;
    temp(1:length(speedMissCueMovement),2) = speedMissCueMovement;
    temp(1:length(speedMIhitCueMovement),3) = speedMIhitCueMovement;
    temp(1:length(speedFACueMovement),4) = speedFACueMovement;
    totalSpeedDimensionCooled{dimension} = temp;
end
%% Plot average speed for each state as a function of cooled and not cooled
figure;hold on
%colors = [0 0.4470 0.7410;0.75 0.75 0.75;0 0.4470 0.7410;190/255 30/255 45/255];
colors = [166/255 14/255 90/255;40/255 153/255 196/255];
% Change X indexing for task variables
% Baseline
dat = cellfun(@(x) horzcat(x(:,1)),totalSpeedDimensionBaseline,'UniformOutput',false);
dat1 = horzcat(dat{:});
errorbar(1:6,nanmean(dat1,1),nanstd(dat1,[],1)./sqrt(size(dat1,1)/50),'color',colors(1,:),'linewidth',2),hold on
% Cooling
dat = cellfun(@(x) horzcat(x(:,1)),totalSpeedDimensionCooled,'UniformOutput',false);
dat2 = horzcat(dat{:});
errorbar(1:6,nanmean(dat2,1),nanstd(dat2,[],1)./sqrt(size(dat2,1)/50),'color',colors(2,:),'linewidth',2),hold on
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([0.5 6.5]),ylim([0 0.1])
ylabel('Trajectory Speed')
xlabel('Dimensions')
%% Do stats
% Extract baseline and cooling data
baselineData = cellfun(@(x) horzcat(x(:,1)),totalSpeedDimensionBaseline,'UniformOutput',false);
baselineData = horzcat(baselineData{:});

coolingData = cellfun(@(x) horzcat(x(:,1)),totalSpeedDimensionCooled,'UniformOutput',false);
coolingData = horzcat(coolingData{:});

% Combine data into a single matrix
coolingDataFix = nan(max([size(baselineData,1),size(coolingData,1)]),6);
coolingDataFix(1:size(coolingData,1),:) = coolingData;
all_data = [baselineData; coolingDataFix];


% Create grouping variables
num_samples = size(baselineData, 1); % Number of rows in baseline
dimensions = repmat(1:6, num_samples * 2, 1); % Dimension grouping (1-6)
conditions = [repmat({'Baseline'}, num_samples, 6); repmat({'Cooling'}, num_samples, 6)]; % Condition grouping

% Reshape data into column vector for ANOVA
all_data_vector = all_data(:);
dimensions_vector = dimensions(:);
conditions_vector = conditions(:);

% Perform two-way ANOVA
[p, tbl, stats] = anovan(all_data_vector, {dimensions_vector, conditions_vector}, ...
    'model', 'interaction', 'varnames', {'Dimension', 'Condition'});

% Display results
disp('ANOVA Table:');
disp(tbl);

% Perform post-hoc analysis if necessary
disp('Post-hoc comparisons:');
multcompare(stats, 'Dimension',1) % compare over neural dimensions
multcompare(stats, 'Dimension',2) % Compare over baseline and cooled

%% Trajectory Deviation
tempCutoff = -12; % Cuttoff of temperature cooling
hitbaseline = [];hitcooled = [];
missbaseline = [];misscooled = [];
FAbaseline = [];FAcooled = [];
dimension = 1;
for n = 1:length(M1DynamicsCooled)
    temperatureId = (dynamics(n).IntanBehaviour.hitTemp<tempCutoff);
    hitbaseline{n} = dynamics(n).neuralDynamics.hit.s(dimension,~temperatureId);
    hitcooled{n} = dynamics(n).neuralDynamics.hit.s(dimension,temperatureId);
    
    temperatureId = (dynamics(n).IntanBehaviour.missTemp<tempCutoff);
    missbaseline{n} = dynamics(n).neuralDynamics.miss.s(dimension,~temperatureId);
    misscooled{n} = dynamics(n).neuralDynamics.miss.s(dimension,temperatureId);
    
    temperatureId = (dynamics(n).IntanBehaviour.FATemp<tempCutoff);
    FAbaseline{n} = dynamics(n).neuralDynamics.MIFA.s(dimension,~temperatureId);
    FAcooled{n} = dynamics(n).neuralDynamics.MIFA.s(dimension,temperatureId);
end
hitbaseline = horzcat(hitbaseline{:});
hitcooled = horzcat(hitcooled{:});
missbaseline = horzcat(missbaseline{:});
misscooled = horzcat(misscooled{:});
FAbaseline = horzcat(FAbaseline{:});
FAcooled = horzcat(FAcooled{:});



temp = nan(max([size(hitbaseline,2),size(hitcooled,2),...
    size(missbaseline,2),size(misscooled,2),...
    size(FAbaseline,2),size(FAcooled,2)]),6);

temp(1:size(hitbaseline,2),1) = hitbaseline;
temp(1:size(hitcooled,2),2) = hitcooled;
temp(1:size(missbaseline,2),3) = missbaseline;
temp(1:size(misscooled,2),4) = misscooled;
temp(1:size(FAbaseline,2),5) = FAbaseline;
temp(1:size(FAcooled,2),6) = FAcooled;
%%
figure,customBarplot(temp);
box off,set(gca,'tickdir','out','fontsize',14),axis square
ylabel('Trajectory Deviation')
ylim([0 20])
%%
[p,t,stats] = anova1(temp)
c = multcompare(stats)


%% Supplementary analysis
figure;
for n = 1:4
    dat = cellfun(@(x) horzcat(x(:,n)),totalSpeedDimension,'UniformOutput',false);
    dat = horzcat(dat{:});
    errorbar(1:6,nanmean(dat,1),nanstd(dat,[],1)./sqrt(size(dat,1)),'k.-','linewidth',2,'color',colors(n,:)),hold on
end
box off,set(gca,'tickdir','out','fontsize',14),axis square,xlim([0.5 6.5]),ylim([0 0.06])



[p, tbl, stats] = anova1(dat);

fprintf('ANOVA p-value: %f\n', p);

if p < 0.05
    c = multcompare(stats, 'Display', 'off');
    disp('Multiple Comparisons:');
    disp(c);
end


%%
figure;
customBarplot(temp,'scatter','off');
ylabel('Neural Trajectory Speed');
%title(['Average Trajectory Speed by Behavioral State for Dimension ' num2str(dimension)]);
box off,set(gca,'tickdir','out','fontsize',14),axis square,ylim([0 0.07])
%ylim([0 0.08])
% Perform statistical comparison (ANOVA)
[p, tbl, stats] = anova1(temp);

fprintf('ANOVA p-value: %f\n', p);

if p < 0.05
    c = multcompare(stats, 'Display', 'off');
    disp('Multiple Comparisons:');
    disp(c);
end

%% functions

function IntanBehaviour = grabTemp(IntanBehaviour)
for n = 1:IntanBehaviour.nCueHit
    IntanBehaviour.hitTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueHitTrace(n).LFPIndex(1));
end
IntanBehaviour.hitTemp = IntanBehaviour.hitTemp-IntanBehaviour.temperature(100);
for n = 1:IntanBehaviour.nCueMiss
    IntanBehaviour.missTemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.cueMissTrace(n).LFPIndex(1));
end
IntanBehaviour.missTemp = IntanBehaviour.missTemp-IntanBehaviour.temperature(100);
for n = 1:length(IntanBehaviour.missTrace)
    IntanBehaviour.FATemp(n,1) = IntanBehaviour.temperature(IntanBehaviour.missTrace(n).LFPIndex(1));
end
IntanBehaviour.FATemp = IntanBehaviour.FATemp-IntanBehaviour.temperature(100);
end

function [p_value, observed_difference, perm_diffs] = permutation_test(data, n_baseline, n_cooled, n_permutations)
    % Inputs:
    %   data: Combined matrix [n_baseline + n_cooled × time_points]
    %   n_baseline: Number of baseline trials
    %   n_cooled: Number of cooled trials
    %   n_permutations: Number of permutations (default 1000)
    
    if nargin < 4
        n_permutations = 1000;
    end
    
    % 1. Compute observed difference
    baseline_mean = mean(data(1:n_baseline, :), 1);
    cooled_mean = mean(data(n_baseline+1:end, :), 1);
    observed_difference = baseline_mean - cooled_mean;
    
    % 2. Initialize permutation results
    perm_diffs = zeros(n_permutations, size(data, 2));
    
    % 3. Permutation loop
    for i = 1:n_permutations
        % Shuffle rows without replacement
        shuffled_idx = randperm(size(data, 1));
        shuffled_data = data(shuffled_idx, :);
        
        % Split into pseudo-groups
        perm_baseline = shuffled_data(1:n_baseline, :);
        perm_cooled = shuffled_data(n_baseline+1:end, :);
        
        % Compute permutation difference
        perm_diffs(i, :) = mean(perm_baseline, 1) - mean(perm_cooled, 1);
    end
    
    % 4. Calculate p-value (two-tailed test)
    abs_observed = (observed_difference);
    abs_permutations = (perm_diffs);
    
    % Count where permuted difference >= observed difference
    extreme_count = sum(abs_permutations >= abs_observed, 1);
    p_value = extreme_count / n_permutations;
end
