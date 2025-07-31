%% Cooling and eOPN lever behaviour: additional statistics
% Additional statistics on perturbation experiments regarding normal lever
% hold and pull time in the same mouse across multiple mice/sessions
% This is to see if the M2->M1 routing is effecting only RT time to
% external stimulus or actually changing the learned motor movement. We
% don't expect a significant difference in lever movement. 

% We only need the lever structure which we can take from the pooled spike
% data in which we saved the quenching analysis. Fastest way to load the
% process data


% Compress traces and calc mean/se
trace = cell2mat(arrayfun(@(x) smoothdata(x.trace,'movmean',100), IntanBehaviour.cueHitTrace, 'UniformOutput', false))';
traceM = mean(trace,1);
traceSE = std(trace,[],1)/sqrt(size(trace,1));
figure(1)
plot(IntanBehaviour.cueHitTrace(1).time,traceM,'k'),hold on
plot(IntanBehaviour.cueHitTrace(1).time,traceM+traceSE,'k')
plot(IntanBehaviour.cueHitTrace(1).time,traceM-traceSE,'k')
box off,set(gca,'tickdir','out','fontsize',16)
%% Pull speed, Pull duration, pull reproducability
% integrated pull speed from cue to movement (RT)
figure(2),hold on
pSpeed = [];
temp = diff(trace,[],2);
for n = 1:size(trace,1)
    pSpeed{n} = (temp(n,1500:(1500+IntanBehaviour.reactionTime(n)*1000)));
    plot(smoothdata(pSpeed{n},'movmean',10),'color',[0 0 0 0.1]),hold on
end
figure(3)
customBoxplot(cellfun(@mean,pSpeed)');
