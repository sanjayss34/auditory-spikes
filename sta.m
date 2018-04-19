function [D] = sta()
load('CohenNeurons.mat');
load('rawStimuli.mat');
trials_stimulus = rawStimCollector;
N = 16;
D = zeros(N, 5000);
load('aud_stream.mat');
fs = Audstim.fs/1000;
for n=1:N
for i=1:numel(CohenNeurons(n).trials)
    contributions = zeros(1, size(D, 2));
    spikes = CohenNeurons(n).trials(i).spikes;
    stimulus = trials_stimulus{i};
    count_spikes = 0;
    for t=1:numel(spikes)
        for j=1:min(size(D, 2), spikes(t)-1)
            if spikes(t) > 0 && spikes(t) <= 2000
                contributions(j) = contributions(j)+stimulus(cast((spikes(t))*fs, 'int32')-j);
                count_spikes = count_spikes+1;
            end
        end
    end
    if count_spikes > 0
        contributions = contributions/count_spikes;
        D(n,:) = D(n,:)+contributions;
    end
end
D(n,:) = D(n,:)/numel(CohenNeurons(n).trials);
end