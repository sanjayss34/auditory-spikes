function [D] = sta()
load('CohenNeurons.mat');
load('rawStimuli.mat');
trials_stimulus = rawStimCollector;
N = 16;
filterlength = 10000;
D = zeros(N, filterlength);
load('aud_stream.mat');
fs = Audstim.fs/1000;
for n=1:N
    for i=1:numel(CohenNeurons(n).trials)
        contributions = zeros(1, size(D, 2));
        spikes = CohenNeurons(n).trials(i).spikes;
        spikes = spikes(spikes>0);
        spikes = round(spikes(spikes<2000));
        count_spikes = numel(spikes);
        stimulus = abs(trials_stimulus{i});
%         index = round((1:2000)*fs);
%         stimulus = stimulus(index);

        for t=1:numel(spikes)
            time = round(spikes(t)*fs);
            if time <= filterlength
                stimulus_before_spike = stimulus(1:time-1);
                contributions(1:time-1) = contributions(1:time-1) + fliplr(stimulus_before_spike);
            end
            if time > filterlength
                contributions = contributions + fliplr(stimulus(time-filterlength:time-1));
            end
        end
        if count_spikes > 0
            D(n,:) = D(n,:) + contributions/count_spikes;
        end
    end
    D(n,:) = D(n,:)/numel(CohenNeurons(n).trials);
end