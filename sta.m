function [D] = sta()
load('CohenNeurons.mat');
load('rawStimuli.mat');
load('aud_stream.mat');
trials_stimulus = rawStimCollector;
N = 16;
spectro = spectrogram(rawStimCollector{1}, 128, 120, 128, Audstim.fs);
filterlength = round(size(spectro, 2)/2);
D = zeros(N, filterlength);
for n=1:N
    for i=1:numel(CohenNeurons(n).trials)
        contributions = zeros(1, size(D, 2));
        spikes = CohenNeurons(n).trials(i).spikes;
        spikes = spikes(spikes>0);
        spikes = round(spikes(spikes<2000));
        count_spikes = numel(spikes);
        stimulus = abs(trials_stimulus{i});
        spectro = abs(spectrogram(stimulus, 128, 120, 128, Audstim.fs, 'yaxis'));
        stimulus = sum(spectro, 1);
%         index = round((1:2000)*fs);
%         stimulus = stimulus(index);

        for t=1:numel(spikes)
            % time = round(spikes(t)*Audstim.fs/1000);
            index = round(spikes(t)*size(D, 2)/2000);
            if index <= filterlength
                stimulus_before_spike = stimulus(1:index-1);
                contributions(1:index-1) = contributions(1:index-1) + fliplr(stimulus_before_spike);
            end
            if index > filterlength
                contributions = contributions + fliplr(stimulus(index-filterlength:index-1));
            end
        end
        if count_spikes > 0
            D(n,:) = D(n,:) + contributions/count_spikes;
        end
    end
    D(n,:) = D(n,:)/numel(CohenNeurons(n).trials);
end