function [D] = sta2d()
load('CohenNeurons.mat');
load('rawStimuli.mat');
trials_stimulus = rawStimCollector;
N = 16;
numfreqs = 65;
filterlength = 10000;
D = cell(N, 1);
load('aud_stream.mat');
fs = Audstim.fs/1000;
% For each neuron
for n=1:N
    % For each trial
    for i=1:numel(CohenNeurons(n).trials)
        % Initialize STA
        tempD = zeros(numfreqs, filterlength);
        % Find spikes > 0 and < trial end
        spikes = CohenNeurons(n).trials(i).spikes;
        spikes = spikes(spikes>0);
        spikes = round(spikes(spikes<2000));
        count_spikes = numel(spikes);
        % Determine stimulus
        stimulus = abs(trials_stimulus{i});

        % For each spike, add pre-spike stimulus to tempD
        for t=1:numel(spikes)
            time = round(spikes(t)*fs);
            if time <= filterlength
                stimulus_before_spike = stimulus(1:numfreqs, 1:time-1);
                tempD(1:numfreqs, 1:time-1) = tempD(1:numfreqs, 1:time-1) + fliplr(stimulus_before_spike);
            end
            if time > filterlength
                tempD = tempD + fliplr(stimulus(1:numfreqs, time-filterlength:time-1));
            end
        end
        % Linear filter is proportional to pre-spike stimulus
        if count_spikes > 0
            D{n} = D{n} + tempD/count_spikes;
        end
    end
    % Average linear filter over trials
    D{n} = D{n}/numel(CohenNeurons(n).trials);
end