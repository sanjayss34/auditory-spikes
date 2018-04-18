function [filtered] = filter_stimulus(binranges)
D = sta();
load('rawStimuli.mat');
load('aud_stream.mat');
load('CohenNeurons.mat');
N = 16;
fs = Audstim.fs/1000;
num_trials = numel(CohenNeurons(1).trials);
filtered = zeros(N, num_trials, numel(binranges));
for n=1:N
    for i=1:num_trials
        stimulus = rawStimCollector{i};
        l = zeros(numel(stimulus));
        for t=1:numel(l)
            for j=1:size(D, 2)
                if t-j*fs >= 1
                    l(t) = l(t)+stimulus(cast(t-j*fs, 'int8'))*D(n, j);
                end
            end
        end
        for t=1:numel(binranges)
            if binranges(t) > 0
                filtered(n, i, t) = l(cast(binranges(t)*fs, 'int8'));
            end
        end
    end
end