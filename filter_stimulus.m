function [filtered] = filter_stimulus(binranges)
tic;
D = sta()
load('rawStimuli.mat');
load('aud_stream.mat');
load('CohenNeurons.mat');
N = 16;
fs = Audstim.fs/1000;
num_trials = numel(CohenNeurons(1).trials);
filtered = zeros(N, num_trials, numel(binranges));
disp('before loop');
toc
for n=1:N
    for i=1:num_trials
        stimulus = rawStimCollector{i};
        lengthD = length(D);
        size(D)
        tempD = [D(n,:) zeros(1, lengthD)];
        l = conv(tempD, stimulus);
        l = l(1:length(l)-lengthD);
        l
        i
        
        for t=1:numel(binranges)
            if binranges(t) > 0
                filtered(n, i, t) = l(cast(binranges(t)*fs, 'int8'));
            end
        end
    end
end