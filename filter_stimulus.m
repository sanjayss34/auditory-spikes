function [filtered] = filter_stimulus(binranges)
tic;
D = sta();
toc
disp('vijay singh')
tic;
load('rawStimuli.mat');
load('aud_stream.mat');
load('CohenNeurons.mat');
N = 16;
fs = Audstim.fs/1000;
num_trials = numel(CohenNeurons(1).trials);
spectro = spectrogram(rawStimCollector{1}, 128, 120, 128, Audstim.fs, 'yaxis');
filtered = zeros(N, num_trials, size(spectro, 2));
for n=1:N
    for i=1:num_trials
        % stimulus = abs(rawStimCollector{i});
        stimulus = sum(abs(spectrogram(rawStimCollector{i}, 128, 120, 128, Audstim.fs, 'yaxis')), 1);
%         numzeros = numel(binranges)-mod(numel(stimulus),numel(binranges));
%         stimulus = [stimulus zeros(1,numzeros)];
%         stimulus = transpose(reshape(stimulus, numel(binranges),[]));
%         stimulus = mean(stimulus);

        tempD = [0 D(n,:)];
        l = conv(tempD, stimulus);
        l = l(1:length(stimulus));
        
        index = round(binranges*fs);
        index(1) = 1;
        % filtered(n, i, :) = l(index);
        filtered(n, i, :) = l;
    end
    n
end
toc
disp('vijay singh 2')
end