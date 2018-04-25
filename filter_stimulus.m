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
filtered = zeros(N, num_trials, numel(binranges));
for n=1:N
    for i=1:num_trials
        stimulus = abs(rawStimCollector{i});
%         numzeros = numel(binranges)-mod(numel(stimulus),numel(binranges));
%         stimulus = [stimulus zeros(1,numzeros)];
%         stimulus = transpose(reshape(stimulus, numel(binranges),[]));
%         stimulus = mean(stimulus);

        tempD = [0 D(n,:)];
        l = conv(tempD, stimulus);
        l = l(1:length(stimulus));
        
        index = round(binranges*fs);
        index(1) = 1;
        filtered(n, i, :) = l(index);
    end
    n
end
toc
disp('vijay singh 2')
end