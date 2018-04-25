function [transform] = stimulus_binning(numbins)    

% Transform is the matrix of the amplitude, extracted from the raw
% stimulus.

    load('rawStimuli.mat')
    totalmax = 0;
    numtrials = numel(rawStimCollector);
    transform = zeros(numtrials, numbins);
    for t = 1%:numtrials
        stim = abs(rawStimCollector{t});
        box = ones(1,1000);
        convolution = conv(box,stim);
        convolution = convolution(1:numel(stim));
        stim = convolution*(max(stim)/max(convolution));
%         plot(stim);
    %     hold on;
        totalmax = max(totalmax, max(stim));
        
        index = round((1:numbins)*numel(stim)/numbins);
        stim = stim(index);
        
%         numzeros = numbins - mod(numel(stim), numbins);
%         stim = [stim zeros(1,numzeros)];
%         stim = mean(transpose(reshape(stim, [], numbins)));
        transform(t,:) = stim;  
        hold on;
        plot(stim)
    end
    
end