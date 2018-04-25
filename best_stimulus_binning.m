function [] = best_stimulus_binning(numbins2)    

% Transform is the matrix of the amplitude, extracted from the raw
% stimulus.
figure;
    numbins = round(numbins2/3);
    load('rawStimuli.mat')
    totalmax = 0;
    numtrials = numel(rawStimCollector);
    transform = zeros(numtrials, numbins2);
    for t = 1:numtrials
        stim = abs(rawStimCollector{t});
        numzeros = numbins - mod(numel(stim), numbins);
        stim = [stim zeros(1,numzeros)];
        stim = transpose(reshape(stim, [], numbins));
        stim = max(stim, [], 2);
        stim = interp1(1:numbins, stim, linspace(1,numbins, 200));
        
%         box = ones(1,7);
%         convolution = conv(box,stim);
%         convolution = convolution(1:numel(stim));
%         stim = convolution*(max(stim)/max(convolution));
        
        transform(t,:) = stim;
        if t == 1 || t == 2
            subplot(2,1,t)
            plot(stim)
        end
        
        
    end
    
end