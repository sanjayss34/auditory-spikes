function [] = intensity_trains()
    load('sorted_trains.mat');
    num_intensities = size(sorted_trains,2);
    for i = 1:num_intensities
        trains = sorted_trains(:,i);
        filename = ['TNR',num2str(i),'_trains.mat'];
        save(filename, 'trains');
    end
end