
Let us say you have recordings over T trials for N neurons and S stimuli.

spike_train_processing.m takes 2 inputs: trial_length and bin_length. It removes everything recorded after the trial has ended, i.e. time > trial_length. It then creates spike trains by binning with bin size = bin_length.

concatenate_trains.m concatenates these spike trains for every tuple (neuron, stimulus). Thus, it returns an NxS array of spike trains, which is saved in 'sorted_trains.mat'

create_neuron_trains.m further concatenates those spike trains across neurons, irrespective of the stimulus. Thus, this is used to create stimulus-independent Ising models and returns an Nx1 array of spike trains, saved in 'neuron_trains.mat'

intensity_trains concatenates the spike trains according to the stimulus. Thus, it returns S distinct arrays, one for each stimulus. These are used to create stimulus-dependent Ising models. Each array is of size Nx1, saved in a file named 'TNRx_trains.mat', where x is an integer representing the stimulus.

All of these functions are called by spike_train_generator.m. Thus, to process and generate all of these spike train files, simply use spike_train_generator and specify (1) trial_length and (2) bin_length.