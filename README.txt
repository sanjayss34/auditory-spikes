tl;dr, Run spike_train_generator(trial_length, bin_length), with the 2 appropriate inputs.

--- SPIKE TRAIN GENERATION ---
Let us say you have recordings over T trials for N neurons and S stimuli. spike_train_generator.m calls the following 4 functions:

spike_train_processing.m trims the firing data and creates spike trains by taking 2 inputs: trial_length and bin_length. It eliminates data after the trial has ended, i.e. time > trial_length. It then creates spike trains by binning with bin size = bin_length.

concatenate_trains.m concatenates these spike trains for every (neuron, stimulus). Thus, it returns an NxS array of spike trains, which is saved in 'sorted_trains.mat'

create_neuron_trains.m further concatenates the spike trains in 'sorted_trains.mat' according to neuron, irrespective of the stimulus. Thus, it returns an Nx1 array of spike trains of neurons over all trials, which can be used for stimulus-independent Ising models and is saved in 'neuron_trains.mat'

intensity_trains.m further concatenates the spike trains in 'sorted_trains.mat' according to stimulus. It returns S distinct arrays, one for each stimulus. These can be used for stimulus-dependent Ising models. Each array is size Nx1 and is saved in a file named 'TNRx_trains.mat', where x is an integer representing the stimulus.
___________________________________________________________________

--- ISING MODEL ESTIMATION ---
After creating the 'neuron_trains.mat' file as steps above, run
[h, J, train_indices] = estimate_ising(ITERS)
to generate the Ising model parameters h and J along with the indices of the data used to train the model. Here, ITERS is the desired number of learning iterations. (1000 is a reasonable number to start)
