%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function generates the input spike trains for the paired and the
% unpaired phase. In the paired phase, one randomly chosen input channel
% receives Poisson input. In the unpaired phase, all but the previously
% paired input channels receive Poisson input.
%
%
% This code is used in the "Biophysical_Model.m" for the manuscript:
%
% Heterosynaptic Plasticity Determines the Set-Point for Cortical Excitatory-
% Inhibitory Balance (2018)
% Rachel Field, James D'amour, Robin Tremblay, Christoph Miehl, Bernardo Rudy, 
% Julijana Gjorgjieva, Robert Froemke
% bioRxiv, doi: https://doi.org/10.1101/282012
%
%
% The code was written by Christoph Miehl (christoph.miehl@brain.mpg.de).
% The concept was developed by Christoph Miehl and Julijana Gjorgjieva (gjorgjieva@brain.mpg.de).
% July 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [presyn_input_mat,paired_channel]=InputGeneration_Biophysical_Model(paired_channel_before,N_exc,N_inh,type_of_input,dt,firing_rate_E,firing_rate_I,bin_size,number_neurons_per_pattern)

rng('shuffle');
bin_size_timesteps=bin_size/dt; % Convert from ms into number of timesteps
presyn_input_mat=zeros(N_exc+N_inh,bin_size_timesteps); % Predefine presynaptic input matrix


if type_of_input==1 % Generates Poisson inputs for all input neurons, besides the previously paired channel. This corresponds to the unpaired phase.
    paired_channel=paired_channel_before;
    
    % Generate Poisson input for excitatory inputs
    for ll2=1:N_exc
        if (ll2<(1+(paired_channel-1)*number_neurons_per_pattern)) || (ll2>paired_channel*number_neurons_per_pattern)
            for ll1=1:bin_size_timesteps              
                if rand <=firing_rate_E/1000*dt
                    presyn_input_mat(ll2,ll1)=1;
                end
            end
        end
    end
    
    % Generate Poisson inputs for inhibitory inputs
    for ll3=1:N_inh
        if (ll3<(N_exc+1+(paired_channel-1)*number_neurons_per_pattern)) || (ll3>(N_exc+paired_channel*number_neurons_per_pattern))
            for ll4=1:bin_size_timesteps              
                if rand <=firing_rate_I/1000*dt
                    presyn_input_mat(N_exc+ll3,ll4)=1;
                end
            end
        end
    end
     
elseif type_of_input==2 % Generates Poisson inputs to one randomly chosen channel only, weights from other channels receive no input. This corresponds to the paired phase.
    
    number_of_possible_patterns=N_exc/number_neurons_per_pattern;
    paired_channel=randi(number_of_possible_patterns);
    
    for ll2=paired_channel*number_neurons_per_pattern-(number_neurons_per_pattern-1):paired_channel*number_neurons_per_pattern
        for ll1=1:bin_size_timesteps
            % For excitatory inputs          
            if rand <=firing_rate_E/1000*dt 
                presyn_input_mat(ll2,ll1)=1;
            end
            % For inhibitory inputs          
            if rand <=firing_rate_I/1000*dt
                presyn_input_mat(N_exc+ll2,ll1)=1;
            end
        end
    end
end








