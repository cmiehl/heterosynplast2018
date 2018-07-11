%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The following code generates excitatory and inhibitory tuning curves where the
% individual weights (channels) are drawn from a uniform distribution. For
% each tuning curve, a randomly-chosen channel is increased/decreased by a 
% fixed amount implementing homosynaptic plasticity and the maximum strength channel 
% is increased/decreased by a fixed amount implementing heterosynaptic plasticity. 
% If the maximum strength channel is the channel which increased/decreased due to
% homosynatic plasticity, the second maximum strength channel is changed by
% heterosynaptic plasticity. This procedure is repeated for a number of
% tuning curve initializations and the E/I correlations are calculated before and
% after plasticity induction. At the end, we calculate the probability that 
% the E/I correlation-after will increase, decrease or stay the same 
% for a certain value of E/I correlation-before.
% 
%
%
% This code is the basis for the "probabilistic model" in the manuscript:
%
% Heterosynaptic Plasticity Determines the Set-Point for Cortical Excitatory-
% Inhibitory Balance (2018)
% Rachel Field, James D'amour, Robin Tremblay, Christoph Miehl, Bernardo Rudy, 
% Julijana Gjorgjieva, Robert Froemke
% bioRxiv, doi: https://doi.org/10.1101/282012
%
%
% This code was used to generate Fig. 3A-C.
% To titrate the plasticity ratio of heterosynaptic to homosynaptic plasticity,
% uncomment lines 55-58, 179-180, 222-226
% i.e. the for-loop with kk2 as an index and the plotting in figure(2).
%
%
% The code was written by Christoph Miehl (christoph.miehl@brain.mpg.de).
% The concept was developed by Christoph Miehl and Julijana Gjorgjieva (gjorgjieva@brain.mpg.de).
% July 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
rng('shuffle');

%% Parameters of the model
channels=12; % Number of input channels (6 or 12)
numb_sims=50000; % Number of randomly generated tuning curves

% Set the absolute homosynaptic/heterosynaptic weight change. The
% corresponding weights are multiplied by this value, therefore if the
% value is >1 the weights will increase and if it is <1 they will decrease.
eLT=1.65; %homosynaptic potentiation of excitatory weights by 65% of the original synaptic strength
iLT=1.65; %homosynaptic potentiation of inhibitory weights by 65% of the original synaptic strength
ehet=0.4; %heterosynaptic depression of excitatory weights by 60% of the original synaptic strength
ihet=0.4; %heterosynaptic depression of inhibitory weights by 60% of the original synaptic strength


% Uncomment this part if you would like to titrate the plasticitity ratio of heterosynaptic vs homosynaptic plasticity
% het_decrease=[0.98:-0.02:0.02];
% for kk2=1:length(het_decrease)
%     ehet=1-het_decrease(kk2);
%     ihet=1-het_decrease(kk2);
    
    %% Initialization of parameters
    prob_matrix=zeros(201,201);
    
    vec_to_check_before=-100:1:100;
    vec_to_check_after=100:-1:-100;
    
    %% Start of simulation
    for kk=1:numb_sims
        
        rand_mat=rand(channels,2); % Generate random initial tuning curves (uniform distribution)
        start_rand_mat=rand_mat;
        
        for jj=3:2:(channels*2-1) % Copy the generated numbers, so the same pairing experiment can be done with all the available channels
            rand_mat(:,jj:jj+1)=rand_mat(:,1:2);
        end
        
        calc_corr_start=corr(rand_mat(:,1),rand_mat(:,2)); % Calculate the initial correlation of the randomly generated tuning curves
        
        counter=1;
        for jj2=1:channels
            
            % Synaptic weight-change based on homosynaptic plasticity
            rand_mat(jj2,counter)=rand_mat(jj2,counter)*eLT;
            rand_mat(jj2,counter+1)=rand_mat(jj2,counter+1)*iLT;
            
            % Synaptic weight-change based on heterosynaptic plasticity
            dummyE=rand_mat(:,counter);
            dummyI=rand_mat(:,counter+1);
            %find the channels that is maximally strong
            [a,maxE]=max(dummyE);
            [a,maxI]=max(dummyI);
            if maxE==jj2 % If the maximum is the paired channel, use the second maximum
                [a,maxE2]=max(dummyE(dummyE<max(dummyE)));
                if maxE2>=jj2
                    maxE=maxE2+1;
                else
                    maxE=maxE2;
                end
            end
            if maxI==jj2
                [a,maxI2]=max(dummyI(dummyI<max(dummyI)));
                if maxI2>=jj2
                    maxI=maxI2+1;
                else
                    maxI=maxI2;
                end
            end
            
            rand_mat(maxE,counter)=rand_mat(maxE,counter)*ehet;
            rand_mat(maxI,counter+1)=rand_mat(maxI,counter+1)*ihet;
            
            
            calc_corr_afterboth(jj2)=corr(rand_mat(:,counter),rand_mat(:,counter+1));
            counter=counter+2;
        end
        
        
        % Puts the correlation values before and after plasticity induction into a matrix,
        % which counts how often a correlation value before will become a certain correlation value after.
        % The result is a matrix where the columns correspond to the correlation value before 
        % (Matrix dimension is 201x201 and therefore the resolution of the correlation values is 0.01, with correlations from -1 to 1)
        % and the rows correspond to the correlation value after the experiment. 
        % The entires in the matrix are the probabilities that a certain E/I correlation value is reached based on a certain
        % E/I correlation before the experiment.
        for jj4=1:length(calc_corr_afterboth)
            check_after_both=round(calc_corr_afterboth(jj4)*100);
            for jj3=1:length(vec_to_check_before)
                if round(calc_corr_start*100)==vec_to_check_before(jj3)
                    index_before=jj3;
                end
                if check_after_both==vec_to_check_after(jj3)
                    index_after_both=jj3;
                end
            end
            prob_matrix(index_after_both,index_before)=prob_matrix(index_after_both,index_before)+1;
        end
        
    end
    
    sum_matrix_both=sum(prob_matrix);
    for ii=1:length(prob_matrix)
        prob_matrix(:,ii)=prob_matrix(:,ii)./sum_matrix_both(ii);
    end
    
    % Calculates the probability that for a certain E/I correlation before plasticity induction, 
    % the E/I correlation after will increase, decrease or stay the same.
    for ii2=1:length(prob_matrix)
        fraction_above(ii2)=sum(prob_matrix(1:(length(prob_matrix)-ii2),ii2));
        fraction_below(ii2)=sum(prob_matrix((length(prob_matrix)+2-ii2):length(prob_matrix),ii2));
        fraction_same(ii2)=prob_matrix(length(prob_matrix)+1-ii2,ii2);
    end
    
    % Calculate the crossing point of the two probability curves
    get_cross=find(fraction_above+fraction_same/2<fraction_below);
    
    if isempty(get_cross)==1
        get_cross=length(fraction_above);
    end
    
    max_value_of_corr=length(prob_matrix(1,:));
    length_we_check=20; % Length we check where the number of max_wrong_ones is not exceeded
    max_wrong_ones=3; % Maximum of wrong values we accept to happen in this interval
    
    for zz=1:length(get_cross)
        if length(get_cross)-zz<length_we_check
            length_we_check=length(get_cross)-zz;
            max_wrong_ones=2;
        elseif length(get_cross)-zz<max_wrong_ones
            zz=length(get_cross);
            break;
        end
        if get_cross(zz+length_we_check)<=(get_cross(zz)+length_we_check+max_wrong_ones)
            break;
        end
    end
    if get_cross(zz)<max_value_of_corr-3
        get_cross=get_cross+2;
    end
    vec_correlations=-1:0.01:1;
    
    crossing_point=vec_correlations(get_cross(zz));

    
% Uncomment this part if you would like to titrate the plasticitity ratio of heterosynaptic vs homosynaptic plasticity
%     save_crossing_points(1,kk2)=crossing_point;
% end

%% Plotting of the simulation results
% Note: comment this text if running Fig 3C, to prevent a figure popping out for each
% tuning curve simulation.
figure(1)
subplot(2,2,1)
hold on
scatter(1:channels,start_rand_mat(:,1),'b')
scatter(1:channels,start_rand_mat(:,2),'r')
plot(1:channels,start_rand_mat(:,1),'b')
plot(1:channels,start_rand_mat(:,2),'r')
xlim([0.9 channels+0.1]);
hold off
xlabel('Stimulus channel')
ylabel('Synaptic strength (a.u.)')
legend('E channels','I channels')
title(['Tuning curve example - before, Corr: ' num2str(calc_corr_start)])

subplot(2,2,3)
hold on
scatter(1:channels,rand_mat(1:channels,1),'b')
scatter(1:channels,rand_mat(1:channels,2),'r')
plot(1:channels,rand_mat(1:channels,1),'b')
plot(1:channels,rand_mat(1:channels,2),'r')
xlim([0.9 channels+0.1]);
hold off
xlabel('Stimulus channel')
ylabel('Synaptic strength (a.u.)')
legend('E channels','I channels')
title(['Tuning curve example - after, Corr: ' num2str(corr(rand_mat(1:channels,1),rand_mat(1:channels,2)))])

subplot(2,2,[2,4])
hold on
scatter(vec_to_check_before./100,fraction_above,'blue')
scatter(vec_to_check_before./100,fraction_below,'green')
scatter(vec_to_check_before./100,fraction_same,'red')
plot([crossing_point,crossing_point],[0,1],'black')
hold off
xlabel('r_{ei}-before')
ylabel('r_{ei} change probability')
ylim([0 1])
xlim([-1 1])
legend('Before<After','Before>After','Before=After')

% figure(2)
% scatter(het_decrease./(eLT-1),save_crossing_points,'black')
% xlabel('Plasticity ratio (het/hom)')
% ylabel('Equilibrium point')
% ylim([0 1])
% legend([num2str(channels) ' stim channels'])