%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements synaptic plasticity according to spike-timing-dependent
% plasticity (STDP) rules. Excitatory synaptic weights change based on a
% classic Hebbian learning window, where a presynaptic spike followed by a
% postsynaptic spike leads to long-term potentiation (LTP) and a
% postsynaptic spike followed by a presynaptic spike leads to long-term
% depression (LTD). Inhibitory weights change based on a symmetrical
% learning window, where weights increase for both, pre-post and post-pre
% pairing.
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

function W=STDP_Biopysical_Model(W,tau_w_E,tau_w_I,r_1,o_1,A_LTD_E,A_LTD_I,A_LTP_E,A_LTP_I,pre_post,E_or_I,bounds)

% Set upper and lower bounds for the weights
W_min_E=bounds(1);
W_max_E=bounds(2);
W_min_I=bounds(3);
W_max_I=bounds(4);

if E_or_I==1 % Excitatory STDP learning rule
    if pre_post==1 % Presynaptic spike
        W_E=W-(o_1.*A_LTD_E)./tau_w_E;
    elseif pre_post==2 % Postsynaptic spike
        W_E=W+(r_1.*A_LTP_E)./tau_w_E;
    end
    
    % Check if weights are out of bound - this almost never happens
    check_out_of_bound=(W_E>W_min_E & W_E<W_max_E); 
    check_out_of_bound2=(W_E<W_max_E); 
    check_out_of_bound3=(W_E>W_min_E); 
    W=W_E.*check_out_of_bound+(ones(length(W_E),1)-check_out_of_bound2).*W_max_E+(ones(length(W_E),1)-check_out_of_bound3).*W_min_E; 
    
elseif E_or_I==2 % Inhibitory STDP learning rule
    if pre_post==1 % Presynaptic spike
        W_I=W+A_LTD_I.*o_1./tau_w_I;
    elseif pre_post==2 % Postsynaptic spike
        W_I=W+A_LTP_I.*r_1./tau_w_I;
    end
     
    % Check if weights are out of bound
    check_out_of_bound=(W_I>W_min_I & W_I<W_max_I);
    check_out_of_bound2=(W_I<W_max_I); 
    check_out_of_bound3=(W_I>W_min_I); 
    W=W_I.*check_out_of_bound+(ones(length(W_I),1)-check_out_of_bound2).*W_max_I+(ones(length(W_I),1)-check_out_of_bound3).*W_min_I;
end
