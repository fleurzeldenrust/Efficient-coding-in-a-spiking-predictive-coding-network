function  [rcc, scc, ncc] = calc_cor_signalnoise(Qmat, neuronvec, tlshow, dt)
% Calculate noise and signal correlations between two spike trains following the method in 
% Bair, W., E. Zohary, and W. T. Newsome. 2001. 
% ?Correlated Firing in Macaque Visual Area MT: Time Scales and Relationship to Behavior.? The Journal of neuroscience 21(5):1676?97.
% Input:
% * Qmat = NtrialxNneuronxNtime
% * neuronvec: which neurons
% * tlshow: time window for plotting
% * dt
% Output:
% * rcc: raw cross-correlogram
% * scc: signal cross-correlogram
% * ncc: noise cross-correlogram
% NB needs av_corr (corrected for number of spikes)

[Ntrial, ~, Ntime] = size(Qmat);

n1  = neuronvec(1);
n2  = neuronvec(2);    

L = 2*Ntime-1;
Lm = ceil(L/2);
tshow = tlshow(1):dt:tlshow(2);
Ntshow = floor(length(tshow)/2);

T = Ntime*dt/1000; % (time in seconds)


%% Raw correlations


% cross-correlations
% disp('Calculate raw cross-correlations')
rcc  = av_corr(Qmat,  n1, n2)*T;

%% Signal correlations = shift predictor
if Ntrial>1
    Qmatav = sum(Qmat)/Ntrial; % PSTH
else
    Qmatav = Qmat;
end
% cross-correlations
% disp('Calculate signal cross-correlations (shift predictor)')
scc  = av_corr(Qmatav,  n1, n2)*T;

%% Noise correlations
% cross-correlations
% disp('Calculate noise cross-correlations')
ncc  = rcc-scc;


%% Other way of computing noise correlations
% Subtract mean from every spike train, and compute cross-correlations
% between corrected spike trains
% gives similar shape, but different y axis

% disp('Other way of computing')
% keyboard
% if Ntrial>1
%     Qmat_corr = Qmat - repmat(Qmatav, Ntrial,1);
% else
%     Qmat_corr = Qmat;
% end
% 
% if Ntrial>1
%     %% Noise correlations
%     % cross-correlations
%     ncc_12_v2  = av_corr(Qmat_corr,  n1, n2)/T;
% end

%% Plot
% figure
% subplot(2,2,1)
% plot(tshow, rcc(Lm-Ntshow:Lm+Ntshow));
% title('Raw Cross-correlations')
% xlabel('lags (ms)')
% 
% subplot(2,2,2)
% hold all
% plot(tshow, scc(Lm-Ntshow:Lm+Ntshow));
% title('Signal correlations: Cross-correlations')
% xlabel('lags (ms)')
% 
% subplot(2,2,3)
% hold all
% plot(tshow, ncc(Lm-Ntshow:Lm+Ntshow));
% title('Noise correlations: Cross-correlations, method 1')
% xlabel('lags (ms)')

% subplot(2,2,4)
% hold all
% plot(tshow, ncc_12_v2(Lm-Ntshow:Lm+Ntshow));
% title('Noise correlations: Cross-correlations, method 2')
% xlabel('lags (ms)')    
end