function  [cc_netneur, sta, tsta] = calc_cor_networkneuron(Qmat, st, dt, kernel, tkernel, plotyn)
%% Correlations between network activity and spike train
% Input:
% * Qmat: zero and one matrix, size Ntrial x Nneuron x Ntime
% * st: either integer (which spiketrain to take from Qmat) or spiketrain
% (size NtrialxNtime)
% * dt
% * kernel, tkernel: for if you want to convolve the spike train with a
% kernel
% * plotyn: 0 or 1, whether you want to plot
% Output
% * cc_netneur: 2*Ntime-1 matrix with the average CCG between the neuron
% and the network activity
% * sta: Ntrial x ... matrix with the STA between the neuron and the
% network for every trial

%% params
if nargin == 3
    kernel  = [];
    tkernel = [];
    plotyn  = 1;
else
    if abs(dt*sum(kernel)-1)<0.00001
        disp('normalized kernel')
        % NB Note that we normalized taking dt into account
    else
        ip = input('non-normalized kernel, do you want to normalize? (y/n)', 's');
        if strcmp(ip, 'y')
            kernel = kernel./sqrt(0.5*dt*sum(kernel.^2));
        end
    end
end

[Ntrial, Nneuron, Ntime] = size(Qmat);
T = Ntime*dt;
Nsp_net = sum(Qmat,3);


la = sum(tkernel<0);
lc = sum(tkernel>0);

[Ntrialst, Ntimest] = size(st);

if Ntimest == 1 && Ntrialst == 1
    % st is neuronnumber from Qmat
    spiketrain = squeeze(Qmat(:,st,:));
elseif ~(Ntrialst == Ntrial) || ~(Ntimest == Ntime)
    error('Give spike train (size NtrialxNtime) or neuron number')
else
    % st is independent spike train
    spiketrain = st;
end
Nsp_neuron = sum(spiketrain,2);

%% Calculate network activity
netact = zeros(Ntrial, Ntime);
for nt = 1:Ntrial
    netact_temp = zeros(1,Ntime);
    for nn=1:Nneuron
        netact_1n = squeeze(Qmat(nt,nn, :));
        if ~isempty(kernel)
            netact_1n = convolve_kernel_acausal(netact_1n', kernel, la, lc);
        end
        netact_temp = netact_temp + netact_1n;
    end
    netact(nt,:) = netact_temp;
end

if ~isempty(kernel)
    neuronact = zeros(Ntrial, Ntime);
    for nt=1:Ntrial
        neuronact(nt,:) = convolve_kernel_acausal(spiketrain(nt,:), kernel, la, lc);
    end
else 
    neuronact = spiketrain;
end 
keyboard
%% Correlate network activity with spike train
L = 2*Ntime-1;
Lm = ceil(L/2);
cc_temp = zeros(L,1);
for nt=1:Ntrial
    cc_temp = cc_temp+xcorr(squeeze(netact(nt,:)), squeeze(neuronact(nt,:)))';
end
N1 = sum(sum(Nsp_net))/(Ntrial);
% You can also choose to use the network activity per neuron (Ntrial ->
% Ntrial*Nneuron) 
N2 = sum(Nsp_neuron)/Ntrial;
cc_netneur = T*cc_temp/(sqrt(N1*N2)*Ntrial*1000);

%% STA neuron - network activity (see Okun, Michael et al. 2015. ?Diverse Coupling of Neurons to Populations in Sensory Cortex.? Nature 521(7553):511?15.)
str = struct('dt',dt,'twindow', [-300 200],'spiketype','spiketrain', 'noplot',1);
tsta = str.twindow(1):dt:str.twindow(2);
sta = zeros(Ntrial, length(tsta));
for nt = 1:Ntrial
    sta(nt,:) = calc_sta_conv(netact(nt,:), spiketrain(nt,:), str);
end

%% Plot
if plotyn == 1
    time = dt*(1:Ntime);
    % Network activity 
    figure
    subplot(3,1,1)
    hold all
    plot(time, netact(1,:), 'k', 'LineWidth',2)
    xlim([10000 12000])
    title('Network activity')
    ylabel('activity (kHz)')

    axes212 = subplot(3,1,2);
    hold all
    for nn=1:Nneuron
        spiketimes = find(Qmat(1,nn,:)==1)*dt;
        plot(spiketimes, nn*ones(size(spiketimes)),  '.k')
    end
    plot(st(1,:), (Nneuron+1)*ones(size(spiketimes)),  '.r')
    plot([0 T],[25.5 25.5],'k','LineWidth',2)
    plot([0 T],[50.5 50.5],'k','LineWidth',2)
    plot([0 T],[75.5 75.5],'k','LineWidth',2)
    set(axes212,'YTickLabel',{'type 1','-type 1','type 2','-type 2'},...
        'YTick',[12.5 37.5 62.5 87.5]);
    xlim([10000 12000])
    xlabel('time (ms)')

    subplot(3,2,5)
    tshow = -150:dt:150;
    Lshow = floor(length(tshow)/2);
    hold all
    plot(tshow, cc_netneur(Lm-Lshow:Lm+Lshow))
    xlabel('lags(ms)')
    ylabel('# coincidences / spike')
    title('Correlation neuron - network activity')

    subplot(3,2,6)
    plot(tsta, mean(sta))
    xlabel('time relative to spike (ms)')
    ylabel('average network activity (kHz)')
    title('Spike-triggered network activity')
end

        

