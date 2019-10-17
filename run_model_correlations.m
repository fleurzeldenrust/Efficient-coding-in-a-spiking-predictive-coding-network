close all
clear variables
clear global

f =  filesep;
addpath(['.' f 'functions'])
%% For testing
test = 1;       % if test == 1: short simulation with report on results for testing

%% Type of experiment
expi = 'nonoise';      % type of experiment: nonoise, indnoise or cornoise

%% Params

% T = 100000;               % total time
T = 20000;
dt = .2;                    % dt
Ntime = T/dt;
time = (1:Ntime)*dt;

% kernels
options.normalizegin = 1;
Tf = 50;                    % time kernels
Nbf = 8;                    % # basis functions for kernels
tfilt = 0:dt:Tf;            % time of filter
delta = 7.5;                % delay time
Ndelta = round(delta/dt);

% network
options.multspikes = 0;     % only one neuron can spike at each dt
options.multspikerand = 1;  % choose which neuron spikes randomly
Nneuron = 100;              % # neurons
net = 't12';                % choose 'homogeneous_t1' or 'heterogeneous' or 't12' or
% 'heterogeneous_matchedt12
seed = 5;

% spike cost
cost = 'abscost';       % abscost or relcost
nu = 1.5;               % relative spike cost (absolute = nu*Th)
aamp = 1.5;             % relative adaptive spike cost 
taua = 60;              % time constant adaptation

%% INPUT SIGNAL
signal = 'exp';             % choose 'exp' 'pink' 'powerlaw'

% exponential signal
tausig = 15;                % time constant signal
Tsigfilt = 5*tausig;        % time filter noisesignal
tsigfilt = 0:dt:Tsigfilt;

% powerlaw signal
% alpha = 1.5;

%% NOISE
noisetype = 'exp';          % choose exp or powerlaw
% exponential noise

% independent noise
taunoisei = 15;              % time constant noise
filtnoisei = exp(-tsigfilt/taunoisei)/sum(exp(-tsigfilt/taunoisei));
% shared noise
taunoises = 15;              % time constant noise
filtnoises = exp(-tsigfilt/taunoises)/sum(exp(-tsigfilt/taunoises));
% powerlaw noise
alphanoisei = 1.5;
alphanoises = 1.5;

% downsampling for saving
ndown = 5;
Ntimedown = Ntime/ndown;
savenr = 10; % save every # of trials

%% Set params according to type of experiment
if strcmp(expi, 'nonoise')
    if strcmp(signal, 'exp') && strcmp(noisetype, 'exp')
        asig = 5.5;
    elseif strcmp(signal, 'powerlaw') && strcmp(noisetype,'powerlaw')
        asig = 5.5;
    end
    Ntrial = 1;
    anoisei = 0;
    anoises = 0;
    noisecop = Nneuron;
elseif strcmp(expi, 'indnoise')
    if strcmp(signal, 'exp') && strcmp(noisetype, 'exp')
        f = 0.25;
        atot = 3;
        asig = sqrt((atot^2)/(1+f));
        anoisei = f*asig;
    elseif strcmp(signal, 'powerlaw') && strcmp(noisetype,'powerlaw')
        f = 0.25;
        atot = 3;
        asig = sqrt((atot^2)/(1+f));
        anoisei = f*asig;
    end
    Ntrial = 300;
    noisecop = Nneuron;
    anoises = 0;
elseif strcmp(expi, 'cornoise')
    % NB remember that standard deviations sum as stdtot =
    % sqrt(std1^2+std2^2); 
    if strcmp(signal, 'exp') && strcmp(noisetype, 'exp')
        % sqrt(2*0.23^2)=0.32!
%         asig = 1;
%         anoisei = 0.23;
%         anoises = anoisei;
        f = 0.25;
        atot = 3;
        asig = sqrt((atot^2)/(1+f));
        anoisetot = f*asig;
        anoisei = anoisetot / sqrt(2);
        anoises = anoisei;
    elseif strcmp(signal, 'powerlaw') && strcmp(noisetype,'powerlaw')
%         asig = 1;
%         anoisei = 0.27;
%         anoises = anoisei;
        f = 0.25;
        atot = 3;
        asig = sqrt((atot^2)/(1+f));
        anoisetot = f*asig;
        anoisei = anoisetot / sqrt(2);
        anoises = anoisei;
    end
    Ntrial = 300;
    noisecop = 10; % how many copies of the noise in network (so # neurons the same noise = Nneuron/noisecop)
end

savestring = ['sim_corr_', net,'_',signal,'signal','_',cost, '_',noisetype,'noise_',expi];
disp(savestring)

if test == 1
    Ntrial = 2;
    T = 2000;
    Ntime = T/dt;
    Ntimedown = Ntime/ndown;
    time = (1:Ntime)*dt;
    savestring = [];
end
    


%% generate representing filters 
kernel = make_kernels_network(tfilt, Nbf, Nneuron, net, seed);

%% generate lateral and input filters and threshold from representing filters
[tg, g, gin, gout, Th ] = generate_filters( tfilt, kernel, [], delta, options);

%% generate signal
if strcmp(signal, 'exp')
    filtsignal = exp(-tsigfilt/tausig)/sum(exp(-tsigfilt/tausig));
    si=randn(1,Ntime);
    si=conv(si,filtsignal,'same');
    si=conv(si,fliplr(filtsignal),'same');
    si = si(1:Ntime);
elseif strcmp(signal,'pink')
    hcn = dsp.ColoredNoise('InverseFrequencyPower',1, 'SamplesPerFrame',Ntime);
    si = step(hcn)';
elseif strcmp(signal, 'powerlaw')
    hcn = dsp.ColoredNoise('InverseFrequencyPower',alpha, 'SamplesPerFrame',Ntime);
    si = step(hcn)';
else
    disp('choose type of signal')
end
si = (si-mean(si))*asig/std(si);     



%% Run
% xestdown = zeros(Ntrial, Ntimedown);
Odown = zeros(Ntrial, Nneuron, Ntimedown);
Vdown = zeros(Ntrial, Nneuron, Ntimedown);
% Thdown = zeros(Ntrial, Nneuron, Ntimedown);
noisemat = zeros(Nneuron, Ntime);

for nt = 1:Ntrial
    %% generate noise 
    
    % independent noise
    disp(['trial number ', num2str(nt)])
    rng('shuffle')
    noisemati = randn(Nneuron, Ntime);
    if strcmp(noisetype, 'exp')
        for nn=1:Nneuron
            noisetempi = conv(noisemati(nn,:),filtnoisei,'same');
            noisetempi = conv(noisetempi,fliplr(filtnoisei),'same');
            noisetempi = noisetempi(1:Ntime);
            noisemati(nn,:) = (noisetempi-mean(noisetempi))/std(noisetempi);
        end
    elseif strcmp(noisetype, 'powerlaw')
        for nn=1:Nneuron
            rng('shuffle')
            hcn = dsp.ColoredNoise('InverseFrequencyPower',alphanoisei, 'SamplesPerFrame',Ntime);
            noisetempi = step(hcn)';
            noisetempi = (noisetempi-mean(noisetempi))/std(noisetempi); 
            noisemati(nn,:) = noisetempi;
        end
    else 
        disp('Choose type of noise')
        keyboard
    end
    noisemati = noisemati*anoisei;

    % correlated noise  
    noisemats = randn(Nneuron, Ntime);
    for nn=1:Nneuron
        if nn>noisecop
            % Take existing copy
            nmod = mod(nn,noisecop);
            if nmod == 0
                nmod = noisecop;
            end
            noisemats(nn,:) = noisemats(nmod,:);
        else
             % Make new copy
            if strcmp(noisetype, 'exp')
                noisetemp = conv(noisemats(nn,:),filtnoises,'same');
                noisetemp = conv(noisetemp,fliplr(filtnoises),'same');
                noisetemp = noisetemp(1:Ntime);
                noisemats(nn,:) = (noisetemp-mean(noisetemp))/std(noisetemp);
            elseif strcmp(noisetype, 'powerlaw')
                rng('shuffle')
                hcn = dsp.ColoredNoise('InverseFrequencyPower',alphanoises, 'SamplesPerFrame',Ntime);
                noisetemp = step(hcn)';
                noisetemp = (noisetemp-mean(noisetemp))/std(noisetemp); 
                noisemats(nn,:) = noisetemp;
            end
        end
    end
    noisemats = noisemats*anoises;

    noisemat = noisemati + noisemats;


    %% run network
    if strcmp(cost, 'abscost')
        if test == 1
            [xest, O, V, Thvec] = run_abscost_noise(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua, noisemat, options);
        else
            [~, O, V, ~] = run_abscost_noise(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua, noisemat, options);
        end
    elseif strcmp(cost, 'relcost')
        if test == 1
            [xest, O, V, Thvec] = run_relcost_noise(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua, noisemat, options);
        else
            [~, O, V, ~] = run_relcost_noise(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua, noisemat, options);
        end
    end

    if ndown>1
%         xestdown(nt, :) = downsample_mean(xest, ndown);
        for nn=1:Nneuron
            Odown(nt,nn,:) = downsample_sum(squeeze(O(nn,:)), ndown);
            Vdown(nt,nn,:) = downsample_mean(squeeze(V(nn,:)), ndown);
%             Thdown(nt,nn,:) = downsample_mean(squeeze(Thvec(nn,:)), ndown);
        end
    end 
    
    %% Plot, display, safe
    if test == 1
        Oreal = [zeros(Nneuron, Ndelta) O];
        Oreal = Oreal(:,1:Ntime);
        
        figure
        subplot(2,1,1)
        plot(time, si, 'k', 'LineWidth',2)
        hold all
        plot(time, xest, 'r','LineWidth',1)
        axes212 = subplot(2,1,2);
        hold all
        for nn=1:Nneuron
            spiketimes = find(Oreal(nn,:)==1)*dt;
            plot(spiketimes, nn*ones(size(spiketimes)),  '.k')
        end
        plot([0 T],[25.5 25.5],'k','LineWidth',2)
        plot([0 T],[50.5 50.5],'k','LineWidth',2)
        plot([0 T],[75.5 75.5],'k','LineWidth',2)
        set(axes212,'YTickLabel',{'type 1','-type 1','type 2','-type 2'},...
            'YTick',[12.5 37.5 62.5 87.5]);
        
        disp(['Mean firing rate = ',num2str(1000*sum(sum(O))/(T*Nneuron))])
        disp(['Mean firing rate = ',num2str(1000*sum(sum(squeeze(Odown(nt,:,:))))/(T*Nneuron))])

        disp(['Mean firing rate type 1 = ',num2str(1000*sum(sum(O(1:round(Nneuron/2),:)))/(T*round(Nneuron/2)))])
        disp(['Mean firing rate type 2 = ',num2str(1000*sum(sum(O(round(Nneuron/2):end,:)))/(T*round(Nneuron/2)))])
        
        if nt == 1
            noisetemp1 = noisemat(1,:);
        end
    end
    if ~(test == 1)
        if mod(nt, savenr)==0
            disp('Save')
            save(savestring, '-v7.3')
        end
    end
    
end

if ~(test == 1)
    disp('Save')
    save(savestring, '-v7.3')
    
    disp('Analyse')
    analyse_correlations(savestring)
else
    figure
    subplot(3,1,1)
    plot(si)
    hold all
    for nt = 1:Ntrial
        plot(xest)
    end

    subplot(3,1,2)
    hold all
    plot(noisetemp1)
    plot(squeeze(noisemat(1,:)))
    title('Noise different trials neuron 1')

    subplot(3,1,3)
    hold all
    nn=1;
    plot(squeeze(noisemat(nn,:)))
    if noisecop<Nneuron
        plot(squeeze(noisemat(nn+noisecop,:)))
        plot(squeeze(noisemat(nn+noisecop+1,:)))
    else 
        plot(squeeze(noisemat(end,:)))
    end
    title('Noise different neurons  trial 1')
end