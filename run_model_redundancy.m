close all
clear all

f =  filesep;
addpath(['.' f 'functions'])
% NB this network runs every signal twice, with the first 500 ms different,
% to test the trial to trial variability. 

%% Parameters

T = 3000;           % total time
dt = .1;             % dt
Ntime = T/dt;
timestart = 500;
Ntimestart = timestart/dt;

% kernels
options.normalizegin = 1;
Tf = 50;           % time kernels
Nbf = 8;            % # basis functions for kernels
tfilt = 0:dt:Tf;    % time of filter
delta = 7.5;         % delay time

% network
options.multspikes = 0;     % only one neuron can spike at each dt
options.multspikerand = 1;  % choose which neuron spikes randomly
Nneuron = 100;              % # neurons
net = 't12';                % choose 'homogeneous_t1' or 'heterogeneous' or 't12' or
% 'heterogeneous_matchedt12
seed = 5;

% spike cost
cost = 'abs';
nu = 1.5;             % relative spike cost (absolute = nu*Th)
aamp = 1.5;           % relative adaptive spike cost 
taua = 60;         % time constant adaptation

% input signal

tausigvec = [1:30 40 50];            % time constant signal
Nt = length(tausigvec);
ampvec = [1:50];               % amplitude signal
Na = length(ampvec);

% Coincidence factor
precision = 2;

Ntrial = 1;

% For saving
% savestring = ['redundancy_',net ]; 

%% Make Network

kernel = make_kernels_network(tfilt, Nbf, Nneuron, net, seed);

[tg, g, gin, gout, Th ] = generate_filters( tfilt, kernel, 1:Nneuron, delta, options);

%% Run
rate =  zeros(Nt, Na, 2*Ntrial);
rate1 =  zeros(Nt, Na, Nneuron);
rate2 =  zeros(Nt, Na, Nneuron);
MSE0 =  zeros(Nt, Na, 2*Ntrial);
MSE =   zeros(Nt, Na, 2*Ntrial);
Gamma = zeros(Nt, Na, 2*Ntrial, Nneuron);
si =    zeros(Ntrial, Ntime);


O = zeros(2, Nneuron, Ntime);
for tn = 1:Nt
    % different autocorrelation times signal
    disp(['tau = ', num2str(tausigvec(tn))])
    %% Make signals
    Tsigfilt = 5*tausigvec(tn);         % time filter noisesignal
    tsigfilt = 0:dt:Tsigfilt;
    filtsignal = exp(-tsigfilt/tausigvec(tn))/sum(exp(-tsigfilt/tausigvec(tn)));    % Filter signal
    for ntr = 1:Ntrial
        disp(['trial = ', num2str(ntr)])
        rng('shuffle')
        si((ntr-1)*2+1,:) = randn(1,Ntime);
        si((ntr-1)*2+1,:)=conv(si((ntr-1)*2+1,:),filtsignal,'same');
        si((ntr-1)*2+1,:)=conv(si((ntr-1)*2+1,:),fliplr(filtsignal),'same');
        
        % change first 500 ms, for test trial to trial variability
        si(ntr*2,:) = si((ntr-1)*2+1,:);
        sistart = randn(1,Ntimestart);
        sistart=conv(sistart,filtsignal,'same');
        sistart=conv(sistart,fliplr(filtsignal),'same');
        si(ntr*2,1:Ntimestart) = sistart;

        for an = 1:Na   
            disp(['amplitude = ', num2str(ampvec(an))])
            % different amplitudes signal
            
            for ntr2 = 1:2
                % Two trials per signal, with different starts
                disp('Make signal')
                sitemp = si((ntr-1)*2+ntr2,:)*ampvec(an)/std(si((ntr-1)*2+ntr2,:));
                MSE0(tn,an,(ntr-1)*2+ntr2) = calc_MSE(sitemp, zeros(size(sitemp)));

                
                disp('Run')
                if strcmp(cost, 'rel')
                    evalc('[xest, O(ntr2,:,:), ~, ~] = run_relcost(dt, sitemp, g, tg, gin, gout, Th, nu, delta, aamp, taua, options)');
                elseif strcmp(cost, 'abs')
                    evalc('[xest, O(ntr2,:,:), ~, ~] = run_abscost(dt, sitemp, g, tg, gin, gout, Th, nu, delta, aamp, taua, options)');
                end
                
                
                
                MSE(tn,an,(ntr-1)*2+ntr2) = calc_MSE(sitemp, xest);
                rate(tn,an,(ntr-1)*2+ntr2) = 1000*sum(sum(O(ntr2,:,:)))/(T*Nneuron);
                
                disp(['Relative MSE = ',num2str(MSE(tn,an,(ntr-1)*2+ntr2)/MSE0(tn,an,(ntr-1)*2+ntr2))])
                disp(['Average firing rate = ',num2str(rate(tn,an,(ntr-1)*2+ntr2))])
                
                
            end

            for nn=1:Nneuron
                eventtimes1 = find(squeeze(O(1, nn, Ntimestart+1:end))==1)*dt;
                eventtimes2 = find(squeeze(O(2, nn, Ntimestart+1:end))==1)*dt;
                rate1(tn, an, nn) = length(eventtimes1)/(T-timestart);
                rate2(tn, an, nn) = length(eventtimes2)/(T-timestart);
            
                [~, ~, ~, Gamma(tn,an,2*ntr-1,nn)] = calccofac_ignoredoublespikes(eventtimes1, eventtimes2, rate2(tn, an, nn), precision,0);
                [~, ~, ~, Gamma(tn,an,2*ntr  ,nn)] = calccofac_ignoredoublespikes(eventtimes2, eventtimes1, rate1(tn, an, nn), precision,0);
            end

            disp(['Average Gamma = ', num2str(nanmean(nanmean(Gamma(tn,an,:,:),4),3))])
            
        end
        

        
    end
    save(savestring)
end



%% Plot
cbarvecmse = 0:0.1:0.1;
% cbarvecmse = 0:0.01:0.05;
cbarvecrate = 0:10:90;
cbarvecgamma = 0:.1:1;

splts = [2 3 1];
splts = splts+6;

% Figure1=figure(1);clf;
% set(Figure1,'defaulttextinterpreter','latex');
subplot(3,3,splts(1))
MSEnorm = squeeze(mean(MSE./MSE0, 3));
% surf(ampvec,tausigvec,MSEnorm,'EdgeColor','none');
% axis xy; axis tight; view(0,90);
h = pcolor_fleur(ampvec,tausigvec,MSEnorm);
set(h, 'EdgeColor','none');
% shading interp
colormap(jet); 
ylim([1 30])
title('$\overline{MSE}$')
ylabel('\tau (ms)')
xlabel('amplitude signal')
c = colorbar('Ticks',0:0.025:0.1);
set(c, 'TickLabels',{'0','.025','.05','.075','>.1'})
%,'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})
caxis([0 0.1])
% c.Label.String = 'MSE';


subplot(3,3,splts(2))
rateav = squeeze(mean(rate, 3));
% surf(ampvec,tausigvec,rateav,'EdgeColor','none');
% axis xy; axis tight; view(0,90);
h = pcolor_fleur(ampvec,tausigvec,rateav);
set(h, 'EdgeColor','none');
% shading interp
colormap(jet); 
ylim([1 30])
title('Activity A (Hz)')
ylabel('\tau (ms)')
xlabel('amplitude signal')
c = colorbar('Ticks',cbarvecrate);
caxis([cbarvecrate(1) cbarvecrate(end)])
% c.Label.String = 'frequency (Hz)';

subplot(3,3,splts(3))
Gammaav = squeeze(nanmean(Gamma, 3));
Gammaav = squeeze(nanmean(Gammaav, 3));
% surf(ampvec,tausigvec,Gammaav,'EdgeColor','none');
% axis xy; axis tight; view(0,90);
h = pcolor_fleur(ampvec,tausigvec,Gammaav);
set(h, 'EdgeColor','none');
% shading interp
colormap(jet); 
ylim([1 30])
title('$\overline{\Gamma}$')
ylabel('\tau (ms)')
xlabel('amplitude signal')
c = colorbar('Ticks',cbarvecgamma);
caxis([cbarvecgamma(1) cbarvecgamma(end)])
% c.Label.String = '$\Gamma$';

% saveas(gcf, savestring)
