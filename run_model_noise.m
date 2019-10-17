close all
clear variables
clear global

f =  filesep;
addpath(['.' f 'functions'])

%% Parameters

T = 3000;               % total time
dt = .1;                % dt
Ntime = T/dt;

% kernels
options.normalizegin = 1;
Tf = 50;                % time kernels
Nbf = 8;                % # basis functions for kernels
tfilt = 0:dt:Tf;        % time of filter
delta = 7.5;            % delay time

% network
options.multspikes = 0;     % only one neuron can spike at each dt
options.multspikerand = 1;  % choose which neuron spikes randomly
Nneuron = 100;              % # neurons
net = 'heterogeneous_matchedt12';   % choose 'homogeneous_t1' or 'heterogeneous' or 't12' or
% 'heterogeneous_matchedt12
seed = 5;

% spike cost
cost = 'abs';
nu = 1.5;             % relative spike cost (absolute = nu*Th)
aamp = 1.5;           % relative adaptive spike cost 
taua = 60;         % time constant adaptation

% input signal
tausig = 15;            % time constant signal
amp = 10;               % amplitude signal
Tsigfilt = 5*tausig;         % time filter noisesignal and added noise
tsigfilt = 0:dt:Tsigfilt;
filtsignal = exp(-tsigfilt/tausig)/sum(exp(-tsigfilt/tausig));    % Filter 

% For saving
savestring = ['noise_',net ];

% noise params
ampvec_noise_rel = 0:0.05:1.2; 
ampvec_noise = ampvec_noise_rel*amp;
Na = length(ampvec_noise);
corvec = [1 2 5 10 20 25 50 100];
Nc = length(corvec);

Ntrial = 3;

%% Make Network

kernel = make_kernels_network(tfilt, Nbf, Nneuron, net, seed);

[tg, g, gin, gout, Th ] = generate_filters( tfilt, kernel, [], delta, options);

%% Run

rate =  zeros(Ntrial, Na, Nc);
MSE0 =  zeros(Ntrial,1);
MSE =   zeros(Ntrial, Na, Nc);
MSEnorm = zeros(Ntrial, Na, Nc);
noisemat = randn(Nneuron, Ntime);

for nt = 1:Ntrial
% for nt = nt:Ntrial
    
    disp(['Trial number ',num2str(nt)'])
    
    %% Make signal
    rng('shuffle')
    si = randn(1,Ntime);
    si=conv(si,filtsignal,'same');
    si=conv(si,fliplr(filtsignal),'same');
    si = si*amp/std(si);
    
    MSE0(nt) = calc_MSE(si, zeros(size(si)));
    
    for na = 1:Na
%     for na = na:Na;
        anoise = ampvec_noise(na);
        disp(['amplitude noise = ', num2str(anoise)])

        for nc = 1:Nc
%         for nc = nc:Nc
            cnoise = corvec(nc);
            disp(['correlation noise = ', num2str(cnoise)])

            %% Make noise
            disp('Make noise-matrix')
            rng('shuffle')

            for nn=1:Nneuron
                if nn>cnoise
                    % Take existing copy
                    nmod = mod(nn,cnoise);
                    if nmod == 0
                        nmod = cnoise;
                    end
                    noisemat(nn,:) = noisemat(nmod,:);
                else
                    % Make new copy
                    rng('shuffle')
                    noisetemp = randn(1,Ntime);
                    noisetemp = conv(noisetemp,filtsignal,'same');
                    noisetemp = conv(noisetemp,fliplr(filtsignal),'same');
                    noisetemp = noisetemp*anoise/std(noisetemp);
                    noisemat(nn,:) = noisetemp;
                end
            end

            %% Run
            disp('Run')
            
            if strcmp(cost, 'abs')
                [xest, O, ~, ~] = run_abscost_noise(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua, noisemat, options);
            elseif strcmp(cost, 'rel')
                [xest, O, ~, ~] = run_relcost_noise(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua, noisemat, options);
            end
            
            disp('Calculate')
            MSE(nt,na,nc) = calc_MSE(si, xest);
            MSEnorm(nt, na, nc) = MSE(nt, na, nc)./MSE0(nt);
            rate(nt,na,nc) = 1000*sum(sum(O))/(T*Nneuron);
                
            disp(['Relative MSE = ',num2str(MSEnorm(nt, na, nc))])
            disp(['Average firing rate = ',num2str(rate(nt,na,nc))])

        end
        save(savestring)
    end
end

%% Plot

% corvecplot = corvec/Nneuron;

Figure1=figure(1);clf;
set(Figure1,'defaulttextinterpreter','latex');
subplot(3,2,1) 
MSEnorm = squeeze(nanmean(MSE./repmat(MSE0,1,Na,Nc), 1));
% MSEnorm = squeeze(MSEnorm(1,:,:));
h = pcolor_fleur(ampvec_noise_rel,1:Nc,MSEnorm');
set(h, 'EdgeColor','none');
set(gca, 'YTickLabel',corvec)
% shading interp
colormap(jet); 
title('$\overline{MSE}$')
ylabel('# noise copies (100 neurons)')
xlabel('relative amplitude noise')
% ylim([0 1])
c = colorbar('Ticks',0:0.05:0.2);
set(c, 'TickLabels',{'0','.05','.1','.15','>.2'})
%,'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})
caxis([0 0.2])
% c.Label.String = 'MSE';


subplot(3,2,2)
rateav = squeeze(nanmean(rate, 1));
% rateav = squeeze(rate(1,:,:));
h = pcolor_fleur(ampvec_noise_rel,1:Nc,rateav');
set(h, 'EdgeColor','none');
set(gca, 'YTickLabel',corvec)
% shading interp
colormap(jet); 
title('Average frequency / neuron')
ylabel('# noise copies (100 neurons)')
xlabel('relative amplitude noise')
cbarvec = 0:10:100;
c = colorbar('Ticks',cbarvec);
caxis([cbarvec(1) cbarvec(end)])
c.Label.String = 'Activity A (Hz)';
% ylim([0 1])

saveas(gcf, savestring)