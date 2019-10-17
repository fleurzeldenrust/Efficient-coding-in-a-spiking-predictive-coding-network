close all
clear all

f =  filesep;
addpath(['.' f 'functions'])

T = 2000;           % total time
dt = .2;             % dt
Ntime = T/dt;
Nneuron = 100;       % # neurons

options.normalizegin = 1;
Tf = 50;             % time kernels
Nbf = 8;             % # basis functions for kernels
tfilt = 0:dt:Tf;     % time of filter
delta = 7.5;         % delay time

cost = 'abs';
nu = 0.5;             % relative spike cost (absolute = nu*Th)
aamp = 0.5;           % relative adaptive spike cost 

options.multspikes = 0;
options.multspikerand = 1;

net = 'homogeneous_t1';  % choose 'homogeneous_t1' or 'heterogeneous' or 't12' or
% 'heterogeneous_matchedt12
seed = 15;


taua = 60;         % time constant adaptation

Tsigfilt = 50;         % time filter noisesignal
tsigfilt = 0:dt:Tsigfilt;
tausig = 20;            % time constant signal
asig = 1;               % amplitude signal



% NB delta = 0 and exponential kernels (1 basis function) corresponds to the 'classical'
% network

% NB Note that here spikes are fired in the past! To get the 'realistic'
% output, shift everything that comes out of lew_run delta into the future




%% generate signal
% noisy
filtsignal = exp(-tsigfilt/tausig)/sum(exp(-tsigfilt/tausig));
% filtsignal = hamming(round(length(tfilt)/8));
% filtsignal = hamming(250);

rng(3)
si=randn(1,Ntime);
si=conv(si,filtsignal,'same');
si = conv(si, fliplr(filtsignal), 'same');
si = si*asig/std(si);

% steps
% si = zeros(1, Ntime);
% si(100/dt:400/dt) = -1;
% si(400/dt:800/dt) = 2;
% si(800/dt:1200/dt) = -2;
% si(1200/dt:1600/dt) = 4;

% oscillator
% t = (1:Ntime)*dt;
% si = sin(2*pi*50*t/1000)+cos(2*pi*30*t/1000)+cos(2*pi*80*t/1000);

% scale
% si = si/3;

% add noise
% si = si + randn(size(si))*std(si)/10;

%% Make kernels network
kernel = make_kernels_network(tfilt, Nbf, Nneuron, net, seed);
[tg, g, gin, gout, Th ] = generate_filters( tfilt, kernel, 1:Nneuron, delta, options);


%% run network

if strcmp(cost, 'rel')
    [xest, O, V, Thvec] = run_relcost(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua, options);
elseif strcmp(cost, 'abs')
    [xest, O, V, Thvec] = run_abscost(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua, options);
end

MSE = calc_MSE(si, xest);
MSE_nospikes = calc_MSE(si, zeros(size(si)));

%% Plot
time = (1:Ntime)*dt;
figure
subplot(3,1,1)
for nn = 1:Nneuron
    hold all
    plot(tfilt, g(nn,:))
end
plot([delta delta],[min(min(g))-0.5 max(max(g))+0.5], 'k', 'LineWidth',2)
xlabel('time')
ylabel('kernel amplitude')
title('Kernels')
grid on

subplot(3,1,2)
plot(time,si)
hold all
plot(time,xest)
xlim([0 T])
xlabel('time')
ylabel('amplitude')
title('Signal/estimate')
legend('signal', 'estimate')
grid on

subplot(3,1,3)
hold all
for nn = 1:Nneuron
    st = O(nn,:);
    plot(time, nn*st, '.')
end
ylim([0.2, Nneuron+0.2])
xlim([0 T])
xlabel('time')
ylabel('neuron')
title(['Ideal network activity (real activity shifted by \Delta)'])
grid on

disp(['total : ' ,num2str(sum(sum(O))), ' spikes, this makes a rate of ', num2str(1000*sum(sum(O))/(Nneuron*T))])
disp(['MSE = ',num2str(MSE), '; relative MSE = ',num2str(MSE/MSE_nospikes)])
disp(['total cost = ', num2str((1000*sum(sum(O))/T)*(MSE/MSE_nospikes))])

