close all
clear all
f =  filesep;
addpath(['.' f 'functions'])

T = 2000;           % total time
Tf = 50;           % time kernels
dt = 1;             % dt
Ntime = T/dt;
Nneuron = 40;       % # neurons
Nbf = 15;            % # basis functions for kernels
tfilt = 0:dt:Tf;    % time of filter
nu = 1;             % relative spike cost (absolute = nu*Th)
aamp = 1;           % relative adaptive spike cost 
taua = 100;         % time constant adaptation

delta = 10;         % delay time

% NB delta = 0 and exponential kernels (1 basis function) corresponds to the 'classical'
% network

% NB Note that here spikes are fired in the past! To get the 'realistic'
% output, shift everything that comes out of run_model delta into the future

%% generate representing filters 

% Make basis funcions

% Gamma functions
% basisf = make_basisfunction(tfilt, 15, 'gamma',[]);

% Exponentials
% params.tau = 2*(1:15);
% params.tau = 10;
% basisf = make_basisfunction(tfilt, Nbf, 'exponential', params);

% Gaussians
params.mu = (5:20);
params.sigma = params.mu/2;
basisf = make_basisfunction(tfilt, 15, 'gaussian', params);

% Make kernels from basis functions
Gamma=randn(Nneuron,Nbf);
kernel=Gamma*basisf;


%% generate signal
si=randn(1,Ntime);
si=conv2(si,exp(-[1:300]/10)/sum(exp(-[1:300]/10)),'same');
si=conv2(si,exp(-[1:300]/10)/sum(exp(-[1:300]/10)),'same');
si = si*10;


%% generate lateral and input filters and threshold from representing filters
[tg, g, gin, gout, Th ] = generate_filters( tfilt, kernel, 1:Nneuron, delta);

%% run network
[xest, O, V, Thvec] = run_model(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua);

%% Plot
figure
subplot(3,1,1)
for nn = 1:Nneuron
    hold all
    plot(tfilt, kernel(nn,:))
end
plot([delta delta],[min(min(kernel))-0.5 max(max(kernel))+0.5], 'k', 'LineWidth',2)
xlabel('time')
ylabel('kernel amplitude')
title('Kernels')
grid on

subplot(3,1,2)
plot(si)
hold all
plot(xest)
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
    plot(nn*st, '.')
end
ylim([0.2, Nneuron+0.2])
xlim([0 T])
xlabel('time')
ylabel('neuron')
title('Ideal network activity (real activity shifted by \Delta)')
grid on
