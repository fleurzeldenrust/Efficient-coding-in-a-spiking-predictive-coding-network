function [disturbvec, dphivec, example, Orealp] = lew_make_prc(tfilt, kernel, nu, delta, aamp, taua, cinput,  plotyn, options )
% Makes a phase-response curve (PRC) of a neuron with predictive field
% kernel (with time relative to the spike tfilt and delay time delta in
% ms). NB absolute cost is used!
% INPUT: Aamp and taua determina spike-frequency adaptation, nu the absolute
% spike cost. The current-input cinput is a struct with fields:
% * amp (background current)
% * Npulse (determines how many pulses to use),  
% * nsp determines which spike (in response to normal pulse) to use
% * pulseamp (amplitude of the pulse)
% * pulselength (length of pulse in ms)
% plotyn determines whether to make a plot
% OUTPUT: dphivec = vector with phase advance or delay (delta period /
% length period), disturbvec = vector with timing of pulses (time since
% last spike / length period)


if nargin == 8
    disp('Allowing multiple spikes')
    disp('Normalize by input kernel')
    options.multspikes = 1; % allow for multiple spikes at each dt (uncoupled neurons)
    options.normalizegin = 1; % normalize filters by input filters
elseif nargin < 8
    error('input missing')
else
    if ~isfield(options, 'multspikes')
        disp('Allowing multiple spikes')
        options.multspikes = 1;
    end
    if ~isfield(options, 'normalizegin')
        disp('Normalize by input kernel')
        options.multspikes = 1;
    end
end

T = 500;
dt = tfilt(2)-tfilt(1);
Ntime = T/dt;
time = (1:Ntime)*dt;
Ndelta = round(delta/dt);

pulsen = round(cinput.pulselength/dt);

%% generate lateral and input filters and threshold from representing filters
[tg, g, gin, gout, Th ] = lew_generate_filters( tfilt, kernel, [], delta, options);


%% Make signal
si = cinput.amp*ones(1, Ntime);
si(1:round(20/dt))=0;

%% Run first
[xest, O, ~, ~] = lew_run_abscost(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua, options);
Oreal = [zeros(1, Ndelta) O];
Oreal = Oreal(1, 1:Ntime);

if plotyn
    fighandle = figure;
    subplot(3,1,1)
    plot(tfilt, kernel)
    hold all
    plot([delta delta], [min(kernel) max(kernel)])

    subplot(3,1,2)
    hold all
    plot(time, si)
    plot(time,xest)
    plot(time, Oreal, '.')
else
    fighandle = [];
end

example.time = time;
example.si = si;
example.xest = xest;
example.Oreal = Oreal;


%% Make starttimes pulses
nspike    = find(Oreal==1);
nsecspike = nspike(cinput.nsp);
tsecspike = nsecspike*dt;
nthispike = nspike(cinput.nsp+1);
tthispike = nthispike*dt;
Tun = tthispike-tsecspike;

pulsestartlength = Tun/((cinput.Npulse-1));
pulsestartn = round(pulsestartlength/dt);
pulsestartvec = (0:cinput.Npulse-1)*pulsestartn;

%% Disturb the system
Op = zeros(cinput.Npulse+1,Ntime);
Orealp = zeros(cinput.Npulse+1,Ntime);
xestv = zeros(cinput.Npulse+1,Ntime);
nthivec = zeros(cinput.Npulse, 1);

Op(1,:) = O;
Orealp(1,:) = Oreal;
xestv(1,:) = xest;

for np = 1:cinput.Npulse
    sitemp = si;
%     sitemp(nsecspike+(np-1)*2*pulsestartn+1:nsecspike+(np-1)*2*pulsestartn+pulsen)=asig+pulseamp;
    sitemp(nsecspike+pulsestartvec(np)+1:nsecspike+pulsestartvec(np)+pulsen)=si(nsecspike+pulsestartvec(np)+1:nsecspike+pulsestartvec(np)+pulsen)+cinput.pulseamp;
    [xestv(np+1,:), Op(np+1,:), ~, ~] = lew_run_abscost(dt, sitemp, g, tg, gin, gout, Th, nu, delta, aamp, taua, options);
   
    Orealptemp = [zeros(1, Ndelta) Op(np+1,:)];
    Orealptemp = Orealptemp(1, 1:Ntime);
    Orealp(np+1,:) = Orealptemp;
    nspike = find(squeeze(Orealp(np+1,:)==1));
    nthivec(np) = nspike(cinput.nsp+1);
    
    if plotyn
        subplot(3,1,2)
        hold all
        plot(time, sitemp)
        plot(time, xestv(np+1,:))
        plot(time, (1+np/cinput.Npulse)*Orealp(np+1,:), '.')
    end
    
end

npervec = nthivec-nsecspike;
Tpervec = npervec*dt;

dphivec = (Tun-Tpervec)/Tun;
disturbvec = (pulsestartvec-pulsestartvec(1))/(Tun/dt);

if plotyn
    subplot(3,1,3)
    plot(disturbvec, dphivec)
    xlim([0 1])
    xlabel('Phase start pulse')
    ylabel('Phase shift')
    title('PRC')
    pause(5)
end


end

