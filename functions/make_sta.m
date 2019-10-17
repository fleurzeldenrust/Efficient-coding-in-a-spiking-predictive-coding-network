function [tsta, sta, MSE, MSE0, rate ] = lew_make_sta( tfilt, kernel, nu, delta, aamp, taua, cinput, stawindow, plotyn, options )
% Makes a spike-triggered average of a neuron with predictive field
% kernel (with time relative to the spike tfilt and delay time delta in
% ms). 
% INPUT: Aamp and taua determina spike-frequency adaptation, nu the absolute
% spike cost. The current-input cinput is a struct with fields:
% * amp (std of current)
% * filt (vec, 2 times filtered with this filter)  
% * T (scalar, length in ms of the used stim)
% plotyn determines whether to make a plot
% OUTPUT: 
% * tsta: time of sta (relative to spike)
% * sta: amplitude of sta
% * MSE: MSE between input and estimate
% * MSE0: MSE between input and 0
% * rate: firing rate (in Hz) of neuron

if nargin == 9
    disp('Allowing multiple spikes')
    disp('Normalize by input kernel')
    options.multspikes = 1; % allow for multiple spikes at each dt (uncoupled neurons)
    options.normalizegin = 1; % normalize filters by input filters
elseif nargin < 9
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

T = cinput.T;           % total time
dt = tfilt(2)-tfilt(1);
Ntime = T/dt;
Ndelta = round(delta/dt);


%% generate lateral and input filters and threshold from representing filters
[tg, g, gin, gout, Th ] = lew_generate_filters( tfilt, kernel, [], delta, options);


%% Make signal
si=randn(1,Ntime);


si=conv(si,cinput.filt,'same');
si=conv(si,fliplr(cinput.filt),'same');
si = si*cinput.amp/std(si);

%% Run
[xest, O, ~, ~] = lew_run_abscost(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua, options);

Oreal = [zeros(1, Ndelta) O];
Oreal = Oreal(1:Ntime);

%% Make STA
str = struct('dt',dt,'twindow', stawindow,'spiketype','spiketrain', 'noplot',1);

sta = calc_sta_conv(si, Oreal, str);

%% Calc MSE
MSE0 = calc_MSE(si, zeros(size(si)));
MSE = calc_MSE(si, xest);
rate = 1000*(sum( O)/T);
tsta = (str.twindow(1):dt:str.twindow(2));

if plotyn
    figure
    hold all
    plot(tsta, sta)
    xlim(str.twindow)
    ylim([min(min(sta)) max(max(sta))])
    xlabel('time relative to spike')
    ylabel('stimulus amplitude')
    title('STA')
    h = get(gca, 'Children');
    set(h, 'LineWidth',2)
    grid on
    box on
end

end

