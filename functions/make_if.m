function [fstart, fend, tdelay, rheo, nspike ] = lew_make_if( tfilt, kernel, nu, delta, aamp, taua, ampvec,  plotyn, options )
% Makes a inputcurrent - outpurfrequency curve of a neuron with predictive field
% kernel (with time relative to the spike tfilt and delay time delta in
% ms). NB absolute cost is used!
% INPUT: Aamp and taua determina spike-frequency adaptation, nu the absolute
% spike cost. The parameter ampvec gives the steps of input current.
% plotyn determines whether to make a plot
% OUTPUT: 
% * fstart = vec with frequency first two spikes
% * fend = vec with frequency last two spikes
% * tdelay = vec with delay to first spike
% * rheo = (scalar) current with first spike (longest delay)

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

T = 2000;
tstart = 200;
T = T+tstart;
dt = tfilt(2)-tfilt(1);
Ntime = T/dt;
time = (1:Ntime)*dt;

%% generate lateral and input filters and threshold from representing filters
[tg, g, gin, gout, Th ] = lew_generate_filters( tfilt, kernel, [], delta, options);


%% Make signal
tdelay = NaN*ones(1,length(ampvec));
fstart = zeros(1,length(ampvec));
fend = zeros(1,length(ampvec));
nspike = zeros(1,length(ampvec));
ii=0;
for a = ampvec
%     disp(['amplitude step = ',num2str(a)])
    ii=ii+1;
    si = zeros(1,length(time));
    si(round(tstart/dt)+1:end) = a;
    if isfield(options, 'noiseamp')
        si = si+options.noiseamp*rand(size(si));
    end
    
    [~, O, ~, ~] = lew_run_abscost(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, taua, options);
    
    spiketimesi = time(find(O==1));
    spiketimes = spiketimesi+delta;
    if ~isempty(spiketimes)
        tdelay(ii) = spiketimes(1)-tstart;
        nspike(ii) = length(spiketimes);
        if length(spiketimes)>1
            fstart(ii) = 1000/(spiketimes(2)-spiketimes(1));
            
            % NB do not take last ones due to boundary effects
            nt = find(abs(spiketimes-1700)==min(abs(spiketimes-1700)));
            
            if length(nt)>1
                nt = nt(end);
            end
            if (length(spiketimes) == 1 || isempty(spiketimes))
                fend(ii) = NaN;
            elseif nt>6
                fend(ii) = 5*1000/(spiketimes(nt)-spiketimes(nt-6));
            else                
                fend(ii) = (nt-1)*1000/(spiketimes(nt)-spiketimes(1));                
            end
            
        end
    end

end


try
    ndelaymax = find(tdelay==max(tdelay));
    rheo = ampvec(ndelaymax(1));
catch
    rheo = NaN;
end

if plotyn
    figure
    subplot(1,2,1)
    plot(tdelay, ampvec)
    xlabel('delay to first spike (ms)')
    ylabel('amplitude input')
    title('Input versus delay to first spike')
   
    
    subplot(1,2,2)
    hold all
    plot(ampvec, fstart)
    plot(ampvec, fend)
    ylabel('output frequency (Hz)')
    xlabel('amplitude input')
    title('Input-Output curve')
    legend('first two spikes','last spikes')
    


end

