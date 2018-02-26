function [tg, g, gin, gout, Th ] = generate_filters( tfilt, filt, plotvec, delta)
% function that gives optimal input, lateral and output kernels and threshold, as a
% function of :
% * tfilt: filter time (vector), note that spike is at 0
% * filt: matrix of the representation filters (Nker x time)
% * plotvec: if 0, no plotting, else: vector of which filters to plot
% * delta (>=0): delay between spike and evaluation, should be multiple of dt

% Gives as output
% * tg that always runs from min(tfilt(1), 0) to max(tfilt(end),delta, 0),
% so t=0 is always included
% * representing filters g (Nker x length tg)
% * input filters gin (Nker x length tg)
% * output filters gout (Nker x Nker x length tg)
% * thresholds Th (Nker)

% NB gout(i,j) sent by j, received by i 
% Note that the output filters work from tspike+dt + tspike+dt+filtertime
% Note that if delta>tfilt(end), the filter will be padded with zeros, if
% delta<tfilt(end), the filter will be cut off

%% Initialize
dt = tfilt(2)-tfilt(1);
tgmin = tfilt(1);
tgmax = tfilt(end);
tg = tfilt;


g = filt;
[Nker, lg] = size(g);

if ~(lg == length(tg))
    error('filter time does not correspond to filter length')
end

%% Set tgmax = delta for large delta (>tgmax)
if nargin == 4
    if delta > tgmax
        % add zeros
        if tgmax>=0
            gtemp = zeros(Nker, lg+round((delta-tgmax)/dt));
            gtemp(:, 1:lg)=g;
            g = gtemp;
            clear gtemp
            tgmax = delta;
            tg = (round(tgmin/dt):round(tgmax/dt))*dt;
        else
            'purely acausal filter'
            gtemp = zeros(Nker,round((abs(tgmin)+delta)/dt)+1);
            gtemp(:, 1:lg)=g;
            g = gtemp;
            clear gtemp
            tgmax = delta;
            tg = (round(tgmin/dt):round(tgmax/dt))*dt;
        end
        [Nker, lg] = size(g);
    end
elseif nargin == 3
    delta = tgmax;
    if tgmax<0
        'purely acausal filter'
        % add zeros
        delta = 0;
        gtemp = zeros(Nker, round(abs(tgmin)/dt)+1);
        gtemp(:, 1:lg)=g;
        g = gtemp;
        clear gtemp
        tgmax = delta;
        tg = (round(tgmin/dt):round(tgmax/dt))*dt;
    end
    [Nker, lg] = size(g);
end

if ~(lg == length(tg))
    error('filter time does not correspond to filter length')
end

if abs(tg(end)-tgmax)>0.5*dt
    keyboard
    error('tg(end) is not equal to tgmax; something went wrong')
end

if abs(tg(1)-tgmin)>0.5*dt
    error('tg(1) is not equal to tgmin; something went wrong')
end

lgc = length(find(tg>0)); % length causal part
lga = length(find(tg<0)); % lenght acausal part

cf = 0;
if ((lgc+lga+1)==lg)
    disp('filter with causal and acausal part')
else
    if lga == 0
        if lgc == lg
            disp('purely causal filters')
            cf=1;
            % add zeros from t=0    
            tg = (0:round(tgmax/dt))*dt;
            tgmin = 0;
            lgtemp = length(tg);
            gtemp = zeros(Nker, lgtemp);
            gtemp(:, end-lg+1:end) = g;
            g = gtemp;
            lg = lgtemp;
            clear gtemp lgtemp
        elseif lgc + 1 == lg
            disp('causal filter including 0')
        else
            error('length causal is not length filter')
        end
    elseif lgc == 0
        if lga+1 == lg
            disp('acausal filters including 0')
        else
            error('length acausal part is not length filter')
        end
    else
        keyboard
        error('length causal+acausal part+1 is not length filter')
    end
end

if ~(lg == length(tg))
    error('filter time does not correspond to filter length (2)')
end

lgc = length(find(tg>0)); % length causal part
lga = length(find(tg<0)); % lenght acausal part

if ((lgc+lga+1)==lg)
    disp('filter with causal and acausal part')
else
    error('filter does not include 0, something went wrong')
end

%% Make input filters
% gin is the same as g, but just with other time labels as long as
% tgmax<delta
% if delta<tgmax: input filter shorter, pad with zeros
    
gtemp = g(:,1:lga+1+round(delta/dt));
tgtemp = tg(1:lga+1+round(delta/dt));
tgin = fliplr(delta-tgtemp);
gin = fliplr(gtemp);
clear gtemp tgtemp

[Nkerin, lgin] = size(gin);

if ~(Nkerin == Nker)
    error('not the same amount of input and representing kernels, something went wrong')
end

if ~(lgin==lg)
    % add zeros
    gtemp = zeros(Nker, lg);
    gtemp(:,1:lgin) = gin;
    gin = gtemp;
    tgin = ((0:lg-1))*dt;
    clear gtemp
end

%% Make output filters and threshold
Th = zeros(1,Nker);
gout = zeros(Nker, Nker, lg);
for ii = 1:Nker
    Th(ii) = 0.5*sum(gin(ii,:).^2)*dt;
    for jj = 1:Nker
        temp = convolve_kernel_acausal([zeros(1,lg) g(jj,:) zeros(1,lg)], gin(ii, :), 0, lg-1);
        % NB t'spike' = lg+lga+1; evaluation time = lg+lga+1+delta/dt; so first
        % time it can influence other neurons is at lg+lga+1+delta/dt+1;
        % NB gout(ii,jj): filter sent by jj, received by ii 
        if ii==jj
            if abs(2*Th(ii)-dt*temp(lg+lga+1+round(delta/dt)))>2*Th(ii)/10000;
                error('starting value output filter is not threshold; something went wrong')
            end
        end
        gout(ii, jj, :) = dt*temp(lg+lga+2+round(delta/dt):2*lg+lga+2+round(delta/dt)-1); % earliest evaluation time = tspike+dt; 'spike' was at lg+lga+1; 
    end
end


%% Plot
if length(plotvec)>0 
    if min(plotvec)>0
        figure
        for plotnr = plotvec
            subplot(2,2,1)
            hold all
            plot(tg, g(plotnr, :))  
            subplot(2,2,2)
            hold all
            plot(tgin, gin(plotnr, :))  
            subplot(2,2,3)
            hold all
            plot((1:lg)*dt, -squeeze(gout(plotnr, plotnr, :))/dt)
            subplot(2,2,4)
            hold all
            for plotnr2 = plotvec
                if plotnr == plotnr2
                else
                    plot((1:lg)*dt, -squeeze(gout(plotnr, plotnr2, :))/dt)
                end
            end
        end
        subplot(2,2,1)
        hold all
        plot([0 0],[-1 1])
        plot([delta delta],[min(min(g))-1 max(max(g))+1], 'k', 'LineWidth',2)
        xlim([tgmin tgmax])
        title('Representing filters g(t)')
        xlabel('time (ms)')
        grid on
        subplot(2,2,2)
        hold all
        plot([delta delta],[min(min(gin))-1 max(max(gin))+1], 'k', 'LineWidth',2)
        plot([0 0],[-1 1])
        xlim([min(tgin) max(tgin)])
        title('Input filters g^{in}(t)')
        xlabel('time (ms)')
        grid on
        subplot(2,2,3)
        xlim([min(tgin) max(tgin)+dt])
        xlabel('t-T+\Delta (ms)')
        title('Output filter g^{out}(t)')
        grid on
        subplot(2,2,4)
        xlim([min(tgin) max(tgin)+dt])
        title('Lateral filter g^{lat}(t)')
        xlabel('t-T+\Delta (ms)')
        grid on
    end
end

