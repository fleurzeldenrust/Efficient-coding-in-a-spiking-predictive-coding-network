close all
clear all

dbstop if error

addpath('./functions/')

factor = 1;
Tf = 50*factor;           % time kernels
dt = .01*factor;             % dt
tfilt = 0:dt:Tf;    % time of filter
Nbf = 8;  
Nneuron = 1;

% Neuron
taua = 60;         % time constant adaptation
mu = 1.5;          % strength adaptation
nu = 1.5;            % spike cost

delta = 7.5*factor;         % delay time
Ndelta = round(delta/dt);

options.multspikes = 0;
options.normalizegin = 1;

% if
ampvec = [0:0.05:0.7, 0.71:0.01:0.89, 0.9:0.05:1.9, 1.91:0.01:2.19, 2.2:0.05:2.5];
% options.noiseamp = 1;
% Ntrialif = 10;
options.noiseamp = 0;
Ntrialif = 1;

% prc
cinput.Npulse = 100;
factorprc = 1/5;
cinput.pulselength = .5*factorprc;



% sta
stawindow = [-50 20];
Tsigfilt = 50;         % time filter noisesignal
tsigfilt = 0:dt:Tsigfilt;
tausig = 1;            % time constant signal
cinput.filt = exp(-tsigfilt/tausig)/sum(exp(-tsigfilt/tausig));
% cinput.T = 50000;
cinput.T = 100000;

%% Make representing filters
basisf = make_basisfunction(tfilt, Nbf, 'gamma',[]);
for type = [1, 2]
    disp(['Type ' num2str(type)])
% for type = 2
    if type == 1
        cinput.amp = 0.9;

        % Pulses for PRC
        cinput.nsp = 1;

        kernel = basisf(3,:);
        color = 'b';
    elseif type == 2
        cinput.amp = 2.2;



        % Pulses for PRC
        cinput.nsp = 1;

        kernel = basisf(3,:).*(0.2.*ones(size(tfilt))-0.8*sin(tfilt*.6/factor));
        color = 'r';
    end

    [tg, g, gin, gout, Th ] = generate_filters( tfilt, kernel, [], delta, options);




    %% Make IF curve with adaptation

    for ntf = 1:Ntrialif
        disp(['IF trial ' num2str(ntf)])
        [fstart(ntf,:), fend(ntf,:), tdelay(ntf,:), rheo(ntf), nspike(ntf, :) ] = make_if( tfilt, g, nu, delta, mu, taua, ampvec,  0, options );
    end
    
%     %% Make IF curve no adaptation
% 
%     for ntf = 1:Ntrialif
%         disp(['IF trial ' num2str(ntf)])
%         [fstart_nadap(ntf,:), fend_nadap(ntf,:), tdelay_nadap(ntf,:), rheo_nadap(ntf), nspike_nadap(ntf, :) ] = make_if( tfilt, g, nu, delta, 0, taua, ampvec,  0, options );
%     end



    %% Make PRC
    cinput.pulseamp = cinput.amp/(3*factorprc);
    [disturbvec, dphivec, example, disturbed_spikes ] = make_prc(tfilt, g, nu, delta, mu, taua, cinput,  1, options );
    
    cinput.nsp = cinput.nsp+1;
    [disturbvec2, dphivec2, example2, disturbed_spikes2 ] = make_prc(tfilt, g, nu, delta, mu, taua, cinput,  1, options );
        
    %% Make STA
    cinput.amp = cinput.amp/2;
    [tsta, sta, MSE, MSE0, rate ] = make_sta( tfilt, kernel, nu, delta, mu, taua, cinput, stawindow, 0, options );


    %% Plot
    figure(1)

    % Filters
    subplot(4,2,1)
    hold all
    plot(tfilt, g, 'Color',color, 'LineWidth',2)
    xlabel('time (ms)')
    ylabel('amplitude filter (a.u.)')
    title('Representing filter g')    
    xlim([0 30])
    grid on
    box on
    
    subplot(4,2,3)
    hold all
    plot(tfilt, squeeze(gout(1,1,:)),  'Color',color, 'LineWidth',2)
    xlabel('time (ms)')
    ylabel('amplitude filter (a.u.)')
    title('Self-filter g^{lat}')    
    xlim([0 30])
    grid on
    box on
    
    %% Example response
%     subplot(2,3,2)
%     hold all
%     if type==1
%         plot(example.time, example.si, 'k', 'LineWidth',2)
%         value1 = 2*example.si(end);
%     end
%     [Ndisturb,~] = size(disturbed_spikes);
%     for np=1:Ndisturb
%         plot(example.time(disturbed_spikes(np,:)==1), (type-0.5+0.1+np/Ndisturb)*value1*ones(1,sum(disturbed_spikes(np,:))),'.k')
%     end
%     if type==2
%         for np=1:Ndisturb
%             plot(example2.time(disturbed_spikes2(np,:)==1), (type+0.5+0.1+np/Ndisturb)*value1*ones(1,sum(disturbed_spikes2(np,:))),'.k')
%         end
%     end
        
%     yticks([value1, 2*value1, 3*value1])
%     yticklabels({'stimulus','type 1','type 2'})
%     xlabel('time (ms)')
%     title('Example respose to step-and-hold input')
%     box on
    
    %% PRC
    subplot(2,2,2)
    hold all
    plot(disturbvec, dphivec, 'Color',color, 'LineWidth',2)
    plot(disturbvec2, dphivec2, 'Color',color, 'LineStyle','--', 'LineWidth',2)
    xlim([0 1])
    legend('type 1 - first two spikes', 'type 1 - spike 2 and 3','type 2 - first two spikes', 'type 2 - spike 2 and 3')
    xlabel('Phase start pulse')
    ylabel('Phase shift')
    title('Phase Response Curve (PRC)')
    grid on
    box on

    %% Nspike
%     subplot(2,3,4)
%     hold all
%     if Ntrialif>1
%         plot(ampvec, mean(nspike), '.k', 'LineWidth',2)
%     else
%         plot(ampvec,nspike, '.k', 'LineWidth',2)
%     end
%     xlabel('amplitude input')
%     ylabel('# spikes')
%     title('Number of spikes in 2000 ms')
%     grid on
%     box on
%     
    % IF
    subplot(4,2,5)
    hold all
    if Ntrialif>1
%         plot(ampvec, mean(fstart), '-k', ampvec, mean(fend), '--k', 'LineWidth',2)
        plot(ampvec, mean(fstart), '.', 'Color',color)
    else
%         plot(ampvec, fstart, '-k', ampvec, fend, '--k', 'LineWidth',2)
        plot(ampvec, fstart, '.', 'Color',color)
    end
    xlabel('amplitude input (a.u.)')
    ylabel('frequency (Hz)')
    title('Frequency first two spikes')
    ylim([0 500])
    grid on
    box on
    
    subplot(4,2,7)
    hold all
    if Ntrialif>1
        plot( ampvec, mean(nspike/2), '.', 'Color',color)
    else
        plot(ampvec, nspike/2, '.k', 'Color',color)
    end
    xlabel('amplitude input (a.u.)')
    ylabel('frequency (Hz)')
    title('Average frequency')
    ylim([0 70])
    grid on
    box on
 
    % STA
    subplot(2,2,4)
    hold all
    plot(tsta, sta, 'Color',color, 'LineWidth', 2)
    xlabel('time relative to spike (ms)')
    ylabel('stimulus amplitude (a.u.)')
    title('STA')
    box on 
    grid on


end

figure(1)
subplot(4,2,1)
hold all
plot([delta, delta],[-2,2], '--k')
legend('type 1', 'type 2', '\Delta')
ylim([-0.6, 1.5])
grid on
box on

subplot(2,2,4)
hold all
plot([0, 0],[-1,2.5], '--k')
xlim([-25 10])
grid on
box on