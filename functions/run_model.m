function [estimate, spikesideal, V, Thvec] = run_model(dt, si, g, tg, gin, gout, Th, nu, delta, aamp, atau)
% Runs the filter version of the predispike model with
% * input signal si in time
% * representing filters g (Nker x length tg)
% * dt (should be the same for filters and signal!)
% * input filters gin (Nker x length tg)
% * output (lateral) filters gout (Nker x Nker x length tg); NB gout(ii,jj): filter sent by jj, received by ii 
% * Thresholds Th (vector size Nker) 
% * spike cost nu (scalar)
% * delay time between spike and decision delta
% * adaptation (increase in threshold) after each spike, with time constant
% atau and relative size aamp (relative to size threshold).

% Gives as output
% * estimate of the network of the signal
% * ideal spike train (= real spike train shifted delta backward in time)
% * membrane potential V of each neuron
% * threshold for each neuron over time: Thvec= Th*(1+nu+aampexp(-t/atau)) 

% NB allows only one spike in the network for every time step: it picks the
% one with maximal V-Th; if there are two it picks the one with the lowest
% index


[Nker, lg] = size(g);
Ntime = length(si);

%% Make adaptation kernel
ta = 0:dt:5*atau;
la = length(ta);
thker = aamp.*exp(-ta/atau);
% NB do not normalize: always starts at 1


%% Run network
V=zeros(Nker,Ntime);                 % 'membrane potential;': sum of effects of all filters on 1 neuron
Thvec=repmat((Th.*(1+nu))',1,Ntime); % Threshold matrix; includes spike cost and adaptation
siconv=zeros(Nker,Ntime);            % input convolved with input filters
oiconv=zeros(Nker,Ntime);            % output convolved with output and lateral filters
O=zeros(Nker,Ntime);                 % spike trains

% input convolved with input kernel
for i=1:Nker
   siconv(i,:)=convolve_kernel_acausal(si, gin(i,:), 0, lg-1)*dt;
end

for tn = round(delta/dt)+1:Ntime

    % subtract output kernels and decide spike
    V(:,tn) = siconv(:,tn)-oiconv(:,tn);
    if max(V(:,tn)-Thvec(:,tn))>0
        % a spike is fired in neuron nn with maximal V-Th at time =
        % tn-delta, so in the past!
        nvec = find(V(:,tn)-Thvec(:,tn) == max(V(:,tn)-Thvec(:,tn)));
        if length(nvec)<1
            keyboard
            error('max V does not exist')
        elseif length(nvec)==1
            nn = nvec;
        else
            disp('spike in multiple neurons; take one with lowest index')
            nn = nvec(1);
        end
        O(nn,tn-round(delta/dt))=1;
        % update oiconv 
        ll = length(tn+1:Ntime);
        for mm = 1:Nker
            if tn+lg<Ntime
                % NB gout(i,j) sent by j, received by i 
                oiconv(mm, tn+1:tn+lg) = oiconv(mm, tn+1:tn+lg) + squeeze(gout(mm, nn ,:))';
            else
                oiconv(mm, tn+1:end) = oiconv(mm, tn+1:end) + squeeze(gout(mm, nn ,1:ll))';
            end
        end
        % update Thvec
        if tn + la < Ntime
            Thvec(nn,tn+1:tn+la) = Thvec(nn,tn+1:tn+la) + Th(nn).*thker;
        else
            Thvec(nn,tn+1:end) = Thvec(nn, tn+1:end) + Th(nn).*thker(1:ll);
        end
        
        

    end
    
end
spikesideal = O;


%% Make estimate
lgc = length(find(tg>0)); % length causal part
lga = length(find(tg<0)); % lenght acausal part

est=zeros(1, Ntime);
for i=1:Nker
    est = est + convolve_kernel_acausal( O(i,:), g(i,:), lga, lgc);
    % is the same as est(nspike-lga:nspike+lgc) = est(nspike-lga:j+lgc) + g(i,:);
end
estimate = est;


end