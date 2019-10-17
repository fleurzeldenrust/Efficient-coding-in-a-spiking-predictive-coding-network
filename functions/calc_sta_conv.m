function sta = calc_sta_conv(inputvector, spiketrain, str)
% Calculate STA with the help of convolution
%
% Input:
% * inputvector
% * spiketrain (zeros and ones of same length as input or spiketimes)
% * options struct:
%       * dt
%       * twindow (time window of STA)
%       * type of spiketrain (spiketimes or spiketrain)
%       * ev. noplot (does not plot if 1)
% Example: str = struct('dt',0.2,'twindow', [-300 200],
% 'spiketype','spiketimes')

% Output sta of size twindow(1):dt:twindow(2)

% varargin = 1) type of spike input (optional) 2) window 
    if nargin == 2
        if length(inputvector) == length(spiketrain)
            disp('Not given spiketrain or spiketimes; assume spiketrain')
            sta = conv(fliplr(spiketrain), inputvector )./sum(spiketrain);
        else
            disp('spiketype is not given; assume spiketimes')
            try 
                spiketrain = make_spiketrain(inputvector, spiketrain, str.dt);
                sta = conv(fliplr(spiketrain), inputvector )./sum(spiketrain);
            catch
                disp('no dt given? cannot compute')
                sta = nan;
            end
        end
    elseif nargin == 3
        try str.spiketype;
            if strcmp(str.spiketype,'spiketrain')
                disp('spikes of spiketrain (0 and 1) form')
                sta = conv(fliplr(spiketrain), inputvector )./sum(spiketrain);
            elseif strcmp(str.spiketype,'spiketimes')
                 disp('spikes of spiketimes form')
                 try
                    spiketrain = make_spiketrain(inputvector, spiketrain, str.dt);
                    sta = conv(fliplr(spiketrain), inputvector )./sum(spiketrain);
                 catch
                     disp('Convertion spiketimes into spiketrain failed')
                 end
            else
                disp('choose spiketype as spiketrain of spiketimes')
            end
        catch
            disp('Not given spiketrain or spiketimes')
            if length(inputvector) == length(spiketrain)
                disp('assume spiketrain')
                sta = conv(fliplr(spiketrain), inputvector )./sum(spiketrain);
            else
                disp('assume spiketimes')
                try 
                    spiketrain = make_spiketrain(inputvector, spiketrain, str.dt);
                    sta = conv(fliplr(spiketrain), inputvector )./sum(spiketrain);
                catch
                    disp('no dt given? cannot compute')
                    sta = nan;
                end
            end
        end
        
        try str.twindow;
            try str.dt;
                twindow = str.twindow;
                tgmin = twindow(1);
                tgmax = twindow(2);
                dt = str.dt;
                lga = tgmin/dt;
                lge = tgmax/dt;
                sta = sta(ceil(length(sta)/2)+lga:ceil(length(sta)/2)+lge);
            catch
                disp('dt not given: give back full sta')
            end
        catch 
            disp('twindow not given: give back full sta')
        end
        
        

        try 
            if str.noplot == 1
                % no plot
                disp('no plot')
                plt = 0;
            elseif str.noplot == 0
                plt = 1;
            end
        catch
            plt = 1;
        end

        try
            if plt == 1
                figure
                plot(tgmin:dt:tgmax, sta)
                hold all
                plot([0 0],[min(sta) max(sta)], 'k')
                xlabel('time relative to spike')
                ylabel('average amplitude')
            end
        catch
            disp('no plot possible')
        end
    end
    
    
end

function st = make_spiketrain(inputvector, spiketimes, dt)
    st = zeros(size(inputvector));
    st(round(spiketimes/dt)) = 1;
end
    
