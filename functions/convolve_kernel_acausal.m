function s = convolve_kernel_acausal( signal, filter, la, lc)
% Convolve signal with kernel with causal and acausal part, give back on same time
% scale.

% Function input:
%   signal(vector) = spike train or input signal
%   filter (vector) = filter
%   la = length acausal part filter (tbegin(negative):0-dt)
%   lc = length causal part filter (dt:tmax/dt)
%   NB lc+la+1 = length(filter) for filters with both causal and acausal
%   parts

% Function output: convolved signal s

cf = 0;
ca = 0;

if ~(la+lc+1==length(filter))
    if la == 0 && length(filter)==lc
        disp('purely causal filter')
        cf = 1;
    elseif lc == 0 && length(filter)==la
        disp('purely acausal filter')
        ca = 1;
    else
        error('give correct lengths of causal and acausal part filter')
    end
end

if cf
    s = conv(signal, [0 filter]);
    s = s(1:length(signal));
elseif ca
    s = conv([signal zeros(1,la)], [filter 0], 'valid');
else
    s = conv([zeros(1,lc) signal zeros(1,la)], filter, 'valid');
end



end

