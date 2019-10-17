function xc = av_corr(ON, neuron1, neuron2)
% Calculate cross-correlogram, normalized to number of spikes and averaged
% over trials.
% NB If you want to normalize by the rates, xc still has to be multiplied
% by the total time in seconds (twice, so without sqrt)
% ON = spikes matrix (Ntrial x Nneuron x Ntime)

if ndims(ON)==3
    [Ntrial, ~, Ntime] = size(ON);
else
    Ntrial = 1;
    [Nneuron, Ntime] = size(ON);
    temp = zeros(1,Nneuron, Ntime);
    temp(1,:,:) = ON;
    ON = temp;
end
L = 2*Ntime-1;

xctemp = zeros(L, 1);
N1=0;
N2=0;
for nt=1:Ntrial
    trace1 = squeeze(ON(nt,neuron1,:));
    trace2 = squeeze(ON(nt,neuron2,:));
    xctemp = xctemp+xcorr(trace1, trace2);
    N1 = N1+sum(trace1);
    N2 = N2+sum(trace2);
end
N1 = N1/Ntrial;
N2 = N2/Ntrial;
xctemp = xctemp/Ntrial;
xc = xctemp/sqrt(N1*N2);