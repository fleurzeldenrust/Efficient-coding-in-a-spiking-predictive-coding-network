function sigd = downsample_mean(signal, nsteps)

[ny, nx] = size(signal);

if ny == 1
    sigd = zeros(1,ceil(nx/nsteps));
elseif nx == 1
    sigd = zeros(ceil(ny/nsteps), nx);
else
    disp('please give one-dimensional signal')
    return
end

for nn=1:length(sigd)-1
    sigd(nn) = mean(signal((nn-1)*nsteps+1:nn*nsteps));
    if isnan(sigd(nn))
        disp('NaN: something went wrong!')
        keyboard
    end
end


if nn*nsteps<length(signal)
    sigd(end) = mean(signal(nn*nsteps+1:end));
end

return