function kernel = lew_make_kernels_network(tfilt, Nbf, Nneuron, net, seed)


if nargin == 4
    rng('shuffle')
elseif nargin>4
    if isempty(seed)
        rng('shuffle')
    else
        rng(seed)
    end
end


%% Make basis funcions

% Gamma functions
basisf = make_basisfunction(tfilt, Nbf, 'gamma',[]);

%% Make kernels from basis functions
if strcmp(net, 'homogeneous_t1')
    % type 1
    kernel = zeros(Nneuron, length(tfilt));
    kernelt1 = basisf(3,:);
    kernel(1:Nneuron/2, :) = repmat(kernelt1, Nneuron/2,1);
    kernel(Nneuron/2+1:end, :) = -repmat(kernelt1, Nneuron/2,1);
elseif strcmp(net, 'homogeneous_t1_nonegative')
    % type 1
    kernel = zeros(Nneuron, length(tfilt));
    kernelt1 = basisf(3,:);
    kernel(1:Nneuron, :) = repmat(kernelt1, Nneuron,1);
elseif strcmp(net, 't12')
    kernel = zeros(Nneuron, length(tfilt));
    % type 1
    kernelt1 = basisf(3,:);
    kernel(1:Nneuron/4, :) = repmat(kernelt1, Nneuron/4,1);
    kernel(Nneuron/4+1:Nneuron/2, :) = -repmat(kernelt1, Nneuron/4,1);
    % type 2
    kernelt2 = basisf(3,:).*(0.2.*ones(size(tfilt))-0.8*sin(tfilt*.6));
    kernel(Nneuron/2+1:3*Nneuron/4, :) = repmat(kernelt2, Nneuron/4,1);
    kernel(3*Nneuron/4+1:end, :) = -repmat(kernelt2, Nneuron/4,1);
elseif strcmp(net, 'heterogeneous')
    freqvec = 1.5*rand(1,Nneuron);
    kernel = zeros(Nneuron, length(tfilt));
    env = basisf(3,:);
    for nn = 1:Nneuron
        if nn<=Nneuron/4
            % sin
            kernel(nn, :) = env.*(0.2.*ones(size(tfilt))+0.8*sin(tfilt*freqvec(nn)));
        elseif Nneuron/4 < nn <= Nneuron/2
            % -sin
            kernel(nn, :) = env.*(0.2.*ones(size(tfilt))-0.8*sin(tfilt*freqvec(nn)));
        elseif Nneuron/2< nn <= 3*Nneuron/4
            % cos
            kernel(nn, :) = env.*(0.2.*ones(size(tfilt))+0.8*cos(tfilt*freqvec(nn)));
        elseif 3*Nneuron/4 < nn 
            % -cos
            kernel(nn, :) = env.*(0.2.*ones(size(tfilt))-0.8*cos(tfilt*freqvec(nn)));
        end
    end
    
elseif strcmp(net, 'heterogeneous_matchedt12')
    % This makes a heterogeneous network, but 
    freqvec = .6*rand(1,Nneuron);
    kernel = zeros(Nneuron, length(tfilt));
    env = basisf(3,:);
    for nn = 1:Nneuron
        if nn<=Nneuron/4
            % sin
            kernel(nn, :) = env.*(0.2.*ones(size(tfilt))+0.8*sin(tfilt*freqvec(nn)));
        elseif Nneuron/4 < nn <= Nneuron/2
            % -sin
            kernel(nn, :) = env.*(0.2.*ones(size(tfilt))-0.8*sin(tfilt*freqvec(nn)));
        elseif Nneuron/2< nn <= 3*Nneuron/4
            % cos
            kernel(nn, :) = env.*(0.2.*ones(size(tfilt))+0.8*cos(tfilt*freqvec(nn)));
        elseif 3*Nneuron/4 < nn 
            % -cos
            kernel(nn, :) = env.*(0.2.*ones(size(tfilt))-0.8*cos(tfilt*freqvec(nn)));
        end
    end

end


