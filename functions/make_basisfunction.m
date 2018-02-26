function basisf = make_basisfunction(tfilt, Nf, type, params)
% Make basis functions for learning filters, given:
% * tfilt: the time of the filter
% * Nf = # basis functions
% * type: type of basis functions ('gamma', or 'exponential', or
% 'gaussian')
% * params: a struct with parameters specific to the type of basis function

% returns Nfxlength(tfilt) matrix with basis functions

if nargin == 3
    % no parameters
end

Nt = length(tfilt);
basisf = zeros(Nf, Nt);

for nf = 1:Nf
    if strcmp(type, 'gamma')
        t = 20*[0:1/Nt:1-1/Nt];
        basisf(nf,:) = (t).^(nf-1).*exp(-t);
    elseif strcmp(type, 'exponential')
        basisf(nf,:) = exp(-tfilt/params.tau(nf));
    elseif strcmp(type, 'gaussian')
        sigma = params.sigma(nf);
        mu = params.mu(nf);
        basisf(nf,:) = (1./(sigma*sqrt(2*pi))).*exp(-((tfilt-mu).^2)./(2*sigma^2));
    end
    basisf(nf,:)=basisf(nf,:)./(sum(basisf(nf,:).^2).^0.5);
end