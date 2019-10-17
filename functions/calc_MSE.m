function MSE = calc_MSE(sig1, sig2)
    MSE = sum((sig1-sig2).^2);
end