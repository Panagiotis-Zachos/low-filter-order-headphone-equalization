function [spf] = spectral_flatness(a)
    pxx = periodogram(a);
    num=geomean(pxx);
    den=mean(pxx);
    spf=num/den ;
end


