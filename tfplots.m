function [logscale,smoothmagn]=tfplots(data, Fs, fract, avg)

%TFPLOTS - Smoothed transfer fucntion plotting
%   [FREQ,MAGN]=TFPLOTS(IMPRESP, Fs, FRACT, AVG)
%   Logarithmic transfer function plot from impluse response IMPRESP. 
%   A half hanning window is applied before an FFT, then the data is colleced
%   into logaritmically spaced bins and the average is computed for
%   each bin (100/octave). Then this is smoothed (convolved) by a hanning window, where
%   FRACT defines the fractional-octave smoothing (default is 3, meaning third-octave).
%   The length of the smoothing hanning window is the double compared to the distance
%   defined by FRACT.
%   The sampling frequency is set by FS (default is 44.1 kHz)
%
%   If the AVG variable is set to 'power' then the power is averaged
%   in the logaritmic bins and during smoothing (this is the default - 
%   on the contrary to the TFPLOT function, where 'comp' is the default),
%   if it is 'abs' then the absolute value, and if to 'comp',
%   it averages the complex transfer function.
%   
%   avg values: ('power', 'complex', 'abs')
%
%   Created by: Balazs Bank, 2007-2010.
%   Edited  by: Panagiotis Zachos

octbin = 100; % how many frequency bins are computed per octave

if nargin<4
    avg = 'power';
end
if nargin<3
    fract = 3;
end
if nargin<2
    Fs = 44100;
end

data = data(:);
FFTSIZE=length(data);

logfact = 2^(1/octbin);
LOGN = floor(log(Fs/2)/log(logfact));
logscale = logfact.^(0:LOGN); % logarithmic scale from 1 Hz to Fs/2

tf = fft(data,FFTSIZE);

compamp = tf(1:FFTSIZE/2);


fstep = Fs/FFTSIZE;
logmagn = zeros(1,LOGN+1);
for k = 0:LOGN
    start = round(logscale(k+1)/sqrt(logfact)/fstep);
    start = max(start,1);
    start = min(start,FFTSIZE/2);
    stop = round(logscale(k+1)*sqrt(logfact)/fstep);
    stop = max(stop,1);
    stop = min(stop,FFTSIZE/2);
    if strcmp(avg, 'complex') % averaging the complex transfer function
        if start ~= stop
            for ii = start:stop
                logmagn(k+1) = logmagn(k+1) + compamp(ii);
            end
            logmagn(k+1) = logmagn(k+1)/(stop-start+1);
        else
            logmagn(k+1) = compamp(start);
        end
    elseif strcmp(avg, 'abs') % averaging absolute value
        logmagn(k+1)=mean(abs(compamp(start:stop))); 
    elseif strcmp(avg, 'power') % averaging power
        logmagn(k+1)=sqrt(mean(abs(compamp(start:stop)).^2)); 
    end
end

%creating hanning window
HL = 2*round(octbin/fract); % fractional octave smoothing
hh = hanning(HL);

L = length(logmagn);
logmagn(L+1:L+HL) = 0;

% Smoothing the log. spaced data by convonvling with the hanning window
if strcmp(avg, 'complex') ||  strcmp(avg, 'abs')
   tmp = fftfilt(hh,logmagn); 
   smoothmagn = tmp(HL/2+1:HL/2+L)/sum(hh);
end
if strcmp(avg, 'power')
   tmp=fftfilt(hh,logmagn.^2); 
   smoothmagn = sqrt(tmp(HL/2+1:HL/2+L)/sum(hh));
end
