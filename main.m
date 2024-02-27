clearvars; clc; close all;

PSTART = 2; % number of pole frequencies to start with
SMOOTH = 6; % fractional octave smoothing, e.g., SMOOTH=6 is 6th octave - 0 is no smoothing
PLOTSECT = 0; % 1: plot the response of the separate sections, 0: don't plot
MAX_POLES = 32;
EQ_START = 40;
EQ_END = 10000;
NFIR = 0; % no FIR part is needed for a mimimum-phase filter
POLE_DIST = 0.2; % This parameter determines how close the algorithm is allowed to place the poles. Smaller value => more poles, Bigger value => Less poles. Very small values may lead to instability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Working in different sampling rates may lead to strange results and
% hasn't been tested. Recommended to resample your inputs and then resample
% again at the end. You may need to also resample Am Bm filter coeffcients
% depending on what you need for your application.
Fs = 44100; 
target_imp = audioread('Harman_Target_minph.wav'); % You can set this to whatever you want. Must be minimum phase.

impresp = [0 1 zeros(1,4094)]; % you can substitute this with any impulse response. Some examples are included in demo_IRs folder. Do not forget to 
[~,impresp] = rceps(impresp);

[~,spH] = tfplots(impresp,Fs,SMOOTH,'complex'); % smoothed input response
[fr,spT] = tfplots(target_imp,Fs,SMOOTH,'complex'); % smoothed target response


w = 2*pi*fr/Fs;
% use frequency range fr(start:stop) in the filter design
% start_search : start frequency (default is 20)
% stop_search : end frequency (default is Fs/2)
start_search = 40;
stop_search = Fs/2;
[~,start]=min(abs(fr-start_search)); % start frequency index
[~,stop]=min(abs(fr-stop_search)); % stop frequency index
clear start_search stop_search

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logarithmic pole positioning as a first guess

fplog=logspace(log10(EQ_START),log10(EQ_END),PSTART);
fp=fplog;
p=freqpoles(fp);

%parallel filter design 
[Bm,Am,FIR] = parfiltidfr(w(start:stop),spH(start:stop),spT(start:stop),p,NFIR);

equalized_imp = parfilt(Bm,Am,FIR,impresp)'; % equalized freq. response
[~,equalized] = tfplots(equalized_imp,Fs,SMOOTH,'complex');

%%%% 
% Intervals
upper_freq_lim = 10000;
[~,pole_stop] = min(abs(fr-upper_freq_lim));
intervals = [start, round((pole_stop+start-1)/2); round((pole_stop+start-1)/2), pole_stop];
errors = [immse(db(spT(intervals(1,1):intervals(1,2))), db(equalized((intervals(1,1):intervals(1,2)))));
          immse(db(spT(intervals(2,1):intervals(2,2))), db(equalized((intervals(2,1):intervals(2,2)))))];

command = 1;
NumPolesPlaced = 0; % Number of poles placed since last pole removal operation 

figure
while command ~= 2
    
    % Magnitude Plots
    plot_responses(p,fp,fr,spH,spT,equalized,w,Am,Bm,FIR,PLOTSECT);
    
%     ginput(1); % Uncomment this line to watch step by step operation of algorithm
     
    zero_indices = [];

    if command == 1 % add new pole
        while 1
            [~, max_ind] = max(errors);
            interval_med = round(sum(intervals(max_ind,:))/2);
            npolepos = fr(interval_med); % find worst MSE interval
            
            tmp = [fp npolepos];
            tmp = sort(tmp);
            if min(log10(tmp(2:end)) - log10(tmp(1:end-1))) > POLE_DIST
                new_intervals = [intervals(max_ind,1), interval_med;
                             interval_med, intervals(max_ind,2)];
                intervals = [intervals(1:max_ind-1,:); new_intervals; intervals(max_ind+1:end,:)];
                fp=[fp npolepos];
                fp=sort(fp);
                NumPolesPlaced = NumPolesPlaced + 1;
                break;
            else
                errors(max_ind) = -1;
                zero_indices = [zero_indices max_ind]; %#ok<*AGROW>
            end
            if max(errors) < 0
                command = 2;
                break;
            end
        end
    end
    if command==3 || NumPolesPlaced==3 % delete poles if the response would better match target
        NumPolesPlaced = 0;
        fp = remove_pole(spH,spT,equalized,fp,w,Fs,start,stop,SMOOTH,NFIR,impresp,MAX_POLES);
    end
    p=freqpoles(fp);

    % parallel filter design 
    [Bm,Am,FIR]=parfiltidfr(w(start:stop),spH(start:stop),spT(start:stop),p,NFIR);
    equalized_imp = parfilt(Bm,Am,FIR,impresp)'; % equalized freq. response
    [~,equalized]=tfplots(equalized_imp,Fs,SMOOTH,'complex');

    errors = zeros(length(intervals),1);
    for k = 1:length(intervals)
        int_start = intervals(k,1);
        fin = intervals(k,2);
        errors(k) = immse(db(spT(int_start:fin)), db(equalized(int_start:fin)));
    end
    for k = 1:length(zero_indices)
        errors(zero_indices(k)) = -1;
    end
    
    spDifAlg = abs(spectral_flatness(db(spT(665:1330))) - spectral_flatness(db(equalized(665:1330))));
end

function plot_responses(p,fp,fr,spH,spT,equalized,w,Am,Bm,FIR,PLOTSECT)

semilogx(fr,db(abs(spH)),'LineWidth',2);
hold on;
semilogx(fr,db(abs(spT)),'LineWidth',2);
semilogx(fr,db(abs(equalized)),'LineWidth',2);
grid on
if PLOTSECT % compute and plot the response of the sections
    Z=exp(-1i*w(:)); %creating z^-1
    Z2=Z.^2; %creating z^-2
    Hsectv=zeros(length(w),size(Am,2));
    for k=1:size(Am,2) %second-order sections
        Hsectv(:,k) = (Bm(1,k)+Bm(2,k)*Z)./(Am(1,k)+Am(2,k)*Z+Am(3,k)*Z2);
    end
    semilogx(fr,db(Hsectv),':'); 
    if NFIR>0 % compute and plot the response of the FIR part
        HFIR=freqz(FIR,1,w);
        semilogx(fr,db(HFIR),'k');
    end
end
plot(fp,0,'ko','LineWidth',2); % plot pole frequencies
hold off;
SIZE=10;
set(gca,'FontName','Times','Fontsize',SIZE);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');

s=sprintf('Total filter order: %i', length(p));
text(25,17,s,'FontName','Times','Fontsize',SIZE);
legend('Original','Target','Equalized','Location','Southwest');

axis([18 22050 -20 20]);
end

function [fp] = remove_pole(spH,spT,equalized,fp,w,Fs,start,stop,SMOOTH,NFIR,impresp,max_poles)
    if length(fp) <= 2
        return
    end
    starting_flatness = abs(spectral_flatness(db(spT(665:1330))) - spectral_flatness(db(equalized(665:1330))));
    starting_frechet = DiscreteFrechetDist(real(db(equalized(665:1330))), real(db(spT(665:1330))));
    starting_frechet = (starting_frechet - 15) / (90 - 15);
    starting_metric = (starting_flatness + starting_frechet) / 2;
    
    metric = zeros(1, length(fp));

    for k = 1:length(fp)
        fp_temp=[fp(1:k-1) fp(k+1:end)];
        p_temp=freqpoles(fp_temp);
        [Bm,Am,FIR]=parfiltidfr(w(start:stop),spH(start:stop),spT(start:stop),p_temp,NFIR);
        equalized_imp = parfilt(Bm,Am,FIR,impresp)'; % equalized freq. response
        [~,equalized]=tfplots(equalized_imp,Fs,SMOOTH,'complex');
        metric(k) = 0.5 * abs(spectral_flatness(db(spT(665:1330))) - spectral_flatness(db(equalized(665:1330))));
        frechet = DiscreteFrechetDist(real(db(equalized(665:1330))), real(db(spT(665:1330))));
        metric(k) = metric(k) + 0.5 * (frechet - 15) / (90 - 15);
    end
    [min_val, ind] = min(metric);
    % Remove pole if smaller metric is achieved without it, or if the desired number of poles (max_poles) has been exceeded
    if min_val <= starting_metric || 2*length(fp) > max_poles
        fp=[fp(1:ind-1) fp(ind+1:end)];
    end
end





