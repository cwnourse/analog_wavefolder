% wavefolder circuit dsp mockup
threshold = .8;  % 0->1
folds=2;
LPF = 0.5;  %0->1


duration = 1;  % second
fs=44100;  % Hz
tt=0:1/fs:duration-1/fs;  % len = duration*fs

L = fs*duration; % samples

f = 1000;  % Hz
A =1;  % say, volts
sig_in = A*sin(2*pi*f*tt);

clipped =zeros(L,folds);
peaks = zeros(L,folds);

cur_peak=sig_in;
for j=1:folds
    thresh_v = threshold*A;
    for i=1:L  % sample wise loop
        if abs(cur_peak(i)) > thresh_v
            clipped(i,j) = sign(cur_peak(i))*thresh_v;
            peaks(i,j) = cur_peak(i)-sign(cur_peak(i))*thresh_v;
        else
            clipped(i,j) = cur_peak(i);
            peaks(i,j) = 0;
        end
    end
    clipped(:,j) = (-1)^(j-1)*clipped(:,j);
    cur_peak = peaks(:,j);  % use this peak data for next folding stage
    peaks(:,j) = (-1)^(j)*peaks(:,j);
    A = A - thresh_v;  % amplitude of peak waveform    
end
% to get folded waveform, sum up clipped and last peak with alternating
% signs

sig_fold = sum(clipped,2)+peaks(:,end);
subplot(2,1,1)
plot(tt,sig_fold)
subplot(2,1,2)
my_fft(sig_fold, L, fs)







% % plot(tt, sig_in)
% my_fft(sig_in, L, fs);
% % comparator system: produce zero-bias peak clippings
%     % how to do this analog with adjustable threshold?
% thresh_v = threshold*A;
% peaks_old = sig_in;
% for j=1:folds
%     for i=1:length(tt)
%         if abs(sig_in(i)) > thresh_v
%             peaks_new(i) = peaks_old(i)-sign(peaks_old(i))*thresh_v;
%         else
%             peaks_new(i) = 0;
%         end
%     end
% end
% 
% % plot(tt, peaks)
% 
% % wavefolding system: fold peaks over onto themselves
% sig_fold = sig_in - 2*peaks;
% % plot(tt, sig_fold)
% my_fft(sig_fold,L, fs);

function p1 = my_fft(sig, L, fs)
    y = fft(sig);
    p2 = abs(y/L);
    p1 = p2(1:L/2+1);
    p1(2:end-1) = 2*p1(2:end-1);
    f = fs*(0:(L/2))/L;
%     figure()
    plot(f, p1)
    title("FFT")
    xlabel("f (Hz)")
    ylabel("|P(f)|")
end

function y = clip_fn(sig, thresh)
    L = length(sig);
    y = zeros(size(sig));
    for i=1:L
        if(abs(sig(i))) > thresh
            y(i)=thresh;
        else
            y(i)=sig(i);  % implicit
        end
    end
end

function y = peaks_fn(sig, thresh)
    L = length(sig);
    y = zeros(size(sig));
    for i=1:L
        if(abs(sig(i))) > thresh
            y(i)=sig(i)-sign(sig(i))*thresh;
        else
            y(i)=0;  % implicit
        end
    end
end

        




