% recursive wavefolder
A = 1;
folds = .1;  % folding thresholds

% generate a test signal
fs=44100;
tt=0:1/fs:1; % 1 second
f=1;  % 1 Hz
sig_in = sin(2*pi*f*tt);

% plot(tt, sig_in)

%fold the wave
sig_out=zeros(1,length(sig_in));
for i=1:length(sig_in)
   cur = sig_in(i);
   if cur >= 0
       if cur > folds
%            cur = folds - (cur - folds);  % = 2*folds - sig_in
           cur = 2*folds - cur;
       end
   else
       fold_neg = -1*folds;
       if -1*cur > folds
%            cur = -1*folds - (cur - -1*folds);  % -2*f - cur
           cur = -2*folds - cur;
       end
   end
   sig_out(i)=cur;
end

plot(tt,sig_out)