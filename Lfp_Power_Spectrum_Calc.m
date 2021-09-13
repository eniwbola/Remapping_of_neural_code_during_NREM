

function [sigFreq, fullSig,findTimes,lfp] = calculateSignalFrequency(sData, skipTime)

findTimes = sData(:,3) >= skipTime;
sData(~findTimes,:) = [];


sData(:,3) = sData(:,3).*20;
uniqueTimes = unique(sData(:,3));
uniqueTimes
histTimes = histc(sData(:,3), uniqueTimes);

spikeVector = zeros(1,max(sData(:,3)));
spikeVector(uniqueTimes) = histTimes;

xVals = -10:0.01:10;
sigma = 2;%0.2%2;
gaussFunc = (1/sqrt(2*pi*sigma^2))*exp(-0.5*(xVals./sigma).^2);
%signal = conv(spikeVector, gaussFunc,'same');
%signal = conv(spikeVector, gaussFunc);
signal = conv(spikeVector, gaussFunc);
% signal=signal(419000:end);%signal;%(380000:end)% be edit

%signal=signal(110000:end);%signal;%(380000:end)% be edit
signal=signal(120000:end);%signal;%(380000:end)% be edit


signal=signal-mean(signal);
lfp=signal;
spikeVector;
freqs = fft(signal);

% figure(77)
% plot(lfp2)
% figure(77)

P2 = abs(freqs/length(signal));
P1 = P2(1:ceil(length(signal)/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = 20000*(0:ceil(length(signal)/2))/ceil(length(signal));
fullSig = [transpose(f(f<=40)) transpose(P1(f<=40))];

findFreqInRange = f >= 5 & f < 20;
f = f(findFreqInRange);
P1 = P1(findFreqInRange);

sigFreq = f(P1 == max(P1));

end
