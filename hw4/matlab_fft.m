clear;
[y, fs] = audioread('gliss.mp4');
y = y(:,1); % Just one channel
N = length(y);

nFFT = 1024; % This can be changed for better/worse resolution
hop = floor(nFFT/4);
nFrames = floor(N/hop) - 1; % Approximately... I guess...
fprintf("%d frames.\n", nFrames);
F = zeros(nFFT, nFrames);
w = hanning(nFFT); % Window choice, again, this can be changed

% w checksum, to verify validity of the hanning function
wchecksum = 0;
for n = 1:numel(w)
    wchecksum = wchecksum + w(n);
end
fprintf("W checksum: %d.\n", wchecksum);

for n = 1:nFrames
    iStart = (n-1) * hop;
    if iStart+nFFT > N, break; end
    F(:,n) = fft(w .* y(iStart+1 : iStart+nFFT));
    G(:,n) = 20*log10(abs(F(1:nFFT/2, n)));
end

% F checksum, to verify validity of the code
Fchecksum = 0;
Frows = size(F, 1); % First row length
Fcols = size(F', 1); % First col length
for i = 1:Frows
    for j = 1:Fcols
        Fchecksum = Fchecksum + F(i, j);
    end
end
fprintf("F checksum: %d.\n", Fchecksum);