function [FF a g est e zf] = LPC_P3(y,Fs,m,zi,x)
L=length(y);
ww = hamming(L,'symmetric');
V=y.*ww;
[a,g]=lpc(V,m);
FF=log(abs(freqz(1,a,1024))); 
[est,zf] = filter([0 -a(2:end)],1,x,zi);
e = y-est;
end 