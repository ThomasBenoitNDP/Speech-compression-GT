clc;
clear all;
load('FTPSig.mat');
len=length(signal);

x=signal;
Fs=32000;
L=length(x);
dt = 1/Fs;                   % seconds per sample
t = (0:dt:(1/20)-dt);     % seconds

figure(1)
Y_mags = fft(signal);
num_bins = length(Y_mags);
plot([0:1/(num_bins -1):1], Y_mags(1:num_bins),'r')
grid on
title('Signal');
xlabel('Normalised frequency (\pi rads/sample)');
ylabel('Magnitude');

%chirp1
A1=1;
Fc1 = 4000;                     % hertz
Fch1 = 40; % hertz
phi1=(2*Fc1*pi*t)+40*cos(2*pi*Fch1*t);
ch1 = A1*cos(phi1);
fch1 = freCh(ch1,Fs);

A2=1;
Fc2 = 8000;                     % hertz
Fch2 = 40;
phi2=(2*pi*Fc2*t)+40*cos(2*pi*Fch2*t);
ch2 = A2*cos(phi2);
fch2 = freCh(ch2,Fs);

A3=1;
Fc3=12000;
Fch3 = 40; % hertz
phi3=(2*Fc3*pi*t)+40*cos((2*pi*Fch3*t));
ch3 = A3*cos(phi3);
fch3 = freCh(ch3,Fs);

%% chirp 2
A1=1;
Fc1 = 4000;                     % hertz
Fch1 = 20; % hertz
phi1=(2*Fc1*pi*t)+40*cos(2*pi*Fch1*t);
ch21 = A1*cos(phi1);
fch21 = freCh(ch21,Fs);

A2=1;
Fc2 = 8000;                     % hertz
Fch2 = 20;
phi2=(2*pi*Fc2*t)+40*cos(2*pi*Fch2*t);
ch22 = A2*cos(phi2);
fch22 = freCh(ch22,Fs);

A3=1;
Fc3=12000;
Fch3 = 20; % hertz
phi3=(2*Fc3*pi*t)+40*cos(2*pi*Fch3*t);
ch23 = A3*cos(phi3);
fch23 = freCh(ch23,Fs);

%%
%chirp3
Fc1=4000;
xch1= -chirp(t,Fc1,1,Fc1+1/Fs,'linear');

Fc2=8000;
xch2= -chirp(t,Fc2,1,Fc2+1/Fs,'linear');

Fc3=12000;
xch3= -chirp(t,Fc3,1,Fc3+1/Fs,'linear');


%%
yband3=[9600 14400];
xb3=bandpass(x,yband3,Fs);
N_Frames = 1600;
N_seg = floor(length(xb3)/N_Frames) ;
Iter=N_seg;
x_frag = zeros(Iter,N_Frames); 
j=(N_Frames);
x_frag(1,:)=xb3(1:N_Frames);

for i=2:Iter
    x_frag(i,:)=xb3(j+1:(j+N_Frames));
    j=j+N_Frames;
end
ych=zeros(1,1600);
YCH3=[];
for i=1:Iter
    x1=x_frag(i,:);
    f3 = freCh(x1,Fs);
    fd31=f3-fch3;
    count1=0;
    for i=1:1600
        if(fd31(i)==0)
            count1=count1+1;
        end 
    end
    if(count1>600)
        ych=ch3;
    else
        fd32=f3-fch23;
        count2=0;
        for i=1:1600
            if(fd32(i)==0)
                count2=count2+1;
            end 
        end
        if(count2>600)
            ych=ch23;
        end
    end
    if(ych(1)==0)
        ych=xch3;
    end
    YCH3=[YCH3 ych];
    ych=zeros(1,1600);
end

%%
yband2=[5600 10400];
xb2=bandpass(signal,yband2,Fs);
N_Frames = 1600;
N_seg = floor(length(xb2)/N_Frames) ;
Iter=N_seg;
x_frag = zeros(Iter,N_Frames); 
j=(N_Frames);
x_frag(1,:)=xb2(1:N_Frames);

for i=2:Iter
    x_frag(i,:)=xb2(j+1:(j+N_Frames));
    j=j+N_Frames;
end
ych=zeros(1,1600);
YCH2=[];
for i=1:Iter
    x1=x_frag(i,:);
    f2 = freCh(x1,Fs);
    fd21=f2-fch2;
    count1=0;
    for i=1:1600
        if(fd21(i)==0)
            count1=count1+1;
        end 
    end
    if(count1>700)
        ych=ch2;
    else
        fd22=f2-fch22;
        count2=0;
        for i=1:1600
            if(fd22(i)==0)
                count2=count2+1;
            end 
        end
        if(count2>700)
            ych=ch22;
        end
    end
    if(ych(1)==0)
        ych=xch2;
    end
    YCH2=[YCH2 ych];
    ych=zeros(1,1600);
end

%%
yband1=[1600 6400];
xb1=bandpass(signal,yband1,Fs);
N_Frames = 1600;
N_seg = floor(length(xb1)/N_Frames) ;
Iter=N_seg;
x_frag = zeros(Iter,N_Frames); 
j=(N_Frames);
x_frag(1,:)=xb1(1:N_Frames);

for i=2:Iter
    x_frag(i,:)=xb1(j+1:(j+N_Frames));
    j=j+N_Frames;
end
ych=zeros(1,1600);
YCH1=[];
for i=1:Iter
    x1=x_frag(i,:);
    f1 = freCh(x1,Fs);
    fd11=f1-fch1;
    count1=0;
    for i=1:1600
        if(fd11(i)==0)
            count1=count1+1;
        end 
    end
    if(count1>1000)
        ych=ch1;
    else
        fd12=f1-fch21;
        count2=0;
        for i=1:1600
            if(fd12(i)==0)
                count2=count2+1;
            end 
        end
        if(count2>1000)
            ych=ych+ch21;
        end
    end
    if(ych(1)==0)
        ych=xch1;
    end
    YCH1=[YCH1 ych];
    ych=zeros(1,1600);
end

%%

YCHb1=bandpass(YCH1,yband1,32000);
YCHb2=bandpass(YCH2,yband2,32000);
YCHb3=bandpass(YCH3,yband3,32000);
YCH=YCHb1+YCHb2+YCHb3;

figure(2)
subplot(4,1,1)
spectrogram(signal(1:32000),256,250,512,32000,'yaxis');
subplot(4,1,2)
spectrogram(YCHb1(1:32000),256,250,512,32000,'yaxis');
subplot(4,1,3)
spectrogram(YCHb2(1:32000),256,250,512,32000,'yaxis');
subplot(4,1,4)
spectrogram(YCHb3(1:32000),256,250,512,32000,'yaxis');

R1=signal-YCH';
[a,g]=lpc(YCH,100);
R=filter(a,1,signal);
figure(3)
subplot(4,1,1)
spectrogram(signal(3000:6000),256,250,512,32000,'yaxis');
title('Signal');
subplot(4,1,2)
spectrogram(YCH(3000:6000),256,250,512,32000,'yaxis');
title('Chirp Signal');
subplot(4,1,3)
spectrogram(R1(3000:6000),256,250,512,32000,'yaxis');
title('Residual Signal1 via subtraction');
subplot(4,1,4)
spectrogram(R(3000:6000),256,250,512,32000,'yaxis');
title('Residual Signal via LPC filtering');

fbp=[300 3500];
x1=bandpass(R,fbp,Fs);

% Removing the noise attached to the speech amplitude

xdm1= demod(R,8000,Fs,'am');
xdm = demod(xdm1,4000,Fs,'am');
audiowrite('Residual1.wav',xdm,32000);
audiowrite('ReconstructedChirps.wav',YCH,32000);

figure(4)
Y_mags = abs(fft(xdm));
num_bins = length(Y_mags);
plot([0:1/(num_bins/2 -1):1], Y_mags(1:num_bins/2),'r')
grid on
title('THE SPEECH');
xlabel('Normalised frequency (\pi rads/sample)');
ylabel('Magnitude');


%% Transfering to G729

G = xdm;

% Firstly, downslamping: Fs=32000 -> Fs=8000 | Downslamped by 4
down_speech = downsample(G,4);

Dw = down_speech/(max(abs(down_speech))*1.001);

% The signal with Fs=8000
audiowrite('INPUT_F8000_16BIT_PCM.wav',Dw,8000);


%% LOADING INITIAL SIGNAL AND SIGNAL GIVEN BY G729

% Loading the input of G729
% The speech part, for Fs=320000
[X_speechInit,Fs_speechInit] = audioread("SPEECH_F320000.wav");

X_speechInitW = X_speechInit/(max(abs(X_speechInit))*1.001);    % Normalized the value to compare effectively


% Loading the input of G729
[X_in,Fs_in]     = audioread('INPUT_F8000_16BIT_PCM.wav');
X_inW = X_in/(max(abs(X_in))*1.001);    % Normalized the value to compare effectively

% Loading the output of G729
fid_out = fopen('OUT.PST', 'r');
X_out = fread(fid_out,inf,"int16");
fclose(fid_out);

D8000 = X_out/(max(abs(X_out))*1.001);  % Normalized the value to compare effectively

% UPSAMPLING THE OUT FILE FROM G729
X_out_3200 = upsample(X_out,4);

D32000 = X_out_3200/(max(abs(X_out_3200))*1.001);   % Normalized the value to compare effectively



%% Computing Measurements

%%% The operating bit rate 
% Note that we have opted to a fixed bit rate (the initial version og
% G7.29) of 8 kbit/s.

%%% Distortion 

% 1) signal to coding noise ratio (power of the signal over power of the difference signal between the original and the decoded one) for a waveform-based and a hybrid type coding system;

%SQNR: G729 INPUT Binary File <-> G729 OUTPUT Binary File
P_DIFF_IO = rms(X_inW - D8000)^2;
P_I = rms(X_inW)^2;

SQNR_IO = 10*log10(P_I/P_DIFF_IO)

%SQNR: F=32000 speech File <-> G729 OUTPUT Binary File upsampled by 4 (8000 -> 32000)
P_DIFF_SO = rms(X_speechInitW - D32000)^2;
P_I = rms(X_speechInitW)^2;

SQNR_SO= 10*log10(P_I/P_DIFF_SO)

% 2) spectral distortion using appropriate measures -> Log spectral distance with Cepstral approximation 

%Log cseptral distance: G729 INPUT Binary File <-> G729 OUTPUT Binary File
LSD_IO = log_spectral_dist(X_inW,D8000,5) 

%Log cseptral distance: F=32000 speech File <-> G729 OUTPUT Binary File upsampled by 4 (8000 -> 32000) 
LSD_SO = log_spectral_dist(D32000,X_speechInitW,5)    

% Likelihood dsitortion ~ Cepstral approximation for "little distortion"

%% GENERATING WAV FILES

audiowrite("F32000_CMPRSS_SPEECH.wav",D32000,32000);
audiowrite("F8000_CMPRSS_SPEECH.wav",D8000,8000);


%% RECONSTRUCTION PART 

Rcstr_overall_signal = YCH' + D32000;
audiowrite("RECONSTRUCTED_OVERALL_SIGNAL.wav",Rcstr_overall_signal,32000);

