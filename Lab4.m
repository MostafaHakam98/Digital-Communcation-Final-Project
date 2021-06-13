%% Clear Workspace

clear all;
clc;

%%  Part(1):

%   Simulation Parameters:
%   Number of Bits/SNR = 1e4 bits
%   Signal To Noise Ratio Range = 0 to 60 db with 4 db steps

N_bits = 1e4;
SNR_range = 0 : 4 : 60;

%%  Part(2):

%   Generate Random Binary Data Vector

% Creating a Sequence of Zeroes of N_bits
% Looping though the whole sequence
% Generating a random number (1, 2) - 1 = (0, 1)
% Adding the random bit to the sequence

bit_seq = zeros(1, N_bits);

for index = 1:N_bits
    temp = randi(2) - 1;
    bit_seq(index) = bit_seq(index) + temp;
end

%%  Part(3.1):

%   Modulating the sequence to: OOK, PRK FSK

% Bit Period = 1/30
% Carrier Frequency = 1

ts = 1/30;
fs = 1/ts;
fc = 1;

tb = 0 : ts : 1-ts;
tt = 0 : ts : N_bits - ts;

% OOK
a1 = cos(2*pi*fc*tb)*1;
a0 = cos(2*pi*fc*tb)*0;

% PRK
p1 = cos(2*pi*fc*tb)*1;
p0 = cos(2*pi*fc*tb)*(-1);

% FSK
f1 = cos(2*pi*fc*tb);
f0 = sin(2*pi*fc*tb);

% Modulation

ook = [];
prk = [];
fsk = [];
for index = 1:N_bits
    if bit_seq(index) == 1
        ook = [ook a1];
        prk = [prk p1];
        fsk = [fsk f1];
    else
        ook = [ook a0];
        prk = [prk p0];
        fsk = [fsk f0];
    end
end

% Plot OOK, PRK, FSK

figure(1)
subplot(411)
stairs(0:10,[bit_seq(1:10) bit_seq(10)],'linewidth',1.5)
axis([0 10 -0.5 1.5])
title('Message Bits');grid on

t = 0 : ts : 10-ts;

subplot(412)
plot(t, ook(1 : 10/ts),'b','linewidth',1.5)
title('OOK Modulation');grid on

subplot(413)
plot(t, prk(1 : 10/ts),'k','linewidth',1.5)
title('PRK Modulation');grid on

subplot(414)
plot(t, fsk(1 : 10/ts),'r','linewidth',1.5)
title('FSK Modulation');grid on
xlabel('Time');

%%  Part(3.2):

%   Converting OOK, PRK, FSK with no Carrier

ts = 1;
fs = 1;

t = 0 : ts : 10;

% OOK
a1 = ones(1, 1/ts);
a0 = zeros(1, 1/ts);

% PRK
p1 = ones(1, 1/ts);
p0 = -1 * ones(1, 1/ts);

% FSK
f1 = ones(1, 1/ts);
f0 = i * ones(1, 1/ts);

% Modulation

ook_ = [];
prk_ = [];
fsk_ = [];
for index = 1:N_bits
    if bit_seq(index) == 1
        ook_ = [ook_ a1];
        prk_ = [prk_ p1];
        fsk_ = [fsk_ f1];
    else
        ook_ = [ook_ a0];
        prk_ = [prk_ p0];
        fsk_ = [fsk_ f0];
    end
end

% Plot OOK, PRK, FSK with no Carrier

figure(2)
subplot(411)
stairs(t,[bit_seq(1:10) bit_seq(10)],'linewidth',1.5)
axis([0 10 -0.5 1.5])
title('Message Bits');grid on

t = 0 : ts : 10;

subplot(412)
stairs(t, ook_(1 : 10/ts + 1),'b','linewidth',1.5)
title('OOK Modulation');grid on

subplot(413)
stairs(t, prk_(1 : 10/ts + 1),'k','linewidth',1.5)
title('PRK Modulation');grid on

subplot(414)
stairs(t, fsk_(1 : 10/ts + 1),'r','linewidth',1.5)
title('FSK Modulation');grid on
xlabel('Time');

%%  Part(4, 5, 6):

%   Part(4):

%   Applying Noise to Bits

OOK_BER_list = [];
PRK_BER_list = [];
FSK_BER_list = [];

OOK_noise_list = [];
PRK_noise_list = [];
FSK_noise_list = [];

for snr = SNR_range
    OOK_noise = awgn(ook_, snr, 'measured');
    PRK_noise = awgn(prk_, snr, 'measured');
    FSK_noise = awgn(fsk_, snr, 'measured');
    
    if snr == 0 | snr == 12 | snr == 60
        OOK_noise_list = [OOK_noise_list OOK_noise];
        PRK_noise_list = [PRK_noise_list PRK_noise];
        FSK_noise_list = [FSK_noise_list FSK_noise];
    end
    
    % Detection in Reciever
    OOK_Reciever = [];
    PRK_Reciever = [];
    FSK_Reciever = [];
    
%   Part(5):
    
    for index = 1:N_bits
        % OOK Detection
        if sum(real(OOK_noise(1+((index-1)*fs) : (index*fs))))/fs > 0.5
            OOK_Reciever = [OOK_Reciever 1];
        else
            OOK_Reciever = [OOK_Reciever 0];
        end
        % PRK Detection
        if sum(real(PRK_noise(1+((index-1)*fs) : (index*fs))))/fs > 0
            PRK_Reciever = [PRK_Reciever 1];
        else
            PRK_Reciever = [PRK_Reciever 0];
        end
        % FSK Detection
        if sum(real(FSK_noise(1+((index-1)*fs) : (index*fs))))/fs > 0.5
            FSK_Reciever = [FSK_Reciever 1];
        else
            FSK_Reciever = [FSK_Reciever 0];
        end
    end
    
%   Part(6):
    
    % Bit Error Rate Calculation
    OOK_BER = biterr(OOK_Reciever, bit_seq);
    PRK_BER = biterr(PRK_Reciever, bit_seq);
    FSK_BER = biterr(FSK_Reciever, bit_seq);
    
    OOK_BER = OOK_BER / N_bits;
    PRK_BER = PRK_BER / N_bits;
    FSK_BER = FSK_BER / N_bits;
    
    OOK_BER_list = [OOK_BER_list OOK_BER];
    PRK_BER_list = [PRK_BER_list PRK_BER];
    FSK_BER_list = [FSK_BER_list FSK_BER];
end

%%  Part(7):

%   Plot OOK, PRK, FSK with Noise at SNR = 0, 12, 60 db

for index = 2:4
    OOK_noise = OOK_noise_list(1 + (index - 2)*N_bits/ts : (index - 1)*N_bits/ts);
    PRK_noise = PRK_noise_list(1 + (index - 2)*N_bits/ts : (index - 1)*N_bits/ts);
    FSK_noise = FSK_noise_list(1 + (index - 2)*N_bits/ts : (index - 1)*N_bits/ts);
    
    if index == 2
        SNR_num = '0';
    elseif index == 3
        SNR_num = '12';
    elseif index == 4
        SNR_num = '60';
    end
    
    figure(index + 1)
    subplot(411)
    stairs(0:10,[bit_seq(1:10) bit_seq(10)],'linewidth',1.5)
    axis([0 10 -0.5 1.5])
    title('Message Bits');grid on

    t = 0 : ts : 10;

    subplot(412)
    stairs(t, OOK_noise(1 : 10/ts + 1),'b','linewidth',1.5)
    title(['OOK with Noise at SNR = ' SNR_num]);grid on

    subplot(413)
    stairs(t, PRK_noise(1 : 10/ts + 1),'k','linewidth',1.5)
    title(['PRK with Noise at SNR = ' SNR_num]);grid on

    subplot(414)
    stairs(t, FSK_noise(1 : 10/ts + 1),'r','linewidth',1.5)
    title(['FSK with Noise at SNR = ' SNR_num]);grid on
    xlabel('Time');
    
end

%%  Part(8):

%   Plotting BER against SNR

figure(6)
semilogy(SNR_range, OOK_BER_list, 'b','linewidth',2)
title('BER Vs SNR')
grid on;
hold on
semilogy(SNR_range, FSK_BER_list,'r','linewidth',2)
semilogy(SNR_range, PRK_BER_list, 'k','linewidth',2)
xlabel('SNR(dB)')
ylabel('BER')
hold off
legend('OOK','FSK','PRK');
