clc; clear; close all;

M = 2; %Constellation Size
N = 10 ^ 6; %No. of bits 
data_points = 0 : M - 1; %Symbols for the modulation scheme
Eb = 1; %Bit Energy
constellation = -sqrt(Eb) * exp(-1i*2*pi*data_points/M); %Constellation Points
bitstream_decoded_rayleigh = [];
bitstream_decoded_awgn = [];

EbNo_dB = [-5 : 0.5 : 15]; %Array of SNR values used for simulation (in dB)
ber_array_rayleigh = [];
ber_array_awgn = [];
EbNo = 10 .^ (EbNo_dB/10); %SNR values in linear scale

bitstream = randi([0 1],N,1); %Random bitstream generation
bitstream_constellation_rep = -sqrt(Eb) * exp(-1i*2*pi*(bitstream)/M); %Converting 
%bitstream representation to constellation points

for j = 1 : length(EbNo)

    N0_vector = sqrt((Eb/EbNo(j))/2) * randn(N,1) + 1i * (sqrt((Eb/EbNo(j))/2) * randn(N,1)); %Generating noise samples from a
    %Circularly Symmetric Gaussian Distribution of variance N0/2 along each dimension
    h = (1/sqrt(2)) * (randn(1,N) + 1i * randn(1,N)); %Rayleigh Flat Fading factor (single tap)
    rcvd_constellation_awgn = bitstream_constellation_rep + N0_vector; %Adding AWGN   
    rcvd_constellation_rayleigh = bitstream_constellation_rep .* h(:) + N0_vector;  
    rcvd_constellation_rayleigh_equalised =  rcvd_constellation_rayleigh./(h(:)); %Zero forcing equaliser
    
    %Minimum Euclidean Distance decoding (for the fading + AWGN corrupted received bits)
    for k = 1 : length(rcvd_constellation_rayleigh_equalised)
        EucD = abs(constellation - rcvd_constellation_rayleigh_equalised(k) * ones(size(constellation)));
        %Computing the Euclidean distance of the received symbol from each contellation point
        %for the given Modulation Scheme
        [~,pos] = min(EucD); %Minimum Euclidean distance computation
        bitstream_decoded_rayleigh(k) = data_points(pos); %Decision based on minimum Euclidean
        %distance
    end

    [~,ber] = biterr(bitstream,bitstream_decoded_rayleigh'); %BER computation  
    ber_array_rayleigh(j) = ber;    

    %Minimum Euclidean Distance decoding (for the AWGN corrupted received bits)
    for k = 1 : length(rcvd_constellation_awgn)
        EucD = abs(constellation - rcvd_constellation_awgn(k) * ones(size(constellation)));
        %Computing the Euclidean distance of the received symbol from each contellation point
        %for the given Modulation Scheme
        [~,pos] = min(EucD); %Minimum Euclidean distance computation
        bitstream_decoded_awgn(k) = data_points(pos); %Decision based on minimum Euclidean
        %distance
    end

    [~,ber] = biterr(bitstream,bitstream_decoded_awgn'); %BER computation  
    ber_array_awgn(j) = ber; 

end

semilogy(EbNo_dB,ber_array_awgn,'-bo','LineWidth',2);hold on;
semilogy(EbNo_dB,ber_array_rayleigh,'-kp','LineWidth',2);
legend('AWGN','Rayleigh fading + AWGN');
xlabel('$\frac{Eb}{N0} (dB)$','Interpreter','latex');
ylabel('BER');
title('Comparison of BER vs. SNR for AWGN and Rayleigh fading channel');
grid on;

save('SISO_BPSK.mat','EbNo_dB','ber_array_awgn','ber_array_rayleigh');