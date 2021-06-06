clc; clear; close all;

M = 2; %Constellation Size
N = 10 ^ 6; %No. of bits 
data_points = 0 : M - 1; %Symbols for the modulation scheme
Eb = 1; %Bit Energy
constellation = -sqrt(Eb/2) * exp(-1i*2*pi*data_points/M); %Constellation Points 
bitstream_decoded_OSTBC_2 = [];

EbNo_dB = [-5 : 0.5 : 15]; %Array of SNR values used for simulation (in dB)
ber_array_OSTBC_2 = [];
EbNo = 10 .^ (EbNo_dB/10); %SNR values in linear scale

bitstream = randi([0 1],N,1); %Random bitstream generation
bitstream_constellation_rep = -sqrt(Eb/2) * exp(-1i*2*pi*(bitstream)/M); %Converting 
%bitstream representation to constellation points (sqrt(1/2)) factor is added to ensure Eb 
%energy is transmitted in every time instant

h = (1/sqrt(2)) * (randn(1,N) + 1i * randn(1,N)); %Rayleigh Flat Fading factor (single tap).
%N/2 coefficients are used for characterising the fading coefficients for Antenna 1 & remaining N/2 for antenna 2
%It is assumed that the channel matrix is constant user two time instances
h_mod = kron(reshape(h,2,N/2),ones(1,2)); %Fading coefficients for each antenna are stored in each row
STBC_coded_data = zeros(2,N);
STBC_coded_data(:,1 : 2 : end) = reshape(bitstream_constellation_rep,2,N/2); %[x1(1:2:N);x2(1:2:N)]
STBC_coded_data(:,2 : 2 : end) = repmat([-1;1],1,N/2) .* flipud(reshape(conj(bitstream_constellation_rep),2,N/2));
%[-x2*(2:2:N);x1*(2:2:N)]

for j = 1 : length(EbNo)

    N0_vector = sqrt((Eb/EbNo(j))/2) * randn(1,N) + 1i * (sqrt((Eb/EbNo(j))/2) * randn(1,N)); %Generating noise samples from a
    %Circularly Symmetric Gaussian Distribution of variance N0/2 along each dimension

    rcvd_data = sum(STBC_coded_data .* h_mod,1) + N0_vector;
    y = reshape(rcvd_data,2,N/2);
    y(2,:) = conj(y(2,:)); % [y1(1 : N/2); y2*(1 : N/2)]

    C = zeros(2,N);
    C(:,[1 : 2 : end]) = reshape(h,2,N/2); 
    C(:,[2 : 2 : end]) = repmat([1;-1],1,N/2) .*  flipud(reshape(h,2,N/2)); 
    C(2,:) = conj(C(2,:));

    C1 = C(:,[1 : 2 : end]);
    C1_norm = sqrt(sum(abs(C1) .^ 2,1));
    C2 = C(:,[2 : 2 : end]);
    C2_norm = sqrt(sum(abs(C2) .^ 2,1));

    rcvd_data_STBC_decoded = zeros(1,N);
    rcvd_data_STBC_decoded(1 : 2 : end) = sum((conj(C1) .* y),1) ./ C1_norm;
    rcvd_data_STBC_decoded(2 : 2 : end) = sum((conj(C2) .* y),1) ./ C2_norm;
   
    %Minimum Euclidean Distance decoding (for the fading + AWGN corrupted received bits)
    for k = 1 : length(rcvd_data_STBC_decoded)
        EucD = abs(constellation - rcvd_data_STBC_decoded(k) * ones(size(constellation)));
        %Computing the Euclidean distance of the received symbol from each contellation point
        %for the given Modulation Scheme
        [~,pos] = min(EucD); %Minimum Euclidean distance computation
        bitstream_decoded_OSTBC_2(k) = data_points(pos); %Decision based on minimum Euclidean
        %distance
    end

    [~,ber] = biterr(bitstream,bitstream_decoded_OSTBC_2'); %BER computation  
    ber_array_OSTBC_2(j) = ber;    

end

semilogy(EbNo_dB,ber_array_OSTBC_2,'-rp','LineWidth',2);
legend('2 x 1 OSTBC with Rayleigh fading + AWGN');
xlabel('$\frac{Eb}{N0} (dB)$','Interpreter','latex');
ylabel('BER');
title('Comparison of BER vs. SNR for (2x1) Transmit diversity for BPSK');
grid on;

save('MISO_BPSK.mat','EbNo_dB','ber_array_OSTBC_2');