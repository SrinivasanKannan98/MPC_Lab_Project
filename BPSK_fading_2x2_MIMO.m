clc; clear; close all;

M = 2; %Constellation Size
N = 10 ^ 6; %No. of Bits
data_points = 0 : M - 1; %Symbols of the constellation
Eb = 1;
constellation = -sqrt(Eb) * exp(-1i*2*pi*data_points/M); %Constellation Points

EbNo_dB = [-5 : 0.5 : 15]; %Array of SNR values used for the simuation (dB)
ber_array_alamouti_2x2 = zeros(size(EbNo_dB));
EbNo = 10 .^ (EbNo_dB/10); %SNR values in the Linear Scale

bitstream = randi([0 1],N,1); %Random bitstream generation
bitstream_constellation_rep = -sqrt(Eb/2) * exp(-1i*2*pi*(bitstream)/M); %Converting the bitstream representation to
%constellation points; sqrt(1/2) factor is added to ensure that Eb energy is transmitted in every time instant

symbolsTxD = reshape(bitstream_constellation_rep,2,N/2); %Every two symbols are encoded using 2x2 Alamouti codes and
%sent from 2 transmitting antennas in 2 time instances

%Alamouti Encoding & Decoding {Every 2 symbols are encoded using 2x2 Alamouti code and corrupted using AWGN and Rayleigh
%fading before decoding}

for ii = 1 : length(EbNo)

	numErrs = 0;

	for iii = 1 : size(symbolsTxD,2)

		s = symbolsTxD(:,iii); %Extract 2 symbols from the bitstream in each iteration
		s_1 = s; %Encoded symbols in the first time instance
		s_2 = alamoutiEncoder2x2MIMO(s.',2);

		S = [s_1 s_2]; %Encoded Symbol Matrix

		h = sqrt(0.5) .* (randn(2,2) + 1i * randn(2,2)); %Channel Fading Matrix for 4x4 links between Tx and Rx, it is expected
		%that the Channel Fading Matrix remains constant over 2 time instances
		N0_matrix = sqrt((Eb/EbNo(ii))) * randn(2,2) + 1i * (sqrt((Eb/EbNo(ii))) * randn(2,2)); %Generating Noise Samples
		%from a Circularly Symmetric Gaussian Distribution of Variance N0/2 along each dimension

		R = h*S + N0_matrix; %Received symbols at each antenna element over 2 time instances

		%Decoding Operations
		Y = [R(:,1); conj(R(:,2))];
		
 		H_comp_1 = h; 
 		H_comp_2 = alamoutiEncoder2x2MIMO(h,2);

 		H_comp = [H_comp_1; H_comp_2]; %Composite Channel Matrix
        C1 = H_comp(:,1);
        C1_norm = sqrt(sum(abs(C1) .^ 2,1));
        C2 = H_comp(:,2);
        C2_norm = sqrt(sum(abs(C2) .^ 2,1));
        
        s_hat_1 = sum((conj(C1) .* Y),1) ./ C1_norm;
        s_hat_2 = sum((conj(C2) .* Y),1) ./ C2_norm;
        s_hat = [s_hat_1; s_hat_2];
         
        s_decoded = zeros(size(s_hat));
		%Minimum Euclidean Distance decoding 
		for k = 1 : length(s_hat)
			EucD = abs(constellation - s_hat(k) * ones(size(constellation)));
	        %Computing the Euclidean distance of the received symbol from each contellation point
	        %for the given Modulation Scheme
	        [~,pos] = min(EucD); %Minimum Euclidean distance computation
	        s_decoded(k) = constellation(pos); %Decision based on Minimum Euclidean Distance
        end 

        s = s > 0;
        s_decoded = s_decoded > 0;
        [Errors,~] = biterr(s,s_decoded);
        numErrs = numErrs + Errors;

	end

	ber_array_alamouti_2x2(ii) = numErrs/N;

end

semilogy(EbNo_dB,ber_array_alamouti_2x2,'-rp','LineWidth',2);
legend('2 x 2 MIMO using Alamouti Codes in Rayleigh Fading + AWGN');
xlabel('$\frac{Eb}{N0} (dB)$','Interpreter','latex');
ylabel('BER');
title('Comparison of BER vs. SNR for (2x2) Transmit diversity for BPSK');
grid on;

save('MIMO_2x2_BPSK.mat','EbNo_dB','ber_array_alamouti_2x2');
