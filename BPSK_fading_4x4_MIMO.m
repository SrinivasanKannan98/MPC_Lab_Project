clc; clear; close all;

M = 2; %Constellation Size
N = 10 ^ 4; %No. of Bits
data_points = 0 : M - 1; %Symbols of the constellation
Eb = 1;
constellation = -sqrt(Eb/4) * exp(-1i*2*pi*data_points/M); %Constellation Points

EbNo_dB = [-5 : 0.5 : 15]; %Array of SNR values used for the simuation (dB)
ber_array_alamouti_4x4 = zeros(size(EbNo_dB));
EbNo = 10 .^ (EbNo_dB/10); %SNR values in the Linear Scale

bitstream = randi([0 1],N,1); %Random bitstream generation
bitstream_constellation_rep = -sqrt(Eb/4) * exp(-1i*2*pi*(bitstream)/M); %Converting the bitstream representation to
%constellation points; sqrt(1/4) factor is added to ensure that Eb energy is transmitted in every time instant

symbolsTxD = reshape(bitstream_constellation_rep,4,N/4); %Every four symbols are encoded using 4x4 Alamouti codes and
%sent from 4 transmitting antennas in 4 time instances

%Alamouti Encoding & Decoding {Every 4 symbols are encoded using 4x4 Alamouti code and corrupted using AWGN and Rayleigh
%fading before decoding}

for ii = 1 : length(EbNo)

	numErrs = 0;

	for iii = 1 : size(symbolsTxD,2)

		s = symbolsTxD(:,iii); %Extract 4 symbols from the bitstream in each iteration
		s_1 = s; %Encoded symbols in the first time instance
		s_2 = alamoutiEncoder4x4MIMO(s.',2);
		s_3 = alamoutiEncoder4x4MIMO(s.',3);
		s_4 = alamoutiEncoder4x4MIMO(s.',4);

		S = [s_1 s_2 s_3 s_4]; %Encoded Symbol Matrix

		h = sqrt(0.5) .* (randn(4,4) + 1i * randn(4,4)); %Channel Fading Matrix for 4x4 links between Tx and Rx, it is expected
		%that the Channel Fading Matrix remains constant over 4 time instances
		N0_matrix = sqrt((Eb/EbNo(ii))/2) * randn(4,4) + 1i * (sqrt((Eb/EbNo(ii))/2) * randn(4,4)); %Generating Noise Samples
		%from a Circularly Symmetric Gaussian Distribution of Variance N0/2 along each dimension

		R = h*S + N0_matrix; %Received symbols at each antenna element over 4 time instances

		%Decoding Operations
		Y = [R(:,1); conj(R(:,2)); conj(R(:,3)); R(:,4)];
		
		H_comp_1 = h; 
		H_comp_2 = alamoutiEncoder4x4MIMO(h,2);
		H_comp_3 = alamoutiEncoder4x4MIMO(h,3);
		H_comp_4 = alamoutiEncoder4x4MIMO(h,4);
		H_comp = [H_comp_1; H_comp_2; H_comp_3; H_comp_4]; %Composite Channel Matrix

		s_hat = pinv(H_comp) * Y; %The estimates of the 4 transmitted symbols is obtained by zero-forming receive combining
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

	ber_array_alamouti_4x4(ii) = numErrs/N;

end

semilogy(EbNo_dB,ber_array_alamouti_4x4,'-rp','LineWidth',2);
legend('4 x 4 MIMO using Alamouti Codes in Rayleigh Fading + AWGN');
xlabel('$\frac{Eb}{N0} (dB)$','Interpreter','latex');
ylabel('BER');
title('Comparison of BER vs. SNR for (4x4) Transmit diversity for BPSK');
grid on;

save('MIMO_4x4_BPSK.mat','EbNo_dB','ber_array_alamouti_4x4');


