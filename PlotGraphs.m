%Script for plotting the different curves

semilogy(EbNo_dB,ber_array_rayleigh,'-rp','LineWidth',2); hold on;
semilogy(EbNo_dB,ber_array_OSTBC_2,'-kd','LineWidth',2);
semilogy(EbNo_dB,ber_array_MRC_1x2,'-co','LineWidth',2);
semilogy(EbNo_dB,ber_array_alamouti_2x2,'-gd','LineWidth',2);
semilogy(EbNo_dB,ber_array_alamouti_4x4,'-m+','LineWidth',2);
semilogy(EbNo_dB,ber_array_awgn,'-b*','LineWidth',2);

xlabel('$\frac{Eb}{N0} (dB)$','Interpreter','latex');
ylabel('BER');
% legend('SISO in Rayleigh Fading + AWGN', '(2 x 1) MISO in Rayleigh Fading + AWGN', ...
%     '(2 x 2) MIMO in Rayleigh Fading + AWGN', '(4 x 4) MIMO in Rayleigh Fading + AWGN','SISO in AWGN');
legend('SISO in Rayleigh Fading + AWGN', '(2 x 1) MISO in Rayleigh Fading + AWGN', ...
    '(1 x 2) SIMO in Rayleigh Fading + AWGN','(2 x 2) MIMO in Rayleigh Fading + AWGN', ...
    '(4 x 4) MIMO in Rayleigh Fading + AWGN','SISO in AWGN');
title('Comparison of BER vs. SNR for BPSK Modulation under different transmit diversity configurations');
grid on;