function [SER,legendtext] = MRC(L)

if nargin < 1
    L = 2;
end

% Find the SER for SIMO case with L channel outputs(antennaes),
% decoding done by Maximal-Ratio Combining(MRC),combining signals from all the receive antennae 
%
% Input:
% L = Row vector with each value tells number of receive antennaes
% eg. L = [1 3] or [1:4]
%
% Output:
% SER = A LxNumOfSnrPoints matrix with each row corresponds to particular L outputs
%
% Plot:
% A plot corresponding to each value in vector L

EbNo_dB = [-5 : 0.5 : 15]; %Array of SNR values used for the simuation (dB)
EbNo = 10.^(EbNo_dB/10);
sd = sqrt(1./EbNo);  % Standard Deviation
var = sd.^2;        % Variance

% Number of Receive antenae to simulate
%L = 2;
% Number of Monte carlo iteration
T = 10 ^ 6;

%% STEP 1 OF MONTE CARLO SIMULATION
%% Generating transmit data
x_bpsk = sign(randn(1, T));

SER = zeros(length(L),length(sd));

for p = 1 : length(L)
    noise = zeros(L(p),T);
    h_channel = zeros(L(p),T);
    y = zeros(L(p),T);
    y_new = zeros(1,T);
    x_hat = zeros(1,T);
    
    for k = 1 : length(sd)
      %% STEP 2 OF M-C SIMULATION
      % Generating COMPLEX RANDOM NOISE(LxT) as L rx chains MEAN 0 & VAR = N0 = 1/SNR
      noise = (sqrt(L)) * ((1/sqrt(2)) * sd(1,k) * randn(L(p),T) + 1i * (1/sqrt(2)) * sd(1,k) * randn(L(p),T));
    
      % Generating COMPLEX RANDOM CHANNEL(LxT) WITH MEAN 0 & VAR = 1
      h_channel = (1/sqrt(2)) * randn(L(p),T) + 1i * (1/sqrt(2)) * randn(L(p),T);
    
      %%STEP 3 OF M-C SIMULATION
      %% RECEIVED SIGNAL ACCORDING TO RELATION y(LxT) = hx + n
      y = h_channel .* repmat(x_bpsk,L(p),1) + noise;
    
      %% STEP 4 OF M-C SIMULATION
      %% Using the detector threshold to do detection and get x_hat
      %% On receiver we do (h/|h|) dot with y, projection, Matched Filter, MRC
      %% Matched Filter/MRC step(Scalar sufficient statistic, dot product considering cols as vectors)
      y_new = dot(h_channel./repmat(sqrt(sum(abs(h_channel).^2,1)),L(p),1),y,1);    
      %y_new = conj(h_channel./abs(h_channel)).*y;
    
      x_hat = sign(real(y_new));    %% since BPSK, so taking real sufficient statistic
            
      %% STEP 5 OF M-C SIMULATION
      %% Evaluate I
      SER(p,k) = mean(x_hat~=x_bpsk);
    end
end

%figure(1)
%plot(abs(y(1,:)))
%plot(real(y(30,:)),imag(y(30,:)),'o')
%%PLOTTING
figure(1);
semilogy(EbNo_dB,SER,'b-s','linewidth',2);
grid on
axis([min(EbNo_dB) max(EbNo_dB) 1e-6 1])
title('BPSK - SIMO(MRC): SER Simulation with varying Receive antennae')
xlabel('signal-to-noise ratio (SNR) [dB]')
ylabel('symbol error rate (SER)')
legendtext ='';
for p = 1:length(L)
    legendtext = [legendtext; sprintf('MRC L = %d ',L(p))];
end
legend(legendtext);

ber_array_MRC_1x2 = SER;
save('SIMO_1x2_BPSK.mat','EbNo_dB','ber_array_MRC_1x2');