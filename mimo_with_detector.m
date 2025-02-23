clear
N = 10^6; 
Eb_N0_dB = [0:25]; 
nTx = 2;
nRx = 2;
for ii = 1:length(Eb_N0_dB)
 
    ip = rand(1,N)>0.5; 
    s = 2*ip-1; 
 
    sMod = kron(s,ones(nRx,1)); 
    sMod = reshape(sMod,[nRx,nTx,N/nTx]); 
 
    h = 1/sqrt(2)*[randn(nRx,nTx,N/nTx) + j*randn(nRx,nTx,N/nTx)]; 
    n = 1/sqrt(2)*[randn(nRx,N/nTx) + j*randn(nRx,N/nTx)]; 
 
    
    y = squeeze(sum(h.*sMod,2)) + 10^(-Eb_N0_dB(ii)/20)*n;
 
    hCof = zeros(2,2,N/nTx)  ; 
    hCof(1,1,:) = sum(h(:,2,:).*conj(h(:,2,:)),1);  
    hCof(2,2,:) = sum(h(:,1,:).*conj(h(:,1,:)),1);  
    hCof(2,1,:) = -sum(h(:,2,:).*conj(h(:,1,:)),1); 
    hCof(1,2,:) = -sum(h(:,1,:).*conj(h(:,2,:)),1); 
    hDen = ((hCof(1,1,:).*hCof(2,2,:)) - (hCof(1,2,:).*hCof(2,1,:)));
    hDen = reshape(kron(reshape(hDen,1,N/nTx),ones(2,2)),2,2,N/nTx);  
    hInv = hCof./hDen; % inv(H^H*H)
 
    hMod =  reshape(conj(h),nRx,N); 
    
    yMod = kron(y,ones(1,2));
    yMod = sum(hMod.*yMod,1); 
    yMod =  kron(reshape(yMod,2,N/nTx),ones(1,2)); 
    yHat = sum(reshape(hInv,2,N).*yMod,1); 
   
    
    ipHat = real(yHat)>0;
 
    
    nErr(ii) = size(find([ip- ipHat]),2);
 
end
 
simBer = nErr/N; % simulated ber
EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbN0Lin).^(-0.5)); 
p = 1/2 - 1/2*(1+1./EbN0Lin).^(-1/2);
theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p)); 
 
close all
figure
semilogy(Eb_N0_dB,theoryBer_nRx1,'bp-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,theoryBerMRC_nRx2,'kd-','LineWidth',2);
semilogy(Eb_N0_dB,simBer,'mo-','LineWidth',2);
axis([0 25 10^-5 0.5])
grid on
legend('theory (nTx=1,nRx=1)', 'theory (nTx=1,nRx=2, MRC)', 'sim (nTx=2, nRx=2, ZF)');
xlabel('Average Eb/No,dB');
ylabel('Bit Error Rate');
title('BER for BPSK modulation with 2x2 MIMO and ZF equalizer (Rayleigh channel)');
