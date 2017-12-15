%this simulation is for IEEE802.11
clear all
%close all

%dur  = 3.2e-6;  
NLOP = 1;        % number of loop
NFFT = 64;          % Number of FFT points
NC   = 48;          % Number of subcarriers
NDS  = 2;           % Number of Data symbol per frame
NS   = NDS*NLOP;    % number of symbols
NP   = 4;           % Number of pilots in symbol –21, –7, 7, and 21
CP   = 16;          % cyclic prefix length
PRE  = 4;           % preamble symbol = 2

%nSUI = 1:1:6;
M = 16;
L = 16;
C = 1*M; %length of computed received samples for Mp
Q = 32;


toff = 31;
tcor = toff+17+3*M;

testcnt = 0;
    
Mp_fail = zeros(50,20);
Mp_miss = zeros(50,20);
Mp_correct = zeros(50,20);

for FOFF = 0,
for nSUI = 0,
for SNR  = 10,
    testcnt = testcnt + 1;
%OFDM TX Create NLOP frames for simulation ================================
    %data
    bit_symbols = round(3*rand(NC, NS));
    %QPSK =================================================================
    QPSK = 2.*mod(bit_symbols,2) - 1 + 1i *(2.*floor(bit_symbols/2)-1);
    QPSK    = QPSK *(1/sqrt(2));   
    dat_mod = QPSK;
    %insert subcarriers & pilots ==========================================
    % pilot ===============================================================
    %Pil = pils(:,1:NS);
%     Pil = zeros(8,NS);
%     symbol = [ zeros(1,NS); QPSK(1  :12, :); ...
%                   Pil(1,:); QPSK(13 :36, :); ...
%                   Pil(2,:); QPSK(37 :60, :); ...
%                   Pil(3,:); QPSK(61 :84, :); ...
%                   Pil(4,:); QPSK(85 :96, :); ...
%                  zeros(NFFT-NC-NP-1,NS); ...
%                             QPSK(97 :108,:); ...    
%                   Pil(5,:); QPSK(109:132,:); ...
%                   Pil(6,:); QPSK(133:156,:); ...
%                   Pil(7,:); QPSK(157:180,:); ...
%                   Pil(8,:); QPSK(181:192,:); ];
    pilots_802_11;
    Pil = repmat(pils(:,1:NDS),1,NLOP);
    symbol = zeros(NFFT,NS);
    symbol(1,:)     = zeros(1,NS);
    symbol(2:7,:)   = dat_mod(1:6, :);
    symbol(8,:)     = Pil(1,NS);
    symbol(9:21,:)  = dat_mod(7:19, :);
    symbol(22,:)    = Pil(2,NS);
    symbol(23:27,:) = dat_mod(20:24, :);
    symbol(39:43,:) = dat_mod(25:29, :);
    symbol(44,:)    = Pil(3,NS);
    symbol(45:57,:) = dat_mod(30:42, :);
    symbol(58,:)    = Pil(4,NS);
    symbol(59:64,:) = dat_mod(43:48, :);
    %IFFT =================================================================
    tx_d =  ifft(symbol, NFFT);
    
    %Add CP ===============================================================
    tx_d = [tx_d(NFFT-CP+1: NFFT,:); tx_d];

    %Add Preamble =========================================================
    
    tx_out = zeros((NFFT+CP), (PRE + NDS)*NLOP);

    preamble_802_11;   
    preamble_nor = [short_pre long_pre]; 
    preamb = reshape(preamble_nor, NFFT+CP, PRE);

    for ii = 0:NLOP -1,
        for jj = 1:PRE,
            tx_out(:,(PRE + NDS)*ii+jj) = preamb(:,jj);
        end
        %tx_out(:,(PRE + NDS)*ii+2) = preamb(:,2);
        if (NDS ~=0 )
            for jj = 1:NDS,
                tx_out(:,(PRE + NDS)*ii+PRE+jj) = tx_d(:,ii*NDS+jj);            
            end
        end
    end
    
    tx_out = reshape(tx_out, 1, (NFFT+CP)*(PRE + NDS)*NLOP);
%==========================================================================   
       
    n=1:(CP+NFFT)*(PRE + NDS);
    freoffs = exp(1i*2*pi*FOFF*(n.'./NFFT));
    %freoffs_cp = [ freoffs(NFFT-CP+1: NFFT); freoffs];
    tx_temp = reshape(tx_out, (CP+NFFT)*(PRE + NDS), NLOP);
    tx_temp = tx_temp .* repmat(freoffs,1,NLOP);   
    tx_out  = reshape(tx_temp,1,length(tx_out));
    
    %SUI channel simulation ===============================================
    if (nSUI ~= 0)
        ['Simulate channel SUI=' num2str(nSUI) '==============================']
        tx_out = reshape(tx_out, (NFFT+CP)*(PRE + NDS),NLOP);
        rx_in  = zeros((NFFT+CP)*(PRE + NDS),NLOP);
        for ii= 1:NLOP,
            channel = channelSUI(nSUI,L/NFFT,20);
            rx_in(:,ii) = SUIChan_SIM(channel(:,1),tx_out(:,ii)')';
        end
        rx_in = reshape(rx_in, 1, (NFFT+CP)*(PRE + NDS)*NLOP);
    else
        ['Simulate channel AWGN with SNR=' num2str(SNR) '=====================']
        rx_in = tx_out;      
    end
    
% Receiver side and simulation of synchronisation==========================
        
    rx_in = reshape(rx_in,  (NFFT+CP)*(PRE + NDS),NLOP);
    % add toff to OFDM frame ==========================================
    rx_in = [zeros(toff,NLOP); rx_in];
    %AWGN channel simulation ==========================================
    rx_in = awgn(rx_in ,SNR,'measured'); 
        
    Flen = toff+ (NFFT+CP)*(PRE + NDS);
    fail_cnt = 0;
    
    Pp = zeros(toff + ((NFFT+CP)*PRE),NLOP);  
    Rp = zeros(toff + ((NFFT+CP)*PRE),NLOP);  
    Mp = zeros(toff + ((NFFT+CP)*PRE),NLOP);

    rx_syn = [zeros(3*M,NLOP); rx_in; zeros(2*M,NLOP)];
    %known_pre = pre64;
    rx_syn = rx_syn ./ max([max(real(rx_syn)) max(imag(rx_syn))]);
    
    ['Calculate Metric']
    abs_pre = abs(short_pre).^2;
    known_coeff = abs_pre;
    %known_coeff = round((abs_pre./max(abs_pre)).*4)./4;
    
    round_coff = 2^Q;
    rx_abs = abs(rx_syn).^2;
    rx_abs_q = fix(rx_abs .* round_coff) ./round_coff;
    
    for d = 1: toff + (NFFT+CP)+M,
        for jj = 1:NLOP,
            Pp(d,jj) = conj(rx_syn(d+(0:C-1),jj)).' * rx_syn(d+M+(0:C-1),jj); 
            Rp(d,jj) = 0;
            for ii = 0:C-1,
                R_q = rx_abs_q(d+M+ii,jj) * known_coeff(1+ii);
                R_q = fix(R_q * round_coff) / round_coff
                Rp(d,jj) = Rp(d,jj) + R_q;
            end
            Pp_Re = abs(real(Pp(d,jj)));
            Pp_Im = abs(imag(Pp(d,jj)));
            if(Pp_Re > Pp_Im),
                Pp_abs = Pp_Re + Pp_Im / 2;
            else
                Pp_abs = Pp_Im + Pp_Re / 2;
            end
            Mp(d,jj) = (Pp_abs / Rp(d,jj));  
        end
    end    
   
    ['Evaluate failure rate']      
    threscnt =0;
    thres_arr = 0.5:0.01:4;
    for thres = 0.5:0.01:4,
        threscnt = threscnt + 1;
        Mp_fail(threscnt,testcnt) = 0;
        Mp_miss(threscnt,testcnt) = 0;
        Mp_correct(threscnt,testcnt) = 0;
        for jj = 1:NLOP,
        mis=1;
        np = 8;
        for ii = 1:tcor +2,
           if ((Mp(ii,jj)> thres)&&(Rp(ii,jj) > 2)),
               np = np -1;
           end
           if(np ==0),
                if (ii > M),
                    avg_peak = Rp(ii:ii+1.25*M-1,jj) + Rp((ii:ii+1.25*M-1)-M,jj);
                else 
                    avg_peak = Rp(ii:ii+1.25*M-1,jj);
                end
                %[val, id] = max(Rp_resh(ii:ii+2*CP-8,jj));
                [val, id] = max(avg_peak);
                if (id == tcor + 1 - ii), Mp_correct(threscnt,testcnt) = Mp_correct(threscnt,testcnt) + 1;
                else           Mp_fail(threscnt,testcnt) = Mp_fail(threscnt,testcnt) + 1;
                end
                mis = 0;
                break;
            end          
        end
        if (mis == 1), Mp_miss(threscnt,testcnt) = Mp_miss(threscnt,testcnt) + 1; end        
        end 
    end
end  
end
end
