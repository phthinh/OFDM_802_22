clear all
close all

NLOP = 2;    % number of loop
NFFT = 2048;        % Number of FFT points
NC   = 1440;        % Number of subcarriers
NDS  = 2;           % Number of Data symbol per frame
NS   = NDS*NLOP;    % number of symbols
NP   = 240;         % Number of pilots in symbol
CP   = (1/4)*NFFT;  % cyclic prefix length
PRE  = 2;           % preamble symbol = 2

N = 2048;
M = N/4;

SNR = 30;
FOFF = 0;
toff = 32;
tcor = toff+33+3*M;

%OFDM TX Create NLOP frames for simulation ================================
%data
bit_symbols = round(3*rand(NC, NS));

%QPSK =================================================================
QPSK    = 2.*mod(bit_symbols,2)-1 + 1i *(2.*floor(bit_symbols/2)-1);
QPSK    = QPSK *(1/sqrt(2));   
dat_mod = QPSK;

% insert subcarriers & pilots =============================================
% pilots in symbol –21, –7, 7, and 21======================================

pilots_802_22;
% pils_mod = pils.*2-1;
pils_mod = repmat([1; -1],120,1);

symbol = zeros(NFFT,NS);

for nn = 0: NLOP-1,
    for ii = 1:NDS,
    dat_cnt = 1;
    pil_cnt = 1;
    for jj =1:NFFT,
%         if (Al_Vec(jj,ii) == 1),      
        if (Al_Vec(jj,1) == 1), % only take the first allocation vector
%             symbol(jj,ii+NDS*nn) = pils_mod(pil_cnt,ii);
            symbol(jj,ii+NDS*nn) = pils_mod(pil_cnt,1);
            pil_cnt       = pil_cnt +1;
%         elseif(Al_Vec(jj,ii) == 2),
        elseif(Al_Vec(jj,1) == 2),
            symbol(jj,ii+NDS*nn) = dat_mod(dat_cnt,ii+NDS*nn);
            dat_cnt       = dat_cnt +1;
        end
    end
    end
end

%IFFT =================================================================
tx_d =  ifft(symbol, NFFT, 1);

%Add CP ===============================================================
tx_d = [tx_d(NFFT-CP+1: NFFT,:); tx_d];

%Add Preamble =========================================================
tx_out = zeros((NFFT+CP), (PRE + NDS)*NLOP);

preamble_802_22;   
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
tx_out = reshape(tx_out, (NFFT+CP)*(PRE + NDS)*NLOP,1);
%==========================================================================   

%frequency offset adding ==============================================
n=0:(CP+NFFT)*(PRE + NDS)-1;
freoffs = exp(1i*2*pi*FOFF*(n.'./NFFT));    
tx_temp = reshape(tx_out, (CP+NFFT)*(PRE + NDS), NLOP);
tx_temp = tx_temp .* repmat(freoffs,1,NLOP);   
tx_out  = reshape(tx_temp,1,length(tx_out));

%AWGN channel simulation ==============================================
%rx_in = tx_out;  
rx_in = reshape(tx_out,(CP+NFFT)*(PRE + NDS), NLOP);
%rx_in = [rx_in(length(rx_in)- toff + 1 : length(rx_in)) rx_in(1:length(rx_in))];
toff_mat = zeros(toff,NLOP);
rx_in = [toff_mat; rx_in];
rx_in = reshape(rx_in,1,((CP+NFFT)*(PRE + NDS) + toff) * NLOP);
rx_in = awgn(rx_in ,SNR,'measured');   
rx_in = 0.5*(rx_in ./ max([max(real(rx_in)) max(imag(rx_in))]));
%rx_in = rx_in .*2;

%write data to file =======================================================
fid = fopen('OFDM_RX_bit_symbols.txt', 'w');
fprintf(fid, '%d ', bit_symbols);
fclose(fid);

Len = length(rx_in);
fid = fopen('OFDM_RX_datin_Re.txt', 'w');
fprintf(fid, '%f ', real(rx_in));
fclose(fid);
fid = fopen('OFDM_RX_datin_Im.txt', 'w');
fprintf(fid, '%f ', imag(rx_in));
fclose(fid);

datin_rtl = rx_in(1:Len) .*(2^15);
datin_Re = typecast(int16(real(datin_rtl)),'uint16');
datin_Im = typecast(int16(imag(datin_rtl)),'uint16');

SNR_w = round(SNR);
if (SNR >15), SNR_w = 15; end
Flen = toff + (NFFT+CP) *(PRE+NDS);
fid = fopen('RTL_OFDM_RX_datin_len.txt', 'w');
fprintf(fid, '%d %d %d %d', NLOP, Flen, SNR_w, toff);
fclose(fid);
fid = fopen('RTL_OFDM_RX_datin_Re.txt', 'w');
fprintf(fid, '%4x ', datin_Re);
fclose(fid);
fid = fopen('RTL_OFDM_RX_datin_Im.txt', 'w');
fprintf(fid, '%4x ', datin_Im);
fclose(fid);

% data cofficient for system synthesis ==================================== 

known_coeff = 2*(imag(short_pre((1:1*M)))<0) + 1*(real(short_pre((1:1*M))<0));
known_coeff_rtl = typecast(int8(known_coeff),'uint8');
fid = fopen('../MY_SOURCES/Synch_known_coeff_802_22.txt', 'w');
fprintf(fid, '%x ', known_coeff_rtl);
fclose(fid);