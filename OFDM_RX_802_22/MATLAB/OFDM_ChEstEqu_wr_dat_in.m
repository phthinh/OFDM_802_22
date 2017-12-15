clear all
close all

NLOP = 1;    % number of loop
NFFT = 2048;        % Number of FFT points
NC   = 1440;        % Number of subcarriers
NDS  = 2;           % Number of Data symbol per frame
NS   = NDS*NLOP;    % number of symbols
NP   = 240;         % Number of pilots in symbol
CP   = (1/4)*NFFT;  % cyclic prefix length
PRE  = 2;           % preamble symbol = 2

SNR = 200;
FOFF = 0;
toff = 0;


%OFDM TX Create NLOP frames for simulation ================================
%data
bit_symbols = round(3*rand(NC, NS));

%QPSK =================================================================
QPSK    = 2.*mod(bit_symbols,2)-1 + 1i *(2.*floor(bit_symbols/2)-1);
QPSK    = QPSK *(1/sqrt(2));   
dat_mod = QPSK;

%insert subcarriers & pilots ==========================================
% pilot ===============================================================
pilots_802_22;
pils_mod = pils.*2-1;
symbol = zeros(NFFT,NS);

for nn = 0: NLOP-1,
    for ii = 1:NDS,
    dat_cnt = 1;
    pil_cnt = 1;
    for jj =1:NFFT,
        if (Al_Vec(jj,ii) == 1),
            symbol(jj,ii+NDS*nn) = pils_mod(pil_cnt,ii);
            pil_cnt       = pil_cnt +1;
        elseif(Al_Vec(jj,ii) == 2),
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
%rx_in = 0.5*(rx_in ./ max([max(real(rx_in)) max(imag(rx_in))]));

%receive and remove CP ==============================================
rx_resh = reshape(rx_in,((CP+NFFT)*(PRE + NDS) + toff), NLOP);
rx = rx_resh(1:(CP+NFFT)*(PRE + NDS),1);
rx = reshape(rx,(CP+NFFT),(PRE + NDS));
rx(1:CP,:) = [];
%fft ======================================================================
rx_sym = fft(rx,NFFT,1);
ch_sym = rx_sym(:,PRE:(PRE + NDS));

%add channel response =====================================================
ch_ph = rand(NFFT,1);
ch_res = exp(1i*2*pi.*ch_ph);
ch_res([2:2:840, 1210:2:2048]) = ch_res([3:2:841, 1209:2:2047]);
ch_res= repmat(ch_res,1,1+NDS);
rx_ch_sym = ch_sym .* ch_res;
% rx_ch_sym = ch_sym;
rx_ch_sym = reshape(rx_ch_sym,NFFT*(1+NDS),1);

%write data to file =======================================================
ch_sym = reshape(ch_sym,NFFT*(1+NDS),1);
fid = fopen('ChEstEqu_sym_Re.txt', 'w');
fprintf(fid, '%f ', real(ch_sym));
fclose(fid);
fid = fopen('ChEstEqu_sym_Im.txt', 'w');
fprintf(fid, '%f ', imag(ch_sym));
fclose(fid);


fid = fopen('ChEstEqu_datin_Re.txt', 'w');
fprintf(fid, '%f ', real(rx_ch_sym));
fclose(fid);
fid = fopen('ChEstEqu_datin_Im.txt', 'w');
fprintf(fid, '%f ', imag(rx_ch_sym));
fclose(fid);


Len = length(rx_ch_sym);
datin_rtl = rx_ch_sym(1:Len) .*(2^11);
datin_Re = typecast(int16(real(datin_rtl)),'uint16');
datin_Im = typecast(int16(imag(datin_rtl)),'uint16');

Flen = NFFT*(1+NDS);
fid = fopen('RTL_ChEstEqu_datin_len.txt', 'w');
fprintf(fid, '%d', Flen);
fclose(fid);
fid = fopen('RTL_ChEstEqu_datin_Re.txt', 'w');
fprintf(fid, '%4x ', datin_Re);
fclose(fid);
fid = fopen('RTL_ChEstEqu_datin_Im.txt', 'w');
fprintf(fid, '%4x ', datin_Im);
fclose(fid);    

% data cofficient for system synthesis ==================================== 

known_coeff = 2*(imag(long_sym([3:2:841, 1209:2:2047]))<0) + 1*(real(long_sym([3:2:841, 1209:2:2047]))<0);
known_coeff_rtl = typecast(int8(known_coeff),'uint8');
fid = fopen('../MY_SOURCES/ChEstEqu_lpre.txt', 'w');
fprintf(fid, '%x ', known_coeff_rtl);
fclose(fid);