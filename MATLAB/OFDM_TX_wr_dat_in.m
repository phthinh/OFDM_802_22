clear all
close all

%dur  = 3.2e-6;  
NLOP = 2;           % number of loop
NFFT = 2048;        % Number of FFT points
NC   = 1440;        % Number of subcarriers
NDS  = 2;           % Number of Data symbol per frame
NS   = NDS*NLOP;    % number of symbols
NP   = 240;         % Number of pilots in symbol
CP   = (1/4)*NFFT;  % cyclic prefix length
PRE  = 1;           % preamble symbol = 1


% data in for TX ==========================================================
bit_symbols = round(63*rand(1, NC*NS));

Len = NC * NDS;
%write data to file =======================================================
fid = fopen('OFDM_TX_bit_symbols_Len.txt', 'w');
fprintf(fid, '%d ', Len);
fprintf(fid, '%d ', NLOP);
fclose(fid);

fid = fopen('OFDM_TX_bit_symbols.txt', 'w');
fprintf(fid, '%d ', bit_symbols);
fclose(fid);

fid = fopen('RTL_OFDM_TX_bit_symbols.txt', 'w');
fprintf(fid, '%x ', bit_symbols);
fclose(fid);
%write Preamble ===========================================================
preamble_802_22;   
%DL_preamble_nor = DL_preamble ./ max(abs(DL_preamble));
DL_preamble_nor = [short_pre];

Preamble_rtl = DL_preamble_nor .*(2^15);
Preamble_Re  = typecast(int16(real(Preamble_rtl)),'uint16');
Preamble_Im  = typecast(int16(imag(Preamble_rtl)),'uint16');

Pre = uint32(Preamble_Im) * (2^16) + uint32(Preamble_Re);
fid = fopen('../MY_SOURCES/Pre.txt', 'w');
fprintf(fid, '%8x ', Pre);
fclose(fid);

pilots_802_22;
Pilot_seq = reshape(pils, 1, 28*240);
Pilot_seq = (Pilot_seq(1:128)<0)*1;
fid = fopen('../MY_SOURCES/Pilot_seq.txt', 'w');
fprintf(fid, '%d ', Pilot_seq);
fclose(fid);

Alloc_seq = [Al_Vec(2:841,1:5); Al_Vec(1209:2048,1:5)];
Alloc_seq = reshape(Alloc_seq,1,1680*5);
fid = fopen('../MY_SOURCES/Al_vec.txt', 'w');
fprintf(fid, '%d ', Alloc_seq);
fclose(fid);