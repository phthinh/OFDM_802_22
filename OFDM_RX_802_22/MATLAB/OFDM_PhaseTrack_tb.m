close all

NLOP = 1;    % number of loop
NFFT = 2048;        % Number of FFT points
NC   = 1440;        % Number of subcarriers
NDS  = 2;           % Number of Data symbol per frame
NS   = NDS*NLOP;    % number of symbols
NP   = 240;         % Number of pilots in symbol
CP   = (1/4)*NFFT;  % cyclic prefix length
PRE  = 2;           % preamble symbol = 2

% Read data in ============================================================
Para_fid = fopen('RTL_PhaseTrack_datin_len.txt', 'r');
Para = fscanf(Para_fid, '%d ');
Flen  = Para(1);
fclose(Para_fid);

% NDS = Flen/(NFFT); %number of Data symbol excluding preamble
NDS = Flen/(NC+NP); %number of Data symbol excluding preamble

datin_fid = fopen('PhaseTrack_Txsym_Re.txt', 'r');
Txsym_Re = fscanf(datin_fid, '%f ');
fclose(datin_fid);
datin_fid = fopen('PhaseTrack_Txsym_Im.txt', 'r');
Txsym_Im = fscanf(datin_fid, '%f ');
fclose(datin_fid);
Txsym = Txsym_Re + 1i*Txsym_Im;
Txsym = Txsym.';
Txsym = reshape(Txsym, NFFT, NDS);


datin_fid = fopen('PhaseTrack_datin_Re.txt', 'r');
dat_Re = fscanf(datin_fid, '%f ');
fclose(datin_fid);

datin_fid = fopen('PhaseTrack_datin_Im.txt', 'r');
dat_Im = fscanf(datin_fid, '%f ');
fclose(datin_fid);

PhaseTrack_datin =  dat_Re + 1i*dat_Im;
PhaseTrack_datin = PhaseTrack_datin.';
% PhaseTrack_datin = reshape(PhaseTrack_datin, NFFT, NDS);
PhaseTrack_datin = reshape(PhaseTrack_datin, (NC+NP), NDS);

% Read data out of RTL ====================================================
datout_fid = fopen('RTL_PhaseTrack_datout_Re.txt', 'r');
PhaseTrack_datout_Re_rtl = fscanf(datout_fid, '%d ');
fclose(datout_fid);
datout_fid = fopen('RTL_PhaseTrack_datout_Im.txt', 'r');
PhaseTrack_datout_Im_rtl = fscanf(datout_fid, '%d ');
fclose(datout_fid);
PhaseTrack_datout_rtl = (PhaseTrack_datout_Re_rtl./2^6) + 1i*(PhaseTrack_datout_Im_rtl./2^6);

% Simulate with data in ===================================================
Txsym = reshape(Txsym,NFFT,NDS);
Txsym(842:1208,:)=[];
Txsym(1,:) = [];
Txsym = reshape(Txsym,1,(NC+NP)*NDS);

pilots_802_22;
% pils_mod = pils.*2-1;
pils_mod = repmat([1; -1],120,1);
Pil_ref = repmat(pils_mod([1:64 177:240] ,1),1,NDS);
% Pil_ind = [(2:7:443) (1601:7:2042)].';
Pil_ind = [((2:7:443)-1) ((1601:7:2042)-368)].';
Rx_Pil = PhaseTrack_datin(Pil_ind,:);
Ph_est = 1/128*(sum(Rx_Pil .* Pil_ref));
% Ph_comp = ones(NFFT,1) * conj(Ph_est);
Ph_comp = ones((NC+NP),1) * conj(Ph_est);
PhaseTrack_datout = PhaseTrack_datin .* Ph_comp;

% PhaseTrack_datin(842:1208,:)=[];
% PhaseTrack_datin(1,:) = [];
PhaseTrack_datin = reshape(PhaseTrack_datin,1,(NC+NP)*NDS);

% PhaseTrack_datout(842:1208,:)=[];
% PhaseTrack_datout(1,:) = [];
PhaseTrack_datout_sim = reshape(PhaseTrack_datout,1,(NC+NP)*NDS);
% Compare Simulation vs RTL ===============================================

figure(1)
hold on 
plot(1:length(PhaseTrack_datout_sim), angle(PhaseTrack_datout_sim),'o-b');
plot(1:length(PhaseTrack_datin), angle(PhaseTrack_datin),'.-r');
title ('PhaseTrack\_datout\_sim vs PhaseTrack\_datout\_in')
legend('PhaseTrack\_datout\_sim','PhaseTrack\_datout\_in')
xlim([1 1000]);
hold off

figure(2)
hold on 
plot(1:length(PhaseTrack_datout_sim), angle(PhaseTrack_datout_sim),'.-r');
plot(1:length(Txsym), angle(Txsym),'o-b');
title ('PhaseTrack\_datout\_sim vs Tx\_sym')
legend('PhaseTrack\_datout\_sim','Tx\_sym')
xlim([1 1000]);
hold off

figure(3)
hold on 
plot(1:length(PhaseTrack_datout_sim), angle(PhaseTrack_datout_sim),'o-b');
plot(1:length(PhaseTrack_datout_rtl), angle(PhaseTrack_datout_rtl),'.-r');
title ('PhaseTrack\_datout\_sim vs PhaseTrack\_datout\_rtl')
legend('PhaseTrack\_datout\_sim','PhaseTrack\_datout\_rtl')
ylim([-4 4]);
hold off
