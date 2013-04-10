S277='C56F36BB65B724B8E5E8D6137C4AF1942307BF5AB264770B41B00';
S488='203805FF2AB99A227875F4D4ECE9163C851F3D4530C410FC15030';
S277_seq =[];
S488_seq =[];
for ii = 1:length(S277),
    S277_seq = [S277_seq dec2bin(hex2dec(S277(ii)),4)];
    S488_seq = [S488_seq dec2bin(hex2dec(S488(ii)),4)];
end
S277_bvec = S277_seq.*1 - 48;
S277_bvec = S277_bvec.*2 -1;
S488_bvec = S488_seq.*1 - 48;
S488_bvec = S488_bvec.*2 -1;

short_sym = zeros(NFFT,1);
k=(4:4:840);
short_sym(1+k) = S488_bvec(1+((k-4)./4));  
short_sym(2048-843+k) = S277_bvec(1+((k-4)./4));  

short_sym = sqrt(1680/420)*short_sym; % multiply with normalized factor.

STS = ifft(short_sym,NFFT);
short_pre = [STS(NFFT-CP+1:NFFT); STS];


S536='F1C4677539900F45F5E42A3418663A12B8F6C1081350487D8D55D344BACF02CD9C9BCD68C4932A67D2AC0473878B1F970A2A938DF';
S115='A877F40C94889D20B91E7FB49616CB714A17845A62EE00A795947CC27EFBBD3E32F5B7E0FE2607056F6669D872C8A0376E8ED764F';
S536_seq =[];
S115_seq =[];
for ii = 1:length(S536),
    S536_seq = [S536_seq dec2bin(hex2dec(S536(ii)),4)];
    S115_seq = [S115_seq dec2bin(hex2dec(S115(ii)),4)];
end
S536_bvec = S536_seq.*1 - 48;
S536_bvec = S536_bvec.*2 -1;
S115_bvec = S115_seq.*1 - 48;
S115_bvec = S115_bvec.*2 -1;

long_sym = zeros(NFFT,1);
k=(2:2:840);
long_sym(1+k) = S536_bvec(1+((k-2)./2));  
long_sym(2048-841+k) = S115_bvec(1+((k-2)./2));  
long_sym = sqrt(1680/840)*long_sym; % multiply with normalized factor.

LTS = ifft(long_sym,NFFT);
long_pre = [LTS(NFFT-CP+1:NFFT); LTS];