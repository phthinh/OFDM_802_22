% Generate allocation vector and Pillot vector for 28 OFDM symbols
P_pattern=[ 1 2 2 2 2 2 2;...
            2 2 2 1 2 2 2;...
            2 2 2 2 2 1 2;...
            2 1 2 2 2 2 2;...
            2 2 2 2 1 2 2;...
            2 2 2 2 2 2 1;...
            2 2 1 2 2 2 2];

Al_Vec = [zeros(28,1) repmat(P_pattern, 4,120) zeros(28,367) repmat(P_pattern,4,120)].';

P_seed = [0 1 1 0 1 1 1 0 0 0 1 0 1 0 1]; % pilot seed, MSB on left
pils = zeros(28,240);
for ii = 1:28,
    for jj = 1:240,
        pils(ii,jj)= xor(P_seed(1),P_seed(2));
        P_seed = [P_seed(2:15) pils(ii,jj)];
    end
end
pils = pils.';