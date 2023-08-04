%% compute intermodulation freq

syms w1 w2;

Vin = cos(w1)+cos(w2);
Vout = Vin + Vin^2 + Vin^3;

factor(Vout)