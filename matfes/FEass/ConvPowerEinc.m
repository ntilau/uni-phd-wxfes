% a = 1.651e-3;
% b = 0.8255e-3;
a = 22.86e-03;
b = 22.86e-03/2;
% a = 19.05e-3;
% b = a/2;

Einc = 5e3:5e3:15e3;
c0 = 299792458;
Z0 = 376.730313461;
f = 10e9; 
fc = c0/(2*a); % TE10
Pinc = Einc.^2 * ( a.*b.*sqrt(1-(fc./f).^2) ) / (4.*Z0);
Einc = sqrt(1./(( a.*b.*sqrt(1-(fc./f).^2) ) / (4.*Z0)))
format longe
Pinc