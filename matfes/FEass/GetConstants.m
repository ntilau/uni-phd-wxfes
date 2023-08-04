function Sys = GetConstants(Sys)
Sys.c0 = 299792458;
Sys.z0 = 120*pi;
Sys.eps0 = 1/(Sys.z0*Sys.c0);
Sys.mu0 = Sys.z0/Sys.c0;
Sys.db = @(x)20*log10(abs(x));
Sys.arg = @(x)unwrap(angle(x))*180/pi;
end
