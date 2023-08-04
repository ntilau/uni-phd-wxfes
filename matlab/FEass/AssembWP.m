function Sys =  AssembWP(Sys, freq)

if ~isfield(Sys,'nnWP')
    Sys.nnWP = 1:Sys.NDOFs;
    remId = [];
    for ibc=1:length(Sys.Dir)
        remId = [remId; Sys.Dir{ibc}];
    end
    for ibc=1:length(Sys.WP)
        remId = [remId; Sys.WP{ibc}];
    end
    Sys.nnWP(remId) = [];
end

k0 = 2*pi*freq/Sys.c0;

nRHS = length(Sys.WP)*Sys.WPnModes;
Sys.A = Sys.S - k0^2*Sys.T;
Sys.B = sparse(length(Sys.nnWP)+nRHS,nRHS);
PP = zeros(length(Sys.WP)*Sys.WPnModes,length(Sys.WP)*Sys.WPnModes);
IP = zeros(length(Sys.nnWP),length(Sys.WP)*Sys.WPnModes);
for ip=1:length(Sys.WP)
    if isfield(Sys,'Einc')
        Sys.WPpowEq = diag( (Sys.Einc).^2 * Sys.Height * ...
            ( sqrt(1-(Sys.WPfc{ip}./freq).^2) ) / (Sys.z0) );
    else
        Sys.WPpowEq = Sys.WPpow * eye(Sys.WPnModes) * 2 / Sys.Height;
    end
    gamma = diag(1i*2*pi*freq/Sys.c0*sqrt(1-(Sys.WPfc{ip}.'/freq).^2));
    gamma = abs(real(gamma)) + 1i*imag(gamma);
    Sys.WPgvec{ip} = Sys.WPvec{ip}*(sqrt(1i*k0*Sys.z0*(gamma\Sys.WPpowEq)));
    PP((1:Sys.WPnModes)+(ip-1)*Sys.WPnModes,(1:Sys.WPnModes)+(ip-1)*Sys.WPnModes) = ...
        Sys.WPgvec{ip}.'* Sys.A(Sys.WP{ip},Sys.WP{ip})*Sys.WPgvec{ip} +...
        1i*k0*Sys.z0*Sys.WPpowEq*eye(Sys.WPnModes);
    IP(:,(1:Sys.WPnModes)+(ip-1)*Sys.WPnModes) = Sys.A(Sys.nnWP,Sys.WP{ip})*...
        Sys.WPgvec{ip};
    Sys.B((1:Sys.WPnModes)+(ip-1)*Sys.WPnModes,...
            (1:Sys.WPnModes)+(ip-1)*Sys.WPnModes) = ...
            1i*k0*Sys.z0*2*Sys.WPpowEq*eye(Sys.WPnModes);
end
PI = IP.';
II = Sys.A(Sys.nnWP,Sys.nnWP);
Sys.A = [PP, PI; IP, II];
end