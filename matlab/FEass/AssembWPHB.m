function Sys =  AssembWPHB(Sys)

if ~isfield(Sys,'nnWP')
    Sys.nnWP = 1:Sys.NDOFs*Sys.nHarms;
    remId = [];
    for ibc=1:length(Sys.Dir)
        remId = [remId; Sys.Dir{ibc}];
    end
    for ibc=1:length(Sys.WP)
        for jh=1:Sys.nHarms
            remId = [remId; Sys.WP{ibc}+Sys.NDOFs*(jh-1)];
        end
    end
    Sys.nnWP(remId) = [];
end

k0 = 2*pi*Sys.freq/Sys.c0;
nRHS = length(Sys.WP)*Sys.WPnModes*Sys.nHarms;
Sys.A = Sys.S - k0^2*Sys.T;
Sys.B = sparse(length(Sys.nnWP)+nRHS,nRHS);
PP = zeros(length(Sys.WP)*Sys.WPnModes*Sys.nHarms,length(Sys.WP)*Sys.WPnModes*Sys.nHarms);
IP = zeros(length(Sys.nnWP),length(Sys.WP)*Sys.WPnModes*Sys.nHarms);
for ip=1:length(Sys.WP)
    for jh=1:Sys.nHarms
        if isfield(Sys,'Einc')
            Sys.WPpowEq = diag( (Sys.Einc).^2 * Sys.Height * ...
                (sqrt(1-(Sys.WPfc{ip}./...
                (Sys.Harms(jh)*Sys.freq)).^2) ) / (Sys.z0));
        else
            Sys.WPpowEq = Sys.WPpow * eye(Sys.WPnModes) * 2 / Sys.Height;
        end
        gamma = diag(1i*2*pi*(Sys.Harms(jh)*Sys.freq)/Sys.c0*...
            sqrt(1-(Sys.WPfc{ip}.'/(Sys.Harms(jh)*Sys.freq)).^2));
        gamma = abs(real(gamma)) + 1i*imag(gamma);
        Sys.WPgvec{ip,jh} = Sys.WPvec{ip}*...
            (sqrt(1i*(Sys.Harms(jh)*k0)*Sys.z0*(gamma\Sys.WPpowEq)));
        PP((1:Sys.WPnModes)+((jh-1)+(ip-1)*Sys.nHarms)*Sys.WPnModes,...
            (1:Sys.WPnModes)+((jh-1)+(ip-1)*Sys.nHarms)*Sys.WPnModes) = ...
            Sys.WPgvec{ip,jh}.'* ...
            Sys.A(Sys.WP{ip}+Sys.NDOFs*(jh-1),Sys.WP{ip}+Sys.NDOFs*(jh-1))*...
            Sys.WPgvec{ip,jh} +...
            1i*Sys.Harms(jh)*k0*Sys.z0*Sys.WPpowEq*eye(Sys.WPnModes);
        IP(:,(1:Sys.WPnModes)+((jh-1)+(ip-1)*Sys.nHarms)*Sys.WPnModes) = ...
            Sys.A(Sys.nnWP,Sys.WP{ip}+Sys.NDOFs*(jh-1))*Sys.WPgvec{ip,jh};
        Sys.B((1:Sys.WPnModes)+((jh-1)+(ip-1)*Sys.nHarms)*Sys.WPnModes,...
            (1:Sys.WPnModes)+((jh-1)+(ip-1)*Sys.nHarms)*Sys.WPnModes) = ...
            1i*Sys.Harms(jh)*k0*Sys.z0*2*Sys.WPpowEq*eye(Sys.WPnModes);
    end
end
PI = IP.';
II = Sys.A(Sys.nnWP,Sys.nnWP);
Sys.A = [PP, PI; IP, II];

SizScat = length(Sys.WP)*Sys.WPnModes*Sys.nHarms;

% Sys.B(1,1) = Sys.B(1,1) * Sys.Pfund;
% Sys.A(1,1) = Sys.A(1,1) * Sys.Pfund;
% Sys.A(SizScat+1:end,1) = Sys.A(SizScat+1:end,1) * sqrt(Sys.Pfund);
% Sys.A(1,SizScat+1:end) = Sys.A(1,SizScat+1:end) * sqrt(Sys.Pfund);
% 
% Sys.B(Sys.WPnum,1) = Sys.B(Sys.WPnum,Sys.WPnum) * Sys.Pitrf;
% Sys.A(Sys.WPnum,Sys.WPnum) = Sys.A(Sys.WPnum,Sys.WPnum) * Sys.Pitrf;
% Sys.A(SizScat+1:end,Sys.WPnum) = Sys.A(SizScat+1:end,Sys.WPnum) * sqrt(Sys.Pitrf);
% Sys.A(Sys.WPnum,SizScat+1:end) = Sys.A(Sys.WPnum,SizScat+1:end) * sqrt(Sys.Pitrf);
end