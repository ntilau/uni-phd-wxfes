function u = ThermalFE(t,u1)
%dydt = [y(2); (1-y(1)^2)*y(2)-y(1)];

global Sys;

A = Sys.T;
b = (Sys.q0 * Sys.f - Sys.K * Sys.S * u1) / Sys.RhoCp;

if isfield(Sys,'Dir')
    for ibc = 1:length(Sys.Dir)
        b = b - A(:,Sys.Dir{ibc}) * ones(length(Sys.Dir{ibc}),1) * Sys.Tedge{ibc} ;
    end
    for ibc = 1:length(Sys.Dir)
        b(Sys.Dir{ibc}) = Sys.Tedge{ibc};
        A(Sys.Dir{ibc},:) = 0; A(:,Sys.Dir{ibc}) = 0;
        A(Sys.Dir{ibc},Sys.Dir{ibc}) = eye(length(Sys.Dir{ibc}));
    end
end
u = A\b;
%u1 = u;
    
