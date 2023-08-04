function Mtrl = GetMtrlParamsFFT(Sys)
if ~isfield(Sys,'OverSampling')
    Sys.OverSampling = 4;
end
Mtrl.nHarms = length(Sys.HBharms);
f1 = Sys.freq*Sys.HBharms(1);
f2 = Sys.freq*Sys.HBharms(2);
df = min([f1 f2 abs(round(f1-f2))]);
Mtrl.bf = gcd(df,f1);
N = ceil(max(Sys.HBharms)*Sys.OverSampling*f1/Mtrl.bf);
Mtrl.N = N + mod(N+1,2); % always set the number of FFT bins to be odd
Mtrl.t = linspace(0,1/Mtrl.bf*(Mtrl.N-1)/Mtrl.N,Mtrl.N);
Mtrl.fnl = Mtrl.bf*[0 1:floor(Mtrl.N/2) -floor(Mtrl.N/2):-1];
Mtrl.posL = uint32(Sys.freq*Sys.HBharms/Mtrl.bf)+1;
Mtrl.posR = Mtrl.N - uint32(Sys.freq*Sys.HBharms/Mtrl.bf)+1;



end