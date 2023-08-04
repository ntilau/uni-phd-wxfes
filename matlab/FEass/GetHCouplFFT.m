function [Mtrl,Field] = GetHCouplFFT(Mtrl, Field)
%%
F(1:2:2*(Mtrl.nHarms),:) = Field.Sin; 
F(2:2:2*(Mtrl.nHarms),:) = Field.Cos;
nlFFT = zeros(length(Mtrl.Harms),size(Field.Sin,2));
for jh=1:length(Mtrl.Harms)
    nlFFT(jh,:) = fft(Mtrl.func(abs(Field.H).', F, Mtrl.Harms(jh)))/Mtrl.N;
end
% nlFFT( :, abs(nlFFT(1,:).') < 1e-12*max(max(abs(nlFFT(1,:))))) = 0;
% semilogy(nlFFT);
couplFFT = fft( ((nlFFT)*(Mtrl.Cos - 1i*Mtrl.Sin)) .* ...
    (F/Mtrl.N) ,[],2);
Mtrl.Dsin = 1i*(couplFFT(:,Mtrl.posL)-couplFFT(:,Mtrl.posR)).';
Mtrl.Dcos = (couplFFT(:,Mtrl.posL)+couplFFT(:,Mtrl.posR)).';
Mtrl.D(1:2:2*(Mtrl.nHarms),:) = Mtrl.Dsin;
Mtrl.D(2:2:2*(Mtrl.nHarms),:) = Mtrl.Dcos;

