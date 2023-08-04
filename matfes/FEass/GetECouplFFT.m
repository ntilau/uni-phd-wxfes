function [Mtrl,Field] = GetECouplFFT(Mtrl, Field)
%%
F(1:2:2*(Mtrl.nHarms),:) = Field.Sin; 
F(2:2:2*(Mtrl.nHarms),:) = Field.Cos;
nlFFT = fft(Mtrl.func(abs(Field.E).'*F))/Mtrl.N;
couplFFT = fft( ones(2*Mtrl.nHarms,1)*((nlFFT)*(Mtrl.Cos - 1i*Mtrl.Sin)) .* ...
    (F/Mtrl.N),[],2);
Mtrl.Dsin = 1i*(couplFFT(:,Mtrl.posL)-couplFFT(:,Mtrl.posR)).';
Mtrl.Dcos = (couplFFT(:,Mtrl.posL)+couplFFT(:,Mtrl.posR)).';
Mtrl.D(1:2:2*(Mtrl.nHarms),:) = Mtrl.Dsin;
Mtrl.D(2:2:2*(Mtrl.nHarms),:) = Mtrl.Dcos;
