function D = GetCouplFFT(Mtrl, Field)
%%
F(1:2:2*(Mtrl.nHarms),:) = Field.Sin; 
F(2:2:2*(Mtrl.nHarms),:) = Field.Cos;
E(1:2:2*(Mtrl.nHarms),:) = Field.E;
E(2:2:2*(Mtrl.nHarms),:) = Field.E;
nlFFT = fft(Mtrl.func(E.'*F))/Mtrl.N;
couplFFT = fft( ones(2*Mtrl.nHarms,1)*(real(nlFFT)*Mtrl.Cos - imag(nlFFT)*Mtrl.Sin) .* ...
    (F/Mtrl.N) ,[],2);
D = (couplFFT(:,Mtrl.posL)-couplFFT(:,Mtrl.posR)).';
% rD = real(D)
D = -imag(D)