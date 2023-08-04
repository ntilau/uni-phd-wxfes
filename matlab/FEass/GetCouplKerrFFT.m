function D = GetCouplKerrFFT(nHarms, E, a, b, MtrlCx, FieldSx, posR, posL)

N = size(MtrlCx,2);
nlFFT = fft((a+b*(E.'*FieldSx).^2))/N;
% disp(nlFFT)
couplFFT = fft( ones(nHarms,1)*(nlFFT*MtrlCx) .* FieldSx./N ,[],2);
D = -imag(couplFFT(:,posL)-couplFFT(:,posR)).';