function D = GetCouplSinFFT(nHarms, E, a, b, MtrlCx, FieldSx, posR, posL)

N = size(MtrlCx,2);
nlFFT = fft((a+b*(E.'*FieldSx).^2))/N;
% figure; plot(1:N,real(nlFFT),'+',1:N,imag(nlFFT),'*');
nl = real(nlFFT(1:size(MtrlCx,1)));
idx = find( abs(nl) > 1e-13*max(max(abs(nl))) );
nl(2:end) = nl(2:end)*2;
Nx1 = ones(nHarms,1)*(nl(idx)*MtrlCx(idx,:));
SxN = FieldSx./N;
NxSx = Nx1 .* SxN;
% plot(NxSx.')
couplFFT = fft( NxSx ,[],2);
% figure; plot(1:N,real(couplFFT(1,:)),'+',1:N,imag(couplFFT(1,:)),'*');
D = 1i*(couplFFT(:,posL)-couplFFT(:,posR)).';
