clear all
warning off
graphics_toolkit gnuplot
%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%
Infty=50;
nBitsPerFrame=2000;
max_nFrame =2000;
fErrLim=100;
snr_db=0:2:20;
Eb_n0=zeros(size(snr_db));  %Initializing Eb/N0
%%%%% SYMBOL CONSTELLATION %%%%%%%%%
tb2=[0;1];
nBitPerSym=1;
M=2^nBitPerSym;
symbolBook=[1 1j];
bitBook=tb2;
%%%%%%%%%%%%%%ERROR MATRICES%%%%%%%%%%%%%%%%%%%%
nSymPerFrame=nBitsPerFrame/nBitPerSym;
errbits=zeros(length(snr_db), 1);
nFrames=zeros(length(snr_db), 1);
fErrs=errbits;
errsyms=zeros(length(snr_db), 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nEN = 1:length(snr_db) % SNR POINTS
    snr_p=snr_db(nEN);
    en = 10^(snr_p/10); % convert SNR from unit db to normal numbers
    Eb_n0(nEN)= 10*log10(en/nBitPerSym); % Calculating Eb/N0 in unit db
    sigma = 1/sqrt(en); % standard deviation of AWGN noise
    nframe = 0;
    while (nframe<max_nFrame) && (fErrs(nEN)<fErrLim)
        errbit_count=0;
        nframe = nframe + 1;
        info_bits=round(rand(1,nBitsPerFrame));
        info_matrix=reshape(info_bits, nBitPerSym, nSymPerFrame);
        sym_vec=ones(1,nSymPerFrame);
        for v=1:nBitPerSym
            sym_vec=sym_vec+info_matrix(v,:).*2^(nBitPerSym-v);
        end
        sym_seq=symbolBook(sym_vec);
        %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
        noise=1/sqrt(2)*[randn(1, nSymPerFrame) + 1j*randn(1,nSymPerFrame)];
        det_seq=zeros(1,nBitsPerFrame);
        rec_sig=sym_seq+sigma*noise;
        %%%%DETECTOR %%%%%%%%%%%%
        CODE_SYMBOLS=repmat(transpose(symbolBook),1,nSymPerFrame);
        REC_SIG=repmat(rec_sig,length(symbolBook),1);
        distance_mat=abs(CODE_SYMBOLS-REC_SIG);
        [~, det_sym_ind]=min(distance_mat,[],1);
        detected_bits=[bitBook(det_sym_ind, :)]';

        err_ind= sym_vec-det_sym_ind;
        errsym=sum(err_ind~=0);
        errsyms(nEN)=errsyms(nEN)+errsym;

        errbit = sum(sum(abs(info_matrix-detected_bits)));
        errbits(nEN)=errbits(nEN)+errbit;
        errbit_count=errbit_count+errbit;
        if errbit_count~=0
            fErrs(nEN)=fErrs(nEN)+1;
        end
    end % End of while loop
    nFrames(nEN)=nframe;

end %end for (SNR points)

BER_BFSK=errbits./nFrames/nBitsPerFrame;
save BFSK.mat BER_BFSK

theorerr=qfunc(sqrt(10.^(snr_db/10) ) );

figure
semilogy(Eb_n0, BER_BFSK ,'-x','Color','r');
hold on
semilogy(Eb_n0, theorerr,'s','Color','r');
hold off
xlabel("Noise Ratio_d_b");
ylabel("Error Amount");
legend('BER vs Eb/N_0');
axis square
grid on
set(gca,'FontSize',14);
ylim([0.00001 1]);

