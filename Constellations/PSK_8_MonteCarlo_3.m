clear all
warning off
%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%
Infty=50;
nBitsPerFrame=6000;
max_nFrame =2000;
fErrLim=100;
snr_db=0:2:20;
Eb_n0=zeros(size(snr_db));  %Initializing Eb/N0
%%%%% SYMBOL CONSTELLATION %%%%%%%%%
tb4=[0 0;0 1; 1 0; 1 1];
tb8=[zeros(4,1) tb4; ones(4,1) tb4];
tb8_gray=[0 0 0; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 1 1; 1 0 1; 1 0 0];
nBitPerSym=3;
M=2^nBitPerSym;
symbolBook=exp(1j*2*pi*[0:(M-1)]/M); 
%bitBook=tb8;%Uniform Case
bitBook=tb8_gray; %Gray Mapping
binary2gray=[1 2 4 3 8 7 5 6];
%%%%%%%%%%%%%%ERROR MATRICES%%%%%%%%%%%%%%%%%%%%
nSymPerFrame=nBitsPerFrame/nBitPerSym;
errbits=zeros(length(snr_db), 1);
nFrames=zeros(length(snr_db), 1);
fErrs=errbits;

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
        sym_vec=binary2gray(sym_vec); %%For gray mapping
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
        
                
        errbit = sum(sum(abs(info_matrix-detected_bits)));
        errbits(nEN)=errbits(nEN)+errbit;
        errbit_count=errbit_count+errbit;
        if errbit_count~=0
            fErrs(nEN)=fErrs(nEN)+1;
        end
    end % End of while loop
    nFrames(nEN)=nframe;
    
end %end for (SNR points)

% BER_8PSK_uniform=errbits./nFrames/nBitsPerFrame;
% save PSK8_uniform.mat BER_8PSK_uniform
BER_8PSK_gray=errbits./nFrames/nBitsPerFrame;
save PSK8_gray.mat BER_8PSK_gray


load PSK8_uniform.mat
load PSK8_gray.mat

therorerr_gray=2*qfunc(sqrt( 2*10.^(snr_db./10) ).*sin(pi/M) )./log2(M);
theorrerr_uniform=therorerr_gray*7/4;

figure
semilogy( snr_db, BER_8PSK_uniform,'-x','Color','r');
hold on
semilogy( snr_db, BER_8PSK_gray,'-x','Color','b');
semilogy( snr_db, theorrerr_uniform,'s','Color','r');
semilogy( snr_db, therorerr_gray,'s','Color','b');
xlabel("SNR_d_b");
ylabel("BER");
legend('Uniform Mapping', 'Gray Mapping'); 
axis square
grid on
set(gca,'FontSize',14);
ylim([0.00001 1]);
