%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ZFE%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=[0.74 -0.514 0.37 0.216 0.062]; %Channel
messagelength=10000;
maxnframe=600;
snr_db= -10:1:15;

tap = 10;
avgerrors=zeros(size(snr_db));

g(1,:)=[h zeros(1,tap-1)];
for i=2:tap
  g(i,:)=circshift(g(i-1,:),1);
  g(i,1)=0;
end
g=g';
for nsnr=1:length(snr_db) %Loop over Snr values

    for nframe=1:maxnframe

    %Creates Message Signal and Noisy Received Signal

    message=2*randi(2,1,messagelength)-3;
    r=conv(message,h,'same');

    snr_p=snr_db(nsnr);
    sn = 10^(snr_p/10);
    sigma = 1/sqrt(sn);
    noise=1/sqrt(2)*randn(1,messagelength);
    r=r+sigma*noise;


    size_=tap+length(h)-1;
    e=zeros(1,size_);
    e(1+size_/2)=1;
    w=(g'*g\g'*e')';

    m_pred3=conv(r,w,'same');
    m_pred3(m_pred3>0)=1;
    m_pred3(m_pred3<=0)=-1;

    err_count=0;
    for i=1:messagelength

       if m_pred3(i)~=message(i)
          err_count=err_count+1;
       end
    end
    avgerrors(nsnr)=avgerrors(nsnr)+err_count/maxnframe;

    end
end

avgerrors=avgerrors/messagelength;

figure(1);
semilogy(snr_db,avgerrors,'-x','Color','r');
hold on
axis square
grid on
set(gca,'FontSize',14);
xlabel("SNR_d_b");
ylabel("BER");
legend("ZFE");
