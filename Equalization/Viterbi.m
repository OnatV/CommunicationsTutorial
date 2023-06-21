clear
clc
close all

h=[0.74 -0.514 0.37 0.216 0.062];
messagelength=10000;
maxnframe=600;
snr_db=-10:1:15;
avgerrors=zeros(size(snr_db));


%Defines States
s1=[-1 -1; -1 1; 1 -1; 1 1];
s2=[-ones(4,1) , s1; ones(4,1),s1];
states=[-ones(8,1) , s2; ones(8,1),s2];



%%Finding Branch Lengths (Required Outputs)
temp=(h(2:end))';
a=-0.74+states*temp; %Output if -1 arrives
b=0.74+states*temp; %Output if 1 arrives
req=zeros(16,16);


for i=1:16
   for j=1:16
       
       if(j==ceil(i/2))
          req(i,j)=a(i);
       elseif (j==ceil(i/2)+8)
           req(i,j)=b(i);      
       else
           req(i,j)=inf;
       end
   end
    
end

for nsnr=1:length(snr_db) %Loop over Snr values
    
    for nframe=1:maxnframe
    
    %Creates Message Signal and Noisy Received Signal

    message=2*randi(2,1,messagelength)-3;
   
    r=zeros(1,messagelength); 
    for i=1:messagelength
       if i==1
          r(i)=0.74*message(i);
          
       elseif i==2
           r(i)=0.74*message(i)-0.514*message(i-1);
       elseif i==3
            r(i)=0.74*message(i)-0.514*message(i-1)+0.37*message(i-2);
       elseif i==4
            r(i)=0.74*message(i)-0.514*message(i-1)+0.37*message(i-2)+0.216*message(i-3);
       else
            r(i)=0.74*message(i)-0.514*message(i-1)+0.37*message(i-2)+0.216*message(i-3)+0.062*message(i-4);
       end
    end
    state_weights=inf*ones(messagelength+1,16); %Length of the path ending at state k at sample i
    state_weights(1,1)=0; %Initial State
    paths=zeros(messagelength+1,16);
    snr_p=snr_db(nsnr);
    sn = 10^(snr_p/10);
    sigma = 1/sqrt(sn);
    noise=1/sqrt(2)*randn(1,messagelength);    
    r=r+sigma*noise;
    
   
    %%Calculates State Weights
    for i=2:messagelength+1
       for k=1:16

           if k<=8 
           cand1=state_weights(i-1,2*k)+(r(i-1)-req(2*k,k))^2;
           cand2=state_weights(i-1,2*k-1)+(r(i-1)-req(2*k-1,k))^2;     
           state_weights(i,k)=min(cand1,cand2);
           
               if min(cand1,cand2)==cand2 %%records previos state
                   paths(i,k)=2*k-1;
               else 
                   paths(i,k)=2*k;
               end
           else 
           cand1=state_weights(i-1,2*(k-8))+(r(i-1)-req(2*(k-8),k))^2;
           cand2=state_weights(i-1,2*(k-8)-1)+(r(i-1)-req(2*(k-8)-1,k))^2;
           state_weights(i,k)=min(cand1,cand2);
               if min(cand1,cand2)==cand2 %%records previos state
                   paths(i,k)=2*(k-8)-1;
               else 
                   paths(i,k)=2*(k-8);
               end
           end



       end

    end

    minstate=find(state_weights(end,:)==min(state_weights(end,:)));

    m_pred=zeros(1,messagelength); %Predicted Message
    
    %Finds the predicted message
    for i=messagelength:-1:1
        m_pred(i)=states(minstate,1);
        minstate=paths(i+1,minstate);

    end
    err_count=0;
    for i=1:messagelength

       if m_pred(i)~=message(i)
          err_count=err_count+1; 
       end
    end
    err_count=err_count/messagelength;
    avgerrors(nsnr)=avgerrors1(nsnr)+err_count;
    end
    avgerrors(nsnr)=avgerrors1(nsnr)/maxnframe;
end



figure(1);
semilogy(snr_db,avgerrors,'-x','Color','r');
hold on
axis square
grid on
set(gca,'FontSize',14);
xlabel("SNR_d_b");
ylabel("BER");
legend("Viterbi");



