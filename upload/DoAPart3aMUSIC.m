                    %DoA MUSIC
MusicDoA()
                    
function MusicDoA()                   
%Known Data                    
M=16; %Number Of elementaries in antenna
N=8;  %Number of Incoming Signals
Pg=1; %Power for our Signals
SNRdb=10; %SNR in Db
SNR=10^(SNRdb/10);   %SNR
Pn=Pg/SNR;     %Power of Noise.From SNR=Pdesirable/Pnoise
theta=50:10:120; %Angle of Signals

 %Calculate A
 A_vector=zeros(M,N); %Initialization of steering incoming vector

 for signal=1:1:N  %Calculate steering incoming vector
       for pos_elem=1:1:M
            A_vector(pos_elem,signal)=exp((1j)*pi*(pos_elem-1)*cosd(theta(1,signal)));
       end
 end
 Rgg=Pg*eye(N,N); %Correletion Matrix of incoming signals  modulation gi
 %Irrelevant with each_other
 Rnn=Pn*eye(M,M); %Correletion Matrix of noise in every element of antenna
 Rxx=A_vector*Rgg*(A_vector')+Rnn;          %Correletion Matrix of incoming signals  Xi in every element of antenna
 [V , ~]=eig(Rxx); %Take eigen values In Diagonical Matrix from Rxx Matrix,Sorted from Minimum to Maximum by default from Matlab
                %Take eigen vectors thats satisfy  A*V = V*D Sorted by
                %default from Matlab V(:,1) goes to eig D(1,1) which is the
                %minimum and goes on
U=V(:,1:(M-N)); %Construct U Noise Vector


%Initialization for plotting Power Diagram
thetav=(0:0.1:180); %Initialization for scanning angle in P
a=zeros(M,size(thetav,2)); %Initialization  variable a_vector
for pos_elem=1:1:M
b=exp((1j)*pi*(pos_elem-1)*cosd(thetav));
a(pos_elem,:)=b;
end

P=zeros(1,size(thetav,2));
for i=1:1:size(thetav,2)
    P(1,i)=1/((a(:,i)')*U*(U')*a(:,i)); %Calculate P Power
end
Pdb=10*log10(abs(P)); %Convert to Db
PdbN=Pdb/max(Pdb);   %Normalize from Pmax(Db)

%Plot P power
figure
plot(thetav,PdbN);
ylabel('P/Pmax');
xlabel('theta');
titleId=('Power Diagram MUSIC DoA');
title(titleId);
legend('P/Pmax(Db)')

  %Calculate Direction Deviations for Signals
 [~,locs] = findpeaks(abs(PdbN)); %take the peaks from |P| and their index in matrix
 Dtheta=zeros(1,N); %Initialization for Deviation for angles of incoming signals
 Peaks=thetav(1,locs(1,1:end));
for i=1:1:N
    Dtheta(1,i)=min(abs(Peaks-theta(i)));
end
 %Display Results
 disp(Dtheta)
 
 %save results in a txt

fileId=fopen('resultsMUSICa.txt','w+');
for i=1:1:size(Dtheta,2)
    fprintf(fileId,' %f',Dtheta(1,i));
    
end
fclose(fileId);
end 