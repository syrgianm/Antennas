                             %DoA MUSIC Find distinctive ability
MusicDistinctive();
function MusicDistinctive ()                             
%Loop for check smaller and smaller delta until we have a same peak for 2
%different signals
for delta=2:(-0.01):0
    
    %Known Data                             
    M=16; %Number Of elementaries in antenna
    N=8;  %Number of Incoming Signals
    Pg=1; %Power for our Signals
    SNRdb=10; %SNR in Db
    SNR=10^(SNRdb/10);   %SNR
    Pn=Pg/SNR;     %Power of Noise.From SNR=Pdesirable/Pnoise



    %Initialize Angles of Incoming Signals
    theta=zeros(1,N);
    for i=1:1:N
        theta(1,i)=50+(i-1)*delta;
    end
  
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


    thetav=(0:0.01:180); %Initialization for scanning angle in P
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
    PdbN=Pdb/max(Pdb);   %Normalize from Pmax
            
        %Plot P power
    %?f we want all the figures just delete % down from %figure
    %figure
    plot(thetav,PdbN);
    ylabel('P/Pmax(Db)');
    xlabel('theta');
    titleId=['Power Diagram Music DoA Delta=' num2str(delta) 'SNR(Db)=' num2str(SNRdb)];
    title(titleId);
    legend('P/Pmax(Db)');
    
    %Check for Peaks 
    pk=findpeaks(PdbN);
    Lowest_Threshold=0.3; %In Normalized Db of radiation Diagram
    sizeOfRealPeaks=0; %Initialize conter
    for i=1:1:size(pk,2)  %Cut Low peaks to keep only from incoming signals
        if pk(i)>=Lowest_Threshold
            sizeOfRealPeaks=sizeOfRealPeaks+1; %Peaks that have a high peak in radiation diagram
        end
    end
   if sizeOfRealPeaks<N  %if we have less peaks than Number of Signals then Display delta and break
        
        %Display Results
         disp(delta);
        
         %save results in a txt

         fileId=fopen('resultsMUSICb.txt','w+');
         fprintf(fileId,'%f',delta);
         fclose(fileId);
          break;
    end
end
end