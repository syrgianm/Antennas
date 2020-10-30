                      % Implementation of MVDR Beamformer


fileId=fopen('resultsMVDR.txt','w+');  %delete if results.txt has data
fclose(fileId);
AoANb=0; %Counter For change AoAdev SINR.txt to save the data

%Calculate Results for different SNRdb,delta
for SNRdb=0:5:20
    for delta=2:2:10
        AoANb=AoANb+1;
        CalcAoA(delta,SNRdb,AoANb);
        CalcRes(AoANb);
    end

end






%function to make AoAdev_SINR.txt 
function CalcAoA(delta,SNRdb,AoANb)

       
N=5; %Number of Incoming Signals
signals_DoA=zeros(1,N); %Initialization for incoming angles
M=16; %Number Of elementaries in antenna
Pg=1; %Power for our Signals
Psd=Pg; %Power of desirable signal
SNR=10^(SNRdb/10);   %SNR
Pn=Pg/SNR;     %Power of Noise.From SNR=Pdesirable/Pnoise
AoANb=int2str(AoANb); %Parameter to take different txt for AoADev
fileName=strcat('AoAdev_SINR',AoANb,'.txt');
fileID=fopen(fileName,'w+');

%Loop for all theta
for polar_angle=30:1:150
        
        
        %Make angles for signals in a vector
        for i=1:1:N
            signals_DoA(1,i)=polar_angle+(i-1)*delta;
        end
        if (polar_angle+4*delta)>150  %Check For Maximum_Angle=150 The worst case in every scenario
            break;                    %If an angle exceed 150 then break the loop
        end
        
        %Take in every loop one desirable and others interference
        for desirable=1:1:size(signals_DoA,2)
              
                theta_0=signals_DoA(desirable); %Desirable Signal
                fprintf(fileID,'%d',theta_0);
                theta_i=zeros(1,N-1); %Initialization for Interference Signals
                counter=0;   
                
                %Seperate desirable angle from the interference angles
                 for i=1:1:N  %Loop for take Interference Signals
                    if i==desirable
                       continue;
                    else
                    counter=counter+1;
                    fprintf(fileID,' %d',signals_DoA(i));
                    theta_i(1,counter)=signals_DoA(i);  %Interference Signal    
                    end
                end
                
                
                %Calculate A,ad vectors
            ad_vector=zeros(M,1); %Initialization of steering vector
            A_vector=zeros(M,N); %Initialization of steering incoming vector
            for signal=0:1:(N-1)  %Calculate Vectors
                for pos_elem=1:1:M
                    if signal==0
           
                       ad_vector(pos_elem,signal+1)=exp((1j)*pi*(pos_elem-1)*cosd(theta_0));
                       A_vector(pos_elem,signal+1)=exp((1j)*pi*(pos_elem-1)*cosd(theta_0));
                    else
                        
                        A_vector(pos_elem,signal+1)=exp((1j)*pi*(pos_elem-1)*cosd(theta_i(1,signal)));
                    end
                    
                end
            end
            
              %Calculate SINR  
            Rgg=Pg*eye(N,N); %Correletion Matrix of incoming signals  modulation gi 
            %Irrelevant with each_other
            Rnn=Pn*eye(M,M); %Correletion Matrix of noise in every element of antenna
            Rxx=A_vector*Rgg*(A_vector')+Rnn;          %Correletion Matrix of incoming signals  Xi in every element of antenna
            
            Wmvdr=(inv(Rxx))*ad_vector; %Weight vector
            Ai=A_vector(:,2:end); %Steering Vector for interference signals
            Rgigi=Rgg(2:end,2:end); %Correletion of incoming interference signals modulation
            SINR=(Psd*(Wmvdr')*ad_vector*(ad_vector')*Wmvdr)/((Wmvdr')*Ai*Rgigi*(Ai')*Wmvdr+(Wmvdr')*Rnn*Wmvdr);
            SINRdb=10*log10(abs(SINR)); %SINR (Db)
            
            %Plot Radiation Diagram
            
            thetav=(0:0.1:180); %Initialization for scanning angle in AF
            a=zeros(M,size(thetav,2)); %Initialization for variable a vector
            for pos_elem=1:1:M
               b=exp((1j)*pi*(pos_elem-1)*cosd(thetav));
               a(pos_elem,:)=b;
            end

            AF=(Wmvdr')*a;  %radiation Diagram
            titleId=['Radiation Diagram MVDR Delta=' num2str(delta) 'SNR(Db)=' num2str(SNRdb)];

            
            %Normalize
            AFN=abs(AF)/max(abs(AF));
            %If we want all the figures just delete % down from %figure
            %figure
            plot(thetav,AFN);
            ylabel('|AF|/|AFmax|');
            xlabel('theta');
            title(titleId);
            legend('|AF|')
            
            Dtheta=zeros(1,N); %Initialization for Deviation for interference angles and |AF| zeros
            
            %Calculate Direction Deviations for Desirable
            [~,locs] = findpeaks(abs(AF)); %take the peaks from |AF| and their index in matrix
            Peaks=thetav(1,locs(1,1:end)); %Angles for the peaks
            Dtheta(1,1)=min(abs(Peaks-theta_0)); %Min Deviation Of Peak and incoming desirable angle
            fprintf(fileID,' %f',Dtheta(1,1));
            
            
             %Calculate Direction Deviations for Interference
            [~,locs] = findpeaks(-abs(AF));
            zeroes=thetav(1,locs(1,1:end)); %Angles for the zeros
            for i=2:1:N
                Dtheta(1,i)=min(abs(zeroes-theta_i(1,(i-1)))); %Min Deviation Of Zero and incoming interference angle
            end 
           for i=2:size(Dtheta,2)
           fprintf(fileID,' %f',Dtheta(1,i));
           end
           fprintf(fileID,' %f\n',SINRdb);
           
        end 
     
end
fclose(fileID);
end
%function to make results.txt
function CalcRes(AoANb)
N=5; %Number of signals

%Read Data From AoAdev_SINR.txt 
AoANb=int2str(AoANb);
fileName=strcat('AoAdev_SINR',AoANb,'.txt');
fileID=fopen(fileName,'r');
sizeA=[11 Inf];
A=fscanf(fileID,'%d %d %d %d %d %f %f %f %f %f %f \n',sizeA);
A=A'; %Transpose array to take the data like in txt
fclose(fileID);

%Calculate min,max,mean std deviation for Dtheta_0(desirable)
vector=A(:,(N+1)); %Take all Dtheta_0 from the array
Dtheta_0max=max(vector);
Dtheta_0min=min(vector);
Dtheta_0mean=mean(vector);
Dtheta_0std=std(vector);

%Calculate min,max,mean std deviation for Dtheta_i(interference signals)

vector=reshape(A(:,(N+2):2*N),[],1);  % convert matrix to column vector every column stuck under the first column
Dtheta_imax=max(vector);
Dtheta_imin=min(vector);
Dtheta_imean=mean(vector);
Dtheta_istd=std(vector);

%Calculate min,max,mean std deviation for SINR(db)
vector=A(:,(2*N+1));
SINRdbmax=max(vector);
SINRdbmin=min(vector);
SINRdbmean=mean(vector);
SINRdbstd=std(vector);

fileId=fopen('resultsMVDR.txt','a+');
fprintf(fileId,'%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f \n',Dtheta_0min,Dtheta_0max,Dtheta_0mean,Dtheta_0std,Dtheta_imin,Dtheta_imax,Dtheta_imean,Dtheta_istd,SINRdbmin,SINRdbmax,SINRdbmean,SINRdbstd);
fclose(fileID);
end