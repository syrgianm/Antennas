            % RLS BeamFormer

  RLSbeamformer();          
            
            
            
            
            
            
            
            
            
            
            
            
function RLSbeamformer()
M=16; %Number Of elementaries in antenna
N=6;  %Number of Incoming Signals
Pg=1; %Power of Incoming Signals
Pn=sqrt(0.1); %Power of Noise
%Incoming signals(polar_angle
theta=50:20:150; %Incoming Signals
 
 
             %Calculate A teering incoming vector
 A_vector=zeros(M,N); %Initialization of steering incoming vector
 for signal=0:1:(N-1)  %Calculate Vectors
       for pos_elem=1:1:M
            A_vector(pos_elem,signal+1)=exp((1j)*pi*(pos_elem-1)*cosd(theta(1,(signal+1))));
       end
 end
         
       
   

        %Initialize Algorithm     
        %q=1
a=0.9; %forgetting factor
Q=50;  %NumberOfSamples
Wrls=zeros(M,1); %Initialize weight-vector
delta=10^6;
invRxx=delta*eye(M,M); %Initialize Rxx^(-1)


    %Repetition Of RLS Algorithm
for q=2:1:Q
    g=Pg*randn(N,1); %Calculate Signals_modulation from normal Distribution
    n=Pn*randn(M,1); %Calculate noise From normal Distribution
    r0=g(1,1);       %Reference Signal = Desirable Signal
    x=A_vector*g+n;  %Calculate x_vector 
    h=(a^(-1)*invRxx*x)/(1+a^(-1)*(x')*invRxx*x); %calculate h
    invRxx=(q/(q-1))*(a^(-1)*invRxx-a^(-1)*h*(x')*invRxx); %Calculate Rxx^(-1)
    Wrls(:,1)=Wrls(:,1)+h*(conj(r0)-(x')*Wrls(:,1));  %Calculate Weight lifting 
    
    
end


%Display Weight Vector after 50 repetitions
disp(Wrls(:,1));

%Plot Radiation Diagram

thetav=(0:0.1:180); %Initialization for scanning angle in AF
a=zeros(M,size(thetav,2)); %Initialization for variable a vector

for pos_elem=1:1:M
b=exp((1j)*pi*(pos_elem-1)*cosd(thetav));
a(pos_elem,:)=b;
end

AF=(Wrls(:,1)')*a;  %radiation Diagram

%Normalize
AFN=abs(AF)/max(abs(AF));
   
%if we want all the figures just delete % down from %figure
%figure
plot(thetav,AFN);
ylabel('|AF|/|AFmax|');
xlabel('theta');
titleId=('Radiation Diagram RLS Q=50 repetitions');
title(titleId);
legend('|AF|')

Dtheta=zeros(1,N); %Initialization for Deviation for angles and |AF| max,zeros

%Calculate Direction Deviations for Desirable
 [~,locs] = findpeaks(abs(AF)); %take the peaks from |AF| and their index in matrix

 Peaks=thetav(1,locs(1,1:end)); %Angles for the peaks
 Dtheta(1,1)=min(abs(Peaks-theta(1,1))); %Min Deviation Of Peak and incoming desirable angle
           
 %Calculate Direction Deviations for Interference
  [~,locs] = findpeaks(-abs(AF));
  zeroes=thetav(1,locs(1,1:end)); %Angles for the zeros
  for i=2:1:N
       Dtheta(1,i)=min(abs(zeroes-theta(1,i))); %Min Deviation Of Zero and incoming interference angle
  end 
 %Display Deviations for every signal
disp(Dtheta);

%save results in a txt

fileId=fopen('resultsRLS.txt','w+');
for i=1:1:size(Dtheta,2)
    fprintf(fileId,' %f',Dtheta(1,i));
    
end
fclose(fileId);
end


