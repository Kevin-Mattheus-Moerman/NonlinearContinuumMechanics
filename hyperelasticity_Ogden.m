%% Define parameters common to all examples

%Define material parameters
N=1; %The Ogden law order
c1=1; %The shear modulus like parameter
m1=12; %The non-linearity parameter
kp=1000; %Bulk modulus like parameter (used for constrained model)
k=kp; %Bulk modulus (used for uncoupled model)

%Derive applied stretch
appliedStretch=1.3; %Set applied stretch
nDataPoints=50; %Number of data points to use for evaluation and graph
lambda_3=linspace(1,appliedStretch,nDataPoints); %The 3 direction stretch

%% The constrained formulation

%Direct stress computation
S3=(c1/m1).*(lambda_3.^m1-lambda_3.^(-m1/2));
S1=zeros(size(S3));
S2=zeros(size(S3));

%Compute Jacobian for plotting
lambda_1=sqrt(1./lambda_3);
lambda_2=lambda_1;
J=lambda_1.*lambda_2.*lambda_3; 

%Visualize stress graphs
figure; hold on; 
title(['Constrained form. Cauchy stress, min: ',num2str(min(S3(:))),...
', max: ',num2str(max(S3(:)))]); %Add title
h1=plot(lambda_3,S1,'r-','LineWidth',20); %The 1 direction principal stress
h2=plot(lambda_3,S2,'g-','LineWidth',15); %The 2 direction principal stress
h3=plot(lambda_3,S3,'b-','LineWidth',10); %The 3 direction principal stress
hl=legend([h1 h2 h3],{'\sigma_1','\sigma_2','\sigma_3'}); %Add legend
set(hl,'FontSize',15,'Location','NorthEastOutside','Box','off'); %Adjust legend
axis tight; axis square; grid on; box on; 
set(gca,'FontSize',15);

%Visualize Jacobian
figure; hold on; 
title(['Constrained form. Jacobian, min: ',num2str(min(J(:))),...
', max: ',num2str(max(J(:)))]); %Add title
h1=plot(lambda_3,J,'k-','LineWidth',20); %The 1 direction principal stress
hl=legend([h1],{'J'}); %Add legend
set(hl,'FontSize',15,'Location','NorthEastOutside','Box','off'); %Adjust legend
axis tight; axis square; grid on; box on; 
set(gca,'FontSize',15);

%% The unconstrained or coupled formulation 

% One approach is to define a function for S1 and to find the J for which it is zero. 
% For this application the fzero function is useful to find J for S1(J)=0.

%Compute Jacobian given boundary conditions S1=S2=0
J=zeros(size(lambda_3)); %Initialize an arry of J values which are all zeros
for q=1:1:nDataPoints %Loop over all data points
    %Create stress function with current lambda
    S1_fun=@(J) kp*(J-1)+(1/J)*(c1/m1)*((sqrt(J/lambda_3(q)).^m1)-1);    
    
    %Find Jacobian for zero stress, use J=1 as initial
    J(q)=fzero(S1_fun,1); %Find root of nonlinear function   
end

%Compute transverse stretches using J values
lambda_1=sqrt(J./lambda_3);
lambda_2=lambda_1; %Due to uniaxial loading

%Compute principal stresses (note, these are not ordered)
S1=kp*(J-1)+(1./J).*(c1/m1).*((lambda_1.^m1)-1);
S2=kp*(J-1)+(1./J).*(c1/m1).*((lambda_2.^m1)-1);
S3=kp*(J-1)+(1./J).*(c1/m1).*((lambda_3.^m1)-1);

%Visualize stress graphs
figure; hold on; 
title(['Unconstrained form. Cauchy stress, min: ',num2str(min(S3(:))),...
', max: ',num2str(max(S3(:)))]); %Add title
h1=plot(lambda_3,S1,'r-','LineWidth',20); %The 1 direction principal stress
h2=plot(lambda_3,S2,'g-','LineWidth',15); %The 2 direction principal stress
h3=plot(lambda_3,S3,'b-','LineWidth',10); %The 3 direction principal stress
hl=legend([h1 h2 h3],{'\sigma_1','\sigma_2','\sigma_3'}); %Add legend
set(hl,'FontSize',15,'Location','NorthEastOutside','Box','off'); %Adjust legend
axis tight; axis square; grid on; box on; 
set(gca,'FontSize',15);

%Visualize Jacobian
figure; hold on; 
title(['Unconstrained form. Jacobian, min: ',num2str(min(J(:))),...
', max: ',num2str(max(J(:)))]); %Add title
h1=plot(lambda_3,J,'k-','LineWidth',20); %The 1 direction principal stress
hl=legend([h1],{'J'}); %Add legend
set(hl,'FontSize',15,'Location','NorthEastOutside','Box','off'); %Adjust legend
axis tight; axis square; grid on; box on; 
set(gca,'FontSize',15);

%% 
% For this approach the function for S1 is evaluate for a range of J which should cover the required J. 
% Where this graph crosses the x-axis S1(J)=0, and this point is approximated using interpolation.

%Compute Jacobian given boundary conditions S1=S2=0
nTestPoints=100; %Set up a number of test values (more=better but slower)
J_test=linspace(0.9,1.1,nTestPoints); %The test J values
J=zeros(size(lambda_3)); %Initialize an arry of J values which are all zeros
for q=1:1:nDataPoints %Loop over all data points
    %Compute test stresses
    S1_test=kp*(J_test-1)+(1./J_test).*(c1/m1).*((sqrt(J_test./lambda_3(q)).^m1)-1);    
    
    %Find Jacobian for S1(J)=0 using interpolation
    % J(q)=interp1(S1_test,J_test,0,'linear'); %linear interpolation   
    J(q)=interp1(S1_test,J_test,0,'pchip'); %piece-wise cubic hermite interpolation   
end

%Compute transverse stretches using J values
lambda_1=sqrt(J./lambda_3);
lambda_2=lambda_1; %Due to uniaxial loading

%Compute principal stresses (note, these are not ordered)
S1=kp*(J-1)+(1./J).*(c1/m1).*((lambda_1.^m1)-1);
S2=kp*(J-1)+(1./J).*(c1/m1).*((lambda_2.^m1)-1);
S3=kp*(J-1)+(1./J).*(c1/m1).*((lambda_3.^m1)-1);

%Visualize stress graphs
figure; hold on; 
title(['Unconstrained form. Cauchy stress, min: ',num2str(min(S3(:))),...
', max: ',num2str(max(S3(:)))]); %Add title
h1=plot(lambda_3,S1,'r-','LineWidth',20); %The 1 direction principal stress
h2=plot(lambda_3,S2,'g-','LineWidth',15); %The 2 direction principal stress
h3=plot(lambda_3,S3,'b-','LineWidth',10); %The 3 direction principal stress
hl=legend([h1 h2 h3],{'\sigma_1','\sigma_2','\sigma_3'}); %Add legend
set(hl,'FontSize',15,'Location','NorthEastOutside','Box','off'); %Adjust legend
axis tight; axis square; grid on; box on; 
set(gca,'FontSize',15);

%Visualize Jacobian
figure; hold on; 
title(['Unconstrained form. Jacobian, min: ',num2str(min(J(:))),...
', max: ',num2str(max(J(:)))]); %Add title
h1=plot(lambda_3,J,'k-','LineWidth',20); %The 1 direction principal stress
hl=legend([h1],{'J'}); %Add legend
set(hl,'FontSize',15,'Location','NorthEastOutside','Box','off'); %Adjust legend
axis tight; axis square; grid on; box on; 
set(gca,'FontSize',15);

%% The uncoupled formulation 
% One approach is to define a function for S1 and to find the J for which it is zero. 
% For this application the fzero function is useful to find J for S1(J)=0.

%Compute Jacobian given boundary conditions S1=S2=0
J=zeros(size(lambda_3)); %Initialize an arry of J values which are all zeros
for q=1:1:nDataPoints %Loop over all data points
    %Create stress function with current lambda
    S1_fun=@(J) k*(log(J)/J)+(1/J)*(c1/m1)*(sqrt(J/lambda_3(q))^m1...
    -((1/3)*(2*(J/lambda_3(q))^(m1/2)+lambda_3(q)^m1)));    
    
    %Find Jacobian for zero stress, use J=1 as initial
    J(q)=fzero(S1_fun,1); %Find root of nonlinear function   
end

%Compute transverse stretches using J values
lambda_1=sqrt(J./lambda_3);
lambda_2=lambda_1; %Due to uniaxial loading

%Compute principal stresses (note, these are not ordered)
S1=k*(log(J)./J)+(1./J).*(c1/m1).*(lambda_1.^m1-((1/3)*(lambda_1.^m1+lambda_2.^m1+lambda_3.^m1)));
S2=k*(log(J)./J)+(1./J).*(c1/m1).*(lambda_2.^m1-((1/3)*(lambda_1.^m1+lambda_2.^m1+lambda_3.^m1)));
S3=k*(log(J)./J)+(1./J).*(c1/m1).*(lambda_3.^m1-((1/3)*(lambda_1.^m1+lambda_2.^m1+lambda_3.^m1)));

%Visualize stress graphs
figure; hold on; 
title(['Uncoupled form. Cauchy stress, min: ',num2str(min(S3(:))),...
', max: ',num2str(max(S3(:)))]); %Add title
h1=plot(lambda_3,S1,'r-','LineWidth',20); %The 1 direction principal stress
h2=plot(lambda_3,S2,'g-','LineWidth',15); %The 2 direction principal stress
h3=plot(lambda_3,S3,'b-','LineWidth',10); %The 3 direction principal stress
hl=legend([h1 h2 h3],{'\sigma_1','\sigma_2','\sigma_3'}); %Add legend
set(hl,'FontSize',15,'Location','NorthEastOutside','Box','off'); %Adjust legend
axis tight; axis square; grid on; box on; 
set(gca,'FontSize',15);

%Visualize Jacobian
figure; hold on; 
title(['Uncoupled form. Jacobian, min: ',num2str(min(J(:))),...
', max: ',num2str(max(J(:)))]); %Add title
h1=plot(lambda_3,J,'k-','LineWidth',20); %The 1 direction principal stress
hl=legend([h1],{'J'}); %Add legend
set(hl,'FontSize',15,'Location','NorthEastOutside','Box','off'); %Adjust legend
axis tight; axis square; grid on; box on; 
set(gca,'FontSize',15);
