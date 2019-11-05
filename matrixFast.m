clear all
close all

Fs = 44100; % Sample rate
k = 1/Fs;   % Time step
T = 1;     % end time
c = 100;     % Wave speed
h = 2*c*k;  % Grid spacing  
L = 1;      % Length of string

Ns = floor(T/k);   % Number of samples
N = floor(L/h);    % N

h = L/N; % Redefine h so it matches with N
b = 1; % Damping coefficient

% Initialising the matrices for system
A = (1 + b*k/2) * eye(N+1,N+1);
C = (-1 + b*k/2) * eye(N+1,N+1);
b1 = 2-(2*c^2*k^2)/(h^2);
b2 = (c^2*k^2)/(h^2);
nOnes = ones(N+1, 1) ;
B = diag(b1*nOnes, 0) + b2*diag(nOnes(1:N), -1) + b2*diag(nOnes(1:N), 1);



% Starting conditions for first two periods
u = zeros(Ns,N+1);
u(1,(round(N/2-5):round(N/2+5))) = hann(1,N+1);
u(2,(round(N/2-5):round(N/2+5))) = hann(1,N+1);
u(1,1) = 0;
u(1,N+1)=0;


% Loop to show the string over time
for i=2:Ns-1
    
    % Fixed end points
    u(i,1) = 0;
    u(i,N+1) = 0;
    
    % Calculates next step in time
    u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)');
    
    % Plots string in real time
    plot(u(i,:));
    axis([0,round(N),-1,1]);
    drawnow
end


 

 
 
 
 
 