clear sound
close all

% This script uses finite difference schemes to realistically simulate a 
% string, with given properties, with an end goal of synthesising sounds 
% from the string simulation

%%% Script Variables %%%

Boundary = 'fixed';    % 'free' both ends of the string free
                       % 'fixed' both ends fixed
                       % 'LfreeRfixed' left side free, right side fixed
                       % 'RfreeLfixed' right side free, left side fixed
                       % 'xylo' fixed like a xylophone
                       
Starting = false;      % 'true' starts the string from a hann function shape
                       % 'false' starts the string from 0 amplitude everywhere 
                        
Force = 'pluck';       % 'off' for nothing
                       % 'pluck'
                       % 'hardpluck'
                       % 'strike'
                       
Effect = 'fdl';        % 'none' for nothing
                       % 'loss' has just one type of loss
                       % 'realistic' has two
                       % 'bar' is a bar
                       % 'stiff' is a soft bar
                       % 'fdl' is frequency dependent loss
                       % 'hamstr' for lossless hammer-string interaction
                       % 'hamlss' for single loss term ham str interaction
                       
Synthesize = 'play';   % Whether you see string ('plot'), hear string ('play'), or see the sound ('plotsound')


%%% Variable Setup %%%

% General variables
Fs = 44100;        % Sample rate
k = 1/Fs;          % Time step
T = 1.5;           % end time
L = 1;             % Length of string
c = 2*L*261;       % Wavespeed
h = c*k;           % Grid spacing  
Ns = floor(T/k);   % Number of samples
N = floor(L/h);    % Number of string chunks
h = L/N;           % Redefine h so it matches with N

% Loss variables
b0 = 2;          % First loss term
b1 = 0.1;        % Second loss term

% Stiffness variables
E = 2;                       % Youngs modulus of material
p = 800;                     % Density
A = 0.2;                     % Cross-sectional area
I = sqrt(A/pi)/2;            % Bar moment of inertia (I = R/2 for cylinder)
kappa = sqrt((E*I)/(p*A));   % Stiffness parameter
theta = 1;                   % 'free parameter'
hampnt=0;                    %Hammer position relative to string 
uh=zeros(Ns,1);              %Hammer position vector   
switch Effect
    case 'stiff'
        Tns= ;              %Tension
        c = sqrt(T/(p*A));  %Need to change T to tension        
    case {'hamstr','hamlss'}
        Tns= ;              %Tension
        ps= ;               %Set linear mass density
        c=Tns/ps;
        M= ;                %Mass of hammer
        alpha= ;            %exp for phi
        m=(k^2)/M+(k^2)/(ps*h);      
        uh(1)= ;,uh(2)= ; %init conds of hammer
        hampnt=round((N+1)/2); 
    case 'hamlss'
        m=(k^2)/(h*(1+(b0*k)/(2*ps)))+k^2/M;
end

% Starting position of the string
u = zeros(Ns,N+1);

% Starting the string from a hann function shape
if Starting 
    u(1,(round(N/2-5):round(N/2+5))) = hann(11);
    u(2,(round(N/2-5):round(N/2+5))) = hann(11);
end

%Hammer distribution
g=zeros(1,N+1);,g(hampnt)=1/h;   


%%% Running Functions For Synthesis %%%

% Initialising the matrices for system
[A, B, C] = effectSwitch(Effect, b0, b1, c, h, k, N, kappa, theta, ps);

% Force switch
[f, a] = forceSwitch(Force, h, N, Ns, g, k, ps);

% Boundary condition switch
B = boundarySwitch(B, N, Boundary); 

% Switch to synthesize the sound
synthesizeSwitch(Synthesize, Effect, A, B, C, u, f, a, c, N, Ns,hampnt, m, kappa, b0);

 