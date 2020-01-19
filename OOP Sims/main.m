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
switch Effect
    case 'stiff'
        c = sqrt(T/(p*A));
end

% Starting position of the string
u = zeros(Ns,N+1);

% Starting the string from a hann function shape
if Starting 
    u(1,(round(N/2-5):round(N/2+5))) = hann(11);
    u(2,(round(N/2-5):round(N/2+5))) = hann(11);
end


%%% Running Functions For Synthesis %%%

% Initialising the matrices for system
[A, B, C] = effectSwitch(Effect, b0, b1, c, h, k, N, kappa, theta);

% Force switch
[f, a] = forceSwitch(Force, h, N, Ns);

% Boundary condition switch
B = boundarySwitch(B, N, Boundary); 

% Switch to synthesize the sound
synthesizeSwitch(Synthesize, A, B, C, u, f, a, N, Ns);

 