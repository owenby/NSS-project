
clear all
close all

% Things to set 
Boundary = 'fixed';    % Boundary condition (free, fixed, LfreeRfixed, RfreeLfixed)
Starting = true;       % Whether the string starts with a starting position (true) or is flat (false)
force = 'off';         % What kind of force applied to the string ('off' for no force)
Effects = 'fdl';       % 'loss' has just one type of loss
                       % 'realistic' has two.
                       % 'bar' is a bar
                       % 'stiff' is a soft bar
                       % 'fdl' is frequency dependent loss
                       
Synthesize = 'play';         % Whether you see string ('plot'), hear string ('play'), or see the sound ('plotsound')




% Variable setup
Fs = 44100;        % Sample rate
k = 1/Fs;          % Time step
T = 1;             % end time
L = 1;             % Length of string
c = 500;           % Wave speed
h = c*k;           % Grid spacing  
Ns = floor(T/k);   % Number of samples
N = floor(L/h);    % N
h = L/N;           % Redefine h so it matches with N
% Stiffness variables
E = 200;             % Youngs modulus of material
p = 5000;             % Density
I = 0.5;             % Bar moment of inertia (I = R/2 for cylinder)
A = 1;             % Cross-sectional area
kappa = sqrt((E*I)/(p*A*L^4));   % Stiffness parameter
theta = 1;         % 'free parameter'


% Starting conditions for string
u = zeros(Ns,N+1);
if Starting 
    u(1,(round(N/2-5):round(N/2+5))) = hann(11);
    u(2,(round(N/2-5):round(N/2+5))) = hann(11);
%   u(1,:) = hann(N+1);
%   u(2,:) = hann(N+1);
    u(1,1) = 0;
    u(1,N+1)=0;
end

% Adding an external force
a = zeros(N+1,1);
r = 4;
a(round(r)) = 1/h;
switch force
    case 'off'
        f = zeros(1,Ns-1);
    case 'pluck'
        f = zeros(Ns,1);
        f(1:round(N/10)) = 0.01*hann(round(N/10));
        f(round(N/20):round(N/10)) = 0;
        f(round((Ns-1)/2):Ns-1) = 0;
    case 'strike'
end







% Initialising the matrices for system
b0 = 1;                % First loss term (damping coefficient)
b1 = 0.05;             % Second loss term (realistic coefficient)

nOnes = ones(N+1, 1) ;

Dxx = (1/h^2)*(-2*diag(nOnes, 0) + diag(nOnes(1:N), -1) + diag(nOnes(1:N), 1));
Dxxxx = (1/h^4)*(diag(nOnes(1:N-1), -2)-4*diag(nOnes(1:N), -1) + 6*diag(nOnes, 0) - 4*diag(nOnes(1:N), 1) + diag(nOnes(1:N-1), 2));
Mx = (1/2) * (diag(nOnes(1:N), -1) + diag(nOnes(1:N), 1));
Dxx = sparse(Dxx);
Dxxxx = sparse(Dxxxx);
Mx = sparse(Mx);

switch Effects
    case 'loss'
        A = (1 + b0*k/2) * eye(N+1,N+1);
        C = (-1 + b0*k/2) * eye(N+1,N+1);
        b1 = 2-(2*c^2*k^2)/(h^2);
        b2 = (c^2*k^2)/(h^2);
        B = diag(b1*nOnes, 0) + b2*diag(nOnes(1:N), -1) + b2*diag(nOnes(1:N), 1);
    case 'realistic'
        A = (1 + b0*k/2) * eye(N+1,N+1) - b1*k/2*Dxx;
        B = 2*eye(N+1,N+1) + c^2*k^2*Dxx;
        C = A - 2*eye(N+1,N+1);
    case 'bar'
        A = theta * eye(N+1,N+1) + (1-theta) * Mx ;
        B = -kappa^2 * k^2 * Dxxxx + 2*A;
        C = -A;
    case 'stiff'
        A = theta * eye(N+1,N+1) + (1-theta) * Mx ;
        B = k^2*c^2*Dxx - kappa^2*k^2*Dxxxx + 2*A;
        C = -A;
    case 'fdl'
        A = (1+b0*k)*eye(N+1,N+1) - b1*k*Dxx;
        B = -2*eye(N+1,N+1) - c^2*k^2*Dxx + kappa^2*k^2*Dxxxx;
        C = (1+b0*k)*eye(N+1,N+1) + b1*k*Dxx;
        A = -A; % It works :D 
end
     




% Making the matrices sparse martrices for optimality 
A = sparse(A);
B = sparse(B);
C = sparse(C);

% Switching the cases of the boundary out of 4 options
switch Boundary
    case 'free'
        B(1,1) = B(1,1)/2;
        B(N+1,N+1) = B(1,1)/2;
    case 'fixed'  
        B(1,:) = 0;
        B(N+1,:) = 0;
    case 'LfreeRfixed'
        B(1,1) = B(1,1)/2;
        B(N+1,:) = 0;
    case 'RfreeLfixed'
        B(1,:) = 0;
        B(N+1,N+1) = B(1,1)/2;
end





% Switch to synthesize the sound

switch Synthesize
    case 'plot'
        % Loop to show the string over time
        for i=2:Ns-1

            % Calculates next step in time
            u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);

            % Plots string in real time
            plot(u(i,:));
            axis([0,N+1,-1,1]);
            drawnow

        end
        
    case 'play'
        out = zeros(1,Ns-1);  % To hear the sound
        % Loop to show the string over time
        for i=2:Ns-1

            % Calculates next step in time
            u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);

            % Vibration of one point on string
            out(i) = u(i,round(N/3));

        end
        soundsc(out)    % To hear sound
       
        
    case 'plotsound'
        out = zeros(1,Ns-1);  % To hear the sound
        % Loop to show the string over time
        for i=2:Ns-1

            % Calculates next step in time
            u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);

            % Vibration of one point on string
            out(i) = u(i,round(N/3));

        end
        plot(out)       % To plot the sound 
end
 
