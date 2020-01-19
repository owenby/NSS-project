function [f, a] = forceSwitch(Force, h, N, Ns)

% This function adds an external force to the string, at a given point on
% the string (a), of different types. 'off' is no force. 'pluck' is pulling
% the string then releasing sharply. 'hardpluck' is a pluck but of a higher
% magnitude and faster pull. 'strike' is like pluck but without the sharp
% release. 

% Initialising the force
f = zeros(Ns,1);
% Where the string is forced
a = zeros(N+1,1);   
a(7) = 1/h;        

% Switching the types of excitation
switch Force
    case 'off'
        % :)
    case 'pluck' 
        f(1:15) = 0.001*hann(15);  % hann function cut short
        f(10:end) = 0;
    case 'hardpluck'
        f(1:5) = 0.02*hann(5);     % Narrower, 'louder' hann functioon cut short
        f(3:5) = 0;
    case 'strike'
        f(5:20) = 0.001*hann(16);  % Excited with a hann functino curve
end

end

