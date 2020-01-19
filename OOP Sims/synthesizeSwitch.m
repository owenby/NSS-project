function out = synthesizeSwitch(Synthesize, A, B, C, u, f, a, N, Ns)

% This functions runs the finite difference scheme, in 3 cases, 'plot'
% where a plot of the moving string is shown. 'play' where the oscillation
% of a point on the string is measured to get a wave to play. 'plotsound'
% plots the graph of the sound produced by 'play'.

switch Synthesize
    case 'plot'
        % Loop to show the string over time
        for i=2:Ns-1
            % Calculates next step in time
            u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);
            % Plots string in real time
            plot(u(i,:), 'LineWidth',2);
            axis([0,N+2,-1,1]);
            drawnow
        end
        
    case 'play'
        out = zeros(1,Ns-1);  % To hear the sound
        % Loop to show the string over time
        for i=2:Ns-1
            % Calculates next step in time
            u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);
            % Vibration of one point on string
            out(i) = u(i,round(N/7));
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

end

