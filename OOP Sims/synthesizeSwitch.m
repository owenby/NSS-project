function out = synthesizeSwitch(Synthesize, Effect, A, B, C, u, uh, g, f, a, c, N, Ns, hampnt, m, kappa, b0, h, k, ps)

% This functions runs the finite difference scheme, in 3 cases, 'plot'
% where a plot of the moving string is shown. 'play' where the oscillation
% of a point on the string is measured to get a wave to play. 'plotsound'
% plots the graph of the sound produced by 'play'.
eta=zeros(Ns,1);
eta(1)=uh(1)-u(1,hampnt);
eps=10^(-6);
G=1;
alpha=1;
x=linspace(0,1,N+1);
for i=2:Ns-1
    switch Effect
        case 'hamstr' %Finds bn and f for the ith timestep            
            eta(i)=uh(i)-u(i,hampnt);
            an=eta(i-1);
            phia=kappa/(alpha+1)*(max(0,an)^(alpha+1));
            dxxb=c*k^2*(u(i,hampnt+1)-2*u(i,hampnt)+u(i,hampnt-1))/h^2;     
            bn=-(2*uh(i)-2*uh(i-1)-2*u(i,hampnt)-dxxb+2*u(i-1,hampnt));
        case 'hamlss'
            eta(i)=uh(i)-u(i,hampnt);
            an=eta(i-1);  
            phia=kappa/(alpha+1)*(max(0,an)^(alpha+1));
            dxxb=c*k^2*(u(i,hampnt+1)-2*u(i,hampnt)+u(i,hampnt-1))/h^2; 
            bn=-(2*uh(i)-2*uh(i-1)-(1/(1+(b0*k/(2*ps))))...
               *((2*u(i,hampnt)+c*k^2*dxxb)+(b0*k/(2*ps)-2)*u(i-1,hampnt)));
        case 'hamstf'
            eta(i)=uh(i)-u(i,hampnt);
            an=eta(i-1); 
            phia=kappa/(alpha+1)*(max(0,an)^(alpha+1));
            dxxb=c*k^2*(u(i,hampnt+1)-2*u(i,hampnt)+u(i,hampnt-1))/h^2;
            dxxxxb=(kappa^2*k^2/ps)*(u(i,hampnt+2)-4*u(i,hampnt+1)+6*u(i,hampnt)...
                   -4*u(i,hampnt-1)+u(i,hampnt-2))/h^4;
            bn=-(2*uh(i)-2*uh(i-1)-(2*u(i,hampnt)+dxxb+dxxxxb)+2*u(i-1,hampnt));            
    end
             
        switch Effect 
          case {'hamstr','hamlss','hamstf'}
               rn=an;
            while abs(G)>eps
                phir=kappa/(alpha+1)*(max(0,rn+an)^(alpha+1));
                G=rn+m*((phir-phia)/rn)+bn;
                Gp=1+(m*kappa*((alpha+1)*rn*max(0,rn+an)^(alpha+1)...
                    -max(0,an)^(alpha+1))/(rn^2*(alpha+1)))...
                    +m*kappa*max(0,an)^(alpha+1)/(rn^2*(alpha+1)^2);
                rn=rn-G/Gp;
            end
            f(i)=(phir-phia)/rn;
        end
        
    
    % Calculates next step in time
    u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a');
    switch Effect
        case {'hamstr','hamlss','hamstf'}
            %Updates hammer position
            uh(i+1)=2*uh(i)-uh(i-1)-(k^2/m)*f(i);
    end
    
    switch Synthesize
        case 'plot'
            % Plots string in real time
            switch Effect
                case {'hamstr','hamstf','hamlss'}
                    gham=zeros(1,N+1);
                    gham(hampnt)=uh(i);
                    plot(x,u(i,:),x,gham,'o');
                otherwise
                    plot(x,u(i,:),'LineWidth',2);
            end
            
            axis([0,1,-1,1]);
            drawnow
        case 'play'
            % Vibration of one point on string
            out(i) = u(i,round(N/7));
        case 'plotsound'
            % Vibration of one point on string
            out(i) = u(i,round(N/7));
    end


switch Synthesize
    case 'play'
        soundsc(out)
    case 'plotsound'
        soundsc(out)
        plot(out)
end

G=1;

end
% switch Synthesize
%     case 'plot'
%         % Loop to show the string over time      
%         for i=2:Ns-1 
%             switch Effect
%                 case 'hamstr' %Finds bn for the ith timestep                
%                  G=1;
%                  eta(i)=uh(i)-u(i,hampnt);
%                  an=eta(i-1);
%                  phia=kappa/(alpha+1)*(max(0,an)^(alpha+1));
%                  dxxb=c*k*(u(i,hampnt+1)-2*u(i,hampnt)+u(i,hampnt-1)); 
%                  bn=-(2*uh(i)-2*uh(i-1)-2*u(i,hampnt)-dxxb+2*u(i-1,hampnt));
%                  rn=1;
%                  while abs(G)>eps
%                     phir=kappa/(alpha+1)*(max(0,rn+an)^(alpha+1));
%                     G=rn+m*((phir-phia)/rn)+bn;
%                     Gp=1+(m*kappa*((alpha+1)*rn*max(0,rn+an)^(alpha+1)...
%                         -max(0,an)^(alpha+1))/(rn^2*(alpha+1)))...
%                         +m*kappa*max(0,an)^(alpha+1)/(rn^2*(alpha+1)^2);
%                     rn=rn-G/Gp;
%                  end
%                  f(i)=(phir-phia)/rn;                 
%             end
%             % Calculates next step in time
%             u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);
%             switch Effect
%                 case 'hamstr'
%                     uh(i+1)=2*uh(i)-uh(i-1)-(k^2/m)f(i);
%             end
%             % Plots string in real time
%             plot(u(i,:), 'LineWidth',2);
%             axis([0,N+2,-1,1]);
%             drawnow
%         end
%        
%     case 'play'
%         out = zeros(1,Ns-1);  % To hear the sound
%         % Loop to show the string over time
%         for i=2:Ns-1
%              switch Effect
%                 case 'hamstr' %Finds bn for the ith timestep                
%                  G=1;
%                  eta(i)=uh(i)-u(i,hampnt);
%                  an=eta(i-1);
%                  phia=kappa/(alpha+1)*(max(0,an)^(alpha+1));
%                  dxxb=c*k*(u(i,hampnt+1)-2*u(i,hampnt)+u(i,hampnt-1)); 
%                  bn=-(2*uh(i)-2*uh(i-1)-2*u(i,hampnt)-dxxb+2*u(i-1,hampnt));
%                  rn=1;
%                  while abs(G)>eps
%                     phir=kappa/(alpha+1)*(max(0,rn+an)^(alpha+1));
%                     G=rn+m*((phir-phia)/rn)+bn;
%                     Gp=1+(m*kappa*((alpha+1)*rn*max(0,rn+an)^(alpha+1)...
%                         -max(0,an)^(alpha+1))/(rn^2*(alpha+1)))...
%                         +m*kappa*max(0,an)^(alpha+1)/(rn^2*(alpha+1)^2);
%                     rn=rn-G/Gp;
%                  end
%                  f(i)=(phir-phia)/rn;                 
%             end
%             % Calculates next step in time
%             u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);
%             switch Effect
%                 case 'hamstr'
%                     uh(i+1)=2*uh(i)-uh(i-1)-(k^2/m)f(i);
%             end
%             % Vibration of one point on string
%             out(i) = u(i,round(N/7));
%         end
%         soundsc(out)    % To hear sound
%        
%     case 'plotsound'
%         out = zeros(1,Ns-1);  % To hear the sound
%         % Loop to show the string over time
%         for i=2:Ns-1
%             switch Effect
%                 case 'hamstr' %Finds bn for the ith timestep                
%                  G=1;
%                  eta(i)=uh(i)-u(i,hampnt);
%                  an=eta(i-1);
%                  phia=kappa/(alpha+1)*(max(0,an)^(alpha+1));
%                  dxxb=c*k*(u(i,hampnt+1)-2*u(i,hampnt)+u(i,hampnt-1)); 
%                  bn=-(2*uh(i)-2*uh(i-1)-2*u(i,hampnt)-dxxb+2*u(i-1,hampnt));
%                  rn=1;
%                  while abs(G)>eps
%                     phir=kappa/(alpha+1)*(max(0,rn+an)^(alpha+1));
%                     G=rn+m*((phir-phia)/rn)+bn;
%                     Gp=1+(m*kappa*((alpha+1)*rn*max(0,rn+an)^(alpha+1)...
%                         -max(0,an)^(alpha+1))/(rn^2*(alpha+1)))...
%                         +m*kappa*max(0,an)^(alpha+1)/(rn^2*(alpha+1)^2);
%                     rn=rn-G/Gp;
%                  end
%             % Calculates next step in time
%             u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);
%             switch Effect
%                 case 'hamstr'
%                     uh(i+1)=2*uh(i)-uh(i-1)-(k^2/m)f(i);
%             end
%             % Vibration of one point on string
%             out(i) = u(i,round(N/3));
%         end
%         plot(out)       % To plot the sound 
% end
% 
% end
% 
