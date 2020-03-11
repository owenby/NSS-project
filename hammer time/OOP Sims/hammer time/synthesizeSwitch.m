function out = synthesizeSwitch(Synthesize, Effect, A, B, C, u, uh, g, f, a, c, N, Ns, hampnt, m, kappa, b0, b1,h, k, ps,Kham, M)

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

switch Effect
    case 'hamlss2'
        gAig=A\g;
        gAig=g'*gAig;        
end

for i=2:Ns-1
    
    switch Effect
        case 'hamlss2'
            gAiB=A\(B*u(i,:)');
            gAiB=g'*gAiB;
            gAiC=A\(C*u(i,:)');
            gAiC=g'*gAiC;
    end
    
    switch Effect
        case 'hamstr' %Finds bn and f for the ith timestep            
            eta(i)=uh(i)-u(i,hampnt);
            an=eta(i-1);
            phia=Kham/(alpha+1)*(max(0,an)^(alpha+1));
            dxxb=c*k^2*(u(i,hampnt+1)-2*u(i,hampnt)+u(i,hampnt-1))/h^2;     
            bn=-(2*uh(i)-2*uh(i-1)-2*u(i,hampnt)-dxxb+2*u(i-1,hampnt));
        case 'hamlss'
            eta(i)=uh(i)-u(i,hampnt);
            an=eta(i-1);  
            phia=Kham/(alpha+1)*(max(0,an)^(alpha+1));
            dxxb=c*k^2*(u(i,hampnt+1)-2*u(i,hampnt)+u(i,hampnt-1))/h^2; 
            bn=-(2*uh(i)-2*uh(i-1)-(1/(1+(b0*k/(2*ps))))...
               *((2*u(i,hampnt)+c*k^2*dxxb)+(b0*k/(2*ps)-2)*u(i-1,hampnt)));
        case 'hamlss2'
            eta(i)=uh(i)-u(i,hampnt);
            an=eta(i-1);
            phia=Kham/(alpha+1)*(max(0,an)^(alpha+1));
            bn=2*uh(i)-2*uh(i-1)-h*gAiB-h*gAiC...
                -u(i-1,hampnt);
            m=-((k^2)/(M)-((h*k^2)/(ps))*gAig);
        case 'hamstf'
            eta(i)=uh(i)-u(i,hampnt);
            an=eta(i-1); 
            phia=Kham/(alpha+1)*(max(0,an)^(alpha+1));
            dxxb=c*k^2*(u(i,hampnt+1)-2*u(i,hampnt)+u(i,hampnt-1))/h^2;
            dxxxxb=(kappa^2*k^2/ps)*(u(i,hampnt+2)-4*u(i,hampnt+1)+6*u(i,hampnt)...
                   -4*u(i,hampnt-1)+u(i,hampnt-2))/h^4;
            bn=-(2*uh(i)-2*uh(i-1)-(2*u(i,hampnt)+dxxb+dxxxxb)+2*u(i-1,hampnt));            
    end
             

    rn=an;
    while abs(G)>eps
        phir=Kham/(alpha+1)*(max(0,rn+an)^(alpha+1));
        G=rn+m*((phir-phia)/rn)+bn;
        Gp=1+(m*Kham*((alpha+1)*rn*max(0,rn+an)^(alpha+1)...
            -max(0,an)^(alpha+1))/(rn^2*(alpha+1)))...
            +m*Kham*max(0,an)^(alpha+1)/(rn^2*(alpha+1)^2);
        rn=rn-G/Gp;
    end
    phir=Kham/(alpha+1)*(max(0,rn+an)^(alpha+1));
    f(i)=(phir-phia)/rn;
     
        
    % Calculates next step in time
    switch Effect 
        case {'hamstr', 'hamlss', 'hamstf'}
            u(i+1,:) = A\(B*u(i,:)'+C*u(i-1,:)' + f(i)*a);
        case 'hamlss2'
            u(i+1,:) = (A\(B*u(i,:)'+C*u(i-1,:)'))+((k^2)/(ps))*(A\g)*f(i);
    end

    %Updates hammer position
    uh(i+1)=2*uh(i)-uh(i-1)-(k^2/m)*f(i);
    
    switch Synthesize
        case 'plot'
            gham=zeros(1,N+1);
            gham(hampnt)=uh(i);
            plot(x,u(i,:),x,gham,'o');

            axis([0,1,-0.1,0.1]);
            drawnow
        case 'play'
            % Vibration of one point on string
            out(i) = u(i,round(N/7));
        case 'plotsound'
            % Vibration of one point on string
            out(i) = u(i,round(N/7));
    end

end

G=1;




switch Synthesize
    case 'play'
        soundsc(out)
    case 'plotsound'
        soundsc(out)
        plot(out)
end
        
end
