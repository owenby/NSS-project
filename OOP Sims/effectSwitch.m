function [A, B, C] = effectSwitch(Effect, b0, b1, c, h, k, N, kappa, theta, ps)

% This function creates the matrices for the finite difference scheme
% solver, each case relates to a desired effect on the string. 'fdl' stands
% for frequency dependent loss. 

% Setting up matrices for solving the finite difference scheme
nOnes = ones(N+1, 1) ;
Dxx = (1/h^2)*(-2*diag(nOnes, 0) + diag(nOnes(1:N), -1) + diag(nOnes(1:N), 1));
Dxxxx = (1/h^4)*(diag(nOnes(1:N-1), -2)-4*diag(nOnes(1:N), -1) + 6*diag(nOnes, 0) - 4*diag(nOnes(1:N), 1) + diag(nOnes(1:N-1), 2));
Mx = (1/2) * (diag(nOnes(1:N), -1) + diag(nOnes(1:N), 1));

% Making the matrices sparse martrices for optimality 
Dxx = sparse(Dxx);
Dxxxx = sparse(Dxxxx);
Mx = sparse(Mx);

switch Effect
    case 'none'
        A = eye(N+1,N+1);
        B = 2*eye(N+1,N+1) + c^2*k^2*Dxx;
        C = -eye(N+1,N+1);
    case 'loss'
        A = (1 + b0*k/2) * eye(N+1,N+1);
        B = 2*eye(N+1,N+1) + c^2*k^2*Dxx;       
        C = (-1 + b0*k/2) * eye(N+1,N+1);
    case 'realistic'        
        A = (1 + b0*k/2) * eye(N+1,N+1) - b1*k/2*Dxx;
        B = 2*eye(N+1,N+1) + c^2*k^2*Dxx;
        C = (-1 + b0*k/2) * eye(N+1,N+1) - b1*k/2*Dxx;   
    case 'bar'
        A = theta * eye(N+1,N+1) + (1-theta) * Mx ;
        B = -kappa^2 * k^2 * Dxxxx + 2*A;
        C = -A;
    case 'stiff'
        A = theta * eye(N+1,N+1) + (1-theta) * Mx ;
        B = k^2*c^2*Dxx - kappa^2*k^2*Dxxxx + 2*A;
        C = -A;
    case 'fdl'
        A = -(1+b0*k)*eye(N+1,N+1) + b1*k*Dxx;
        B = -2*eye(N+1,N+1) - c^2*k^2*Dxx + kappa^2*k^2*Dxxxx;
        C = (1+b0*k)*eye(N+1,N+1) + b1*k*Dxx;
    case 'hamstr'
        A = eye(N+1,N+1);
        B = 2*eye(N+1,N+1) + c*k^2*Dxx;
        C = -A;
    case 'hamlss'
        A = (1+(b0*k)/(2*ps))*eye(N+1,N+1);
        B = 2*eye(N+1,N+1) + c*k^2*Dxx;
        C = ((b0*k)/(2*ps)-1)*eye(N+1,N+1);
    case 'hamstf'
        A = eye(N+1,N+1);
        B = 2*A+c*k^2*Dxx+((kappa^2*k^2)/ps)*Dxxxx;
        C = -A;
end

% Making the matrices sparse martrices for optimality 
A = sparse(A);
B = sparse(B);
C = sparse(C);

end

