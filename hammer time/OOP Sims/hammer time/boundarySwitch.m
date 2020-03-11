function B = boundarySwitch(B, N, Boundary)

% This function sets the boudnary conditions for the string. 'xylo' means
% that it is fixed at some point in the middle, much like a xylophone is,
% with the fixed points having a fixed point either side of it.

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
    case 'xylo'  % This is bad
        B(5,:) = 0;
        B(N+1-4,:) = 0;
        B(4,:) = 0;
        B(N+1-5,:) = 0;
        B(6,:) = 0;
        B(N+1-3,:) = 0;
end

end
