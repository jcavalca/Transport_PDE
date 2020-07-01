function [x, t, w] = advectionEquation1DLW(c, f, xdom, tdom, M, N)

% x, t, w
% Discretize domain in time and space

x = linspace(xdom(1), xdom(2), M)';
t = linspace(tdom(1), tdom(2), N)';

% Define step sizes in space and time

h = (xdom(2) - xdom(1))/(M-1);
k = (tdom(2) - tdom(1))/(N-1);

sigma = c*k/h;

% Initialize temperature matrix

w = zeros(N, M);


%Set the initial condition


w(1, :) = f(x);

A = zeros(M, M);

for i = 1:1:M
    if i == 1
        A(1, 2) = sigma*(-1 + sigma)/2;
        A(1, end - 1) = sigma*(1 + sigma)/2;
    elseif i == M
        A(end,2) = sigma*(-1 + sigma)/2;
        A(end, end - 1) = sigma*(1 + sigma)/2;
    else 
        A(i, i-1) = sigma*(1 + sigma)/2;
        A(i, i+1) = sigma*(-1 + sigma)/2;
    end
    A(i, i) = (1-sigma^2);
    
end

for j = 1:1:N-1

    wj = w(j, 1:M)';
    
    w(j+1,1:M) = A*wj;
end

end