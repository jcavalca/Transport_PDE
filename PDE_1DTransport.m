clc; close all; clear;

%% Problem 3 - Joao Cavalcanti

c = 0.5;
L = 2;
f = @(x) exp(-100*(x-1).^2);
h = 0.01;
k = 0.005;

xdom = [0 2];
tdom = [0 5];

M = (xdom(2) - xdom(1))/h + 1;
N = (tdom(2) - tdom(1))/k + 1;

[x_1,t_1, w_1] = advectionEquation1DFTCS(c, f, xdom, tdom, M, N);

figure(1)
plot(x_1, w_1(end, :))
title("@ t = 5 FTCS")


% The plot shows that for a large time (t = 5), the solution is not
% accurate because the graph is unconditionally stable. 

[x_2,t_2, w_2] = advectionEquation1DLF(c, f, xdom, tdom, M, N);


figure(2)

plot(x_2, w_2(end, :))
title("@ t = 5 Lax-Friedrichs")
snapnow

figure(3)
subplot(2, 2, 1)
plot(x_2, w_2(1, :))
title("@ t = 0 Lax-Friedrichs")

subplot(2, 2, 2)
plot(x_2, w_2(101, :))
title("@ t = 0.5 Lax-Friedrichs")

subplot(2, 2, 3)
plot(x_2, w_2(201, :))
title("@ t = 1 Lax-Friedrichs")

subplot(2, 2, 4)
plot(x_2, w_2(301, :))
title("@ t = 1.5 Lax-Friedrichs")
snapnow
% As t increases, the amplitude decreases and the wave curve flattens, 
% which shows the diffusion. 

[x_3,t_3, w_3] = advectionEquation1DLW(c, f, xdom, tdom, M, N);


figure(4)

plot(x_3, w_3(end, :))
title("@ t = 5 Lax-Wendroff")
snapnow

figure(5)
subplot(2, 2, 1)
plot(x_3, w_3(1, :))
title("@ t = 0 Lax-Wendroff")

subplot(2, 2, 2)
plot(x_3, w_3(101, :))
title("@ t = 0.5 Lax-Wendroff")

subplot(2, 2, 3)
plot(x_3, w_3(201, :))
title("@ t = 1 Lax-Wendroff")

subplot(2, 2, 4)
plot(x_3, w_3(301, :))
title("@ t = 1.5 Lax-Wendroff")
snapnow
% The solution does mitigate the diffusion problem, since our amplitude
% doesn't flatten out as t increases. 