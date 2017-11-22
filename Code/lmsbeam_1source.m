format compact; clear; clc; close all;

L = 16;

w0 = 2*pi/2048;
mu = 0.01;
dth = 20;

w = ones(L,1)*0.01;
w2= w;
N = 1000; % Number of look directions
dt = 1/N; % phase unwrap in sin(theta) domain
u = (-N:N)'; % increment in sin(theta) domain
theta=asin(u*dt); % theta vector [deg]
a = exp(1j*pi*dt*u*(0:(L-1)));
e = zeros(1,1e5);

for K = [64]
    figure('units','normalized','outerposition',[0 0 1 1])
    d1 = sign(rand(K,1)-0.5);
    d2 = sign(rand(K,1)-0.5);
    for n = 1:1e5
        th1 = 30+dth*sin(w0*n);
        th2 = 30-dth*sin(w0*n);
        psi1 = pi*sind(th1);
        psi2 = pi*sind(th2);
        [S1,~,~] = ArrayMeasurementPlusNoiseGenerator(10,psi1,L);
        [S2,~,~] = ArrayMeasurementPlusNoiseGenerator(10,psi2,L);
        
        S = S1*d1(mod(n+K-1,K)+1) + S2*d2(mod(n+K-1,K)+1);
        
        %% Target 1
        y = w(:,n)'*S;
        e(n) = d1(mod(n+K-1,K)+1)-y;
        g = -2*S*e(n)';
        w(:,n+1) = w(:,n)-mu*g;
        look = (abs(a*w(:,n)));
        look = look./max(look);
        
        %% Target 2
        y2 = w2(:,n)'*S;
        e2(n) = d2(mod(n+K-1,K)+1)-y2;
        g2 = -2*S*e2(n)';
        w2(:,n+1) = w2(:,n)-mu*g2;
        look2 = (abs(a*w2(:,n)));
        look2 = look2./max(look2);
        
        %% Plotting
        if mod(n-1,10)==0
            subplot(121)
            polarplot(th1*pi/180,1,'*','linewidth',2)
            hold on
            polarplot(th2*pi/180,1,'*','linewidth',2)
            polarplot(theta,look,'linewidth',2)
            polarplot(theta,look2,'linewidth',2)
            hold off
            
            subplot(122)
            loglog((1:n)-1,abs(e(1:n)./max(abs(e(1:n)))),'-*','linewidth',2)
            hold on
            loglog((1:n)-1,abs(e2(1:n)./max(abs(e2(1:n)))),'-*','linewidth',2)
            grid on
            hold off
            drawnow
        end
    end
end