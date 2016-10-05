% Octave file to plot HG(omega)

close all ; clear all ; clc

%data = dlmread('./hgRe500.dat','',2,0);
data = dlmread('./step2_Re1000.dat','',2,0);

figure
%plot(data(:,1),data(:,2:end),'k')
semilogy(data(:,1)/2,data(:,2:end).^0.5,'k')
axis([0 2 1e-0 1e4])
%plot(data(:,1)/pi,data(:,2:end).^0.5,'k')
grid on
xlabel('St')      % xlabel('\omega')
ylabel('G_k')

figure
%plot(data(:,1),data(:,2:end),'k')
plot(data(:,1)/2,data(:,2:end).^0.5,'k')
axis([0 2 0 1e4])
%plot(data(:,1)/pi,data(:,2:end).^0.5,'k')
grid on
xlabel('St')      % xlabel('\omega')
ylabel('G_k')
