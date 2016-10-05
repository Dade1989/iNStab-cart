close all
clear
clc


% Questo semplice script calcola il numero di Reynolds da inserire come input
% in PASTA (che è definito in base ai parametri relativi al solo tubo interno)
% a partire dalla nuova definizione di Reynolds (calcolato utilizzando la bulk
% velocity complessiva dei due tubi e il diametro del tubo interno).


% -> Dati:

Re_new = 2500 % Reynolds desiderato, basato sulla bulk velocity complessiva

u0i = 1; % massima velocità tubo interno [m/s], fissata!

ru = 1 % rapporto di velocità desiderato (ru = u0o/u0i)

Ri = 0.5; % raggio tubo interno [m]:
Ro1 = 0.6; % raggio inferiore tubo esterno [m]:
Ro2 = 1; % raggio superiore tubo esterno [m]

kb = 0.861379644; % coefficiente che lega la bulk velocity del profilo tanh
                  % con la velocità massima del profilo stesso, cioè:
                  % u_bulk = kb*u0;
                  % calcolato con Wolfram a partire dall'integrale analitico


% -> Calcolo velocità massima tubo esterno:

u0o = ru*u0i;


% -> Calcolo bulk velocity per i due tubi:

ubi = kb*u0i; % interno
ubo = kb*u0o; % esterno


% -> Calcolo bulk velocity complessiva:

ubt = ( 2*Ri*ubi + (Ro2-Ro1)*ubo ) / ( 2*Ri + (Ro2-Ro1) );


% -> Rapporto tra Reynolds vecchio (basato su u0i) e nuovo (basato su ubt):
%
%    Re_old = rho*u0i*Di/mu ; Re_new = rho*ubt*Di/mu

rr = u0i/ubt;


% -> Calcolo Re vecchio a partire da Re nuovo:

Re_old = rr*Re_new

% NB: Re_old, ru da utilizzare come input per PASTA!
