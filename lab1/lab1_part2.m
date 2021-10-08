clc
clear all 
close all 


beta = 0.2;
gamma = 0.1;
n = 10; 
tau = 25;
tend = 1505; 

x = mackeyGlass(beta,gamma,n,tau,tend);
t = 301:1500;
plot(t,x(t))
title('mackeyGlass')
xlabel('t')
ylabel('x(t)')

input = [x(t-20); x(t-15); x(t-10); x(t-5); x(t)];
output = x(t+5);

ntrain = 700;
nvali = 1000;
ntest = 1200;

patterns = input(:,1:ntrain);
targets = output(:,1:ntrain);

valiD = input(:,ntrain:nvali);
valiT = output(:,ntrain:nvali);

testD = input(:,nvali:ntest);
testT = output(:,nvali:ntest);


options =  trainingOptions('rmsprop',)

trainNetwork(blbla, layers, options )






