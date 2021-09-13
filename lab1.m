clc
clear all 
close all 

%% 1.1

n = 100;
mA = [ 2.0, 0.5]; sigmaA = 0.5;
mB = [-2.0, 0.0]; sigmaB = 0.5;
classA(1,:) = randn(1,n) .* sigmaA + mA(1);
classA(2,:) = randn(1,n) .* sigmaA + mA(2);
classB(1,:) = randn(1,n) .* sigmaB + mB(1);
classB(2,:) = randn(1,n) .* sigmaB + mB(2);


%plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')


%% 1.2

patterns = [classA,classB];
targets = [ones(1,n),-ones(1,n)];

X = [patterns; ones(1,2*n)];
T = targets;

W = randn(1,3)*0.01;

eta = 0.001;
epochs = 20;
plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on 
for i = 1:epochs
    deltaW = -eta*((W*X)-T)*X';
    W = W + deltaW;
    
    x = linspace(-4,4,10);
    y = -W(1)*(x*W(3)/W(2));
    plot(x,y)
    title(i)
    drawnow
end













