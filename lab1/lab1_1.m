clc
clear all 
close all 

rng(0)

% 1.1

n = 100;
mA = [ 2.0, 0.5]; sigmaA = 0.5;
mB = [-2.0, 0.0]; sigmaB = 0.5;
classA(1,:) = randn(1,n) .* sigmaA + mA(1);
classA(2,:) = randn(1,n) .* sigmaA + mA(2);
classB(1,:) = randn(1,n) .* sigmaB + mB(1);
classB(2,:) = randn(1,n) .* sigmaB + mB(2);


%plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')


%% 1.2
% delta rule 

patterns1 = [classA,classB];
targets1 = [ones(1,n),-ones(1,n)];

X1 = [patterns1; ones(1,2*n)];
T1 = targets1;

W = randn(1,3)*0.001;

eta = 0.001;
epochs = 20;
subplot(2,2,1)
plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on 

x = linspace(-4,4,10);
y = -(W(1)*x+W(3))/W(2);
axis([-4 4 -4 4])
plot(x,y)
drawnow

error = zeros(epochs+1,1);
error(1) = mean(abs(sign(W*X1)-T1)/2);

for i = 1:epochs
    deltaW = -eta*((W*X1)-T1)*X1';
    W = W + deltaW;
    error(i+1) = mean(abs(sign(W*X1)-T1)/2);

    
    y = -(W(1)*x+W(3))/W(2);
    axis([-4 4 -4 4])
    plot(x,y)
    title(i)
    drawnow
end

subplot(2,2,2)
plot(error)
axis([0 epochs -1 2])

title('Error Delta rule')

% WHY WEIGHTS A LINE? 


%% preception 


patterns2 = [classA,classB];
targets2 = [ones(1,n),-ones(1,n)];

X2 = [patterns2; ones(1,2*n)];
T2 = targets2;

W = randn(1,3)*0.01;
x = linspace(-4,4,10);

subplot(2,2,3)
eta = 0.001;
epochs = 20;
plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on

x = linspace(-4,4,10);
y = -(W(1)*x+W(3))/W(2);
axis([-4 4 -4 4])
plot(x,y)
drawnow

error = zeros(epochs+1,1);
error(1) = mean(abs(sign(W*X2)-T2)/2);

for i = 1:epochs
    deltaW = -eta*(sign(W*X2)-T2)*X2'/2; % changed here from delta rule
    error(i+1) = mean(abs(sign(W*X2)-T2)/2);
    
    W = W + deltaW;
    
    y = -(W(1)*x+W(3))/W(2);
    %subplot(2,2,3)
    axis([-4 4 -4 4])
    plot(x,y)
    title(i)
    drawnow
end

subplot(2,2,4)
plot(error)
axis([0 epochs -1 2])
title('Error preception')


%% Sequential learning
% delta rule 

patterns3 = [classA,classB];
targets3 = [ones(1,n),-ones(1,n)];

X3 = [patterns3; ones(1,2*n)];
T3 = targets3;

W = randn(1,3)*0.01;

eta = 0.001;
epochs = 5;
figure(2)
subplot(1,2,1)
plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on 

x = linspace(-4,4,10);
y = -(W(1)*x+W(3))/W(2);
axis([-4 4 -4 4])
plot(x,y)
drawnow

error = zeros(epochs+1,1);
error(1) = mean(abs(sign(W*X3)-T3)/2);

for i = 1:epochs
    
    for j = 1:length(patterns3)
        
        deltaW = -eta*((W*X3(:,j))-T3(j))*X3(:,j)';
        
        W = W + deltaW;
 
        y = -(W(1)*x+W(3))/W(2);
        axis([-4 4 -4 4])
        plot(x,y)
        drawnow
        
    end
    error(i+1) = mean(abs(sign(W*X3)-T3)/2);
    
end

subplot(1,2,2)
plot(error)
axis([0 epochs -1 2])

title('Error Delta rule')




%% 1.2
% delta rule, no bias

patterns4 = [classA,classB];
targets4 = [ones(1,n),-ones(1,n)];

X4 = patterns4;%[patterns1; ones(1,2*n)];
T4 = targets4;

W = randn(1,2)*0.001;

eta = 0.01;
epochs = 20;
subplot(1,2,1)
plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on 

x = linspace(-4,4,10);
y = -(W(1)*x)/W(2);
axis([-4 4 -4 4])
plot(x,y)
drawnow

error = zeros(epochs+1,1);
error(1) = mean(abs(sign(W*X4)-T4)/2);

for i = 1:epochs
    deltaW = -eta*((W*X4)-T4)*X4';
    W = W + deltaW;
    error(i+1) = mean(abs(sign(W*X4)-T4)/2);

    
    y = -(W(1)*x)/W(2);
    axis([-4 4 -4 4])
    plot(x,y)
    title(i)
    drawnow
end

subplot(1,2,2)
plot(error)
axis([0 epochs -1 2])
title('Error Delta rule')






