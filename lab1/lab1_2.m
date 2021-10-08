clc
clear all 
close all 

rng(0)
ndata = 100;
mA = [ 1.0, 0.3]; sigmaA = 0.2;
mB = [ 0.0, -0.1]; sigmaB = 0.3;
classA(1,:) = [ randn(1,round(0.5*ndata)) .* sigmaA - mA(1), randn(1,round(0.5*ndata)) .* sigmaA + mA(1)];
classA(2,:) = randn(1,ndata) .* sigmaA + mA(2);
classB(1,:) = randn(1,ndata) .* sigmaB + mB(1);
classB(2,:) = randn(1,ndata) .* sigmaB + mB(2);

plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')


%% 25 % of all data 
% delta rule 

n = round(ndata*0.75);

rem = randperm(ndata,ndata*0.25);

classA(:,rem) = [];
classB(:,rem) = [];

patterns1 = [classA(:,1:n),classB(:,1:n)];
targets1 = [ones(1,n),-ones(1,n)];


X1 = [patterns1; ones(1,2*n)];
T1 = targets1;

W = randn(1,3)*0.001;

eta = 0.01;
epochs = 20;
subplot(1,2,1)
plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on 

x = linspace(-4,4,10);
y = -(W(1)*x+W(3))/W(2);
axis([-2 2 -2 2])
plot(x,y)
drawnow

error = zeros(epochs+1,1);
error(1) = mean(abs(sign(W*X1)-T1)/2);

for i = 1:epochs
    deltaW = -eta*((W*X1)-T1)*X1';
    W = W + deltaW;
    error(i+1) = mean(abs(sign(W*X1)-T1)/2);

    
    y = -(W(1)*x+W(3))/W(2);
    plot(x,y)
    title(i)
    drawnow
end

subplot(1,2,2)
plot(error)
axis([0 epochs -0.1 1])

title('Error Delta rule')

%% 50 % of class A 


n = ndata;

rem = randperm(ndata,ndata*0.5);

classA(:,rem) = [];

patterns1 = [classA,classB];
targets1 = [ones(1,length(classA)),-ones(1,n)];


X1 = [patterns1; ones(1,(length(classA)+length(classB)))];
T1 = targets1;

W = randn(1,3)*0.001;

eta = 0.01;
epochs = 20;
subplot(1,2,1)
plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on 

x = linspace(-4,4,10);
y = -(W(1)*x+W(3))/W(2);
axis([-2 2 -2 2])
plot(x,y)
drawnow

error = zeros(epochs+1,1);
error(1) = mean(abs(sign(W*X1)-T1)/2);

for i = 1:epochs
    deltaW = -eta*((W*X1)-T1)*X1';
    W = W + deltaW;
    error(i+1) = mean(abs(sign(W*X1)-T1)/2);

    
    y = -(W(1)*x+W(3))/W(2);
    plot(x,y)
    title(i)
    drawnow
end

subplot(1,2,2)
plot(error)
axis([0 epochs -0.1 1])

title('Error Delta rule')

%%  50% class B 


n = ndata;

rem = randperm(ndata,ndata*0.5);

classB(:,rem) = [];

patterns1 = [classA,classB];
targets1 = [ones(1,n),-ones(1,length(classB))];


X1 = [patterns1; ones(1,(length(classA)+length(classB)))];
T1 = targets1;

W = randn(1,3)*0.001;

eta = 0.01;
epochs = 20;
subplot(1,2,1)
plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on 

x = linspace(-4,4,10);
y = -(W(1)*x+W(3))/W(2);
axis([-2 2 -2 2])
plot(x,y)
drawnow

error = zeros(epochs+1,1);
error(1) = mean(abs(sign(W*X1)-T1)/2);

for i = 1:epochs
    deltaW = -eta*((W*X1)-T1)*X1';
    W = W + deltaW;
    error(i+1) = mean(abs(sign(W*X1)-T1)/2);

    
    y = -(W(1)*x+W(3))/W(2);
    plot(x,y)
    title(i)
    drawnow
end

subplot(1,2,2)
plot(error)
axis([0 epochs -0.1 1])

title('Error Delta rule')



%% 20 % 80 % of class A


n = ndata;

rem = [randperm(ndata*0.5,ndata*0.5*0.20), randperm(ndata*0.5,ndata*0.5*0.80)+50];

classA(:,rem) = [];

patterns1 = [classA,classB];
targets1 = [ones(1,length(classA)),-ones(1,length(classB))];


X1 = [patterns1; ones(1,(length(classA)+length(classB)))];
T1 = targets1;

W = randn(1,3)*0.001;

eta = 0.01;
epochs = 20;
figure(4)
subplot(1,2,1)
plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on 

x = linspace(-4,4,10);
y = -(W(1)*x+W(3))/W(2);
axis([-2 2 -2 2])
plot(x,y)
drawnow

error = zeros(epochs+1,1);
error(1) = mean(abs(sign(W*X1)-T1)/2);

for i = 1:epochs
    deltaW = -eta*((W*X1)-T1)*X1';
    W = W + deltaW;
    error(i+1) = mean(abs(sign(W*X1)-T1)/2);

    
    y = -(W(1)*x+W(3))/W(2);
    plot(x,y)
    title(i)
    drawnow
end

subplot(1,2,2)
plot(error)
axis([0 epochs -0.1 1])

title('Error Delta rule')





















