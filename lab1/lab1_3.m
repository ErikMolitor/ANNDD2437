clc
clear all 
close all


rng(0)

n = 100;
mA = [ 0.5, 0.5]; sigmaA = 0.5;
mB = [-1.0, 0.0]; sigmaB = 0.5;
classA(1,:) = randn(1,n) .* sigmaA + mA(1);
classA(2,:) = randn(1,n) .* sigmaA + mA(2);
classB(1,:) = randn(1,n) .* sigmaB + mB(1);
classB(2,:) = randn(1,n) .* sigmaB + mB(2);


%% 3.2.1 uppgifr 1

ndata = length(classA) +length(classB);

patterns = [classA,classB];
targets = [ones(1,length(classA)),-ones(1,length(classB))];

epochs = 20; 

[nin,~] = size(patterns);
[nout,~] = size(targets);
nhidden = 10;

w = randn(nhidden,nin+1);
v = randn(nout,nhidden+1);

alpha = 0.9;
dw = 0;%zeros(nhidden,nin+1);
dv = 0; %zeros(nout,nhidden+1);
eta = 0.001;


for i=1:epochs
    
    hin = w * [patterns ; ones(1,ndata)];
    hout = [2 ./ (1+exp(-hin)) - 1 ; ones(1,ndata)];
    oin = v * hout;
    out = 2 ./ (1+exp(-oin)) - 1;

    delta_o = (out - targets) .* ((1 + out) .* (1 - out)) * 0.5;
    delta_h = (v' * delta_o) .* ((1 + hout) .* (1 - hout)) * 0.5;
    delta_h = delta_h(1:nhidden, :);
    
    dw = (dw .* alpha) - (delta_h * [patterns; ones(1,ndata)]') .* (1-alpha);
    dv = (dv .* alpha) - (delta_o * hout') .* (1-alpha);
    w = w + dw .* eta;
    v = v + dv .* eta;    
    
    
    
    
end

plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on

x = linspace(-4,4,10);

for i = 1:nhidden
    
    y = -(w(i,1)*x+w(i,3))/w(i,2);
    axis([-4 4 -4 4])
    plot(x,y)
    hold on 
    drawnow
end
%{
for i = 1:nout
    
    y = -(v(i,1)*x+v(i,3))/v(i,2);
    axis([-4 4 -4 4])
    plot(x,y)
    hold on 
    drawnow
end
%}



%% 3.2.1 uppgf 2 80 % traning 20% test


%n = round(ndata*0.80);
ndata = length(classA) +length(classB);
rem = randperm(length(classA),ndata*0.20);

valiA = classA(:,rem);
valiB = classB(:,rem);
nvali = length(valiA)+length(valiB);
tvali = [ones(1,length(valiA)),-ones(1,length(valiB))];

classA(:,rem) = [];
classB(:,rem) = [];

ndata = length(classA) +length(classB);

patterns = [classA,classB];
targets = [ones(1,length(classA)),-ones(1,length(classB))];
epochs = 200; 

[nin,~] = size(patterns);
[nout,~] = size(targets);
nhidden = 2;

w = randn(nhidden,nin+1);
v = randn(nout,nhidden+1);

alpha = 0.9;
dw = 0; %zeros(nhidden,nin+1);
dv = 0; %zeros(nout,nhidden+1);
eta = 0.001;


error = zeros(epochs,1);

for i=1:epochs
    
    %traning porcess
    hin = w * [patterns ; ones(1,ndata)];
    hout = [2 ./ (1+exp(-hin)) - 1 ; ones(1,ndata)];
    oin = v * hout;
    out = 2 ./ (1+exp(-oin)) - 1;

    delta_o = (out - targets) .* ((1 + out) .* (1 - out)) * 0.5;
    delta_h = (v' * delta_o) .* ((1 + hout) .* (1 - hout)) * 0.5;
    delta_h = delta_h(1:nhidden, :);
    
    dw = (dw .* alpha) - (delta_h * [patterns; ones(1,ndata)]') .* (1-alpha);
    dv = (dv .* alpha) - (delta_o * hout') .* (1-alpha);
    w = w + dw .* eta;
    v = v + dv .* eta;
    
    % validation process
    
    hin = w * [valiA, valiB ; ones(1,nvali)];
    hout = [2 ./ (1+exp(-hin)) - 1 ; ones(1,nvali)];
    oin = v * hout;
    out = 2 ./ (1+exp(-oin)) - 1;
    
    error(i) = mean((out-tvali).^2);
    
    
    
    
end

%{
figure(2)
plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on

x = linspace(-4,4,10);

for i = 1:nhidden
    
    y = -(w(i,1)*x+w(i,3))/w(i,2);
    axis([-4 4 -4 4])
    plot(x,y)
    hold on 
    drawnow
end

for i = 1:nout
    
    y = -(v(i,1)*x+v(i,3))/v(i,2);
    axis([-4 4 -4 4])
    plot(x,y)
    hold on 
    drawnow
end

%}

figure(3)
plot(error)
axis([0 epochs -1 2])
title('Error')
xlabel('epochs')
ylabel('error')



%% 3.2.1 uppgf 2 20 % traning 80% test


%n = round(ndata*0.80);
rem = randperm(length(classA),length(classA)*0.80);

valiA = classA(:,rem);
valiB = classB(:,rem);
nvali = length(valiA)+length(valiB);
tvali = [ones(1,length(valiA)),-ones(1,length(valiB))];

classA(:,rem) = [];
classB(:,rem) = [];

ndata = length(classA) +length(classB);

patterns = [classA,classB];
targets = [ones(1,length(classA)),-ones(1,length(classB))];
epochs = 200; 

[nin,~] = size(patterns);
[nout,~] = size(targets);
nhidden = 2;

w = randn(nhidden,nin+1);
v = randn(nout,nhidden+1);

alpha = 0.9;
dw = 0; %zeros(nhidden,nin+1);
dv = 0; %zeros(nout,nhidden+1);
eta = 0.001;


error = zeros(epochs,1);

for i=1:epochs
    
    %traning porcess
    hin = w * [patterns ; ones(1,ndata)];
    hout = [2 ./ (1+exp(-hin)) - 1 ; ones(1,ndata)];
    oin = v * hout;
    out = 2 ./ (1+exp(-oin)) - 1;

    delta_o = (out - targets) .* ((1 + out) .* (1 - out)) * 0.5;
    delta_h = (v' * delta_o) .* ((1 + hout) .* (1 - hout)) * 0.5;
    delta_h = delta_h(1:nhidden, :);
    
    dw = (dw .* alpha) - (delta_h * [patterns; ones(1,ndata)]') .* (1-alpha);
    dv = (dv .* alpha) - (delta_o * hout') .* (1-alpha);
    w = w + dw .* eta;
    v = v + dv .* eta;
    
    % validation process
    
    hin = w * [valiA, valiB ; ones(1,nvali)];
    hout = [2 ./ (1+exp(-hin)) - 1 ; ones(1,nvali)];
    oin = v * hout;
    out = 2 ./ (1+exp(-oin)) - 1;
    
    error(i) = mean((out-tvali).^2);
    
    
    
    
end

%{
figure(2)
plot(classA(1,:),classA(2,:), 'o',classB(1,:),classB(2,:), 'o')
hold on

x = linspace(-4,4,10);

for i = 1:nhidden
    
    y = -(W(i,1)*x+W(i,3))/W(i,2);
    axis([-4 4 -4 4])
    plot(x,y)
    hold on 
    drawnow
end

for i = 1:nout
    
    y = -(V(i,1)*x+V(i,3))/V(i,2);
    axis([-4 4 -4 4])
    plot(x,y)
    hold on 
    drawnow
end

%}

figure(3)
plot(error)
axis([0 epochs -1 2])
title('Error')
xlabel('epochs')
ylabel('error')



