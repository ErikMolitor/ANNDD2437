clc
clear all 
close all

%{
rng(0)
n = 100;
mA = [ 2.0, 0.5]; sigmaA = 0.5;
mB = [-2.0, 0.0]; sigmaB = 0.5;
classA(1,:) = randn(1,n) .* sigmaA + mA(1);
classA(2,:) = randn(1,n) .* sigmaA + mA(2);
classB(1,:) = randn(1,n) .* sigmaB + mB(1);
classB(2,:) = randn(1,n) .* sigmaB + mB(2);
%}

x=[-5:0.5:5]';
y=[-5:0.5:5]';
z=exp(-x.*x*0.1) * exp(-y.*y*0.1)' - 0.5;
figure(1)
mesh(x, y, z);
title('Target')

ndata = length(z)*length(z);

targets = reshape (z, 1, ndata);
[xx, yy] = meshgrid (x, y);
patterns = [reshape(xx, 1, ndata); reshape(yy, 1, ndata)];


epochs = 100; 

[nin,~] = size(patterns);
[nout,~] = size(targets);
nhidden = 10;

w = randn(nhidden,nin+1);
v = randn(nout,nhidden+1);

alpha = 0.9;
dw = 0;%zeros(nhidden,nin+1);
dv = 0; %zeros(nout,nhidden+1);
eta = 0.005;
gridsize = length(x);

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
    figure(2)

    
    zz = reshape(out, gridsize, gridsize);
    mesh(x,y,zz);
    axis([-5 5 -5 5 -0.7 0.7]);
    title(['100 epochs. 10 Nodes.'])
    drawnow;
    
    
end







