clc
close all
clearvars
f = 500;
fs = 1000;
t = linspace(0,1,fs);


x = 10*cos(t*f);
x1 = cos((t*f)+pi/4);

q = x1.*x;
q1 = lowpass(q, 1,fs);

fig1 = figure(1);
plot(t,x)
hold on;grid on
plot(t,q1)

fig2 = figure(2);

plot((1:1000),abs(fft(q1)));
