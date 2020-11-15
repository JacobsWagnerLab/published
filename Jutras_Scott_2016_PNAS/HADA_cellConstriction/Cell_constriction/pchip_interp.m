function y_interp = pchip_interp(signal)
x = (1:length(signal)) / length(signal) * 100;
xhat = (1/3):(1/3):100;
y_interp = pchip(x, signal, xhat);