close all
clear

arr_3d = preProcess("NHMS.jpeg");

bit_stream = bitStreamConvert(arr_3d, 4);

alpha = 0.5;
K = 6;
n = 32; % Half Sine Samples
A1 = 1;
A2 = 1;
T = 1;
t1  = linspace(0, 1, n+1);
t1 = t1(1:end-1);
t2 = linspace(-1*K*T, K*T, (2*K*n)+1);
t2 = t2(1:end-1);

half_sine = sin(pi*t1/T);
half_sine = half_sine/norm(half_sine);

x_num = sin(pi*t2*(1-alpha)/T)+(4*alpha*t2.*cos(pi*t2*(1+alpha)/T)/T);
x_den = pi*t2*(1/T).*(1-(4*alpha*t2/T).^2);
x_t = x_num ./ x_den;

f = 0:0.01:2*((1+alpha)/(2*T));
x_f = zeros(size(f));
x_f(0<= abs(f) & abs(f)<= (1-alpha)/(2*T)) = sqrt(T);
x_f((1-alpha)/(2*T)<=abs(f) & abs(f)<=(1+alpha)/(2*T)) = sqrt(T)*cos((pi*T/(2*alpha))*(abs(f((1-alpha)/(2*T)<=abs(f) & abs(f)<=(1+alpha)/(2*T)))-((1-alpha)/(2*T))));


% x_t(0) = 1-alpha + (4 * alpha / pi);
% x_t(abs(t2) == T / (4 * alpha)) = (alpha/sqrt(2))*((1+(2/pi)*sin(pi/(4*alpha)))+(1-(2/pi))*cos(pi/(4*alpha)));

% srrc = conv(sym_up, A2 * x_t);

srrc = rcosdesign(alpha, 2*K, n);
srrc = srrc(1:end-1);
srrc = srrc/norm(srrc);

% Question 1
Nfft = 2^nextpow2(4 * max(length(half_sine), length(srrc)));
f_axis = (-Nfft/2 : Nfft/2-1) / Nfft;

half_sine_freq = fftshift(fft(half_sine, Nfft));
srrc_freq = fftshift(fft(srrc, Nfft));

figure()
set(gcf, "Theme", "light");
sgtitle("Pulse Shaping Functions")
subplot(2, 1, 1)
plot(t1, half_sine);
title("Half-Sine")
xlabel("Time (s)")
grid on

subplot(2, 1, 2)
plot(t2, srrc)
title("SRRC")
xlabel("Time (s)")
grid on

figure()
subplot(2, 1, 1)
set(gcf, "Theme", "light");
plot(f_axis, 20*log10(abs(half_sine_freq)));
ylabel("Magnitude (dB)")
xlabel("Frequency (Hz)")
grid
subplot(2, 1, 2)
plot(f_axis, angle(half_sine_freq))
ylabel("Phase Response")
xlabel("Frequency (Hz)")
% xlabel("Frequency (Hz)")
grid on
sgtitle("Frequency Response of Half-Sine")

figure()
set(gcf, "Theme", "light");
subplot(2, 1, 1)
plot(f_axis, 20*log10(abs(srrc_freq)))
grid
ylabel("Magnitude (dB)")
xlabel("Frequency (Hz)")
subplot(2, 1, 2);
plot(f_axis, angle(srrc_freq))
ylabel("Phase Response")
xlabel("Frequency (Hz)")
sgtitle("Frequency Response of SRRC")
% xlabel("Frequency (Hz)")
grid on

% Question 2

bits = randi([0, 1], 1, 10);

symbols = 2*bits - 1;
sym_up = zeros(1, length(bits) * n);
sym_up(1:n:end) = symbols;

half_sine_mod = conv(sym_up, half_sine);
srrc_mod = conv(sym_up, srrc);

figure()
set(gcf, "Theme", "light")
subplot(2, 1, 1)
plot(half_sine_mod);
title("Half-Sine")
xlabel("Samples")
grid
subplot(2, 1, 2)
plot(srrc_mod)
title("SRRC")
xlabel("Samples")
grid
sgtitle("Signals Modulated for 10 bits")

% Question 3
Nfft = 2^nextpow2(4 * max(length(half_sine_mod), length(srrc_mod)));
f_axis = (-Nfft/2 : Nfft/2-1) / Nfft;

half_sine_freq = fftshift(fft(half_sine_mod, Nfft));
srrc_freq = fftshift(fft(srrc_mod, Nfft));

figure()
subplot(2, 1, 1)
set(gcf, "Theme", "light");
plot(f_axis, 20*log10(abs(half_sine_freq)));
ylabel("Magnitude (dB)")
xlabel("Frequency (Hz)")
grid
subplot(2, 1, 2)
plot(f_axis, angle(half_sine_freq))
ylabel("Phase Response")
xlabel("Frequency (Hz)")
grid on
sgtitle("Frequency Response of Modulated Half-Sine")

figure()
set(gcf, "Theme", "light");
subplot(2, 1, 1)
plot(f_axis, 20*log10(abs(srrc_freq)))
grid
ylabel("Magnitude (dB)")
xlabel("Frequency (Hz)")
subplot(2, 1, 2);
plot(f_axis, angle(srrc_freq))
ylabel("Phase Response")
xlabel("Frequency (Hz)")
sgtitle("Frequency Response of Modulated SRRC")
grid on

% Question 4
figure()
set(gcf, "Theme", "light")
eyediagram(half_sine_mod(17:end), n);

figure()
set(gcf, "Theme", "light")
srrc_eye = eyediagram(srrc_mod, n);

% Question 5

h = [1 0.5 0.75 (-2/7)];
h_upsample = upsample(h, n);

half_sine_channel = conv(h_upsample, half_sine_mod);
srrc_channel = conv(h_upsample, srrc_mod);

figure()
set(gcf, "Theme", "light")
freqz(h)

figure()
eyediagram(half_sine_channel(17:end), n)

figure()
eyediagram(srrc_channel, n)

% Question 7
sigma = 1;

noise = sigma * randn(size(half_sine_channel));
eyediagram(noise(17:end) + half_sine_channel(17:end), n);

sigma = 0.1;

noise = sigma * randn(size(half_sine_channel));
eyediagram(noise(17:end) + half_sine_channel(17:end), n);

sigma = 0.01;

noise = sigma * randn(size(half_sine_channel));
eyediagram(noise(17:end) + half_sine_channel(17:end), n);

sigma = 1;

noise = sigma * randn(size(srrc_channel));
eyediagram(noise + srrc_channel, n);

sigma = 0.1;

noise = sigma * randn(size(srrc_channel));
eyediagram(noise + srrc_channel, n);

sigma = 0.01;

noise = sigma * randn(size(srrc_channel));
eyediagram(noise + srrc_channel, n);


function arr_3d = preProcess(imageName)
    image = imread(imageName);
    image_arr = im2double(im2gray(image));
    
    fun = @(block_struct) dct2(block_struct.data);
    image_dct = blockproc(image_arr, [8 8], fun);
    
    image_max = max(image_dct, [], 'all');
    image_min = min(image_dct, [], 'all');
    
    scaled_arr = (image_dct - image_min) / (image_max - image_min);

    [m, n] = size(scaled_arr);
    arr_3d = reshape(scaled_arr, [8 8 (m*n)/64]);
end

function bit_stream = bitStreamConvert(arr_3d, N)
    col_group = reshape(arr_3d(:, :, 1:N), [64*N 1]);
    col_int = round(col_group *  255);
    bit_stream = reshape(int2bit(col_int, 8), [1 64*8*N]);

end