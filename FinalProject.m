close all
clear

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
h = [1 0.5 0.75 (-2/7)];

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

k_start_half_sine = n+1;
k_start_srrc = 2*K*n+1;

[arr_3d, image_max, image_min, image_size] = preProcess("NHMS.jpeg");

bit_stream = bitStreamConvert(arr_3d, size(arr_3d, 3));
num_bits = length(bit_stream);

mod_sig = signal_mod(bit_stream, n, half_sine);

channel_out = transmit_channel(h, n, mod_sig, 0);

% --- Receiver: matched filter + equalizer ---
mf_out = matched_filter(channel_out, half_sine);
eq_out = qz(mf_out, h);   % or call mmse(...) instead

% --- Sample and detect ---
% k_start: pulse delay + matched filter delay + channel/equalizer delay
% For half-sine + ZF (filter-based): k_start = length(half_sine) = 32
% but you need to verify empirically and add equalizer delay
k_start = length(half_sine);   % adjust based on your delay analysis
detected_bits = detect_bits(eq_out, k_start, n, num_bits);

% --- Reconstruct image ---
image_recovered = convert_to_image(detected_bits, image_max, image_min, image_size);


figure()
imshow(image_recovered)
title(sprintf('Recovered image, \\sigma = %.3f', sigma))

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

bits = randi([0, 1], 1, 1000);

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

% figure()
% eyediagram(half_sine_channel(17:end), n)
% 
% figure()
% eyediagram(srrc_channel, n)
% 
% % Question 7
% sigma = 1;
% 
% noise = sigma * randn(size(half_sine_channel));
% eyediagram(noise(17:end) + half_sine_channel(17:end), n);
% 
% sigma = 0.1;
% 
% noise = sigma * randn(size(half_sine_channel));
% eyediagram(noise(17:end) + half_sine_channel(17:end), n);
% 
% sigma = 0.01;
% 
% noise = sigma * randn(size(half_sine_channel));
% eyediagram(noise(17:end) + half_sine_channel(17:end), n);
% 
% sigma = 1;
% 
% noise = sigma * randn(size(srrc_channel));
% eyediagram(noise + srrc_channel, n);
% 
% sigma = 0.1;
% 
% noise = sigma * randn(size(srrc_channel));
% eyediagram(noise + srrc_channel, n);
% 
% sigma = 0.01;
% 
% noise = sigma * randn(size(srrc_channel));
% eyediagram(noise + srrc_channel, n);

% Question 8

half_sine_mf = flip(half_sine);
srrc_mf = flip(srrc);

x_impulse_response = [1, zeros(1, 500)];

h_half_sine = filter(half_sine, 1, x_impulse_response);
h_srrc = filter(srrc, 1, x_impulse_response);

figure()
subplot(2, 1, 1)
plot(h_half_sine)
title("Half Sine MF Impulse Response")
subplot(2, 1, 2)
plot(h_srrc)
title("SRRC MF Impulse Response")


figure()
freqz(half_sine_mf)
title("Half Sine MF")

figure()
freqz(srrc_mf)
title("SRRC MF")

% Question 9

sigma = 0.01;
noise = sigma * randn(size(half_sine_channel));

eyediagram(conv(half_sine_channel + noise, half_sine_mf), n)
title("Half Sine Matched Filter Output : 1 Bit")
eyediagram(conv(half_sine_channel + noise, half_sine_mf), 2*n)
title("Half Sine Matched Filter Output : 2 Bit")

matched_out_half_sine = conv(half_sine_channel + noise, half_sine_mf);

noise = sigma * randn(size(srrc_channel));
eyediagram(conv(srrc_channel + noise, srrc_mf), n)
title("SRRC Matched Filter Output : 1 Bit")
eyediagram(conv(srrc_channel + noise, srrc_mf), 2*n)
title("SRRC Matched Filter Output : 2 Bit")

matched_out_srrc = conv(srrc_channel + noise, srrc_mf);

% Question 10
q_zf = filter(1, h, [1, zeros(1, 10000)]);

figure()
plot(q_zf)
title("Impulse Response q_{zf}")

[Q_zf, w] = freqz(1, h);

figure()
freqz(1, h)
title("Frequency Response Q_{zf}")

% Question 11
% THIS IS WRONG, GO BACK AND REPLACE CHANNEL WITH THE NOISEY CHANNEL PASSED
% SIGNALS
sigma = 0;
noise = sigma * randn(size(half_sine_channel));
half_sine_mf_out_0 = conv(half_sine_channel + noise, half_sine_mf);
eq_sig = conv(half_sine_mf_out_0, q_zf);
eyediagram(eq_sig, n)
title("Half Sine EQ Sigma = 0")

sigma = 0.1;
noise = sigma * randn(size(half_sine_channel));
half_sine_mf_out_005 = conv(half_sine_channel + noise, half_sine_mf);
eq_sig = conv(half_sine_mf_out_005, q_zf);
eyediagram(eq_sig, n)
title("Half Sine EQ Sigma = 0.05")

sigma = 1;
noise = sigma * randn(size(half_sine_channel));
half_sine_mf_out_1 = conv(half_sine_channel + noise, half_sine_mf);
eq_sig = conv(half_sine_mf_out_1, q_zf);
eyediagram(eq_sig, n)
title("Half Sine EQ Sigma = 1")

sigma = 0;
noise = sigma * randn(size(srrc_channel));
srrc_mf_out_0 = conv(srrc_channel + noise, srrc_mf);
eq_sig = conv(srrc_mf_out_0, q_zf);
eyediagram(eq_sig, n)
title("SRRC EQ Sigma = 0")

sigma = 0.1;
noise = sigma * randn(size(srrc_channel));
srrc_mf_out_005 = conv(srrc_channel + noise, srrc_mf);
eq_sig = conv(srrc_mf_out_005, q_zf);
eyediagram(eq_sig, n)
title("SRRC EQ Sigma = 0.05")

sigma = 1;
noise = sigma * randn(size(srrc_channel));
srrc_mf_out_1 = conv(srrc_channel + noise, srrc_mf);
eq_sig = conv(srrc_mf_out_1, q_zf);
eyediagram(eq_sig, n)
title("SRRC EQ Sigma = 1")

% Question 12

[H, w] = freqz(h_upsample, 1, n, 'whole');

Eb = 1;
Q_mmse =  conj(H)./((abs(H).^2) + (sigma.^2) / Eb);

q_mmse = ifft(Q_mmse);

figure()
freqz(Q_mmse, 1)
title("Frequency Response of Q(f)_{MMSE}")

% Question 13 OFFSET IS WRONG MIGHT NEED TRUCATION FROM Q9 WE NEVER
% IMPLEMENTED
sigma = 0;
noise = sigma * randn(size(half_sine_channel));
half_sine_mf_out_0 = conv(half_sine_channel + noise, half_sine_mf);
eq_sig = conv(half_sine_mf_out_0, q_mmse);
eyediagram(eq_sig, n)
title("Half Sine EQ Sigma = 0")

sigma = 0.1;
noise = sigma * randn(size(half_sine_channel));
half_sine_mf_out_005 = conv(half_sine_channel + noise, half_sine_mf);
eq_sig = conv(half_sine_mf_out_005, q_mmse);
eyediagram(eq_sig, n)
title("Half Sine EQ Sigma = 0.05")

sigma = 1;
noise = sigma * randn(size(half_sine_channel));
half_sine_mf_out_1 = conv(half_sine_channel + noise, half_sine_mf);
eq_sig = conv(half_sine_mf_out_1, q_mmse);
eyediagram(eq_sig, n)
title("Half Sine EQ Sigma = 1")

sigma = 0;
noise = sigma * randn(size(srrc_channel));
srrc_mf_out_0 = conv(srrc_channel + noise, srrc_mf);
eq_sig = conv(srrc_mf_out_0, q_mmse);
eyediagram(eq_sig, n)
title("SRRC EQ Sigma = 0")

sigma = 0.1;
noise = sigma * randn(size(srrc_channel));
srrc_mf_out_005 = conv(srrc_channel + noise, srrc_mf);
eq_sig = conv(srrc_mf_out_005, q_mmse);
eyediagram(eq_sig, n)
title("SRRC EQ Sigma = 0.05")

sigma = 1;
noise = sigma * randn(size(srrc_channel));
srrc_mf_out_1 = conv(srrc_channel + noise, srrc_mf);
eq_sig = conv(srrc_mf_out_1, q_mmse);
eyediagram(eq_sig, n)
title("SRRC EQ Sigma = 1")

% Question 14
function mod_out = signal_mod(bits, n, pulse)
    symbols = 2*bits - 1;
    sym_up = zeros(1, length(bits) * n);
    sym_up(1:n:end) = symbols;
    
    mod_out = conv(sym_up, pulse);
end


function sent_image = convert_to_image(detected_bits, image_max, image_min, image_size)
    m = image_size(1);
    n = image_size(2);
    num_pixels = m * n;
    
    bit_matrix = reshape(detected_bits, 8, num_pixels);
    pixel_ints = bit2int(bit_matrix, 8);
    
    scaled_vals = double(pixel_ints) / 255;
    
    dct_vals = scaled_vals * (image_max - image_min) + image_min;
    
    num_blocks = num_pixels / 64;
    dct_blocks_3d = reshape(dct_vals, [8 8 num_blocks]);
    dct_image = reshape(dct_blocks_3d, [m n]);
    
    fun = @(block_struct) idct2(block_struct.data);
    sent_image = blockproc(dct_image, [8 8], fun);
end

function detected_bits = detect_bits(eq_out, k_start, n, num_bits)
    sample_indices = k_start : n : k_start + (num_bits - 1)*n;
    samples = eq_out(sample_indices);
    
    detected_bits = double(samples > 0);
end

function channel_out = transmit_channel(h, n, signal_in, sigma)
    h_upsample = upsample(h, n);

    channel = conv(h_upsample, signal_in);

    noise = sigma * randn(size(channel));

    channel_out = channel + noise;
end

function mf_out = matched_filter(channel_out, pulse_shape)
    mf = flip(pulse_shape);

    mf_out = conv(channel_out, mf);
end

function qz_out = qz(mf_out, h)

    q_zf = filter(1, h, [1, zeros(1, 10000)]);
    
    qz_out = conv(mf_out, q_zf);

end

function mmse_out = mmse(mf_out, n)

    [H, w] = freqz(h_upsample, 1, n, 'whole');
    
    Eb = 1;
    Q_mmse =  conj(H)./((abs(H).^2) + (sigma.^2) / Eb);
    
    q_mmse = ifft(Q_mmse);

    mmse_out = conv(mf_out, q_mmse);


end

function [arr_3d, image_max, image_min, image_size] = preProcess(imageName)
    image = imread(imageName);
    image = imresize(image, 0.1);
    image_arr = im2double(im2gray(image));

    [m, n] = size(image_arr);
    m = floor(m/8)*8;
    n = floor(n/8)*8;
    image_arr = image_arr(1:m, 1:n);
    
    
    fun = @(block_struct) dct2(block_struct.data);
    image_dct = blockproc(image_arr, [8 8], fun);
    
    image_max = max(image_dct, [], 'all');
    image_min = min(image_dct, [], 'all');
    
    scaled_arr = (image_dct - image_min) / (image_max - image_min);

    image_size = size(scaled_arr);
    arr_3d = reshape(scaled_arr, [8 8 (m*n)/64]);
end

function bit_stream = bitStreamConvert(arr_3d, N)
    col_group = reshape(arr_3d(:, :, 1:N), [64*N 1]);
    col_int = round(col_group * 255);
    bit_stream = reshape(int2bit(col_int, 8), [1, 64*8*N]);
end
