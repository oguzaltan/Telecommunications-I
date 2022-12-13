clear all; close all;

%random bit generate
num_bits = 1e4;

rand_bit = round(rand(1,num_bits));
signal = [];

for i = 1:num_bits
    if (rand_bit(i) == 1)
        signal = [signal triangleSignal(20)];
    else
        signal = [signal (-1)*triangleSignal(20)];
    end
end

figure;
plot(0:length(signal)-1,signal);
title("Transmitted Signal");
grid on;

% Generate Noise
% Choose an SNR value of 5 dB
SNR_db = 5;
SNR = 10^(SNR_db/10);
No = 1/SNR;
n = sqrt(No/2).* randn(1,length(signal));

r = signal + n;

figure;
plot(0:length(r)-1,r);
title("Received Signal r");
grid on;

filter = triangleSignal(20);

figure;
plot(0:length(filter)-1,filter);
title("Matched Filter");
grid on;

%Part 1
disp("Part 1");

error_vec1  = zeros(1,100);
for i = 1:100
    n = sqrt(No/2).* randn(1,length(signal));
    r = signal + n;
    conv_signal = conv(r,filter);
    
    for j = 1:length(rand_bit)
        if (conv_signal(20*j)>0)
            check_bits(j) = 1;
        else
            check_bits(j) = 0;
        end
    end
    
    error = sum(abs(check_bits-rand_bit))/length(rand_bit)*100;
    error_vec(i) = error;
end

error_av = mean(error_vec);
disp(['Average Error Rate for Sampling Rate 20: ', num2str(error_av), '%']);

%Part 2
disp("Part 2");

signal2 = [];

for i = 1:num_bits
    if (rand_bit(i) == 1)
        signal2 = [signal2 triangleSignal(50)];
    else
        signal2 = [signal2 (-1)*triangleSignal(50)];
    end
end

filter2 = triangleSignal(50);

error_vec2 = zeros(1,100);
for i = 1:100
    n2 = sqrt(No/2).* randn(1,length(signal2));
    r2 = signal2 + n2;
    conv_signal2 = conv(r2,filter2);
    
    for j = 1:length(rand_bit)
        if (conv_signal2(50*j)>0)
            check_bits2(j) = 1;
        else
            check_bits2(j) = 0;
        end
    end
    
    error2 = sum(abs(check_bits2-rand_bit))/length(rand_bit)*100;
    error_vec2(i) = error2;
end

error_av2 = mean(error_vec2);
disp(['Average Error Rate for Sampling Rate 50: ', num2str(error_av2), '%']);

%Part 3
disp("Part 3");

t3 = 0:1/20:1-1/20;
filter3 = ones(1,length(t3));
% plot(t3-1,filter3);

error_vec3  = zeros(1,100);

for i = 1:100
    n = sqrt(No/2).* randn(1,length(signal));
    r = signal + n;
    conv_signal3 = conv(r,filter3);
    
    for j = 1:length(rand_bit)
        if (conv_signal3(20*j)>0)
            check_bits3(j) = 1;
        else
            check_bits3(j) = 0;
        end
    end
    
    error3 = sum(abs(check_bits3-rand_bit))/length(rand_bit)*100;
    error_vec3(i) = error3;
end

error_av3 = mean(error_vec3);
disp(['Average Error Rate for Not Ideal Filter: ', num2str(error_av3), '%']);

%Part 4
disp("Part 4");

deltaT_vec = [-10 -5 -3 3 5 10];
error_av = zeros(1,6);

for ts = 1:length(deltaT_vec)
    
    for i = 1:100
        n = sqrt(No/2).* randn(1,length(signal));
        r = signal + n;
        conv_signal = conv(r,filter);
        
        for j = 1:length(rand_bit)
            if (conv_signal(20*j + deltaT_vec(ts)) > 0)
                check_bits(j) = 1;
            else
                check_bits(j) = 0;
            end
        end
        
        error = sum(abs(check_bits-rand_bit))/length(rand_bit)*100;
        error_vec(i) = error;
    end
    
    error_av(ts) = mean(error_vec);
end

disp(['Error for T-10deltaT ', num2str(error_av(1)), '%'])
disp(['Error for T-5deltaT ', num2str(error_av(2)), '%'])
disp(['Error for T-3deltaT ', num2str(error_av(3)), '%'])
disp(['Error for T+3deltaT ', num2str(error_av(4)), '%'])
disp(['Error for T+5deltaT ', num2str(error_av(5)), '%'])
disp(['Error for T+10deltaT ', num2str(error_av(6)), '%'])

%Part 5
disp("Part 5");

T = 20;
deltaT_vec = [-10 -5 -3 3 5 10];
error_q5 = zeros(1,6);

for i = 1:length(deltaT_vec)
    check_bits_q5 = zeros(1,num_bits);
    
    %for previous message
    if (deltaT_vec(i) < 0)
        first = T + 1;
        last = 2*T;
        prev_and_current = [zeros(1,T) r];
        
        for j = 1:num_bits
            conv_prev = conv(prev_and_current(first:last),filter);
            prev_signal = prev_and_current(first-T:last-T);
            temp = conv(prev_signal,filter);
            temp = temp(T:length(temp));
            conv_prev(1:T) = conv_prev(1:T) + temp;
            
            if conv_prev(deltaT_vec(i) + T) > 0
                check_bits_q5(j) = 1;
            else
                check_bits_q5(j) = 0;
            end
            first = first + T;
            last = last + T;
        end
        
        %for next message
    else
        first = 1;
        last = T;
        current_and_next = [r zeros(1,T)];
        
        for j = 1:num_bits
            conv_next = conv(current_and_next(first:last),filter);
            next_signal = current_and_next(first+T:last+T);
            temp = conv(next_signal,filter);
            temp = temp(1:T);
            conv_next(T:length(conv_next)) = conv_next(T:length(conv_next)) + temp;
            
            if conv_next(deltaT_vec(i) + T) > 0
                check_bits_q5(j) = 1;
            else
                check_bits_q5(j) = 0;
            end
            first = first + T;
            last = last + T;
        end
    end
    
    error_q5(i) = sum(abs(check_bits_q5-rand_bit))/length(rand_bit)*100; 
end

disp(['Error for T-10deltaT ', num2str(error_q5(1)), '%'])
disp(['Error for T-5deltaT ', num2str(error_q5(2)), '%'])
disp(['Error for T-3deltaT ', num2str(error_q5(3)), '%'])
disp(['Error for T+3deltaT ', num2str(error_q5(4)), '%'])
disp(['Error for T+5deltaT ', num2str(error_q5(5)), '%'])
disp(['Error for T+10deltaT ', num2str(error_q5(6)), '%'])

%Part 6
disp("Part 6");

signal_fsk = [];
check_bits_fsk = zeros(1,length(rand_bit));

Ts = 20;
f1 = 1;
f2 = 4;
t_fsk = 0:1/Ts:1-1/Ts;

for i = 1:num_bits
    if rand_bit(i) == 1
        signal_fsk = [signal_fsk sqrt(2/Ts)*cos(2*pi*f1*t_fsk)];
    else
        signal_fsk = [signal_fsk sqrt(2/Ts)*cos(2*pi*f2*t_fsk)];
    end
end

figure;
one_signal_fsk = [signal_fsk(1:Ts) 0.316227766016838];
plot(0:length(one_signal_fsk)-1,one_signal_fsk);
grid on;
title("Sinusoidal Signal");
xlabel("Time");
ylabel("Amplitude");

n_fsk = sqrt(No/2).* randn(1,length(signal_fsk));
r_fsk = signal_fsk + n_fsk;

figure;
one_signal_fsk_r = r_fsk(1:Ts);
plot(1:length(one_signal_fsk_r),one_signal_fsk_r);
grid on;
title("Received Signal");
xlabel("Time");
ylabel("Amplitude");

basis1_fsk = sqrt(2/Ts)*cos(2*pi*f1*t_fsk);
basis2_fsk = sqrt(2/Ts)*cos(2*pi*f2*t_fsk);

for i = 1:length(rand_bit)
    
    c1 = sum((r_fsk((i-1)*Ts+1:i*Ts).*basis1_fsk));
    c2 = sum((r_fsk((i-1)*Ts+1:i*Ts).*basis2_fsk));
    
    if c1 >= c2
        check_bits_fsk(i) = 1;
    else
        check_bits_fsk(i) = 0;
    end
end

error_fsk = 100*sum(abs(check_bits_fsk-rand_bit))/length(rand_bit);
disp(['Estimated Error Rate for FSK: ', num2str(error_fsk), '%']);

%Part 7
disp("Part 7");

basis11_nc_fsk = sqrt(2/Ts)*cos(2*pi*f1*t_fsk);
basis12_nc_fsk = sqrt(2/Ts)*sin(2*pi*f1*t_fsk);

basis21_nc_fsk = sqrt(2/Ts)*cos(2*pi*f2*t_fsk);
basis22_nc_fsk = sqrt(2/Ts)*sin(2*pi*f2*t_fsk);

check_bits_nc_fsk = zeros(1,length(rand_bit));

for i = 1:length(rand_bit)
    
    c11 = sum(r_fsk((i-1)*Ts+1:i*Ts).*basis11_nc_fsk);
    c12 = sum(r_fsk((i-1)*Ts+1:i*Ts).*basis12_nc_fsk);
    c21 = sum(r_fsk((i-1)*Ts+1:i*Ts).*basis21_nc_fsk);
    c22 = sum(r_fsk((i-1)*Ts+1:i*Ts).*basis22_nc_fsk);
    
    if (c11^2 + c12^2) > (c21^2 + c22^2)
        check_bits_nc_fsk(i) = 1;
    else
        check_bits_nc_fsk(i) = 0;
    end
end

error_nc_fsk = 100*sum(abs(check_bits_nc_fsk-rand_bit))/length(rand_bit);
disp(['Estimated Error Rate for Non-Coherent FSK: ', num2str(error_nc_fsk), '%']);

%generate triangular signal
function signal = triangleSignal(sample_rate)

f = 1; %fundamental freq
T = 1/f;
Ts = T/sample_rate; %sampling period
fs = 1/Ts; %sample rate

t = 0:1/fs:T - 1/fs; %calculate t
signal = (1+sawtooth(2*pi*f*t,1/2))/2; %constructing triangle using sawtooth

%normalization
E = sum(signal.*signal); %energy of the signal
signal = signal/sqrt(E);

end