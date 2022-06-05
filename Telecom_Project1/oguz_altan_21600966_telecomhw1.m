clear all
clc
close all

%preprocessing of the audio data
nBits = 8;
nchannels = 1;
nID = -1;
Fs = 200000;
Fc = 4000;
time = 5;

% recObj = audiorecorder(Fs,nBits,nchannels,nID);
% recObj.StartFcn = 'disp(''Start speaking.'')';
% recObj.StopFcn = 'disp(''End of recording.'')';
% recordblocking(recObj,time);
% x = getaudiodata(recObj);
% y = lowpass(x,Fc,Fs);

load ('x.mat');
x1 = x;
x1_hist = x1;

figure;
plot(1:length(x1),x1);
title('Speech Audio Data');
xlabel('Data Points');
ylabel('Audio Signal Amplitude');
grid on;

x1 = nonzeros(x1);
figure;
plot(1:length(x1),x1)
title('Optimized Speech Audio Data');
xlabel('Data Points');
ylabel('Audio Signal Amplitude');
grid on;

y1 = lowpass(x1,Fc,Fs);
y1 = normalize(x1,'range',[-1,1]);
figure;
h1 = histogram(x1,'normalization','probability');
title('Normalized Histogram to Estimate PDF');
grid on;

pdf_x1 = h1.Values;
x1 = linspace(h1.BinLimits(1),h1.BinLimits(2),h1.NumBins);
pol1 = polyfit(x1,pdf_x1,10);
polx1 = linspace(h1.BinLimits(1),h1.BinLimits(2));

figure;
plot(polx1,polyval(pol1,polx1));
title('Fitted and Normalized Polynomial');
xlabel('x1');
ylabel('f_x(x1)');
grid on

for i = 1:11
    if(pol1(i) < 0)
        pol1(i) = 0;
    end
end

%Part I Question 2
N = [16 64 128]; %number of quantization levels
sqnr1_uni_vec = zeros(1,3);
sqnr1_nonuni_vec = zeros(1,3);
distortion1_uni_vec = zeros(1,3);
distortion1_nonuni_vec = zeros(1,3);
sqnr1_uni_vec_abs = zeros(1,3);
sqnr1_nonuni_abs_vec = zeros(1,3);
distortion1_nonuni_abs_vec = zeros(1,3);
distortion1_uni_abs_vec = zeros(1,3);

for nn = 1:3
    
    %UNIFORM QUANTIZATION
    a = zeros(1,N(nn)+1);
    a(1) = -1;
    a(N(nn)+1) = 1;
    
    log_vec = zeros(1,length(a));
    
    uniqua_error_sum1 = 0;
    uniqua_error_sum1_abs = 0;
    MSE_y1 = 0;
    
    for i = 1:N(nn)+1
        a(i) = a(1) + (a(N(nn)+1)-a(1))/N(nn)*(i-1);
    end
    
    level_dif1 = a(2)-a(1);
    
    for j = 1:N(nn)
        x_quan1(j) = (a(j) + a(j+1))/2;
    end
    
    for ally = 1:length(x1)
        log_vec = y1(ally)- a < level_dif1;
        MSE_y1 = MSE_y1 + y1(ally).^2;
        uniqua_error_sum1 = uniqua_error_sum1 + (y1(ally) - x_quan1(first_one(log_vec))).^2;
        uniqua_error_sum1_abs = uniqua_error_sum1_abs + abs(y1(ally) - x_quan1(first_one(log_vec)));
    end
    
    %error calculation for uniform quantizer
    D1 = uniqua_error_sum1/length(x1);
    D1_abs = uniqua_error_sum1_abs/length(x1);
    P1 = MSE_y1/length(x1);
    distortion1_uni_vec(nn) = D1;
    sqnr1_uni = P1/D1;
    sqnr1_uni_abs = P1/D1_abs;
    distortion1_uni_abs_vec(nn) = D1_abs;
    
    sqnr1_uni_vec(nn) = sqnr1_uni;
    sqnr1_uni_vec_abs(nn) = sqnr1_uni_abs;
    
    %NONUNIFORM QUANTIZATION
    
    a1_nonqua = a; %initializing with uniform quantizer ai's
    x1_nonqua = zeros(N(nn),1);
    y1_nonqua = zeros(length(y1),1);
    
    %integral of x times probability density function of x
    numer_pdf1= @(x_poly1) pol1(1)*(x_poly1.^10)+pol1(2)*(x_poly1.^9)+pol1(3)*(x_poly1.^8)+pol1(4)*(x_poly1.^7)+pol1(5)*(x_poly1.^6)+pol1(6)*(x_poly1.^5)+pol1(7)*(x_poly1.^4)+pol1(8)*(x_poly1.^3)+pol1(9)*(x_poly1.^2)+pol1(10)*(x_poly1)+pol1(11);
    
    %integral of probability density function of x
    denom_x_pdf1 = @(x_poly1) (x_poly1.*(pol1(1)*(x_poly1.^10)+pol1(2)*(x_poly1.^9)+pol1(3)*(x_poly1.^8)+pol1(4)*(x_poly1.^7)+pol1(5)*(x_poly1.^6)+pol1(6)*(x_poly1.^5)+pol1(7)*(x_poly1.^4)+pol1(8)*(x_poly1.^3)+pol1(9)*(x_poly1.^2)+pol1(10)*(x_poly1)+pol1(11)));
    
    D1_nonqua = 0.8;
    
    %calculates xi's using the ratio of integrals over two consecutive ai levels
    while D1_nonqua > 0.002 %loops until distortion changes more than 0.001 threshold value. Stops if distortion does not change more than 0.001
        
        %calculates xi's using the ratio of x*pdf integral and pdf integral
        for i= 1:N(nn)
            x1_nonqua(i) = integral(denom_x_pdf1,a1_nonqua(i),a1_nonqua(i+1))/integral(numer_pdf1,a1_nonqua(i),a1_nonqua(i+1));
        end
        
        %ai's are midpoints of two consecutive xi's
        for j = 2:N(nn)
            a1_nonqua(j)=(x1_nonqua(j)+x1_nonqua(j-1))/2;
        end
        
        for j = 1:length(y1)
            for lev = 1:N(nn)
                if y1(j) > a1_nonqua(lev) && y1(j) < a1_nonqua(lev+1)
                    y1_nonqua(j) = x1_nonqua(lev);
                end
            end
        end
        
        D1_nonqua = mean((y1-y1_nonqua).^2);
        P1_nonqua = mean((y1).^2);
        sqnr1_nonuni_vec(nn) = P1_nonqua/D1_nonqua;
        distortion1_nonuni_vec(nn) = D1_nonqua;
    end
    
    %NONUNIFORM QUANTIZATION WITH ABSOLUTE VALUE
    
    a1_init_abs = a; %initializing with uniform quantizer ai's
    x1_i_hat = zeros([1,N(nn)]);
    y1_nonqua_abs = zeros(length(y1),1);
    
    abs1_integral = @(x1_i_hat, amax, amin) integral(@(x_poly1) abs(x_poly1-x1_i_hat).*(pol1(1)*(x_poly1.^10)+pol1(2)*(x_poly1.^9)+pol1(3)*(x_poly1.^8)+pol1(4)*(x_poly1.^7)+pol1(5)*(x_poly1.^6)+pol1(6)*(x_poly1.^5)+pol1(7)*(x_poly1.^4)+pol1(8)*(x_poly1.^3)+pol1(9)*(x_poly1.^2)+pol1(10)*(x_poly1)+pol1(11)), amax, amin);
    
    quantizer_level = 0;
    
    %divides the quantization level intervals to small points and
    %tries to find the optimal point that gives smallest integral
    %(quantization error)
    
    while quantizer_level < N(nn)
        for j = 1:length(x1_i_hat)
            
            sub_interv = a1_init_abs(j):0.01:a1_init_abs(j+1);
            integ = zeros([1,length(sub_interv)]);
            
            for k = 1:length(sub_interv)
                integ(k)= abs1_integral(sub_interv(k), a1_init_abs(j), a1_init_abs(j+1));
            end
            
            %finds the index of the optimal quantization level that
            %minimizes the integral and locates that level
            [Y,I] = min(integ);
            x1_i_hat(j) = sub_interv(I);
        end
        
        %ai's are middle point of optimized quantization level xi's
        for i = 2:N(nn)
            a1_init_abs(i) = (x1_i_hat(i-1) + x1_i_hat(i))*0.5;
        end
        quantizer_level = quantizer_level + 1;
    end
    
    for j = 1:length(y1)
        for lev = 1:N(nn)
            if y1(j) > a1_init_abs(lev) && y1(j) < a1_init_abs(lev+1)
                y1_nonqua_abs(j) = x1_i_hat(lev);
            end
        end
    end
    
    D1_nonqua_abs = mean(abs(y1 - y1_nonqua_abs));
    P1_nonqua_abs = mean((y1).^2);
    sqnr1_nonuni_abs_vec(nn) = P1_nonqua_abs/D1_nonqua_abs;
    distortion1_nonuni_abs_vec(nn) = D1_nonqua_abs;
end

%second speech signal

load voice.mat
x = [x;x;x(1:199998)];
x2 = x;
x2_hist = x2;

nBits = 8;
nchannels = 1;
nID = -1;
Fs = 200000;
Fc = 4000;
time = 5;

% recObj = audiorecorder(Fs,nBits,nchannels,nID);
% recObj.StartFcn = 'disp(''Start speaking.'')';
% recObj.StopFcn = 'disp(''End of recording.'')';
% recordblocking(recObj,time);
% x2 = getaudiodata(recObj);
% y2 = lowpass(x2,Fc,Fs);

figure;
plot(1:length(x2),x2);
title('Speech Audio Data of Second Speech');
xlabel('Data Points');
ylabel('Audio Signal Amplitude');
grid on;
x2 = nonzeros(x2);

figure;
plot(1:length(x2),x2)
title('Optimized Speech Audio Data of Second Speech');
xlabel('Data Points');
ylabel('Audio Signal Amplitude');
grid on;

y2 = lowpass(x2,Fc,Fs);
y2 = normalize(y2,'range',[-1,1]);
figure;
h2 = histogram(x2,'normalization','probability');
title('Normalized Histogram to Estimate PDF of Second Speech');
grid on;

pdf_x2 = h2.Values;
x2 = linspace(h2.BinLimits(1),h2.BinLimits(2),h2.NumBins);
pol2 = polyfit(x2,pdf_x2,10);
polx2 = linspace(h2.BinLimits(1),h2.BinLimits(2));

figure;
plot(polx2,polyval(pol2,polx2));
title('Fitted and Normalized Polynomial of Second Speech');
xlabel('x2');
ylabel('f_x(x2)');
grid on

for i = 1:11
    if(pol2(i) < 0)
        pol2(i) = 0;
    end
end

%Part I Question 2
N = [16 64 128]; %number of quantization levels
sqnr2_uni_vec = zeros(1,3);
sqnr2_nonuni_vec = zeros(1,3);
sqnr2_nonuni_abs_vec = zeros(1,3);
distortion2_uni_vec = zeros(1,3);
distortion2_nonuni_vec = zeros(1,3);
distortion2_nonuni_abs_vec = zeros(1,3);
mse2_nonuni_vec = zeros(1,3);

for nn = 1:3
    
    %UNIFORM QUANTIZATION
    a = zeros(1,N(nn)+1);
    a(1) = -1;
    a(N(nn)+1) = 1;
    
    log_vec = zeros(1,length(a));
    
    uniqua_error_sum2 = 0;
    uniqua_error_sum2_abs = 0;
    nonuniqua_error_sum2 = 0;
    nonuniqua_error_sum2_abs = 0;
    MSE_y2 = 0;
    MSE_y2_nonuni = 0;
    
    for i = 1:N(nn)+1
        a(i) = a(1) + (a(N(nn)+1)-a(1))/N(nn)*(i-1);
    end
    
    level_dif2 = a(2)-a(1);
    
    for j = 1:N(nn)
        x_quan2(j) = (a(j) + a(j+1))/2;
    end
    
    for ally = 1:length(x2)
        log_vec = y2(ally)- a < level_dif2;
        MSE_y2 = MSE_y2 + y2(ally).^2;
        uniqua_error_sum2 = uniqua_error_sum2 + (y2(ally) - x_quan2(first_one(log_vec))).^2;
        uniqua_error_sum2_abs = uniqua_error_sum2_abs + abs(y2(ally) - x_quan2(first_one(log_vec)));
    end
    
    %error calculation for uniform quantizer
    D2 = uniqua_error_sum2/length(x2);
    D2_abs = uniqua_error_sum2_abs/length(x2);
    P2 = MSE_y2/length(x2);
    distortion2_uni_vec(nn) = D2;
    sqnr2_uni = P2/D2;
    sqnr2_uni_abs = P2/D2_abs;
    
    sqnr2_uni_vec(nn) = sqnr2_uni;
    sqnr2_uni_vec_abs(nn) = sqnr2_uni_abs;
    
    %NONUNIFORM QUANTIZATION
    
    a2_nonqua = a; %initializing with uniform quantizer ai's
    x2_nonqua = zeros(N(nn),1);
    y2_nonqua = zeros(length(y2),1);
    
    D2_nonqua = 0.8; %initial distortion values
    
    %integral of x times probability density function of x
    denom_pdf2= @(x_poly2) pol2(1)*(x_poly2.^10)+pol2(2)*(x_poly2.^9)+pol2(3)*(x_poly2.^8)+pol2(4)*(x_poly2.^7)+pol2(5)*(x_poly2.^6)+pol2(6)*(x_poly2.^5)+pol2(7)*(x_poly2.^4)+pol2(8)*(x_poly2.^3)+pol2(9)*(x_poly2.^2)+pol2(10)*(x_poly2)+pol2(11);
    
    %integral of probability density function of x
    numer_x_pdf2 = @(x_poly2) (x_poly2.*(pol2(1)*(x_poly2.^10)+pol2(2)*(x_poly2.^9)+pol2(3)*(x_poly2.^8)+pol2(4)*(x_poly2.^7)+pol2(5)*(x_poly2.^6)+pol2(6)*(x_poly2.^5)+pol2(7)*(x_poly2.^4)+pol2(8)*(x_poly2.^3)+pol2(9)*(x_poly2.^2)+pol2(10)*(x_poly2)+pol2(11)));
    
    %calculates xi's using the ratio of integrals over two consecutive ai levels
    while D2_nonqua > 0.002 %loops until distortion changes more than 0.002 threshold value. Stops if distortion does not change more than 0.001
        
        %calculates xi's using the ratio of x*pdf integral and pdf integral
        %condiional expectation
        for i= 1:N(nn)
            x2_nonqua(i) = integral(numer_x_pdf2,a2_nonqua(i),a2_nonqua(i+1))/integral(denom_pdf2,a2_nonqua(i),a2_nonqua(i+1));
        end
        
        %ai's are midpoints of two consecutive xi's
        for j = 2:N(nn)
            a2_nonqua(j)=(x2_nonqua(j)+x2_nonqua(j-1))/2;
        end
        
        %locates the quantization level
        for j = 1:length(y2)
            for lev = 1:N(nn)
                if y2(j) > a2_nonqua(lev) && y2(j) < a2_nonqua(lev+1)
                    y2_nonqua(j) = x2_nonqua(lev);
                end
            end
        end
        
        D2_nonqua = mean((y2-y2_nonqua).^2);
        P2_nonqua = mean((y2).^2);
        sqnr2_nonuni_vec(nn) = P2_nonqua/D2_nonqua;
        distortion2_nonuni_vec(nn) = D2_nonqua;
        
        %NONUNIFORM QUANTIZATION WITH ABSOLUTE VALUE
        
        a2_init_abs = a; %initializing with uniform quantizer ai's
        x2_i_hat = zeros([1,N(nn)]);
        y2_nonqua_abs = zeros(length(y2),1);
        
        %the integral for the mean absolute error
        abs2_integral = @(x2_i_hat, amax, amin) integral(@(x_poly2) abs(x_poly2-x2_i_hat).*(pol2(1)*(x_poly2.^10)+pol2(2)*(x_poly2.^9)+pol2(3)*(x_poly2.^8)+pol2(4)*(x_poly2.^7)+pol2(5)*(x_poly2.^6)+pol2(6)*(x_poly2.^5)+pol2(7)*(x_poly2.^4)+pol2(8)*(x_poly2.^3)+pol2(9)*(x_poly2.^2)+pol2(10)*(x_poly2)+pol2(11)), amax, amin);
        
        quantizer_level = 0;
        
        %divides the quantization level intervals to small points and
        %tries to find the optimal point that gives smallest integral
        %(quantization error)
        
        while quantizer_level < N(nn)
            for j = 1:length(x2_i_hat)
                
                %divides the intervals to small points
                sub_interv = a2_init_abs(j):0.01:a2_init_abs(j+1);
                
                %calculates the integral (quantization error)
                integ = zeros([1,length(sub_interv)]);
                for k = 1:length(sub_interv)
                    integ(k)= abs2_integral(sub_interv(k), a2_init_abs(j), a2_init_abs(j+1));
                end
                
                %finds the index of the optimal quantization level that
                %minimizes the integral and locates that level
                [Y,I] = min(integ);
                x2_i_hat(j) = sub_interv(I);
            end
            
            %ai's are middle point of optimized quantization level xi's
            for lev = 2:N(nn)
                a2_init_abs(lev) = (x2_i_hat(lev-1) + x2_i_hat(lev))*0.5;
            end
            quantizer_level = quantizer_level + 1;
        end
        
        for j = 1:length(y2)
            for a = 1:N(nn)
                if y2(j) > a2_init_abs(a) && y2(j) < a2_init_abs(a+1)
                    y2_nonqua_abs(j) = x2_i_hat(a);
                end
            end
        end
        
        D2_nonqua_abs = mean(abs(y2 - y2_nonqua_abs));
        P2_nonqua_abs = mean((y2).^2);
        sqnr2_nonuni_abs_vec(nn) = P2_nonqua_abs/D2_nonqua_abs;
        distortion2_nonuni_abs_vec(nn) = D2_nonqua_abs;
        
        %plots the 2d joinlty pdf of the consecutive source samples
        h3 = histogram2(x1_hist,x2_hist,"Normalization","Probability")
        
    end
end

function one_index = first_one(test_vec)
for i = 1:length(test_vec)
    if test_vec(i) == 1;
        one_index = i;
        break
    end
end
end