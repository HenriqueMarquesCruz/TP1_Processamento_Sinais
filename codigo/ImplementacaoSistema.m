clear; close all; clc;


% 1. Carregamento de dados e apresentação de características básicas 

%% ------- Parâmetros -------
base_dir = fullfile('..','material_fornecido');

audio_file = fullfile(base_dir, 'audio_corrompido.wav');
num_file   = fullfile(base_dir, 'coefs_num.mat');
den_file   = fullfile(base_dir, 'coefs_den.mat');
nfft_spec  = 16384;    % pontos para FFT. Poderia ser outra potência de 2, mas 16384 é um compromisso:
%suficientemente grande para ter boa resolução em Hz, mas sem deixar a FFT lenta
n_impulse  = 1000;     % número de amostras para resposta ao impulso

%% ------- 1.1 Carregamento do áudio e reprodução -------
if ~exist(audio_file, 'file')
    error('Arquivo de audio não encontrado: %s', audio_file);
end

[x, fs] = audioread(audio_file);      % lê áudio
if size(x,2) > 1
    x = mean(x,2);                    % converte para mono (média dos canais)
end
N = length(x);
t = (0:N-1)/fs;

%fprintf('Amostras: %d, fs = %d Hz, duração = %.3f s\n', N, fs, N/fs);
%Extra apenas

%% Perguntar ao usuário se deseja ouvir
fprintf('\n=== Seção 1: Carregamento de dados e apresentação de características básicas  ===\n');

choice = input('Deseja ouvir o áudio? (s/n): ', 's');

if lower(choice) == 's'
    try
        fprintf('Reproduzindo áudio... (aguarde)\n\n');
        sound(x, fs);
        pause(N/fs); % toca até o fim
    catch ME
        warning(ME.identifier, 'Não foi possível reproduzir áudio: %s', ME.message);
    end
else
    fprintf('Ok, o áudio não será reproduzido.\n\n');
end

%% ------- 1.2 Gráfico da forma de onda (tempo em s) -------
figure('Name','1. Carregamento de dados e apresentação de características básicas',...
       'NumberTitle','off','Position',[125 100 1050 600]);
subplot(5,2,[1 2]); % pega colunas 1 e 2 da primeira linha
plot(t, x);
xlabel('Tempo (s)');
ylabel('Amplitude normalizada');
title('Forma de onda do áudio corrompido');
grid on;
print(gcf,'fig_forma_de_onda.png','-dpng','-r150');

%% ------- 1.3 Espectros de amplitude |X(e^{j\omega})|×f (kHz) e fase θ(ω)×f (kHz) -------
Nfft = max(nfft_spec, 2^nextpow2(N));
X = fft(x, Nfft);
Xs = fftshift(X);

freqs = (-Nfft/2 : Nfft/2-1) * (fs / Nfft);   % Hz
freqs_khz = freqs / 1000;                     % kHz

amp = abs(Xs);
phase = angle(Xs);   % fase embrulhada (–π..π)

% --- Plot amplitude ---
subplot(5,2,3);
plot(freqs_khz, amp, 'b');
xlabel('Frequência (kHz)');
ylabel('|X(e^{j\omega})|');
title('Espectro de amplitude do áudio corrompido');
grid on;
xlim([min(freqs_khz) max(freqs_khz)]);

% --- Plot fase ---
subplot(5,2,4);
plot(freqs_khz, phase, 'b');
xlabel('Frequência (kHz)');
ylabel('Fase (rad)');
title('Espectro de fase do áudio corrompido');
grid on;
xlim([min(freqs_khz) max(freqs_khz)]);


%% ------- 1.4 Carregamento dos coeficientes do filtro -------
if ~exist(num_file,'file') || ~exist(den_file,'file')
    error('Arquivos de coeficientes não encontrados. Esperados: %s e %s', num_file, den_file);
end
s_num = load(num_file);    % coeficientes numerador
s_den = load(den_file);    % coeficientes denominador 

num = s_num.num; % conteúdo do campo chamado num dentro da struct s_num.
den = s_den.den;

%% EXTRA: Impressão da função de transferência H(z)
%fprintf('\n=== Função de transferência H(z) ===\n');
%Hz = tf(num, den, -1)   % -1 força o domínio z^-1

% Avaliar no ponto z = 1
%z = 1;
%H1 = polyval(num, z) / polyval(den, z);

%fprintf('H(z=1) = %.6e\n', H1);
% ver tamanho/estrutura
%whos num den

%% ------- 1.5 Respostas de magnitude e fase |H(e^{jω})| e θ(ω) -------

Nfft = 16384;

% ----------- a) Escala linear simétrica (-fs/2..+fs/2) ----------
[H_whole, w_whole] = freqz(num, den, Nfft, 'whole', fs); 
Hshift = fftshift(H_whole);
wshift = (-fs/2 : fs/Nfft : fs/2 - fs/Nfft);
wshift_khz = wshift/1000;

% Magnitude linear
subplot(5,2,5);
plot(wshift_khz, abs(Hshift));
xlabel('f (kHz)');
ylabel('|H(e^{j\omega})|');
title('Magnitude (linear)');
grid on;

% Fase linear
subplot(5,2,7);
plot(wshift_khz, angle(Hshift));
xlabel('f (kHz)');
ylabel('\theta(\omega) (rad)');
title('Fase');
grid on;

% ----------- b) Escala log (dB) 0..fs/2 ----------
[H, w] = freqz(num, den, Nfft, fs);
w_khz = w/1000;

% Magnitude em dB
subplot(5,2,6);
plot(w_khz, 20*log10(abs(H)));
xlabel('Frequência (kHz)');
ylabel('Magnitude (dB)');
title('Magnitude (dB)');
grid on;
ylim([-120 5]); % para destacar

% Fase em graus (unwrap)
subplot(5,2,8);
plot(w_khz, unwrap(angle(H))*180/pi);
xlabel('Frequência (kHz)');
ylabel('Fase (graus)');
title('Fase');
grid on;


%% ------- 1.6 Resposta ao impulso do filtro -------
n_samples = 1000;   % quantidade de amostras da resposta ao impulso
[h, n] = impz(num, den, n_samples);

% Plotar a resposta ao impulso
subplot(5,2,[9 10]);
stem(n, h, 'filled');
xlabel('n');
ylabel('h[n]');
title('Resposta ao impulso do filtro');
grid on;


%% ===============================================================
% 2. Implementação de funções para filtragem
%% ===============================================================

%% ------- 2.1 Filtragem pela equação de diferenças -------
function y = filtragemPorEqDif(x, num, den)
    % Normalização para garantir que den(1) = 1
    if den(1) ~= 1
        num = num / den(1);
        den = den / den(1);
    end
    
    Nx = length(x);
    M = length(num);
    N = length(den);
    y = zeros(size(x));
    
    % Laço amostra a amostra
    for n = 1:Nx
        acc = 0;
        
        % Parte do numerador (entradas x[n-k])
        for k = 1:M
            if (n-k+1) > 0
                acc = acc + num(k) * x(n-k+1);
            end
        end
        
        % Parte do denominador (saídas passadas y[n-k])
        for k = 2:N
            if (n-k+1) > 0
                acc = acc - den(k) * y(n-k+1);
            end
        end
        
        y(n) = acc;
    end
end


%% ------- 2.2 Filtragem pela convolução com resposta ao impulso -------
function [y, h_trunc] = filtragemPorConv(x, h)
    %% (a) Truncagem da resposta ao impulso
    % Critério: manter todas as amostras até o último índice cujo valor
    % seja maior ou igual a 1% do valor de pico da resposta.
    limiar = 0.01 * max(abs(h));
    idx = find(abs(h) >= limiar, 1, 'last');
    h_trunc = h(1:idx);

    %% (b) Apresentação da resposta truncada
    % Número de amostras após truncagem
    Nh = length(h_trunc);
    %fprintf('Nh: %d', Nh); Apenas para conferir: Nh = 277

    %% (c) Filtragem por convolução circular (aqui usamos conv -> linear)
    % A convolução linear equivale à convolução circular com zero-padding
    % de tamanho adequado (Nx+Nh-1).
    y = conv(x, h_trunc);
end

%% ------- 2.3 Filtragem pela multiplicação da FFT -------
function y = filtragemPorFFT(x, h_trunc)
    % Comprimentos
    Nx = length(x);
    Nh = length(h_trunc);

    % Tamanho da FFT (potência de 2 >= Nx+Nh-1)
    Nfft = 2^(nextpow2(Nx+Nh-1));  

    % FFT do sinal e da resposta truncada
    X = fft(x, Nfft);
    H = fft(h_trunc, Nfft);

    % Multiplicação no domínio da frequência
    Y = X .* H;

    % IFFT e truncagem para tamanho correto
    y = real(ifft(Y));
    y = y(1:Nx+Nh-1);
end


%% ===============================================================
% 3. Filtragem do sinal
%% ===============================================================
fprintf('\n=== Seção 3: Filtragem do sinal  ===\n');

% --- 3.1 Filtragem pelas três formas ---

% Eq. Diferenças
y_eqdif = filtragemPorEqDif(x, num, den);

% Convolução com resposta truncada
[y_conv, h_trunc] = filtragemPorConv(x, h);

% FFT
y_fft = filtragemPorFFT(x, h_trunc);

%% --- 3.2 Apresentação dos sinais filtrados no tempo e na frequência ---
Nfft = 16384;  % garantir mesmo tamanho para todos
freqs = (-Nfft/2 : Nfft/2-1) * (fs / Nfft);
freqs_khz = freqs / 1000;

figure('Name','3. Filtragem do sinal','NumberTitle','off','Position',[125 100 1050 600]);

% --- Eq. Diferenças ---
subplot(3,3,1);
plot(t, y_eqdif);
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Eq. Diferenças (tempo)'); grid on;

subplot(3,3,4);
Yeq = fftshift(fft(y_eqdif, Nfft));
plot(freqs_khz, 20*log10(abs(Yeq)));
xlabel('f (kHz)'); ylabel('|Y(f)| dB');
title('Eq. Diferenças (freq)'); grid on;

subplot(3,3,7);
plot(freqs_khz, angle(Yeq));
xlabel('f (kHz)'); ylabel('Fase (rad)');
title('Eq. Diferenças (fase)'); grid on;

% --- Convolução ---
subplot(3,3,2);
plot(t, y_conv(1:N)); 
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Convolução (tempo)'); grid on;

subplot(3,3,5);
Yc = fftshift(fft(y_conv, Nfft));
plot(freqs_khz, 20*log10(abs(Yc)));
xlabel('f (kHz)'); ylabel('|Y(f)| dB');
title('Convolução (freq)'); grid on;

subplot(3,3,8);
plot(freqs_khz, angle(Yc));
xlabel('f (kHz)'); ylabel('Fase (rad)');
title('Convolução (fase)'); grid on;

% --- FFT ---
subplot(3,3,3);
plot(t, y_fft(1:N));
xlabel('Tempo (s)'); ylabel('Amplitude');
title('FFT (tempo)'); grid on;

subplot(3,3,6);
Yf = fftshift(fft(y_fft, Nfft));
plot(freqs_khz, 20*log10(abs(Yf)));
xlabel('f (kHz)'); ylabel('|Y(f)| dB');
title('FFT (freq)'); grid on;

subplot(3,3,9);
plot(freqs_khz, angle(Yf));
xlabel('f (kHz)'); ylabel('Fase (rad)');
title('FFT (fase)'); grid on;

%% --- 3.3 Reprodução dos sinais filtrados ---

tocou = false; % flag para saber se algum áudio foi reproduzido

choice = input('a) Deseja ouvir o áudio após filtragem por Equação de Diferenças? (s/n): ','s');
if lower(choice) == 's'
    fprintf('Reproduzindo áudio filtrado (Eq. Diferenças)...\n');
    sound(y_eqdif, fs); pause(N/fs);
    tocou = true;
end

choice = input('b) Deseja ouvir o áudio após filtragem por Convolução? (s/n): ','s');
if lower(choice) == 's'
    fprintf('Reproduzindo áudio filtrado (Convolução)...\n');
    sound(y_conv, fs); pause(N/fs);
    tocou = true;
end

choice = input('c) Deseja ouvir o áudio após filtragem por FFT? (s/n): ','s');
if lower(choice) == 's'
    fprintf('Reproduzindo áudio filtrado (FFT)...\n');
    sound(y_fft, fs); pause(N/fs);
    tocou = true;
end

if ~tocou
    fprintf('Ok, nenhum áudio será reproduzido.\n');
end