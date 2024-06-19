clc;
close all
clearvars;

% RS encoder/ decoder parameters
code_word_length = 255;                                                    %codeword
message_length = 239;                                                      %message length
short_message_length = 188;                                                
generator_polynomial = rsgenpoly(code_word_length, message_length, [], 0); %generates the generator polynomial
number_of_words = 200 ;                                                    %Number of words to encode

%RS with Interleaver
rsEncoder = comm.RSEncoder(code_word_length, message_length, generator_polynomial, short_message_length, 'BitInput', true);
rsDecoder_hard = comm.RSDecoder(code_word_length, message_length, generator_polynomial, short_message_length, 'BitInput', true);
rsDecoder_soft = comm.RSDecoder(code_word_length, message_length, generator_polynomial, short_message_length, 'BitInput', true);

%RS without Interleaver
%rsEncodernointr = comm.RSEncoder(code_word_length, message_length, generator_polynomial, short_message_length, 'BitInput', true);
rsDecoder_hardnointr = comm.RSDecoder(code_word_length, message_length, generator_polynomial, short_message_length, 'BitInput', true);
rsDecoder_softnointr = comm.RSDecoder(code_word_length, message_length, generator_polynomial, short_message_length, 'BitInput', true);

% trellis parameters
constraint_length = 7;
generator_polynomial = [171, 133];
puncturing_matrix = [1; 1; 0; 1];
tb_depth = 96;

% simulation parameters
num_of_bits = log2(code_word_length + 1) * short_message_length * number_of_words;
snr_db = 0:1:7;

% filter design parameters
filter_length = 6;
over_sampling_factor = 8;
roll_off_factor = 0.35;

% creating a ber array
ber_hard_values = zeros(1, numel(snr_db));
ber_soft_values = zeros(1, numel(snr_db));
ber_theory = berawgn(snr_db, 'psk', 4, 'nondiff');     %Comparing with uncoded BER


% generating data
bits = randi([0 1], num_of_bits, 1);

for each_iteration = 1:numel(snr_db)
    snr_value = 10^(0.1 * snr_db(each_iteration));

    rs_encoded_bits = rsEncoder(bits);                                      %RS ENCODER
    interleaved_bits = convintrlv(rs_encoded_bits, 1, 1);                   %Interleaver
    trellis = poly2trellis(constraint_length, generator_polynomial);
    encoded_data = [convenc(interleaved_bits, trellis, puncturing_matrix)];          %Viterbi encoding
    symbols = complex(1 - 2 * encoded_data(1:2:end), 1 - 2 * encoded_data(2:2:end)); %Symbols with Interleaver

    encoded_datanointr = [convenc(rs_encoded_bits, trellis, puncturing_matrix)];                       %RS without ENCODER
    symbolsnointr = complex(1 - 2 * encoded_datanointr(1:2:end), 1 - 2 * encoded_datanointr(2:2:end)); %Symbols without Interleaver


     channel_filter = rcosdesign(roll_off_factor, filter_length, over_sampling_factor);     
     tx_signal = upfirdn(symbols, channel_filter, over_sampling_factor);                %Tx with interleaver
     tx_signalnointr = upfirdn(symbolsnointr, channel_filter, over_sampling_factor);    %Tx without interleaver
     
     %Rx with interleaver
    signal_power = mean(abs(tx_signal).^2);
    noise_power = 2 * (signal_power * over_sampling_factor) / (2 * 2 * snr_value);
    awgn_noise = sqrt(noise_power) * complex(randn(length(tx_signal), 1), randn(length(tx_signal), 1));
    rx_signal = upfirdn(tx_signal + awgn_noise, channel_filter, 1, over_sampling_factor);

    rx_symbols = zeros(numel(encoded_data), 1);
    rx_symbols(1:2:end) = real(rx_signal(filter_length+1:end-filter_length));
    rx_symbols(2:2:end) = imag(rx_signal(filter_length+1:end-filter_length));

    hard_decoded_bits_viterbi = vitdec(rx_symbols < 0, trellis, tb_depth, 'trunc', 'hard', puncturing_matrix);  %Viterbi encoder with interleaver
    soft_decoded_bits_viterbi = vitdec(rx_symbols, trellis, tb_depth, 'trunc', 'unquant', puncturing_matrix);

    %Rx without interleaver
    signal_powernointr = mean(abs(tx_signalnointr).^2);
    noise_powernointr = 2 * (signal_powernointr * over_sampling_factor) / (2 * 2 * snr_value);
    awgn_noisenointr = sqrt(noise_powernointr) * complex(randn(length(tx_signalnointr), 1), randn(length(tx_signalnointr), 1));
    rx_signalnointr = upfirdn(tx_signalnointr + awgn_noisenointr, channel_filter, 1, over_sampling_factor);

    rx_symbolsnointr = zeros(numel(encoded_datanointr), 1);
    rx_symbolsnointr(1:2:end) = real(rx_signalnointr(filter_length+1:end-filter_length));
    rx_symbolsnointr(2:2:end) = imag(rx_signalnointr(filter_length+1:end-filter_length));

    hard_decoded_bits_viterbnointr = vitdec(rx_symbolsnointr < 0, trellis, tb_depth, 'trunc', 'hard', puncturing_matrix);  %Viterbi encoder without interleaver
    soft_decoded_bits_viterbinointr = vitdec(rx_symbolsnointr, trellis, tb_depth, 'trunc', 'unquant', puncturing_matrix);

    %Deinterleaver
    hard_deinterleaved_bits = convdeintrlv(double(hard_decoded_bits_viterbi), 1, 1);
    soft_deinterleaved_bits = convdeintrlv(double(soft_decoded_bits_viterbi), 1, 1);

    %RS Decoder with interleaver
    hard_decoded_bits = rsDecoder_hard(hard_deinterleaved_bits(1:end));
    soft_decoded_bits = rsDecoder_soft(soft_deinterleaved_bits(1:end));
    fprintf("Hard Decoding BER @SNR: %f = %f\n", snr_db(each_iteration), length(find(abs(hard_decoded_bits - bits) > 0)) / numel(bits));
    %fprintf("Soft Decoding BER @SNR: %f = %f\n", snr_db(each_iteration), length(find(abs(soft_decoded_bits - bits) > 0)) / numel(bits));

    %RS Decoder without interleaver
    hard_decoded_bitsnointr = rsDecoder_hardnointr(hard_decoded_bits_viterbnointr(1:end));
    soft_decoded_bitsnointr = rsDecoder_softnointr(soft_decoded_bits_viterbinointr(1:end));
    fprintf("Hard Decoding BER without interleaver@SNR: %f = %f\n", snr_db(each_iteration), length(find(abs(hard_decoded_bitsnointr - bits) > 0)) / numel(bits));
    % fprintf("Soft Decoding BER @SNR: %f = %f\n", snr_db(each_iteration), length(find(abs(soft_decoded_bitsnointr - bits) > 0)) / numel(bits));
    
    %BER with interleaver
    ber_hard_values(each_iteration) = length(find(abs(hard_decoded_bits - bits) > 0)) / numel(bits);
    ber_soft_values(each_iteration) = length(find(abs(soft_decoded_bits - bits) > 0)) / numel(bits);

    %BER without interleaver
    ber_hard_valuesnointr(each_iteration) = length(find(abs(hard_decoded_bitsnointr - bits) > 0)) / numel(bits);
    ber_soft_valuesnointr(each_iteration) = length(find(abs(soft_decoded_bitsnointr - bits) > 0)) / numel(bits);
end

figure;
semilogy(snr_db, ber_theory, 'LineWidth', 2, 'color', 'red', 'LineStyle', ':', 'Marker', '>');
hold on;
grid on;
semilogy(snr_db, ber_hard_values, 'LineWidth', 2, 'color', 'green', 'LineStyle', '-.', 'Marker', '^');
semilogy(snr_db, ber_soft_values, 'LineWidth', 2, 'color', 'black', 'LineStyle', '--', 'Marker', '*');
semilogy(snr_db, ber_hard_valuesnointr, 'LineWidth', 2, 'color', 'blue', 'LineStyle', '-.', 'Marker', 'o');
semilogy(snr_db, ber_soft_valuesnointr, 'LineWidth', 2, 'color', 'magenta', 'LineStyle', '--', 'Marker', '+');
legend("Theoretical BER", "Hard Decoding BER with Interleaver", "Soft Decoding BER with Interleaver","Hard Decoding BER without Interleaver", "Soft Decoding BER without Interleaver");
title('BER PLOT')
xlabel('Eb/No')
ylabel('BER')