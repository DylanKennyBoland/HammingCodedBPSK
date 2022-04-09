%% Digital Communications (EEEN40060) Assignment 2:
% This program attempts to simulate a Hamming-coded Binary Phase-Shift
% Keying (BPSK) system. It also makes a comparison with a non-encoded BPSK
% system. The channel model used is AWGN (Additive white Gaussian noise).
%
% Author: Dylan Boland (Student)

% Boolean variables below allow the user to decide whether to use a Hex
% code to specify the colour, or an RGB triplet:
plotBERcurves = true;
plotSimulatedBERcurve = true;
% Boolean variable to control how the colour is specified
useHexColourCode = true;

if (useHexColourCode == true)
    simulatedBERcolour = "#66ff33";
    theoreticalBERcolour = "#6600ff";
else % in case the user wants to use RGB triplets instead...
    simulatedBERcolour = [0.3 0.9 0];
    theoreticalBERcolour = [0.3 0 1];
end

numDataBits = 48000; % should be a multiple of 4 for convenience
numDataBits = numDataBits - mod(numDataBits, 4); % help make divisible by 4
txBits = randi([0 1], 1, numDataBits); % our bitstream...

% Define the generator matrix for the Hamming encoding scheme. The Hamming
% code is a (7,4) block code. It encodes "blocks" of bits. 4 come in, 7 go
% out. The code rate is therefore: r = 4/7
G = [1 1 0 1 0 0 0; 0 1 1 0 1 0 0; 1 1 1 0 0 1 0; 1 0 1 0 0 0 1];
P = G(:, 1:size(G, 2)-size(G, 1)); % the P matrix
H = [eye(size(G, 2)-size(G, 1)) P']; % the parity check matrix

% Computing these once so that we do not need to keep calling the size()
% function:
msgWordLen = size(G, 1);
codeWordLen = size(G, 2);

% Pass the bit stream to the Hamming encoder. In order to allow for
% efficient matrix multiplication, let us reshape the txBit vector into a
% [4 X numDataBits/4] matrix:
numBlocks = numDataBits/msgWordLen;
txBlock = reshape(txBits, msgWordLen, numBlocks);
encodedTxBlock = zeros(codeWordLen, numBlocks); % encoded bits go in here

% Encoding stage:
for block = 1:numBlocks % block-by-block encoding
    % mod() function, with "2" as the second argument will allow us to do
    % modulo-2 arithmetic, needed for doing binary arithmetic:
    encodedTxBlock(:, block) = mod(sum(txBlock(:, block).*G, 1), 2);
end

% Modulation stage:
% I would like all the 0s to map to "-1", and all the 1s to map to "1". A
% general rule in order to map each bit to the correct BPSK symbol is:
% e -> e + (e - 1)
%
% For example, if I have a "0". Then the symbol will be: 0 + (0-1) = -1.
%
% And if I have a "1", then the symbol will be: 1 + (1-1) = 1. Apply the
% above idea to form the modulated transmit stream:
txStream = encodedTxBlock(:) + (encodedTxBlock(:) - 1);
Es = mean(abs(txStream).^2); % average symbol energy
Eb = Es; % average energy in a transmit bit
r = msgWordLen/codeWordLen; % the code rate

EbNo = 10.^(-1:0.2:3); % a vector of Eb/No values
No = Eb./(r*EbNo); % the corresponding vector of No values based on the value of Eb...

simulatedBER = zeros(length(No), 1); % a vector to hold the BER values
rxBlock = zeros(codeWordLen, numBlocks); % matrix to hold received codewords

% Now to simulate the BER for the various values of the noise power density
% (No):
for j = 1:length(No)
    sigma = sqrt(No(j)/2);
    n = sigma*randn(length(txStream), 1); % the noise vector
    rxStream = txStream + n; % the received stream
    % Demodulation stage. The decision boundary is half way between 1 and
    % -1 - i.e., 0.
    rxBits = double(rxStream > 0); % this line will return a vector of 0s and 1s. If a
    % given element of "rxStream" is to the right of "0", the corresponding
    % element in "rxBits" will be 1, as we would want. Likewise, if a given
    % element in "rxStream" is to the left of 0, then after this comparison
    % the corresponding element in "rxBits" will be 0, as we would also
    % like. This is because during the modulation stage a 1 mapped to "1"
    % and a 0 mapped to a "-1":
    rxBlock = reshape(rxBits, codeWordLen, numBlocks);
    for block = 1:numBlocks
        s = mod(sum(rxBlock(:, block).*H', 1), 2); % compute the Syndrome (s)
        % Instead of calling the mod() and sum() function on each iteration
        % of the current for loop, we will wait until we exit, and then
        % call the mod() function once on the rxBlock matrix:
        switch num2str(s)
            case "0  0  1"
                rxBlock(:, block) = rxBlock(:, block) + [0, 0, 1, 0, 0, 0, 0]';
            case "0  1  0"
                rxBlock(:, block) = rxBlock(:, block) + [0, 1, 0, 0, 0, 0, 0]';
            case "0  1  1"
                rxBlock(:, block) = rxBlock(:, block) + [0, 0, 0, 0, 1, 0, 0]';
            case "1  0  0"
                rxBlock(:, block) = rxBlock(:, block) + [1, 0, 0, 0, 0, 0, 0]';
            case "1  0  1"
                rxBlock(:, block) = rxBlock(:, block) + [0, 0, 0, 0, 0, 0, 1]';
            case "1  1  0"
                rxBlock(:, block) = rxBlock(:, block) + [0, 0, 0, 1, 0, 0, 0]';
            case "1  1  1"
                rxBlock(:, block) = rxBlock(:, block) + [0, 0, 0, 0, 0, 1, 0]';
        end   
    end
    % Now we modulo division by 2, once, on the received block...
    rxBlock = mod(rxBlock, 2);
    % We used a systematic encoder. This means that after a message word is
    % encoded, it appears either at the beginning or end of the codeword.
    % In our case, it appears at the end of the codeword:
    msgRxBlock = rxBlock(msgWordLen:end, :); % extract the message words
    numErrors = numDataBits - sum(msgRxBlock(:) == txBits(:));
    simulatedBER(j) = numErrors/numDataBits; % add BER value to vector
end

theoreticalBER = qfunc(sqrt(2*EbNo)); % compute the theoretical BER values

if (plotBERcurves == true)
    figure(1)
    semilogy(10*log10(EbNo), simulatedBER, '*', 'MarkerSize', 10, 'MarkerFaceColor', simulatedBERcolour, 'LineStyle', '--', 'color', simulatedBERcolour)
    title("\fontsize{14}\fontname{Georgia}Bit Error Rate (BER) Vs E_{b}/N_{0} (Number of data bits transmitted = " + numDataBits +")");
    xlabel('\fontname{Georgia}\bf E_{b}/N_{0} (dB)');
    ylabel('\fontname{Georgia}\bf BER');
    set(gca,'Fontname', 'Georgia');
    ylim([10^-6 10^0]); % setting the limits, so that the graph is clearer...
    hold on
    semilogy(10*log10(EbNo), theoreticalBER, '*', 'MarkerSize', 10, 'MarkerFaceColor', theoreticalBERcolour, 'LineStyle', '--', 'color', theoreticalBERcolour)
    legend("\fontsize{14}\fontname{Georgia}\bfSimulated BER (Hamming encoding)", "\fontsize{14}\fontname{Georgia}\bfTheoretical SER (No encoding)");
    hold off
end

if (plotSimulatedBERcurve == true)
    figure(2)
    semilogy(10*log10(EbNo), simulatedBER, '*', 'MarkerSize', 10, 'MarkerFaceColor', simulatedBERcolour, 'LineStyle', '--', 'color', simulatedBERcolour)
    title("\fontsize{14}\fontname{Georgia}Bit Error Rate (BER) Vs E_{b}/N_{0} (Number of data bits transmitted = " + numDataBits +")");
    xlabel('\fontname{Georgia}\bf E_{b}/N_{0} (dB)');
    ylabel('\fontname{Georgia}\bf BER');
    set(gca,'Fontname', 'Georgia');
    legend("\fontsize{14}\fontname{Georgia}\bfSimulated BER (Hamming encoding)");
    ylim([10^-6 10^0]); % setting the limits, so that the graph is clearer...
end