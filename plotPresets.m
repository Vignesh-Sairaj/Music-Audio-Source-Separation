load('./presetscores.mat');

figure(); plot([nmf(:, 1) snpa(:, 1) greedy(:, 1)]); legend('nmf', 'snpa', 'greedy'); title('Min SNR'); xlabel('iter'); ylabel('snr');
figure(); plot([nmf(:, 2) snpa(:, 2) greedy(:, 2)]); legend('nmf', 'snpa', 'greedy'); title('Mean SNR'); xlabel('iter'); ylabel('snr');
figure(); plot([nmf(:, 3) snpa(:, 3) greedy(:, 3)]); legend('nmf', 'snpa', 'greedy'); title('Max SNR'); xlabel('iter'); ylabel('snr');
figure(); plot([nmf(:, 4) snpa(:, 4) greedy(:, 4)]); legend('nmf', 'snpa', 'greedy'); title('SNR Reconstr'); xlabel('iter'); ylabel('snr');
