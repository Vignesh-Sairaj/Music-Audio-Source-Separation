clearvars;

%% Takes output of STFT(Magnitude)
% Stores in file facts.mat (Look below)


sample = load('./stftMag.mat');



M = double(sample.stftMag);
r = sample.numComps; % # of basis columns
[m, n] = size(M);



%% ==================           SNPA         ==================================================
disp('Running snpa...');
[J_snpa, H_snpa_1] = SNPA(M, r, 1);
H_snpa= nnlsHALSupdt(M, M(:,J_snpa));
W_snpa = M(:, J_snpa);


%% ================== If Greedy too slow, temporarily ========================================================

% J_greedy = J_snpa; H_greedy = H_snpa; W_greedy = W_snpa;



%% ==================        Greedy SNPA         ==============================================
disp('Running greedy...');
[J_greedy, H_greedy_1] = greedy2passSNPA(M, r, 1);
H_greedy = nnlsHALSupdt(M, M(:,J_greedy));



%% ==================         Errors         ==================================================
disp('Computing errors...');

scale = diag(1./sqrt(sum(M.^2)));

% l2_err2_snpa = sqrt(sum(sum(M-M(:, J_snpa)*H_snpa).^2))/sqrt(sum(sum(M.^2)))
% l2_err3_snpa = sum(sum(M-M(:, J_snpa)*H_snpa).^2)/sum(sum(M.^2))

l2_err_snpa = sqrt(sum(sum(((M-M(:, J_snpa)*H_snpa)*scale).^2)))/sqrt(sum(sum((M*scale).^2)));
l2_err_greedy = sqrt(sum(sum(((M-M(:, J_greedy)*H_greedy)*scale).^2)))/sqrt(sum(sum((M*scale).^2)));


l1_err_snpa = sum(sum(abs(M-M(:, J_snpa)*H_snpa)))/sum(sum(abs(M)));
l1_err_greedy = sum(sum(abs(M-M(:, J_greedy)*H_greedy)))/sum(sum(abs(M)));

%% ============================================================================================

save('facts.mat', 'M', 'r', 'J_greedy', 'H_greedy', 'J_snpa', 'H_snpa');
save('errors.mat', 'l2_err_snpa', 'l2_err_greedy', 'l1_err_snpa', 'l2_err_greedy');

