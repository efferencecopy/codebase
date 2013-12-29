%%
fnames = findfile(fnamesFromTxt2('N:\NexFiles\nexfilelists\Zack\SednaRelMotTraining.txt'));
stros = cellfun(@(x) {nex2stro(x)}, fnames);
data = cellfun(@(x) {RelMotTrainingAnalysis(x)}, stros);

guess_prob_per_target = cellfun(@(x) x(4), data);
figure; plot(cell2mat(guess_prob_per_target), '-', 'linewidth', 3);
legend('I (RF)', 'II', 'III', 'IV');
xlabel('file idx'); ylabel('how likely is she guessing?'); title('Sedna RelMotOddManOut training data');

correct_responses_per_target = cellfun(@(x) x(3), data);
figure; plot(100*cell2mat(correct_responses_per_target), '-', 'linewidth', 3);
legend('I (RF)', 'II', 'III', 'IV');
xlabel('file idx'); ylabel('% correct'); title('Sedna RelMotOddManOut training data');