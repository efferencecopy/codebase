% Here's a function that tests all permutations of inputs to TranslateToColourModeMex
function test_ttcmmex() %#ok<*NASGU>
rand8bit = randi(2^8, [randi([50 300], 1, 2) 3])-1;
rand8bit2d = randi(2^8, randi([50 300], 1, 2))-1;
rand16bit = randi(2^16, [randi([50 300], 1, 2) 3])-1;
rand16bit2d = randi(2^16, randi([50 300], 1, 2))-1;

% Notice there aren't tests for empty input arguments: those are easily verified manually.
input_names = {'rand8bit', 'rand8bit2d', 'rand16bit', 'rand16bit2d', ...
    'rand16bit/65535', 'rand16bit2d/65535', ...
    'uint8(rand8bit)', 'uint8(rand8bit2d)', ...
    'uint16(rand8bit)', 'uint16(rand8bit2d)', 'uint16(rand16bit)', 'uint16(rand16bit2d)'
    };
useoldstyle16bit = {0 1};
ccmodes = {-1 0 1 2};

func_args = allcombs(input_names, useoldstyle16bit, ccmodes);
% The order of the tests depends on the order of func_args given by `allcombs.`
% `allcombs` gives the Cartesian product of n sets in a depth-first manner:
% 1st element of set A1, ..., 1st element of set A(n-1), product with all elements of set An,
% 1st element of set A1, ..., 2nd element of set A(n-1), product with all elements of set An, and so on.
tests = {
'isequal(output(:,1:2:end,:),rand8bit) & isequal(output(:,2:2:end,:),zeros(size(rand8bit)))';
'isequal(output,rand8bit/255)';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit)))';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),rand8bit/255)';
'isequal(output(:,1:2:end,:),rand8bit) & isequal(output(:,2:2:end,:),zeros(size(rand8bit)))';
'isequal(output,rand8bit/255)';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit)))';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),rand8bit/255)';
'isequal(output(:,1:2:end,:),rand8bit2d) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d)))';
'isequal(output,rand8bit2d/255)';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d)))';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),rand8bit2d/255)';
'isequal(output(:,1:2:end,:),rand8bit2d) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d)))';
'isequal(output,rand8bit2d/255)';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d)))';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),rand8bit2d/255)';
'isequal(output(:,1:2:end,:),floor(rand16bit/256)) & isequal(output(:,2:2:end,:),rem(rand16bit,256))';
'isequal(output,rand16bit/65535)';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit)))';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,2:2:end,:),rand16bit/65535)';
'isequal(output(:,1:2:end,:),floor(rand16bit/256)) & isequal(output(:,2:2:end,:),rem(rand16bit,256))';
'isequal(output,rand16bit/65535)';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit)))';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,2:2:end,:),rand16bit/65535)';
'isequal(output(:,1:2:end,:),floor(rand16bit2d/256)) & isequal(output(:,2:2:end,:),rem(rand16bit2d,256))';
'isequal(output,rand16bit2d/65535)';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit2d)))';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,2:2:end,:),rand16bit2d/65535)';
'isequal(output(:,1:2:end,:),floor(rand16bit2d/256)) & isequal(output(:,2:2:end,:),rem(rand16bit2d,256))';
'isequal(output,rand16bit2d/65535)';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit2d)))';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,2:2:end,:),rand16bit2d/65535)';
'isequal(output,rand16bit/65535)';
'isequal(output,rand16bit/65535)';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit)))';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,1:2:end,:),rand16bit/65535)';
'isequal(output,rand16bit/65535)';
'isequal(output,rand16bit/65535)';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit)))';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,1:2:end,:),rand16bit/65535)';
'isequal(output,rand16bit2d/65535)';
'isequal(output,rand16bit2d/65535)';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit2d)))';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,1:2:end,:),rand16bit2d/65535)';
'isequal(output,rand16bit2d/65535)';
'isequal(output,rand16bit2d/65535)';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit2d)))';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,1:2:end,:),rand16bit2d/65535)';
'isequal(output(:,1:2:end,:),uint8(rand8bit)) & isequal(output(:,2:2:end,:),zeros(size(rand8bit),''uint8''))';
'isequal(output,rand8bit/255)';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit)))';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),rand8bit/255)';
'isequal(output(:,1:2:end,:),uint8(rand8bit)) & isequal(output(:,2:2:end,:),zeros(size(rand8bit),''uint8''))';
'isequal(output,rand8bit/255)';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit)))';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),rand8bit/255)';
'isequal(output(:,1:2:end,:),uint8(rand8bit2d)) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d),''uint8''))';
'isequal(output,rand8bit2d/255)';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d)))';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),rand8bit2d/255)';
'isequal(output(:,1:2:end,:),uint8(rand8bit2d)) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d),''uint8''))';
'isequal(output,rand8bit2d/255)';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d)))';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),rand8bit2d/255)';
'isequal(output(:,1:2:end,:),uint16(rand8bit)) & isequal(output(:,2:2:end,:),zeros(size(rand8bit),''uint16''))';
'isequal(output,rand8bit/255)';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit)))';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),rand8bit/255)';
'isequal(output(:,1:2:end,:),uint16(rand8bit)) & isequal(output(:,2:2:end,:),zeros(size(rand8bit),''uint16''))';
'isequal(output,rand8bit/255)';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit)))';
'isequal(output(:,1:2:end,:),rand8bit/255) & isequal(output(:,2:2:end,:),rand8bit/255)';
'isequal(output(:,1:2:end,:),uint16(rand8bit2d)) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d),''uint16''))';
'isequal(output,rand8bit2d/255)';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d)))';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),rand8bit2d/255)';
'isequal(output(:,1:2:end,:),uint16(rand8bit2d)) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d),''uint16''))';
'isequal(output,rand8bit2d/255)';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),zeros(size(rand8bit2d)))';
'isequal(output(:,1:2:end,:),rand8bit2d/255) & isequal(output(:,2:2:end,:),rand8bit2d/255)';
'isequal(output(:,1:2:end,:),uint16(floor(rand16bit/256))) & isequal(output(:,2:2:end,:),uint16(rem(rand16bit,256)))';
'isequal(output,rand16bit/65535)';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit)))';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,2:2:end,:),rand16bit/65535)';
'isequal(output(:,1:2:end,:),uint16(floor(rand16bit/256))) & isequal(output(:,2:2:end,:),uint16(rem(rand16bit,256)))';
'isequal(output,rand16bit/65535)';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit)))';
'isequal(output(:,1:2:end,:),rand16bit/65535) & isequal(output(:,2:2:end,:),rand16bit/65535)';
'isequal(output(:,1:2:end,:),uint16(floor(rand16bit2d/256))) & isequal(output(:,2:2:end,:),uint16(rem(rand16bit2d,256)))';
'isequal(output,rand16bit2d/65535)';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit2d)))';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,2:2:end,:),rand16bit2d/65535)';
'isequal(output(:,1:2:end,:),uint16(floor(rand16bit2d/256))) & isequal(output(:,2:2:end,:),uint16(rem(rand16bit2d,256)))';
'isequal(output,rand16bit2d/65535)';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,2:2:end,:),zeros(size(rand16bit2d)))';
'isequal(output(:,1:2:end,:),rand16bit2d/65535) & isequal(output(:,2:2:end,:),rand16bit2d/65535)';
};

assert(length(tests) == length(func_args), ...
    'The number of tests don''t match the number of function argument sets');

test_results = zeros(length(tests),1);
for ii = 1:length(tests)
    output = eval(sprintf('TranslateToColourModeMex(%s,%d,%d)', func_args{ii,:}));
    test_results(ii) = eval(tests{ii});
end

failures = find(~test_results)';
if isempty(failures)
    fprintf('******** no failed tests ********\n');
else
    fprintf('inputs that failed tests:\n');
    for failure = failures
        fprintf('\t%s\t%d\t%d\n', func_args{failure,:});
    end
end
