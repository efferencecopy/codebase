function gitStatusReport()

presdir = pwd;

% grab the status of the git repository
codeDir = which('fin');
fsep = regexp(codeDir, filesep);
codeDir = codeDir(1:fsep(end-1)-1);


cd(codeDir)
[~, stdOut] = system('git status');
lineBreaks = regexpi(stdOut, '\n');

% get the branch name
branch = stdOut(1:lineBreaks(1)-1);
namestart = regexpi(branch, 'On branch ', 'end');
branch = branch(namestart+1:end);

% what's the status of the local repository
gitDirectoryClean = any(regexpi(stdOut, 'working directory clean'));

% alert the user
switch lower(branch)
    case 'master'
        if gitDirectoryClean
            fprintf(stdOut)
        else
            fprintf('\n ########## GIT STATUS REPORT ########## \n\n')
            warning('git local repository ''MASTER'' needs attention')
        end
    otherwise
        if gitDirectoryClean
            fprintf('\n ########## GIT STATUS REPORT ########## \n\n')
            warning(' on git  repository %s', branch)
        else
            fprintf('\n ########## GIT STATUS REPORT ########## \n\n')
            warning(' git repository %s need attention', upper(branch))
        end
end

cd(presdir)

