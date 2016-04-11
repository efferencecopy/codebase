%% SPECIFY WHICH EXPERIMENTS SHOULD CONTRIBUTE, LOAD THE DATA

fin

% load in the workbook that contains all the experimental data
wb_path = [GL_DOCUPATH, 'Other_workbooks', filesep, 'wholeCellSTPCellList.xlsx'];
[~,~,wb_expt] = xlsread(wb_path, 1);

% generate the header index info
for i_atrib = 1:size(wb_expt,2)
    fldname = wb_expt{1,i_atrib};
    fldname(isspace(fldname)) = [];
    hidx.(fldname) = i_atrib;
end
    
% now that the header is formed, delete the first row.
wb_expt(1,:) = [];

% make each row it's own cell array so that it can be passed as a single
% argument to the function that does the major unpacking
attributes = {};
for i_ex = 1:size(wb_expt,1)
    attributes{i_ex,1} = wb_expt(i_ex,:);
end


%
% LOAD THE DATA FILES
%

dat = {};
Nexpts = numel(attributes);
for i_ex = 1:Nexpts
    dat{i_ex} = wcstp_compile_data(attributes{i_ex}, hidx);
end





