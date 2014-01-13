% Code to plot pre- and post-habituation detection data. I'm trying to do
% this in a way that's comparable to Krauskopf et., al.

controlFiles = {'C032710001.nex',...
                'C032710002.nex',...
                'C032710003.nex'};
habitFiles = {'C032810001.nex',...
              'C032810002.nex',...
              'C032810003.nex'};
          
          
for a = 1:length(controlFiles)
    [t{a}, clrCont{a}, sfs{a}] = habitUnpack(nex2stro(findfile(controlFiles{a})));
end

%the folowing manipulation only works when there is a single sf and the
%colors for each exp are identical.
thresh = horzcat(t{:}); %dims = nColors x nExpts
thresh = permute(thresh, [1,3,2]);
log_Tcontrol = log10(thresh);
avg_log_Tcontrol = mean(log_Tcontrol,3);

for a = 1:length(habitFiles)
    [t{a}, clrHabit{a}, sfs{a}] = habitUnpack(nex2stro(findfile(habitFiles{a})));
end
%the folowing manipulation only works when there is a single sf and the
%colors for each exp are identical.
thresh = horzcat(t{:}); %dims = nColors x nExpts
thresh = permute(thresh, [1,3,2]);
log_Thabit = log10(thresh);
avg_log_Thabit = mean(log_Thabit,3);

