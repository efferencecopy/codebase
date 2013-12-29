%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return a 1 if the user is pressing the ESC key, 0 otherwise.
function stopnow = CheckForESCKey()
   [keyisdown, secs, keycode] = KbCheck();
   stopnow = keycode(27);  % keycode(27) = escape key
   % GDLH Added a means for debugging when things break.
 %  if (stopnow)
 %      keyboard;
 %  end
end