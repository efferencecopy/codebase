function sendToRex(udpCom, toRex, dataType, msg)

    %      SEND REX A PTB CALIBRATION STRUCTURE
    % 
    %  EXAMPLE: sendToRex(udpCom, toRex.gammaTable, 'double', message); ==> Sends the gammaTable
    %           sendToRex(udpCom, toRex.numFrames, 'integer', message); ==> initiates transfer
    %
    % This script is designed to feed Rex the pertinent components of a
    % calibration structure. These components include the GammaTable, the
    % monitor emission spectra, and the M matrix. The message is sent as a
    % hexadecimal character array with elements of the cal structure delimited
    % by spaces. This function requires udpCom to be a structure array that
    % specifies the udp socket, port and rexip. toRex should be an element of
    % the calibration structure that the calling function (to sendCals) has
    % access to.
    %
    % CAH 11/20/07

    %******** CONSTANTS *******%
    REXBUFF = 8000;       %The size of the UDP read buffer on the Rex machine in bytes

    %**************************%
    
    %Turn the numbers that rex requests into a long string of space
    %delimited hexadecimal numbers. Do this differently for integers and
    %doubles
    if strcmpi(dataType, 'double')
        calElement = sprintf(' %bx', toRex); %front pad with a space!
%     elseif ~isempty(strmatch(dataType, strvcat('int', 'integer'), 'exact'));
    elseif any(strcmpi(dataType, {'int' 'integer' 'long'}))
        calElement = mydec2hex(toRex);
    else
        error('Unspecified or inadequate data type specifier');
    end
    
    %The first udp message will be from Rex's main state set requesting a
    %number or list of numbers. Echo this request (i.e the EXACT same
    %string that Rex sent) so that Rex can call the appropriate subfunction
    %from with in the udp "listening" state set. If the requested data will
    %fit with this "echo" go ahead and send them together. This will force
    %the rex to use the small-transfer protocol instead of the large-batch
    %transfer protocal.
    maxMsgLength = REXBUFF-1; %buffer size minus \0
    numElements = mydec2hex(numel(toRex));
    totalLengthInChar = length(calElement) + length(msg) + numElements + 5; %four for both >>'s one for \0
    if (totalLengthInChar < maxMsgLength); 
        pnet(udpCom.sock, 'write', [msg, '>>', numElements, '>>' calElement]);
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        return %exit here for small-batches
    else
        pnet(udpCom.sock, 'write', [msg, '>>', numElements, '>>']);
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        pause(0.005); %pause to let Rex catch up.
    end
    
    %Chop up the cal struct into the appropriate number of segments, then
    %send it off. The last segment will cause an out of bounds error unless
    %length(calElement) is an integer multiple of maxMsgLength. Avoid this
    %error by catching the last chunk and including just the last few
    %elements in the snippet.
    for a = 1:maxMsgLength:length(calElement);
        try
            chunk = calElement(a:(a+maxMsgLength-1));
        catch
            chunk = calElement(a:length(calElement));
        end
      
        pnet(udpCom.sock, 'write', chunk)
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        pause(0.005); %avoid Rex buffer overflow issues (necessary?)
    end

    %if the length of calElement is greater than the size of the Rex read
    %buffer Rex will have to concatinate the transmitted set of hex numbers
    %into a character array in a while loop. In order to break out of this
    %while loop send a message that indicates that the transmission has
    %been sent in it's entirety.
    pnet(udpCom.sock, 'write', 'transferComplete'); %must be identical to that used by the transmission functions on rex
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    
    
    %questions  **************************
    %how long of a pause is necessary (if any at all)
    
    %remove the header convention. The transmission is uniquely defined by
    %the "message" echo so just use "message" as the headder. Delimit by
    %>>. Make a function on the Rex side that takes off the header from UDP
    %transmissions and returns the header as well as the rest of the
    %string. evaluate the string based off the header.
    
