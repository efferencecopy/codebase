function mac_add = getmacaddr()

% this should be a multi-platform way of obtaining the mac address so that
% I can easily identify which computer I'm at. I stole the code from this
% blog:
%
% http://undocumentedmatlab.com/blog/unique-computer-id/

switch computer('arch')
    case {'maci','maci64'}
        
        [~,a]=system('ifconfig');
        c=strfind(a,'en0');
        if ~isempty(c)
            a=a(c:end);
        end
        c=strfind(a,'en1');
        if ~isempty(c)
            a=a(1:c-1);
        end
        
        % find the mac address
        b=strfind(a,'ether');
        mac_add=a(1,b(1)+6:b(1)+22);
        
    case {'win32','win64'}
        
        [~,a]=system('getmac');
        b=strfind(a,'=');
        mac_add=a(b(end):b(end)+19);
        
        % might need to be this line:
        mac_add = a(b(end)+1:b(end)+19);
        
    case {'glnx86','glnxa64'}
        [~,a]=system('ifconfig');
        b=strfind(a,'Ether');
        mac_add=a(1,b(1)+17:b(1)+33);
        
    otherwise
        mac_add=[];
end