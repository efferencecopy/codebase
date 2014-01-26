function ID = whoami


switch getmacaddr
    case '345'
        ID = 'glick_rig1';
        
    case 'c4:2c:03:2a:29:3e'
        ID = 'hass_mbp';
        
    case '123'
        
    case '00-25-90-C1-B6-C2'
        ID = 'nuke';
        
    otherwise
        ID = '';
end