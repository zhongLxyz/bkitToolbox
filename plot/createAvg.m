% find % of biases for each participant in each condition.
for i = 1:length(cw12)
    
    cw = cw12{i};
    facing = facing12{i};
    above = above12{i};
    
    biases = zeros(3,3);
    
    biases(1,1) = mean(findCCW(cw.AngularVelocity,cw.CameraAzimuth,cw.Response,cw.CameraElevation));
    biases(1,2) = mean(findFTV(cw.AngularVelocity,cw.CameraAzimuth,cw.Response,cw.CameraElevation));
    biases(1,3) = mean(findVFA(cw.AngularVelocity,cw.CameraAzimuth,cw.Response,cw.CameraElevation));
    
    biases(2,1) = mean(findCCW(facing.AngularVelocity,facing.CameraAzimuth,facing.Response,facing.CameraElevation,'ftv'));
    biases(2,2) = mean(findFTV(facing.AngularVelocity,facing.CameraAzimuth,facing.Response,facing.CameraElevation,'ftv'));
    biases(2,3) = mean(findVFA(facing.AngularVelocity,facing.CameraAzimuth,facing.Response,facing.CameraElevation,'ftv'));
    
    biases(3,1) = mean(findCCW(above.AngularVelocity,above.CameraAzimuth,above.Response,above.CameraElevation,'vfa'));
    biases(3,2) = mean(findFTV(above.AngularVelocity,above.CameraAzimuth,above.Response,above.CameraElevation,'vfa'));
    biases(3,3) = mean(findVFA(above.AngularVelocity,above.CameraAzimuth,above.Response,above.CameraElevation,'vfa'));
    
    output(i) = {biases};
end