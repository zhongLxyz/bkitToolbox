% returns the 4 quantities needed for converting the values into the corresponding
% perception in each dataset
function [angularVel, azimuth, response, elevation] = returnValues(dataset)

elevation = dataset.CameraElevation;
angularVel = dataset.AngularVelocity;
azimuth = dataset.CameraAzimuth;
response = dataset.Response;

end