function converted = convertAzimuth(origAzimuth)
% this function convert the azimuth back to -30, 0 and 30 in Zhong's 3rd
% experiment.
% input: origAzimuth - original column of azimuth for one participant
% output: converted - the converted column of azimuth

origAzimuth(origAzimuth==-40 | origAzimuth==-20) = -30;
origAzimuth(origAzimuth==-10| origAzimuth==10) = 0;
origAzimuth(origAzimuth==20| origAzimuth==40) = 30;
origAzimuth(origAzimuth==-140| origAzimuth==-160) = -150;
origAzimuth(origAzimuth==-170| origAzimuth==170) = 180;
origAzimuth(origAzimuth==160| origAzimuth==140) = 150;
converted = origAzimuth;