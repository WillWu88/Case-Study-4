% simulating light proporgation through optical lens on a 2-d plane

clear;
close all;

% all units are in milimeters
focalLength = 100;
z1 = 200;
% assuming the lens is positioned at 0, has thickness of zero

% testing model, at double focal lens using lens property
z2 = 1 * z1;

Mz1 = [1, z1; 0, 1];
Mf = [1, 0; -1 / focalLength, 1];
Mz2 = [1, z2; 0, 1];
M = Mz2 * Mf * Mz1;

% assumptions:
% lens is 2mm tall, centered at 0
rays_ypos = [1 1 1];
rays_angle = [0 (1/(2*focalLength)) (1/focalLength)];

[rays_out_ypos, rays_out_angle] = simRayProp(M,rays_ypos,rays_angle);


z1Range = [1:0.1:1000];
z2 = zeros(length(z1Range),1);
theoreticalZ2 = zeros(length(z1Range),1);
for z = 1:length(z1Range)
    z2(z) = solveZ2(z1Range(z),focalLength);
    theoreticalZ2(z) = 1/((1/focalLength) - (1/z1Range(z)));
end

figure('Name','Z1 vs Z2 from 1mm to 1m, theo overlay');
plot(z1Range,z2);
axis([0,1000, -500, 500]);
pause(3);
hold on;
plot(z1Range,theoreticalZ2);
hold off;
legend('Real','Theoretical');

%{
theoreticalmag = zeros(length(z1Range),1);

magCoeff = zeros(length(z1Range),1);

for s = 1:length(z1Range);
    mz1 = [1,z1Range(s);0,1];
    mf = [1,0;(-1/focalLength),1];
    mz2 = [1,z2(s);0,1];
    
    updateM = mz2*mf*mz1;
    
    [y2,~] = simRayProp(updateM,1,0);
    magCoeff(s) = y2/1;
    theoreticalmag(s) = 1/(1-(z1Range(s)/focalLength));
end

figure('Name','m vs. z1');
plot(z1Range,magCoeff);
hold on;
plot(z1Range,theoreticalmag);
hold off;
legend('Real','Theoretical');

%}

function plotRayDiagram(rayIn, thetaIn, m, focalLength, z1, z2)
%generate a plot for the ray Diagram 
    lensPos = (z1+z2)/2;
    [rayOut, thetaOut] = simRayProp(m, rayIn, thetaIn);
    lensHeight = [rayIn(1),-rayIn(1)];
    % fill in graphing options
    plot();
    plot();
end

function z2 = solveZ2(z1, f)
%solveZ2 - given a ray in, simulate 3 rays in, solve for z2
%
% Syntax: z2 = solveZ2(z1, rayIn, mz1, mf)
%
    z2 = z1/((z1/f)-1);
end
