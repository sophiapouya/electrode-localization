clear; clc; close all

[electrode_file,electrode_dir]=uigetfile('.fcsv','Select the file with the electrode trajectory (.fcsv)'); % select the pointlist file created in 3D slicer
merShift=input('What are the MER displacements [ML,AP,SI]? '); % enter a matrix of values (mm) indicating medial-lateral,anterir-posterior, and superior-inferior (depth) position of the microelectode relatve to the DBS electrode
[x,y,z] = importDBSlead([electrode_dir,electrode_file]); % import the coodrinates

minno=find(z==min(z)); % find the deepest (minimum) point on the electrode
maxo=find(z==max(z));  % find the highest (maximum) point on the electrode
coord=[-x,-y,z]; % invert the x and y coordiates (3D slicer exports as LPS space but we want RAS space)
[coeff,score,roots] = pca(coord); % perform principal component analysis on the electrode coordiantes, the first component will define the trajectory of the DBS electrode


[n,p] = size(coord); % get the number of points in the electrode space
meanX = mean(coord,1); % calculate the center point of the electrode array
Xfit1 = repmat(meanX,n,1) + score(:,1)*coeff(:,1)'; % project the [x y z] handmarked electrodes onto the new axis of the DBS electrode calculated from PCA
origin=Xfit1(minno,:); % the origin is the tip of the electrode
endpt=Xfit1(maxo,:)'; % the last point marked on the electrode as it exits the brain
dirVect=endpt-origin'; % the direction vector is the vector pointing from the tip to the last point on the elctrode
%shifted = [meanX + (min(score(:,1)))*dirVect'];

plot3([origin(1),endpt(1)],[origin(2),endpt(2)],[origin(3),endpt(3)],'k-'); hold on; % plot the axis of the DBS electrode
X1 = [coord(:,1) Xfit1(:,1) nan*ones(n,1)];
X2 = [coord(:,2) Xfit1(:,2) nan*ones(n,1)];
X3 = [coord(:,3) Xfit1(:,3) nan*ones(n,1)];
hold on
plot3(X1',X2',X3','b-', coord(:,1),coord(:,2),coord(:,3),'bo'); % plot the marked electrode points with line projecting them onto the DBS axis

% newVec=origin+dirVect';  
% plot3([origin(1),newVec(1)],[origin(2),newVec(2)],[origin(3),newVec(3)],'-m'); hold on;
xlabel('x');
ylabel('y');
zlabel('z');


% calculate a basis that is orthogonal to the electrode vector, this will
% be used for calculating medial-lateral and anterior-posterio translations
normV=dirVect./norm(dirVect); 
x=1; y=0; 
newZ= (-(normV(1)*(x)+normV(2)*(y)))/(normV(3));
ax2=[x y newZ];
ax2=ax2/norm(ax2);
newBasis=[ax2',cross(ax2,normV)', normV];
if newBasis(1,1)<0;
    newBasis(:,1)=-newBasis(:,1);
end
if newBasis(2,2)<0;
    newBasis(:,2)=-newBasis(:,2);
end
newOrigin=origin;
newCoord(:,1)=newBasis*[20 0 0]'+origin';
newCoord(:,2)=newBasis*[0 20 0]'+origin';
newCoord(:,3)=newBasis*[0 0 20]'+origin';

% plot the vectors that form the basis
cols=[1 0 0; 0 1 0; 0 0 1];
for i=1:3;
    plot3([origin(1),newCoord(1,i)],[origin(2),newCoord(2,i)],[origin(3),newCoord(3,i)],'-or','Color',cols(i,:));
    hold on;
end

% calculate the position of the microelectrode in this basis
newRel=newBasis*(merShift')+origin'; newRel=newRel';

% plot the position of the microelectrode tip
plot3([origin(1),newRel(1,1)], [origin(2),newRel(1,2)], [origin(3),newRel(1,3)],':ok','MarkerSize',5','MarkerFaceColor','k');
xlabel('x');
ylabel('y');
zlabel('z');

% round the micro-electrode location to the neareast 0.001 mmm
newRel=round(newRel*1000)/1000;
title(['Adjusted MER locations: [',num2str(newRel(1)),',',num2str(newRel(2)),',',num2str(newRel(3)),']']);
newRel

