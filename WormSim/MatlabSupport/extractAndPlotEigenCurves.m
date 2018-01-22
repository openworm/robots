%extract the eigenvectors
%This assumes that the CSV file you are using has curves made up of angles
%in it.  If your curves are parameterized x(s), y(s) curves, this file will
%not work.  Enter the file name in the csvread's field
[vecs,vals] = familyEigenVectorsFromAngles(csvread('../foodTracking.txt'));

%plot the first four in a pretty way
%plot the first four eigenworms
figure(1)
subplot(2,2,1)
plot(linspace(0,1,length(vecs(:,end))),vecs(:,end));
xlabel('s');
ylabel('\theta');
title('Eigenworm 1');
subplot(2,2,2)
plot(linspace(0,1,length(vecs(:,end))),vecs(:,end-1));
xlabel('s');
ylabel('\theta');
title('Eigenworm 2');
subplot(2,2,3)
plot(linspace(0,1,length(vecs(:,end))),vecs(:,end-2));
xlabel('s');
ylabel('\theta');
title('Eigenworm 3');
subplot(2,2,4)
plot(linspace(0,1,length(vecs(:,end))),vecs(:,end-3));
xlabel('s');
ylabel('\theta');
title('Eigenworm 4');
