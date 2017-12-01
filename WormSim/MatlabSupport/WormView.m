% Program to visualise the output of WormSim 

%Import worm data
data = importdata('../Model/simdata.csv', ',');

%Check if objects.csv exists and if so, import object data
%NOTE: if not using objects, you must delete any objects.csv file generated
%during previous runs.
usingObjects = exist('../Model/objects.csv','file');
if usingObjects 
    Objects = importdata('../Model/objects.csv');
    Sz = size(Objects);
    Nobj = Sz(1);
end

Sz = size(data);
Nt = round(Sz(1));              
Dt = (data(2,1) - data(1,1))
Duration = data(end,1);

Nbar = (Sz(2)-1)/3;
NSEG = Nbar-1;
D = 80e-6;

CoM = zeros(Nt,Nbar,3);
CoMplot = zeros(Nt,2);
Dorsal = zeros(Nbar,2);
Ventral = zeros(Nbar,2);
Midline = zeros(Nbar,2);
act_D = zeros(Nt,Nbar-1);
act_V = zeros(Nt,Nbar-1);
L_D = zeros(Nt,Nbar-1);
L_V = zeros(Nt,Nbar-1);
X = zeros(Nt,25);
Y = zeros(Nt,25);

%NOTE: Create figure then check box to fill axes!
figure('Position',[1 1 824 588])
XYratio = 1.3333;   %This is the ratio of Xrange/Yrange required to give equal axes

R = D/2.0*abs(sin(acos(((0:Nbar-1)-NSEG./2.0)./(NSEG/2.0 + 0.2))));

for i = 1:Nt    
   frame = i;
    for j = 1:Nbar        
        CoM(i,j,1) = data(frame,1 + (j-1)*3 + 1);
        CoM(i,j,2) = data(frame,1 + (j-1)*3 + 2); 
        CoM(i,j,3) = data(frame,1 + (j-1)*3 + 3);
    end           
end

% Do visualization

%Would you like to output frames for making a movie? (1=yes, 0=no)
MakeMovie = 0;

simMaxX = max(max(CoM(:,:,1))) + 0.1e-3;
simMinX = min(min(CoM(:,:,1))) - 0.1e-3;
simMaxY = max(max(CoM(:,:,2)));
simMinY = min(min(CoM(:,:,2)));
simXrange = simMaxX-simMinX;
simYrange = simMaxY-simMinY;

if simXrange >= simYrange*XYratio
    plotMaxX = simMaxX;
    plotMinX = simMinX;
    Yerr = simXrange/XYratio - simYrange;
    plotMaxY = simMaxY + Yerr/2;
    plotMinY = simMinY - Yerr/2;
else
    plotMaxY = simMaxY;
    plotMinY = simMinY;
    Xerr = simYrange*XYratio - simXrange;
    plotMaxX = simMaxX + Xerr/2;
    plotMinX = simMinX - Xerr/2;
end

for i = 1:Nt
    for j = 1:Nbar
        dX = R(j)*cos(CoM(i,j,3));
        dY = R(j)*sin(CoM(i,j,3));
        Dorsal(j,1) = CoM(i,j,1) + dX;
        Dorsal(j,2) = CoM(i,j,2) + dY;    
        Ventral(j,1) = CoM(i,j,1) - dX;    
        Ventral(j,2) = CoM(i,j,2) - dY; 
    end
    
    %Plot objects if, if they are present
    if usingObjects
        angles = (0:0.05:1).*(2*pi);
        outline = zeros(length(angles),2);
        for j = 1:Nobj
            outline(:,1) = Objects(j,1) + cos(angles).*Objects(j,3);
            outline(:,2) = Objects(j,2) + sin(angles).*Objects(j,3);
            plot(outline(:,1),outline(:,2),'b','linewidth',2)
            hold on
        end  
    end
    
    %Plot worm    
    plot(Dorsal(:,1),Dorsal(:,2),'k','linewidth',4)
    hold on
    plot(Ventral(:,1),Ventral(:,2),'k','linewidth',4)
    plot([Ventral(1,1) Dorsal(1,1)],[Ventral(1,2) Dorsal(1,2)],'k','linewidth',4)
    plot([Ventral(end,1) Dorsal(end,1)],[Ventral(end,2) Dorsal(end,2)],'k','linewidth',4)     
    set(gca,'xticklabel','','yticklabel','','ytick',[],'xtick',[])    
    xlim([plotMinX plotMaxX])
    ylim([plotMinY plotMaxY])    
    set(gcf,'paperpositionmode','auto')   
    DataTitle = ['Time: ',num2str(i*Dt),'s of: ',num2str(Nt*Dt),'s'];
    title(DataTitle);
    pause(Dt)
    hold off
    
    if MakeMovie
        %Frame output     
        if i < 10
            name = strcat('TempFrames/000',num2str(i),'.jpg');
        elseif i < 100
            name = strcat('TempFrames/00',num2str(i),'.jpg');
        elseif i < 1000
            name = strcat('TempFrames/0',num2str(i),'.jpg');
        else
            name = strcat('TempFrames/',num2str(i),'.jpg');
        end    
        export_fig(name,'-jpg')
    end    
end





