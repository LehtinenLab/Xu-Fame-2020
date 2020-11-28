%% THIS IS THE ONLY SECTION YOU SHOULD HAVE TO MAKE CHANGES IN
clear all;
%put copy the folder path that contains the image you want to analyze
fdir = 'Z:\Fred\Ryann\Mitos\9928';
%copy in the base name. 
%E.g. for file 'Ryann 1144 P0_1003_mito.roi.zip', 
%ID = 'Ryann 1144 P0_1003'.
ID = '9928019';

%set to true if you want to save image with borders as PDF
save_output = true;


%% set file paths and load the image and ROIs
impath = dir(strcat(fdir,filesep,ID,'*.tif'));
im = imread(strcat(fdir,filesep,impath.name));

ROIpath = dir(strcat(fdir,filesep,ID,'*.zip'));
ROIs = ReadImageJROI(strcat(fdir,filesep,ROIpath.name));
%% draw borders!
%checks to see if you've done this before and if so, load old drawings
try
    load(strcat(fdir,filesep,ID,'_variables.mat'));
catch
    %basal first
    F = figure;
    imshow(im,[]);
    %draw the basal border
    b = drawpolyline('Color','b');
    %make a mask of the basal border
    basal_edge = createMask(b);

    %apical next
    %draw the basal border
    a = drawpolyline('Color', 'r');
    %draw the apical border
    apical_edge = createMask(a);

    if save_output
        print(strcat(fdir,filesep,ID,'_apical_basal_border.pdf'),'-dpdf','-fillpage');

    end

    close(F);
end
%% extract the masks and distance transforms
% if ~exist('basal_edge')
%     load(strcat(fdir,filesep,ID,'_variables.mat'));
% end

w = waitbar(0,'extracting masks....');
for i = 1:numel(ROIs)
    %get ROI coordinates
    C = ROIs{1,i}.mnCoordinates;
    %make binary mask from the coordinates
    mask = poly2mask(C(:,1),C(:,2),size(im,1),size(im,2));
    %do distance transformation of binary mask
    D = bwdist(mask);
    %get centroid of mask
    stats = regionprops(mask);
    mask_centroid(i,:) = stats.Centroid;

    %distances along basal edge
    D_b = D(basal_edge);
    %find minimum basal dist
    [Basal_D(i), Basal_I(i)] = min(D_b);
    
    %distances along apical edge
    D_a = D(apical_edge);
    %find minimum apical dist
    [Apical_D(i), Apical_I(i)] = min(D_a);
    waitbar(i/numel(ROIs));
end
close(w);
%% get min XY locations
[XX, YY] = meshgrid(1:size(im,2),1:size(im,1));

BasalLoc = [XX(basal_edge), YY(basal_edge)];
ApicalLoc = [XX(apical_edge), YY(apical_edge)];

BasalMinLoc = BasalLoc(Basal_I,:);
ApicalMinLoc = ApicalLoc(Apical_I,:);    
%% plot distruibution

dist_fig = figure;
set(gcf, 'Position',  [500, 500, 1500, 400])

%plot the shortest distance to the basal membrane for each mito
subplot(1,3,1);
scatter(1:numel(Basal_D),sort(Basal_D));
xlabel('mito #');
ylabel('distance (px)');
title('minimum distance to basal edge');
axis square;

%plot the shortest distance to the apical membrane for each mito
subplot(1,3,2);
scatter(1:numel(Apical_D),sort(Apical_D));
xlabel('mito #');
ylabel('distance (px)');
title('minimum distance to apical edge');
axis square;

%plot the ratio of shortest basal to shortest apical distance
%n.b. 0 = on the basal membrane, 1 = on the apical membrane
AB_ratio = Basal_D ./ (Basal_D + Apical_D);
subplot(1,3,3);
scatter(1:numel(AB_ratio),sort(AB_ratio));
xlabel('mito #');
ylabel('relative distance');
title('distribution of basal to apical (B = 0; A = 1)');
axis square;

%save distribution figure
if save_output
    print(strcat(fdir,filesep,ID,'_distribution.pdf'),'-dpdf','-fillpage');
end

%% save excel sheet with all ofthe basal/apical distances, mito locations, and distances
ROI_index = 1:numel(ROIs);
DistanceTable = table(ROI_index',mask_centroid, Basal_D', BasalMinLoc, Apical_D', ApicalMinLoc, AB_ratio', 'VariableNames',{'ROIidx','MitoCentroid','BasalDistance','BasalMinLoc','ApicalDistance', 'ApicalMinLoc','ABRatio'});
if save_output
    saveexcel = strcat(fdir,filesep,ID,'_distances.xlsx');
    writetable(DistanceTable,saveexcel);
end
%% save all of the variables/locations/borders for later
if save_output
    savevar = strcat(fdir,filesep,ID,'_variables.mat');
    save(savevar,'basal_edge','apical_edge','DistanceTable','im','ID');
end

%% plot the image and the nearest distances

map_fig = figure;
imshow(im,[]);
hold on;

basal_bound = bwboundaries(basal_edge);
p1 = plot(basal_bound{1}(:,2),basal_bound{1}(:,1),'b','LineWidth',7.5);
p1.Color(4) = 0.5;

apical_bound = bwboundaries(apical_edge);
p2 = plot(apical_bound{1}(:,2),apical_bound{1}(:,1),'r','LineWidth',7.5);
p2.Color(4) = 0.5;

for idx = 1:numel(ROIs)
    C = ROIs{1,idx}.mnCoordinates;
    %make binary mask from the coordinates
    mask = poly2mask(C(:,1),C(:,2),size(im,1),size(im,2));
    mask_boundaries = bwboundaries(mask);
    p1 = plot([mask_centroid(idx,1),BasalMinLoc(idx,1)],[mask_centroid(idx,2),BasalMinLoc(idx,2)],'Color',[0.5,0.5,1],'LineWidth',2);
    p1.Color(4) = 0.5;
    p2 = plot([mask_centroid(idx,1),ApicalMinLoc(idx,1)],[mask_centroid(idx,2),ApicalMinLoc(idx,2)],'Color',[1,0.5,0.5],'LineWidth',2);
    p2.Color(4) = 0.5;
    for j = 1:length(mask_boundaries)
        boundary = mask_boundaries{j};
        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
    end
end
hold off;

if save_output
    print(strcat(fdir,filesep,ID,'_map.pdf'),'-dpdf','-fillpage');
end

%%

