Data=[5,1,4,1,9;3,7,8,4,1;8,3,6,6,2;6,10,9,5,9];
[W2,pc2,latent,tsquared,explained,mu] = pca(Data);

pcclusters2 = clusterdata(pc2(:,1:2),'maxclust',2,'linkage','av'); 
labels2 = cellstr( num2str(([1:4]')) ); 
gscatter(pc2(:,1),pc2(:,2),pcclusters2) 
text(pc2(:,1), pc2(:,2), labels2, ...
    'VerticalAlignment','bottom', ...                              
    'HorizontalAlignment','right') 
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Principal Component Scatter Plot Students');


