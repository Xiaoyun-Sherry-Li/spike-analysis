all_dir = 'Z:\Sherry\ephys_acquisition\SLV123\7turn_240731_133751';
n_clusters = 3;
[gm, idx,labels] = GMM_make(all_dir, true, n_clusters)