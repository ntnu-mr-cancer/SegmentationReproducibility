function im3d_bfc = n4BiasFieldCorrection(im3d)

% make tmporary folder on scratch
mkdir('/home/BFC');

% convert images to StrDatax and write to meta IO file
addpath('.../ElastixFromMatlab') % add elxelastic tool
StrDatax = elxIm3dToStrDatax(im3d);
elxStrDataxToMetaIOFile(StrDatax,'/home/BFC/tmpIn.mhd',0);

% do bias field correction from command line
system('/home/c3d/bin/c3d /home/BFC/tmpIn.mhd -biascorr -o /home/BFC/tmpOut.mhd');

% read in corrected image
StrDatax_bfc = elxMetaIOFileToStrDatax('/home/BFC/tmpOut.mhd',0);
im3d_bfc = elxStrDataxToIm3d(StrDatax_bfc);

% clean up
rmdir('/home/BFC','s');

end
