#! /bin/bash
# recipe conditions
#
#
LAMBDA=0.07
START=0
#
##
## run 2756
#AN="IBAnalyzerTrackCount"
AN="IBAnalyzerEM"
#AN="IBAnalyzerPoca"
RUN=2756
FILE=/mnt/mutom-gluster/data/ispra/lnl/pattrec/Radmufit_r${RUN}_0_100000000ev_PR_1Mev.root
MIN=280
VOX=1.5
#
./build/bin/misureLNL ${FILE}  /mnt/mutom-gluster/data/ispra/lnl/vtk/r${RUN}_rototranslat_img_vox${VOX}_${MIN}min  ${MIN} ${START} ${VOX} --detector.transform=true --detector.findMatrix=true --detector.file1="LNLfiles/run2756_ave5.txt" --detector.file2="LNLfiles/run2755_ave5.txt" --iterations='1 1000 5000'
#
MIN=2520
./build/bin/misureLNL ${FILE}  /mnt/mutom-gluster/data/ispra/lnl/vtk/r${RUN}_rototranslat_img_vox${VOX}_${MIN}min  ${MIN} ${START} ${VOX} --detector.transform=true --detector.findMatrix=true --detector.file1="LNLfiles/run2756_ave5.txt" --detector.file2="LNLfiles/run2755_ave5.txt" --iterations='1 1000 5000'
#
### test transformation goodness
#./build/bin/misureLNL ${FILE}  test_r${RUN}_img_vox${VOX}_${MIN}min  ${MIN} ${START} ${VOX}  --analyzers=${AN} --detector.dump="true" --detector.transform=true --detector.findMatrix=true --detector.voxelReferenceSystem=true --detector.file1="LNLfiles/r2756_block_vox_vertex.txt" --detector.file2="LNLfiles/r2755_block_vox_vertex.txt" --iterations='1 50 50'
#./build/bin/misureLNL ${FILE}  test_r${RUN}_img_vox${VOX}_${MIN}min  ${MIN} ${START} ${VOX}  --analyzers=${AN} --detector.dump="true" --detector.transform=true --detector.findMatrix=true --detector.voxelReferenceSystem=false --detector.file1="LNLfiles/run2756_ave5.txt" --detector.file2="LNLfiles/run2755_ave5.txt" --iterations='1 50 50'
#./build/bin/misureLNL ${FILE}  test_r${RUN}_img_vox${VOX}_${MIN}min  ${MIN} ${START} ${VOX}  --analyzers=${AN} --detector.dump="true" --detector.transform=true --detector.findMatrix=true --detector.voxelReferenceSystem=false --detector.file1="LNLfiles/run2756_ave7.txt" --detector.file2="LNLfiles/run2755_ave7.txt" --iterations='1 50 50'
#./build/bin/misureLNL ${FILE}  test_r${RUN}_img_vox${VOX}_${MIN}min  ${MIN} ${START} ${VOX}  --analyzers=${AN} --detector.dump="true" --detector.transform=true --detector.findMatrix=true --detector.voxelReferenceSystem=false --detector.file1="LNLfiles/run2756_full5.txt" --detector.file2="LNLfiles/run2755_full5.txt" --iterations='1 50 50'
#./build/bin/misureLNL ${FILE}  test_r${RUN}_img_vox${VOX}_${MIN}min  ${MIN} ${START} ${VOX}  --analyzers=${AN} --detector.dump="true" --detector.transform=true --detector.findMatrix=true --detector.voxelReferenceSystem=false --detector.file1="LNLfiles/run2756_full7.txt" --detector.file2="LNLfiles/run2755_full7.txt" --iterations='1 50 50'
##
### run 2755
#RUN=2755
#FILE=/mnt/HDD/ext/data/experiments/radmu/ispra/lnl/pattrec/Radmufit_r2755_0_100000000ev_PR_16-20Mev.root
#./build/bin/misureLNL ${FILE}  test_r${RUN}_img_vox${VOX}_${MIN}min  ${MIN} ${START} ${VOX} --detector.transform=false --detector.findMatrix=true --detector.voxelReferenceSystem=false --detector.file1="LNLfiles/run2756_full5.txt" --detector.file2="LNLfiles/run2755_full5.txt" --iterations='1 50 0'
##



