



magick montage Sigmap1.png Sigmap2.png Sigmap3.png Sigmap4.png Sigmap5.png Sigmap6.png Sigmap7.png Sigmap8.png Sigmap9.png Sigmap10.png Sigmap11.png Sigmap12.png Sigmap13.png -tile 3x5 -geometry +0+50 test.png

magick montage MNIcorr_* -tile 3x3 -geometry +0+0 FigS1.png

magick montage MouseCoordCorr_* -tile 3x3 -geometry +0+0 FigS9.png

magick montage MouseVsHumanGene_* -tile 3x3 -geometry +0+0 FigS13.png


magick montage mPC1* -tile 3x3 -geometry +0+0 all_mPC1.png
magick montage mPC2* -tile 3x3 -geometry +0+0 all_mPC2.png
magick montage mPC3* -tile 3x3 -geometry +0+0 all_mPC3.png



magick montage Sigmap* -tile 4x5 -geometry +0+50 PC2_neuromaps.png


magick montage Sigmap* -tile 3x4 -geometry +0+50 PC3_neuromaps.png

magick montage Sigmap* -tile 3x5 -geometry +0+50 PC1_neuromaps.png