>>simulate a density map at 4A
$IMPFAST python ../../bin/sample_density_map.py data/helix27.pdb 4 data/helix27_4.mrc

>>fit a GMM to that target density map
$IMPFAST python ~/imp/modules/isd_emxl/bin/create_gmm.py data/helix27_4.mrc 6 target_gmm/helix27_c6.txt -m target_gmm/helix27_c6.mrc

>>run the modeling script
$IMPFAST python em_rigd.py
