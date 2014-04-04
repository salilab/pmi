>>simulate a density map at 4A
$IMPFAST python ../../bin/sample_density_map.py data/helix27.pdb 4 data/helix27_4.mrc

$IMPFAST python ../../bin/sample_density_map.py data/helix3.pdb 4 data/helix3_4.mrc


>>fit a GMM to that target density map
$IMPFAST python ~/imp/modules/isd_emxl/bin/create_gmm.py data/helix27_4.mrc 6 target_gmm/helix27_c6.txt -m target_gmm/helix27_c6.mrc

$IMPFAST python ~/imp/modules/isd_emxl/bin/create_gmm.py data/helix3_4.mrc 20 target_gmm/helix3_c20.txt -m target_gmm/helix3_c20.mrc

$IMPFAST python ~/imp/modules/isd_emxl/bin/create_gmm.py data/helix3_4.mrc 10 target_gmm/helix3_c10.txt -m target_gmm/helix3_c10.mrc

$IMPFAST python ~/imp/modules/isd_emxl/bin/create_gmm.py data/helix3_4.mrc 5 target_gmm/helix3_c5.txt -m target_gmm/helix3_c5.mrc

$IMPFAST python ~/imp/modules/isd_emxl/bin/create_gmm.py data/helix3_4.mrc 3 target_gmm/helix3_c3.txt -m target_gmm/helix3_c3.mrc

>>run the modeling script
$IMPFAST python em_rigid.py
