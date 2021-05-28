# foldingathome_adaptive
Scripts to run adaptive sampling on Folding@home

To run an adaptive sampling project on Folding@home:
1. Set up your project in the usual way (note the systems MUST be the same for adaptive sampling, note tautomers have to be the same in particular!)
2. Additionally, in the project folder: create a `SEEDS` folder and copy in the `next_seed.py` script, changing the path in it. Add `python $home/next_state.py $results` to the `<next-gen-command>` section of the `project.xml` as the first line (before `mv -f $results/checkpointState.xml.bz2 $jobdir/state.xml.bz2`).
3. The above have created a project where if any new seeds are present in the `SEEDS` folder, each new gen will randomly draw one of them instead of continuing the trajectory. If `SEEDS` is empty, vanilla (i.e. continue to elongate the trajectories) mode of F@h operation automatically resumes (i.e. this framework is 'very' asynchronous, production can always continue now matter how fast and if seeds are being produced).
4. After some initial data has been collectd, continuously run one of the adaptive sampling scripts, using multiple threads (change number in the script, also change two paths in the main()). The folder in which the script is run must contain:

* system.xml (copied from one of the RUNs, and bz2-decompressed, remember all systems must be the same!)
* integrator.xml (as above)
* `SEEDS` folder where new seed state.xml's will be written and copied to the `SEEDS` folder in the F@h project folder
* `seeds_pdb` folder where pdbs of the new seeds will be written for visual inspection if desired
* `temp` folder which is used as temporary copying location for the features.npy due to needed sudo for file writing in the F@h data folder

The script will output (in addition to the seeds as described above) on every iteration:
* `current_dtrajs.npy` - discrete trajectories from k-means clustering of the tICA projection
* `current_trajs_glob.npy` - list of trajectory paths used in the last iteration, in the same order as the discrete trajectories
* `current_goal.txt` - if using exploit-including methods, the furthest along the sought direction of the reaction coordinate the simulation has progressed

Warnings and TODO: 
* This setup has only been verified to speed up a slow event (DFG flip in a kinase) using very exploit-heavy settings. The default FAST settings as set now have NOT been verified that way, this should be an urgent research item to complete, perhaps on the DFG flip or protein folding.
* The testing was also done using a small (50) number of GPUs, same number as number of seeds. Scaling up to the full F@h will require more testing, there is a potential problem that RPW is thinking about now with velocity randomization --- as written currently if there are more GPUs then seeds being produced, multiple trajectories with same seed and same starting velocities will be wastefully produced. Until this is resolved, use only in internal and BETA modes where only a few tens of GPUs will be available, and adjust the number of seeds accordingly.
