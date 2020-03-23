# Selection Variables and Plotting

Hey there! Welcome to my (Davio's) selection.
This should work pretty simply and I think I commented everything quite nicely, but here's some added help:


1) fist thing you're gonna want to run is the MakeFinalVertexFiles-prime script for each of your samples to combine all of the separate variable files into one very convenient package. This is also where I define and store all of the reco variables we use, so take a look in here if you want to play with variables.
2) next, you'll want to use my PickleFiles notebook to take all of the constituent parts used for preselection (which should be outlined fairly clearly in the notebook) and line everything up in convenient dataframes!
3) The BDTTraining notebooks load these pickles up and can train BDTs on whatever variables you'd like! It'll output the weights for you.
4) The Selection notebooks load these pickles and your BDT weights and can perform selection for you!

# FNAL Workflow

Documentation by Taritree.

The DLMerged file (from wire-cell chain) is the assumed input.
For each file, we need to run:
* Shower Reco
* MPID
* Selection variable maker


## Example running the command

Right now, this is just for Taritree to remember how to call the test command in the `test_dir` folder (not in this repo).


### non-lar executable

```
python ../bin/run_ubdlana_lar.py -d merged_dlreco_3e2d592c-e729-4ad2-9844-e517dd0a90b6.root -c ../CalibrationMaps_MCC9.root -m out/multipid_out_0_WC.root -s out/out_showerreco.json -t Overlay -o out/ -g test --ismc -pmt -oh opHitCalib
```

### `lar.py` call

```
lar.py -c ubdlana.fcl -s test_dir/merged_dlreco_3e2d592c-e729-4ad2-9844-e517dd0a90b6.root`
```