# UBDLANA

This branch is a significant deviation from Davio's analysis scripts.
This includes the code used by the MicroBooNE DL working group for the Gen-1 Low energy electron search.
This repository contains the code used to run 'dlana' and 'dlfilter' jobs for production.
One can run the code using the `lar.py` framework used by the MicroBooNE production group.
You can also run it independly of the `lar.py` framework through the `bin/run_ubdlana_lar.py` script.
The dlfilter jobs use the `run_ubdl_multifilter_lar.py` script.
These scripts are used to run on MC files on the Tufts cluster.

## Versions

* v1_1_4: used in Gen-1 analysis. Ensemble BDT setup for 1m1p and 1e1p selection.
* v1_1_3: used in Gen-1 analysis. single BDT setup for 1m1p and 1e1p selection.

See READMEs in `bdt_models` and `showercnn_models` folders for information on the models used.


## UBDLANA

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

### Testing development environment

This is for setting up a situation where you are developing this code AND `dllee_unified`.
This might be the case if you are developing the shower reco code or any other module.

Steps:

* build `dllee_unified` on FNAL or setup your existing copy.  Refer to the README in the `dllee_unified` folder. 
  But, as a reminder, there are scripts to setup and build `dllee_unified`. Also, you want the `mcc9_prod1_29e_rc` branch.
  If you have access to CVMFS, you may also setup a copy of dllee_unified via UPS.
* go into `configure.sh` in this repo and comment out 

      setup dllee_unified develop -q e17:prof

  This line sets up `dllee_unified` using UPS. you don't want to do this, if you have changed that dependency.
* now run `configure.sh`

      source configure.sh


* Now you can use the `lar.py` call in the subsection above.