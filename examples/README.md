# Some Examples

## Quick Example

A really quick example that will introduce you to some of the AraProc features
can be found in `run_example.py`. Run it like:
```
python run_example.py
```
It should generate a panel of 16 waveforms and a skymap.

## Process Full Run

The script that demonstrates how you might use the framework in mass processing
is in `process_run.py`. You can run it like:

```
python process_run.py \
    --input_file /data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/event005626.root \
    --is_simulation 0 \
    --station 5 \
    --output_file station_5_run_005626.root
```

It takes as arguments information about the run, pedestal files, etc.
It also generates a ROOT output file containing the computed variables.
