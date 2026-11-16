# Measured miscorrection counts for "Lab: Measuring Miscorrection"

`miscorrection_counts.json` holds the outcome counts that the chapter "Lab: Measuring Miscorrection" of *Error Correction Coding* quotes and charts. Every count in it was measured on the SystemVerilog decoder in `hdl/secded72/`, simulated by Verilator under cocotb, by `tb/miscorrection/run_miscorrection.py` with its default settings, seed 31 and 100,000 error patterns of weight 3 and 100,000 of weight 4 for each code. The recorded run used these tools, stored under the key `tools`:

```
Python 3.14.9
cocotb 2.1.0
Verilator 5.054 2026-11-01 rev vUNKNOWN-built20261101
```

## How it was measured

The runner builds `hdl/secded72/secded72_bench.sv` as it stands, once for each of the three (72,64) matrices of "SECDED and Hsiao Codes" (the balanced Hsiao matrix, the Hsiao matrix with its weight-5 columns taken in dictionary order, and extended Hamming), with the uncorrectable flag in its `nomatch` form. The test in `tb/miscorrection/test_miscorrection.py` draws every pattern with `draw()` from `python/ch31_lab_measuring_miscorrection.py`, through `random.Random(31)` and its `random()` method alone: two draws make the 64-bit data word, and each further draw, scaled to 72 and truncated, names one codeword position, drawn again whenever it repeats one already chosen. All weight-3 patterns come first and then all weight-4 patterns, from the one generator, so the three codes see the same patterns on the same data words.

For each case the test writes the data word and the error pattern into the bench, reads `data_out`, `corrected` and `uncorrectable`, and checks that they equal what `tb/secded72/secded72_model.py` decodes for the same received word, and that the wrong data bits and the two flags equal what the error pattern alone predicts. It then sorts the case into one of four outcomes from the outputs alone: flagged when `uncorrectable` is 1, correct when `data_out` equals the data word, miscorrected when `corrected` is 1 and the data are still wrong, and undetected when the data are wrong with neither flag set. Last, it writes a data word drawn from `random.Random(32)` under each of the 70 windows of three adjacent cells in three layouts of the word, the natural order, stride 11 and the order found by the layout search of the chapter's script. No search runs for extended Hamming, since no layout takes it below 46 miscorrected windows, so its searched layout repeats stride 11 and it sees 140 distinct windows rather than 210.

The runner repeats the campaign in software from the same seed and exits non-zero at the first count that differs, so a file written by it always agrees with the software model.

## The fields

- `seed`, `draws` and `weights` give the campaign's settings, `adjacent_seed` the seed of the data words under the adjacent windows, and `checkpoints` the draw counts at which running totals were taken.
- `codes` holds one record for each of `hsiao`, `hsiao_dictionary` and `hamming`.
- `bench_code` is the value of the bench's `CODE` parameter, and `columns` the 72 columns of the check matrix in `tb/secded72/secded72_model.py`, the model the bench's SystemVerilog package is generated from, as whole numbers whose bit $i$ is row $i$.
- `campaign` gives, for `weight_3` and `weight_4`, the number of patterns with each of the four outcomes.
- `running` gives, for each weight, the number of patterns miscorrected (weight 3) or undetected (weight 4) among the first 100, 200, 500 and so on up to all 100,000, which the chapter's convergence chart plots.
- `adjacent` gives the four outcome counts of the 70 adjacent windows in each of the three layouts.
- `cases` is the number of cases driven into the bench, and `model_mismatches` and `prediction_mismatches` the cases on which the hardware differed from the model or from the prediction, both 0.
- `paired` gives the weight-3 patterns miscorrected by the balanced Hsiao matrix alone and by the dictionary matrix alone, from which the chapter measures the difference between the two on the same patterns.

## Checking and regenerating

```bash
python3 python/ch31_lab_measuring_miscorrection.py
python tb/miscorrection/run_miscorrection.py
```

The first reads the file and checks every count against its own software campaign and every column against the matrices "SECDED and Hsiao Codes" builds. The second, run in the Python that holds cocotb, measures everything again and rewrites the file. The counts follow from the seed, the matrices and the decoder alone, so a correct decoder gives the same file on any machine, with only the tool versions changed. Run with any other seed, number of draws or weights, the runner prints its counts and intervals and leaves the file alone.

## The run beside it

Over two runs of the default campaign each code's 200,210 cases took 3.7 to 6.1 seconds, and the whole run 11.7 to 18.7 seconds of wall time, with the three builds already in place. A first build from an empty `sim_build/` folder took about 2.5 seconds for each code. Those times are not stored in the file, because they depend on the machine.
