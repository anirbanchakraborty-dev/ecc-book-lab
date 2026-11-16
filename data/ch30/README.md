# Gate counts for "Project: A SEC-DED Encoder and Decoder"

`secded72_gates.json` holds the gate counts and logic depths that the chapter "Project: A SEC-DED Encoder and Decoder" of *Error Correction Coding* quotes, charts and compares with the hand counts of "SECDED and Hsiao Codes". Every number in it was measured by `synth/synth_secded72.py` from the SystemVerilog in `hdl/secded72/`, with this Yosys build:

```
Yosys 0.71+post (git sha1 5628c13240f670b97aaaca3d47940c4afb8ff572, Release, AppleClang clang++ 21.0.0.21000401)
```

The same string is stored under the key `yosys`, so a reader can tell which release produced the file.

## How it was measured

Yosys reads the design and hands the logic to ABC through `abc -g AND,OR,XOR`, which maps onto two-input AND, OR and XOR gates plus the NOT gate Yosys always adds, with the script `strash; balance -x; map`. In the words of ABC's own help, the first step transforms the logic into an AIG, a graph of two-input ANDs and inverters, the second makes it well balanced, minimizing levels and balancing multi-input XORs, and the third performs standard-cell mapping, with area recovery on, onto the library Yosys writes, in which every cell, the inverter included, has the same delay. Depth is the number of two-input gates on the longest path to an output, read from the netlist Yosys writes. An inverter adds no level here and is not counted as a gate, as in the hand count, so the chapter quotes `two_input`. ABC writes an XNOR as an XOR followed by an inverter, so XNOR never appears in the counts.

Each design is measured in two ways. As **separate trees**, each syndrome bit and each check bit of the encoder is synthesized on its own, with every other output removed, so no gate can serve two outputs. That is the model the hand count uses, and the script checks that each tree has exactly $W - 1$ gates and $\lceil \log_2 W \rceil$ levels for a row of $W$ ones. As a **whole block**, the syndrome stage, the encoder and the full decoder are each synthesized in one piece, and ABC may share any gate between outputs, as it would in a real netlist. Whole-block counts come from a heuristic, so they are measurements of what this tool chain built, not proven minimums.

## The fields

The file holds one record per code under `codes`, for `hsiao` and `hamming` (the two (72,64) matrices of "SECDED and Hsiao Codes"), `hsiao_dictionary` (the Hsiao columns taken in dictionary order) and `hsiao_39_32` (the (39,32) Hsiao code of an exercise).

- `n`, `k` and `r` are the code length, the data bits and the check bits.
- `row_weights` lists the ones in each row of the check matrix, row 0 first, and `check_fan_in` the data bits in each check bit's equation.
- `syndrome_trees` and `encoder_trees` list each separate tree as `{bit, gates, xor, depth}`.
- `syndrome_whole` and `encoder_whole` give the whole-block result: `cells` by type (`and`, `or`, `xor`, `not`), `total` cells, `two_input` cells (the total less the inverters) and `depth`, which is a list with one entry per syndrome bit for the syndrome stage and a single number for the encoder.
- `decoder_whole_flag0` and `decoder_whole_flag1` give the whole decoder with the no-match flag and with the parity flag: `cells`, `total` and `two_input` as above, and the depth to the data outputs (`data_depth`), to `corrected` (`corrected_depth`) and to `uncorrectable` (`uncorrectable_depth`).

## Checking and regenerating

From the repository root:

```bash
python3 synth/synth_secded72.py
python3 synth/synth_secded72.py --write
```

The first synthesizes every code again and compares each measurement with this file. It exits 0 when everything matches, 1 when anything differs and 2 when no `yosys` is on the `PATH`. The separate-tree counts follow from the matrices alone and should match on any release, but a different Yosys or ABC release can change the whole-block counts, and then the second command rewrites this file from the new run. The chapter's companion script, `python/ch30_project_secded_encoder_and_decoder.py`, reads the file and checks each number the chapter quotes, so a rewritten file that moves a quoted number makes that script fail until the text is brought up to date.

## The runs beside it

The same chapter times the simulation bench and the formal proofs. On the recorded run, `tb/secded72/run_secded72.py` swept 2,802,514 cases for each (72,64) code, in 44.1 seconds for Hsiao and 45.2 seconds for extended Hamming, 95.45 seconds of wall time in all, and `formal/secded72/run_formal72.py` finished all nine SymbiYosys tasks in 489.6 seconds, with the Hsiao proof taking 18.1 seconds and the extended Hamming proof 189.4 seconds. Those times are not stored in this file, because they depend on the machine.
