<div align="center">

<a href="https://dranirbanchakraborty.com/learn/error-correction-coding"><picture><source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/anirbanchakraborty-dev/anirbanchakraborty-dev/output/repo-ecc-book-lab-banner-dark.svg"><img alt="ecc-book-lab, the companion code for Error Correction Coding. Beside the description, seven data bits gain an even-parity bit, one bit flips in the channel, and the XOR of all eight received bits comes out 1, so the checker reports an error." src="https://raw.githubusercontent.com/anirbanchakraborty-dev/anirbanchakraborty-dev/output/repo-ecc-book-lab-banner-light.svg" width="100%"></picture></a>

[The book](https://dranirbanchakraborty.com/learn/error-correction-coding) · [Running it](#running-it) · [Layout](#layout) · [The parity checker](hdl/parity)

</div>

Companion repository for *Error Correction Coding*, the online book at [dranirbanchakraborty.com/learn/error-correction-coding](https://dranirbanchakraborty.com/learn/error-correction-coding). Every number the book prints and every chart it draws is computed by a program in this repository, and every hardware listing the book supplies lives here in finished form, with the testbench that exercises it.

The book does not teach SystemVerilog or cocotb, and it teaches only as much Python as reading the scripts needs. It hands over each listing complete and explains what the hardware does. For the languages themselves, the *RTL Design and Verification* book, in preparation, and the HDL material in *Computer Architecture* on the same site start from nothing.

Each kind of file feeds the book its own way.

```mermaid
flowchart LR
  numbers["python/ scripts of numbers"] -->|"checked against the chapter"| text["every worked number in the book"]
  charts["python/ scripts of charts"] -->|"written as SVG"| figures["every chart in the book"]
  hdl["hdl/ SystemVerilog"] -->|"driven by cocotb under Verilator"| tb["tb/ testbenches"]
  hdl -->|"synthesised by Yosys"| synth["synth/ gate count and depth"]
```

## Running it

The chapter scripts need Python 3.10 or newer and nothing outside the standard library.

```bash
python3 python/ch01_why_codes_exist.py
```

That reproduces every worked number in Why Codes Exist, the first chapter of Part I, prints each one beside the value the chapter states, and stops with an error the moment any of them disagrees. A run that ends with `All chapter 1 numbers reproduced.` is a passing run. Script numbers count from that chapter, so they differ from the chapter numbers the website shows, and the two front-matter chapters before it have scripts named `chf1` and `chf2`.

```bash
python3 python/ch01_charts.py --out figures
```

That writes the chapter's two charts as SVG files into `figures/`. The book's copies are made the same way, so a chart is corrected by fixing the script and running it again, never by editing the SVG.

The hardware chapters also need Verilator 5.036 or newer and cocotb 2.1, and Yosys for the gate counts. The chapter "Lab: A First Hardware Checker" gives the install steps for macOS, Linux and Windows. From the repository root:

```bash
python tb/parity/run_parity.py
python synth/synth_parity.py
```

The first builds the parity bench three times under Verilator and runs its cocotb tests, and the second synthesizes the parity generator and checker with Yosys and checks their gate counts and depths against the formulas the chapter proves. Simulation output, the waveform included, lands in `tb/parity/sim_build/`, which git ignores.

## Layout

The folders arrive as the chapters need them.

- `python/` holds one script of numbers and, where a chapter has charts, one script of charts per chapter, named after the chapter.
- `hdl/` holds the SystemVerilog, one folder per design. `hdl/parity/` is the parity generator and checker of "Lab: A First Hardware Checker", built from one parameterized XOR tree, with the bench that joins them through a simulated channel.
- `tb/` holds the cocotb testbenches that drive those designs under Verilator, one folder per design, each with a runner script that builds and runs it. `tb/parity/` also holds a Surfer command file that opens the waveform with the bench's signals in view.
- `synth/` holds the open-source synthesis scripts that report each design's gate count and depth. `synth/synth_parity.py` is the one for "Lab: A First Hardware Checker".
- `formal/` will hold the properties that prove a design correct for every input rather than for the inputs a test happened to try.
- `data/` will hold the measured data sets a chapter plots, each with a note on where it came from.

## License

MIT. See [LICENSE](LICENSE).
