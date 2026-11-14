# ecc-book-lab

Companion repository for *Error Correction Coding*, the online book at [dranirbanchakraborty.com/learn/error-correction-coding](https://dranirbanchakraborty.com/learn/error-correction-coding). Every number the book prints and every chart it draws is computed by a program in this repository, and every hardware listing the book supplies lives here in finished form, with the testbench that exercises it.

The book does not teach SystemVerilog, cocotb or Python. It hands over each listing complete and explains what the hardware does. For the languages themselves, the *RTL Design and Verification* book and the HDL material in *Computer Architecture* on the same site start from nothing.

## Running it

The chapter scripts need Python 3.10 or newer and nothing outside the standard library.

```bash
python3 python/ch01_why_codes_exist.py
```

That reproduces every worked number in the first chapter, prints each one beside the value the chapter states, and stops with an error the moment any of them disagrees. A run that ends with `All chapter 1 numbers reproduced.` is a passing run.

```bash
python3 python/ch01_charts.py --out figures
```

That writes the chapter's two charts as SVG files into `figures/`. The book's copies are made the same way, so a chart is corrected by fixing the script and running it again, never by editing the SVG.

## Layout

The folders arrive as the chapters need them.

- `python/` holds one script of numbers and, where a chapter has charts, one script of charts per chapter, named after the chapter.
- `hdl/` will hold the SystemVerilog for each encoder and decoder the book builds.
- `tb/` will hold the cocotb testbenches that drive those designs under Verilator.
- `formal/` will hold the properties that prove a design correct for every input rather than for the inputs a test happened to try.
- `synth/` will hold the open-source synthesis scripts that report each design's gate count and depth.
- `data/` will hold the measured data sets a chapter plots, each with a note on where it came from.

## License

MIT. See [LICENSE](LICENSE).
