// secded72_syndrome.sv
//
// The syndrome stage of the decoder in secded72_dec.sv, kept as a module of its
// own so that synth/synth_secded72.py can count its gates apart from the rest.
// Syndrome bit i is one XOR tree over the received bits whose column holds a 1
// in row i, the row the package stores in its ROWS table. For the balanced Hsiao
// (72,64) code every row holds 27 positions, so the stage is eight trees of 27
// inputs. A codeword gives the syndrome 0, and a received word that differs from
// a codeword by an error pattern e gives the XOR of the columns e marks.
//
// From "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
// Needs secded72_pkg.sv. Synthesizable. Purely combinational.

module secded72_syndrome #(
    parameter int CODE = secded72_pkg::HSIAO,       // which code of secded72_pkg.sv
    parameter int N    = secded72_pkg::n_of(CODE),  // derived from CODE; never set by hand
    parameter int R    = secded72_pkg::r_of(CODE)   // derived from CODE; never set by hand
) (
    input  logic [N-1:0] received,                  // bit j is codeword position j
    output logic [R-1:0] syndrome                   // bit i is syndrome row i
);

    localparam int MAX_N = secded72_pkg::MAX_N;
    localparam logic [secded72_pkg::MAX_R*MAX_N-1:0] ROWS = secded72_pkg::rows_of(CODE);

    for (genvar i = 0; i < R; i++) begin : g_row
        localparam logic [N-1:0] ROW = ROWS[MAX_N*i +: N];
        assign syndrome[i] = ^(received & ROW);
    end

endmodule
