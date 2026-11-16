// secded72_enc.sv
//
// Encoder for the SEC-DED codes of secded72_pkg.sv, the balanced Hsiao (72,64)
// code by default and the extended Hamming (72,64) code with CODE set to
// secded72_pkg::HAMMING. It places each data bit at its codeword position and
// computes each check bit as one XOR tree over the data bits the package lists
// for it, which makes every row of the check matrix sum to 0 over a codeword.
// For the Hsiao code that is eight trees of 26 inputs each, and the codeword is
// the 64 data bits in positions 0 to 63 followed by the 8 check bits.
//
// The extended Hamming code keeps the bit order of "SECDED and Hsiao Codes":
// its check bits sit at positions 0, 1, 2, 4, 8, 16, 32 and 64 and its data
// bits fill the other 64 positions in increasing order.
//
// From "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
// Needs secded72_pkg.sv. Synthesizable. Purely combinational.

module secded72_enc #(
    parameter int CODE = secded72_pkg::HSIAO,       // which code of secded72_pkg.sv
    parameter int N    = secded72_pkg::n_of(CODE),  // derived from CODE; never set by hand
    parameter int K    = secded72_pkg::k_of(CODE)   // derived from CODE; never set by hand
) (
    input  logic [K-1:0] data,                      // the data bits
    output logic [N-1:0] codeword                   // bit j is codeword position j
);

    localparam int R = secded72_pkg::r_of(CODE);
    localparam int MAX_K = secded72_pkg::MAX_K;
    localparam logic [MAX_K*8-1:0] DATA_POS = secded72_pkg::data_pos_of(CODE);
    localparam logic [secded72_pkg::MAX_R*8-1:0] CHECK_POS = secded72_pkg::check_pos_of(CODE);
    localparam logic [secded72_pkg::MAX_R*MAX_K-1:0] CHECK_MASKS = secded72_pkg::check_masks_of(CODE);

    for (genvar b = 0; b < K; b++) begin : g_data
        localparam int P = int'(DATA_POS[8*b +: 8]);
        assign codeword[P] = data[b];
    end

    for (genvar t = 0; t < R; t++) begin : g_check
        localparam int P = int'(CHECK_POS[8*t +: 8]);
        localparam logic [K-1:0] MASK = CHECK_MASKS[MAX_K*t +: K];
        assign codeword[P] = ^(data & MASK);
    end

endmodule
