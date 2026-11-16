// secded72_dec.sv
//
// Decoder for the SEC-DED codes of secded72_pkg.sv, the balanced Hsiao (72,64)
// code by default and the extended Hamming (72,64) code with CODE set to
// secded72_pkg::HAMMING. It takes the received codeword and returns the data
// bits with two flags, in three stages:
//
//     syndrome     eight XOR trees, one per row of the check matrix, in
//                  secded72_syndrome.sv
//     column match one 8-bit comparator per codeword position, 72 of them,
//                  each 1 when the syndrome equals that position's column
//     correction   one XOR per data bit, flipping the bit whose comparator
//                  fired; the comparators of the check positions flip nothing
//
// and two flags:
//
//     corrected       some comparator fired, so a single error was found and
//                     corrected, in a data bit or in a check bit
//     uncorrectable   with FLAG = 0, the syndrome is nonzero and no comparator
//                     fired, which every double error produces; the data bits
//                     then pass through exactly as received
//
// With FLAG = 1 the uncorrectable flag takes its parity form instead: for a
// Hsiao code, a nonzero syndrome with an even number of ones, and for the
// extended Hamming code, nonzero Hamming bits with the overall parity bit 0.
// "Project: A SEC-DED Encoder and Decoder" proves that the two forms agree on
// every error pattern of weight 2 or less. They differ only beyond it, where
// FLAG = 1 lets some uncorrectable words through unflagged.
//
// Every output depends on the received word only through the syndrome and the
// received data bits, which is the shape the theorem in "When One Codeword
// Stands for All" asks for.
//
// From "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
// Needs secded72_pkg.sv and secded72_syndrome.sv. Synthesizable. Purely
// combinational: no register, so the decoder adds no clock cycle to a read.

module secded72_dec #(
    parameter int CODE = secded72_pkg::HSIAO,       // which code of secded72_pkg.sv
    parameter int FLAG = 0,                         // 0: no-match form, 1: parity form
    parameter int N    = secded72_pkg::n_of(CODE),  // derived from CODE; never set by hand
    parameter int K    = secded72_pkg::k_of(CODE)   // derived from CODE; never set by hand
) (
    input  logic [N-1:0] received,                  // bit j is codeword position j
    output logic [K-1:0] data,                      // the corrected data bits
    output logic         corrected,                 // a single error was corrected
    output logic         uncorrectable              // the word holds an error it cannot correct
);

    localparam int R = secded72_pkg::r_of(CODE);
    localparam int MAX_R = secded72_pkg::MAX_R;
    localparam int PARITY_ROW = secded72_pkg::parity_row_of(CODE);
    localparam logic [secded72_pkg::MAX_N*MAX_R-1:0] COLUMNS = secded72_pkg::columns_of(CODE);
    localparam logic [secded72_pkg::MAX_K*8-1:0] DATA_POS = secded72_pkg::data_pos_of(CODE);

    // ---- syndrome stage
    logic [R-1:0] syndrome;

    secded72_syndrome #(.CODE(CODE)) u_syndrome (
        .received (received),
        .syndrome (syndrome)
    );

    // ---- column-match stage: match[j] is 1 when the syndrome equals column j
    logic [N-1:0] match;

    for (genvar j = 0; j < N; j++) begin : g_match
        localparam logic [R-1:0] COLUMN = COLUMNS[MAX_R*j +: R];
        assign match[j] = (syndrome == COLUMN);
    end

    // ---- correction stage: flip the data bit whose column matched
    for (genvar b = 0; b < K; b++) begin : g_correct
        localparam int P = int'(DATA_POS[8*b +: 8]);
        assign data[b] = received[P] ^ match[P];
    end

    // ---- flags
    assign corrected = |match;

    if (FLAG == 0) begin : g_nomatch
        assign uncorrectable = (|syndrome) & ~corrected;
    end else if (PARITY_ROW < 0) begin : g_even
        assign uncorrectable = (|syndrome) & ~(^syndrome);
    end else begin : g_parity_bit
        assign uncorrectable = (|syndrome[PARITY_ROW-1:0]) & ~syndrome[PARITY_ROW];
    end

endmodule
