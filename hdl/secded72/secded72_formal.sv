// secded72_formal.sv
//
// The formal wrapper "Project: A SEC-DED Encoder and Decoder" hands to a solver.
// It is not a circuit to build: it states what the (72,64) encoder and decoder
// must guarantee, and the solver searches every input for a case that breaks it.
//
//     d --> secded72_enc --> codeword --XOR--> received --> secded72_dec --> data_out,
//                                        ^                                 corrected,
//                                        e                                 uncorrectable
//
// The data word d is a free input, so the solver may pick any of its 2^64
// values. The error pattern e is built from two optional positions: p1 and p2
// are free position numbers, v1 and v2 say whether each is used, and e holds a 1
// at every used position. That makes e range over exactly the patterns of weight
// 0, 1 and 2, and wt counts its ones. Four assertions state the guarantees:
//
//     wt = 0   data_out = d, neither flag raised
//     wt = 1   data_out = d, corrected raised, uncorrectable low
//     wt = 2   data_out = the received data bits, unchanged; corrected low,
//              uncorrectable raised
//     any      data_out = d unless uncorrectable is raised, so the decoder never
//              hands out wrong data without saying so
//
// Three cover statements ask the solver to show one input of each weight. They
// guard against a vacuous proof: when the assumptions admit no input of some
// weight, every assertion about that weight holds trivially, and only a cover
// that cannot be reached reveals it. VACUOUS = 1 plants that mistake on purpose,
// an assumption that leaves the second position unused, so the runner can show
// the weight-2 cover failing.
//
// From "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
// Needs secded72_pkg.sv, secded72_enc.sv, secded72_syndrome.sv and
// secded72_dec.sv. For formal tools only: it uses assume, assert and cover.

module secded72_formal #(
    parameter int CODE    = secded72_pkg::HSIAO,       // which code of secded72_pkg.sv
    parameter int FLAG    = 0,                         // passed to the decoder
    parameter int VACUOUS = 0,                         // 0 is the intended assumption; 1 plants a mistake
    parameter int N       = secded72_pkg::n_of(CODE),  // derived from CODE; never set by hand
    parameter int K       = secded72_pkg::k_of(CODE)   // derived from CODE; never set by hand
) (
    input logic [K-1:0] d,                             // any data word
    input logic [6:0]   p1,                            // first error position, if v1
    input logic [6:0]   p2,                            // second error position, if v2
    input logic         v1,                            // the first position holds an error
    input logic         v2                             // the second position holds an error
);

    localparam logic [secded72_pkg::MAX_K*8-1:0] DATA_POS = secded72_pkg::data_pos_of(CODE);

    logic [N-1:0] e, codeword, received;
    logic [K-1:0] data_out, received_data;
    logic         corrected, uncorrectable;
    logic [1:0]   wt;

    assign e = (v1 ? (N'(1) << p1) : '0) | (v2 ? (N'(1) << p2) : '0);
    assign wt = (v1 && v2 && p1 != p2) ? 2'd2 : (v1 || v2) ? 2'd1 : 2'd0;

    secded72_enc #(.CODE(CODE)) u_enc (
        .data     (d),
        .codeword (codeword)
    );

    assign received = codeword ^ e;

    secded72_dec #(.CODE(CODE), .FLAG(FLAG)) u_dec (
        .received      (received),
        .data          (data_out),
        .corrected     (corrected),
        .uncorrectable (uncorrectable)
    );

    for (genvar b = 0; b < K; b++) begin : g_received_data
        localparam int P = int'(DATA_POS[8*b +: 8]);
        assign received_data[b] = received[P];
    end

    always_comb begin
        assume (int'(p1) < N);
        assume (int'(p2) < N);
        if (VACUOUS == 1) assume (!v2);

        if (wt == 2'd0) assert (data_out == d && !corrected && !uncorrectable);
        if (wt == 2'd1) assert (data_out == d && corrected && !uncorrectable);
        if (wt == 2'd2) assert (data_out == received_data && !corrected && uncorrectable);
        assert (data_out == d || uncorrectable);

        cover (wt == 2'd0);
        cover (wt == 2'd1);
        cover (wt == 2'd2);
    end

endmodule
