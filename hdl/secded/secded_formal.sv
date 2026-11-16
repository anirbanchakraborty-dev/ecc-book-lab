// secded_formal.sv
//
// The formal wrapper "Verifying ECC Hardware" hands to a solver. It is not a
// circuit to build: it states what the Hsiao (13,8) encoder and decoder must
// guarantee, and the solver searches every input for a case that breaks it.
//
//     d --> secded_enc --> codeword --XOR--> received --> secded_dec --> data_out,
//                                     ^                                corrected,
//                                     e                                uncorrectable
//
// The data word d and the error pattern e are free inputs, so the solver may
// pick any of their 2^8 * 2^13 values. One assumption keeps it inside the design
// envelope, wt(e) <= 2, and four assertions state the guarantees:
//
//     wt(e) = 0   data_out = d, neither flag raised
//     wt(e) = 1   data_out = d, corrected raised, uncorrectable low
//     wt(e) = 2   data_out = d ^ e[7:0] (data passed through untouched),
//                 corrected low, uncorrectable raised
//     any         data_out = d unless uncorrectable is raised, so the decoder
//                 never hands out wrong data without saying so
//
// Three cover statements ask the solver to show one input of each weight. They
// guard against a vacuous proof: when the assumptions admit no input of some
// weight, every assertion about that weight holds trivially, and only a cover
// that cannot be reached reveals it.
//
// MUTANT is passed to the decoder. VACUOUS swaps the envelope assumption for
// one of two plausible mistakes, so the chapter can show what each one hides:
//
//     0   the intended assumption, wt(e) <= 2
//     1   an off-by-one, wt(e) < 2: the inputs that remain satisfy it, so the
//         solver runs, but no double error is ever tried, the weight-2
//         assertion is never exercised and only the weight-2 cover fails
//     2   one assumption per weight class, written where an implication was
//         meant: together they demand wt(e) = 0, 1 and 2 at once, no input
//         satisfies them, and every assertion would hold trivially
//
// From "Verifying ECC Hardware" in Error Correction Coding.
// For formal tools only: it uses assume, assert and cover.

module secded_formal #(
    parameter int MUTANT  = 0,          // passed to the decoder; 0 is the correct design
    parameter int VACUOUS = 0           // 0 is the intended assumption; 1 and 2 plant a mistake
) (
    input logic [7:0]  d,               // any data word
    input logic [12:0] e                // any error pattern
);

    logic [12:0] codeword, received;
    logic [7:0]  data_out;
    logic        corrected, uncorrectable;

    secded_enc u_enc (
        .data     (d),
        .codeword (codeword)
    );

    assign received = codeword ^ e;

    secded_dec #(.MUTANT(MUTANT)) u_dec (
        .received      (received),
        .data          (data_out),
        .corrected     (corrected),
        .uncorrectable (uncorrectable)
    );

    logic [3:0] wt;                     // the weight of e, 0 to 13
    always_comb begin
        wt = 4'd0;
        for (int i = 0; i < 13; i++) wt = wt + {3'd0, e[i]};
    end

    always_comb begin
        if (VACUOUS == 0) begin
            assume (wt <= 4'd2);
        end else if (VACUOUS == 1) begin
            assume (wt < 4'd2);
        end else begin
            assume (wt == 4'd0);
            assume (wt == 4'd1);
            assume (wt == 4'd2);
        end

        if (wt == 4'd0) assert (data_out == d && !corrected && !uncorrectable);
        if (wt == 4'd1) assert (data_out == d && corrected && !uncorrectable);
        if (wt == 4'd2) assert (data_out == (d ^ e[7:0]) && !corrected && uncorrectable);
        assert (data_out == d || uncorrectable);

        cover (wt == 4'd0);
        cover (wt == 4'd1);
        cover (wt == 4'd2);
    end

endmodule
