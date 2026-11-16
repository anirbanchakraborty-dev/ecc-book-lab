// secded_bench.sv
//
// The bench "Verifying ECC Hardware" simulates: the Hsiao (13,8) encoder, a
// channel that flips the codeword bits wherever the error pattern holds a 1,
// and the decoder.
//
//     data --> secded_enc --> codeword --XOR--> received --> secded_dec --> data_out,
//                                         ^                               corrected,
//                                       error                             uncorrectable
//
// The testbench drives `data` and `error` and compares the decoder's three
// outputs with a reference model written in Python. The channel lives here, in
// the bench, so the encoder and the decoder stay synthesizable blocks a real
// design could use as they are.
//
// From "Verifying ECC Hardware" in Error Correction Coding.

module secded_bench #(
    parameter int MUTANT = 0            // passed to the decoder; 0 is the correct design
) (
    input  logic [7:0]  data,           // data word driven by the testbench
    input  logic [12:0] error,          // 1 where the channel flips a bit
    output logic [12:0] codeword,       // what the encoder sends
    output logic [12:0] received,       // what the decoder sees
    output logic [7:0]  data_out,       // the decoder's data bits
    output logic        corrected,      // the decoder's corrected flag
    output logic        uncorrectable   // the decoder's uncorrectable flag
);

    secded_enc u_enc (
        .data     (data),
        .codeword (codeword)
    );

    assign received = codeword ^ error;

    secded_dec #(.MUTANT(MUTANT)) u_dec (
        .received      (received),
        .data          (data_out),
        .corrected     (corrected),
        .uncorrectable (uncorrectable)
    );

endmodule
