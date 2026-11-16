// secded72_bench.sv
//
// The bench "Project: A SEC-DED Encoder and Decoder" simulates: the (72,64)
// encoder, a channel that flips the codeword bits wherever the error pattern
// holds a 1, and the decoder, for the code the CODE parameter names.
//
//     data --> secded72_enc --> codeword --XOR--> received --> secded72_dec --> data_out,
//                                           ^                                 corrected,
//                                         error                               uncorrectable
//
// The testbench drives `data` and `error` and compares the codeword and the
// decoder's three outputs with the model in tb/secded72/secded72_model.py. The
// channel lives here, in the bench, so the encoder and the decoder stay
// synthesizable blocks a real design could use as they are.
//
// From "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
// Needs secded72_pkg.sv, secded72_enc.sv, secded72_syndrome.sv and
// secded72_dec.sv. For simulation; it adds nothing a design would keep.

module secded72_bench #(
    parameter int CODE = secded72_pkg::HSIAO,       // which code of secded72_pkg.sv
    parameter int FLAG = 0,                         // passed to the decoder
    parameter int N    = secded72_pkg::n_of(CODE),  // derived from CODE; never set by hand
    parameter int K    = secded72_pkg::k_of(CODE)   // derived from CODE; never set by hand
) (
    input  logic [K-1:0] data,                      // data word driven by the testbench
    input  logic [N-1:0] error,                     // 1 where the channel flips a bit
    output logic [N-1:0] codeword,                  // what the encoder sends
    output logic [N-1:0] received,                  // what the decoder sees
    output logic [K-1:0] data_out,                  // the decoder's data bits
    output logic         corrected,                 // the decoder's corrected flag
    output logic         uncorrectable              // the decoder's uncorrectable flag
);

    secded72_enc #(.CODE(CODE)) u_enc (
        .data     (data),
        .codeword (codeword)
    );

    assign received = codeword ^ error;

    secded72_dec #(.CODE(CODE), .FLAG(FLAG)) u_dec (
        .received      (received),
        .data          (data_out),
        .corrected     (corrected),
        .uncorrectable (uncorrectable)
    );

endmodule
