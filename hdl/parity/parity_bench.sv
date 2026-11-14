// parity_bench.sv
//
// The bench chapter 12 simulates: a parity generator, a channel that flips
// the codeword bits wherever the error mask holds a 1, and a parity checker.
//
//     data --> parity_gen --> codeword --XOR--> received --> parity_check --> error
//                                         ^
//                                    error_mask
//
// The testbench drives `data` and `error_mask` and reads the other three
// signals. A 1 in bit i of the mask flips bit i of the codeword, which is how
// a binary symmetric channel acts on a word. The channel lives here, in the
// bench, and not in the generator or the checker, which stay synthesizable
// blocks a real design could use as they are.
//
// Chapter 12 of Error Correction Coding, "Lab: A First Hardware Checker".

module parity_bench #(
    parameter int K = 7                 // data bits; codewords have K+1
) (
    input  logic [K-1:0] data,          // message driven by the testbench
    input  logic [K:0]   error_mask,    // 1 where the channel flips a bit
    output logic [K:0]   codeword,      // what the generator sends
    output logic [K:0]   received,      // what the checker sees
    output logic         error          // the checker's verdict
);

    parity_gen #(.K(K)) u_gen (
        .data     (data),
        .codeword (codeword)
    );

    // The channel: flip every bit the mask marks.
    assign received = codeword ^ error_mask;

    parity_check #(.K(K)) u_check (
        .received (received),
        .error    (error)
    );

endmodule
