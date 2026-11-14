// parity_check.sv
//
// Checker for the single parity check code. It takes a received word of
// n = K+1 bits and raises `error` when the word holds an odd number of 1s,
// which no codeword does. The flag is the exclusive OR of all n bits, made by
// an xor_tree of K gates. A word with an even number of flipped bits passes
// unflagged, and that is the code's blind spot, not a fault in the checker.
//
// Chapter 12 of Error Correction Coding, "Lab: A First Hardware Checker".
// Synthesizable. Purely combinational.

module parity_check #(
    parameter int K = 7                 // data bits; the word has K+1
) (
    input  logic [K:0] received,        // the word as it arrived
    output logic       error            // 1 when the word cannot be a codeword
);

    xor_tree #(.W(K + 1)) u_tree (
        .in_bits (received),
        .out_bit (error)
    );

endmodule
