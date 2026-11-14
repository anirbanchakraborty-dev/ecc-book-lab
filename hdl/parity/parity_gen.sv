// parity_gen.sv
//
// Even-parity generator for the single parity check code. It takes K data
// bits and returns the codeword of n = K+1 bits: the K data bits unchanged,
// with one parity bit on top chosen so that the codeword holds an even number
// of 1s. The parity bit is the exclusive OR of the data bits, made by an
// xor_tree of K-1 gates.
//
// Chapter 12 of Error Correction Coding, "Lab: A First Hardware Checker".
// Synthesizable. Purely combinational.

module parity_gen #(
    parameter int K = 7                 // data bits; the codeword has K+1
) (
    input  logic [K-1:0] data,          // the message
    output logic [K:0]   codeword       // {parity, data}, always of even weight
);

    logic parity;

    xor_tree #(.W(K)) u_tree (
        .in_bits (data),
        .out_bit (parity)
    );

    // The parity bit goes in the top position, bit K.
    assign codeword = {parity, data};

endmodule
