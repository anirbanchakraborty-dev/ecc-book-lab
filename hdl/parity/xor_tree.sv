// xor_tree.sv
//
// The exclusive OR of W input bits, built as a balanced tree of two-input
// XOR gates. Chapter 12 of Error Correction Coding, "Lab: A First Hardware
// Checker", uses it twice: once to make the parity bit of a data word, and
// once to check a received word.
//
// The tree is built by splitting. One bit needs no gate: the bit is its own
// exclusive OR. Wider inputs split into a low half of floor(W/2) bits and a
// high half of the rest, each half gets a tree of its own, and one more gate
// combines the two results. A tree over W bits therefore holds W-1 gates, and
// its deepest path passes through ceil(log2(W)) of them. The chapter proves
// both counts.
//
// Synthesizable. Purely combinational: no clock, no reset, no state.

module xor_tree #(
    parameter int W = 8                 // number of input bits, at least 1
) (
    input  logic [W-1:0] in_bits,       // the bits to combine
    output logic         out_bit        // their exclusive OR, 1 when an odd number are 1
);

    if (W == 1) begin : g_leaf
        // A single bit is its own exclusive OR. No gate.
        assign out_bit = in_bits[0];
    end else begin : g_split
        localparam int L = W / 2;       // bits in the low half, rounded down
        localparam int H = W - L;       // bits in the high half
        logic low_bit, high_bit;

        xor_tree #(.W(L)) u_low (
            .in_bits (in_bits[L-1:0]),
            .out_bit (low_bit)
        );

        xor_tree #(.W(H)) u_high (
            .in_bits (in_bits[W-1:L]),
            .out_bit (high_bit)
        );

        // The one gate this level adds.
        assign out_bit = low_bit ^ high_bit;
    end

endmodule
