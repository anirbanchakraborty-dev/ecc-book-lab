// secded_dec.sv
//
// Decoder for the Hsiao (13,8) SEC-DED code that secded_enc.sv encodes. It takes
// the 13 received bits and returns the 8 data bits with two flags:
//
//     corrected       the syndrome equals one column of the check matrix, so a
//                     single error was found; when the column is a data column
//                     that data bit is flipped back, and when it is a check
//                     column the data bits already are right and pass unchanged
//     uncorrectable   the syndrome is nonzero and equals no column, which is
//                     what every double error produces; the data bits then pass
//                     through exactly as received, with no bit flipped
//
// With no error the syndrome is 0, the data passes unchanged and neither flag
// is raised. The syndrome is declared [0:4], so a literal such as 5'b01101
// reads row 0 first, the order in which "SECDED and Hsiao Codes" prints every
// column and syndrome.
//
// MUTANT plants one small bug, so that "Verifying ECC Hardware" can measure how
// strong each check is. 0 builds the correct decoder and is the only value a
// real design would use. The others:
//
//     1   wiring slip: the XOR tree of syndrome bit 1 leaves out received bit 6
//     2   flag slip: uncorrectable is raised only for a syndrome with exactly
//         two 1s, so double errors whose syndrome has four 1s go unflagged
//     3   steering slip: the comparator for d5's column flips d4, not d5
//     4   data-dependent slip: d0's correction is gated by received bit 7, so it
//         fires only while that data bit is 0
//
// From "Verifying ECC Hardware" in Error Correction Coding.
// Synthesizable. Purely combinational.

module secded_dec #(
    parameter int MUTANT = 0            // 0 is the correct decoder; 1 to 4 plant a bug
) (
    input  logic [12:0] received,       // {c4 .. c0, d7 .. d0} as it arrived
    output logic [7:0]  data,           // the data bits, corrected when possible
    output logic        corrected,      // a single error was found and handled
    output logic        uncorrectable   // an error was found that cannot be corrected
);

    // Syndrome bit i is check bit c_i XORed with the data bits of row i.
    localparam logic D6_IN_S1 = (MUTANT != 1);

    /* verilator lint_off ASCRANGE */
    logic [0:4] syndrome;               // ascending on purpose, so literals read row 0 first
    /* verilator lint_on ASCRANGE */
    assign syndrome[0] = received[8]  ^ received[0] ^ received[1] ^ received[2] ^ received[3];
    assign syndrome[1] = received[9]  ^ received[0] ^ received[1] ^ received[4] ^ received[5]
                       ^ (received[6] & D6_IN_S1);
    assign syndrome[2] = received[10] ^ received[2] ^ received[3] ^ received[4] ^ received[5]
                       ^ received[7];
    assign syndrome[3] = received[11] ^ received[0] ^ received[2] ^ received[4] ^ received[6]
                       ^ received[7];
    assign syndrome[4] = received[12] ^ received[1] ^ received[3] ^ received[5] ^ received[6]
                       ^ received[7];

    // One comparator per column of the check matrix, row 0 first.
    logic [7:0] data_hit;               // data_hit[j]: the syndrome is d_j's column
    logic [4:0] check_hit;              // check_hit[i]: the syndrome is c_i's column
    assign data_hit[0] = (syndrome == 5'b11010);
    assign data_hit[1] = (syndrome == 5'b11001);
    assign data_hit[2] = (syndrome == 5'b10110);
    assign data_hit[3] = (syndrome == 5'b10101);
    assign data_hit[4] = (syndrome == 5'b01110);
    assign data_hit[5] = (syndrome == 5'b01101);
    assign data_hit[6] = (syndrome == 5'b01011);
    assign data_hit[7] = (syndrome == 5'b00111);
    assign check_hit[0] = (syndrome == 5'b10000);
    assign check_hit[1] = (syndrome == 5'b01000);
    assign check_hit[2] = (syndrome == 5'b00100);
    assign check_hit[3] = (syndrome == 5'b00010);
    assign check_hit[4] = (syndrome == 5'b00001);

    // The correction mask: which data bits to flip back.
    logic [7:0] flip;
    always_comb begin
        flip = data_hit;
        if (MUTANT == 3) begin
            flip[4] = data_hit[4] | data_hit[5];
            flip[5] = 1'b0;
        end
        if (MUTANT == 4) begin
            flip[0] = data_hit[0] & ~received[7];
        end
    end

    // Exactly two 1s in the syndrome, which only mutant 2 uses.
    logic two_ones;
    always_comb begin
        int count;
        count = 0;
        for (int i = 0; i < 5; i++) count = count + int'(syndrome[i]);
        two_ones = (count == 2);
    end

    logic nonzero;
    assign nonzero = |syndrome;

    assign data      = received[7:0] ^ flip;
    assign corrected = |{data_hit, check_hit};
    assign uncorrectable = (MUTANT == 2) ? (nonzero & ~corrected & two_ones)
                                         : (nonzero & ~corrected);

endmodule
