// secded_enc.sv
//
// Encoder for the Hsiao (13,8) SEC-DED code. It takes 8 data bits and appends
// the 5 check bits of the code's check matrix, column for column the matrix
// that "SECDED and Hsiao Codes" draws. Codeword bits 0 to 7 carry the data
// unchanged and bits 8 to 12 carry check bits c0 to c4. Check bit c_i is the
// exclusive OR of the data bits whose column holds a 1 in row i, which makes
// every row of the check matrix sum to 0 over a codeword:
//
//     c0 = d0 ^ d1 ^ d2 ^ d3
//     c1 = d0 ^ d1 ^ d4 ^ d5 ^ d6
//     c2 = d2 ^ d3 ^ d4 ^ d5 ^ d7
//     c3 = d0 ^ d2 ^ d4 ^ d6 ^ d7
//     c4 = d1 ^ d3 ^ d5 ^ d6 ^ d7
//
// From "Verifying ECC Hardware" in Error Correction Coding.
// Synthesizable. Purely combinational.

module secded_enc (
    input  logic [7:0]  data,           // the 8 data bits d0 to d7
    output logic [12:0] codeword        // {c4 .. c0, d7 .. d0}
);

    logic [4:0] check;                  // check[i] is c_i

    assign check[0] = data[0] ^ data[1] ^ data[2] ^ data[3];
    assign check[1] = data[0] ^ data[1] ^ data[4] ^ data[5] ^ data[6];
    assign check[2] = data[2] ^ data[3] ^ data[4] ^ data[5] ^ data[7];
    assign check[3] = data[0] ^ data[2] ^ data[4] ^ data[6] ^ data[7];
    assign check[4] = data[1] ^ data[3] ^ data[5] ^ data[6] ^ data[7];

    assign codeword = {check, data};

endmodule
