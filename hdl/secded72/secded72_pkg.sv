// secded72_pkg.sv
//
// Check matrices of the SEC-DED encoder and decoder in hdl/secded72/. Written
// by tb/secded72/gen_secded72_pkg.py from the column lists in
// tb/secded72/secded72_model.py; edit the model and run the generator, never
// this file, and tb/secded72/check_secded72.py fails while the two disagree.
//
// A module picks its code with a CODE parameter set to one of the numbers
// below. Every table is laid out for the widest code, 72 positions of 8 rows,
// and a smaller code leaves the slots above its own size at 0. A column is
// written row 0 in the low bit, so 8'b00000111 has ones in rows 0, 1 and 2.
//
// From "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
// Synthesizable: constants and constant functions only.

package secded72_pkg;

    localparam int HSIAO = 0;
    localparam int HAMMING = 1;
    localparam int HSIAO_DICT = 2;
    localparam int HSIAO_39_32 = 3;

    localparam int MAX_N = 72;
    localparam int MAX_K = 64;
    localparam int MAX_R = 8;

    // ---- HSIAO: balanced Hsiao (72,64), hsiao_72_64() of python/ch20_secded_and_hsiao_codes.py

    localparam int HSIAO_N = 72;
    localparam int HSIAO_K = 64;
    localparam int HSIAO_R = 8;
    localparam int HSIAO_PARITY_ROW = -1;

    // column of the check matrix at each codeword position, row 0 in bit 0
    localparam logic [72*8-1:0] HSIAO_COLUMNS = {
        8'b10000000, 8'b01000000, 8'b00100000, 8'b00010000, 8'b00001000, 8'b00000100, 8'b00000010, 8'b00000001,  // 71 .. 64
        8'b10001111, 8'b11000111, 8'b11100011, 8'b11110001, 8'b11111000, 8'b01111100, 8'b00111110, 8'b00011111,  // 63 .. 56
        8'b11100000, 8'b11010000, 8'b10110000, 8'b01110000, 8'b11001000, 8'b10101000, 8'b01101000, 8'b10011000,  // 55 .. 48
        8'b01011000, 8'b00111000, 8'b11000100, 8'b10100100, 8'b01100100, 8'b10010100, 8'b01010100, 8'b00110100,  // 47 .. 40
        8'b10001100, 8'b01001100, 8'b00101100, 8'b00011100, 8'b11000010, 8'b10100010, 8'b01100010, 8'b10010010,  // 39 .. 32
        8'b01010010, 8'b00110010, 8'b10001010, 8'b01001010, 8'b00101010, 8'b00011010, 8'b10000110, 8'b01000110,  // 31 .. 24
        8'b00100110, 8'b00010110, 8'b00001110, 8'b11000001, 8'b10100001, 8'b01100001, 8'b10010001, 8'b01010001,  // 23 .. 16
        8'b00110001, 8'b10001001, 8'b01001001, 8'b00101001, 8'b00011001, 8'b10000101, 8'b01000101, 8'b00100101,  // 15 .. 8
        8'b00010101, 8'b00001101, 8'b10000011, 8'b01000011, 8'b00100011, 8'b00010011, 8'b00001011, 8'b00000111  // 7 .. 0
    };

    // the same matrix by rows: bit j of row i is set when position j is in row i
    localparam logic [8*72-1:0] HSIAO_ROWS = {
        72'h80f8ed348d221a4420,  // 7
        72'h407cdaaa4a91152210,  // 6
        72'h203eb65926488c9108,  // 5
        72'h101f71c711c4438884,  // 4
        72'h088f0fc0f03c207842,  // 3
        72'h04c7003ff003e007c1,  // 2
        72'h02e300000fffe0003f,  // 1
        72'h01f1000000001fffff  // 0
    };

    // codeword position of each data bit
    localparam logic [64*8-1:0] HSIAO_DATA_POS = {
        8'd63, 8'd62, 8'd61, 8'd60, 8'd59, 8'd58, 8'd57, 8'd56, 8'd55, 8'd54, 8'd53, 8'd52, 8'd51, 8'd50, 8'd49, 8'd48,  // 63 .. 48
        8'd47, 8'd46, 8'd45, 8'd44, 8'd43, 8'd42, 8'd41, 8'd40, 8'd39, 8'd38, 8'd37, 8'd36, 8'd35, 8'd34, 8'd33, 8'd32,  // 47 .. 32
        8'd31, 8'd30, 8'd29, 8'd28, 8'd27, 8'd26, 8'd25, 8'd24, 8'd23, 8'd22, 8'd21, 8'd20, 8'd19, 8'd18, 8'd17, 8'd16,  // 31 .. 16
        8'd15, 8'd14, 8'd13, 8'd12, 8'd11, 8'd10, 8'd9, 8'd8, 8'd7, 8'd6, 8'd5, 8'd4, 8'd3, 8'd2, 8'd1, 8'd0  // 15 .. 0
    };

    // codeword position of each check bit
    localparam logic [8*8-1:0] HSIAO_CHECK_POS = {
        8'd71, 8'd70, 8'd69, 8'd68, 8'd67, 8'd66, 8'd65, 8'd64  // 7 .. 0
    };

    // data bits whose XOR is each check bit
    localparam logic [8*64-1:0] HSIAO_CHECK_MASKS = {
        64'hf8ed348d221a4420,  // 7
        64'h7cdaaa4a91152210,  // 6
        64'h3eb65926488c9108,  // 5
        64'h1f71c711c4438884,  // 4
        64'h8f0fc0f03c207842,  // 3
        64'hc7003ff003e007c1,  // 2
        64'he300000fffe0003f,  // 1
        64'hf1000000001fffff  // 0
    };
    // ---- HAMMING: extended Hamming (72,64), extended_hamming(72, 7) of the same script

    localparam int HAMMING_N = 72;
    localparam int HAMMING_K = 64;
    localparam int HAMMING_R = 8;
    localparam int HAMMING_PARITY_ROW = 7;

    // column of the check matrix at each codeword position, row 0 in bit 0
    localparam logic [72*8-1:0] HAMMING_COLUMNS = {
        8'b11000111, 8'b11000110, 8'b11000101, 8'b11000100, 8'b11000011, 8'b11000010, 8'b11000001, 8'b11000000,  // 71 .. 64
        8'b10111111, 8'b10111110, 8'b10111101, 8'b10111100, 8'b10111011, 8'b10111010, 8'b10111001, 8'b10111000,  // 63 .. 56
        8'b10110111, 8'b10110110, 8'b10110101, 8'b10110100, 8'b10110011, 8'b10110010, 8'b10110001, 8'b10110000,  // 55 .. 48
        8'b10101111, 8'b10101110, 8'b10101101, 8'b10101100, 8'b10101011, 8'b10101010, 8'b10101001, 8'b10101000,  // 47 .. 40
        8'b10100111, 8'b10100110, 8'b10100101, 8'b10100100, 8'b10100011, 8'b10100010, 8'b10100001, 8'b10100000,  // 39 .. 32
        8'b10011111, 8'b10011110, 8'b10011101, 8'b10011100, 8'b10011011, 8'b10011010, 8'b10011001, 8'b10011000,  // 31 .. 24
        8'b10010111, 8'b10010110, 8'b10010101, 8'b10010100, 8'b10010011, 8'b10010010, 8'b10010001, 8'b10010000,  // 23 .. 16
        8'b10001111, 8'b10001110, 8'b10001101, 8'b10001100, 8'b10001011, 8'b10001010, 8'b10001001, 8'b10001000,  // 15 .. 8
        8'b10000111, 8'b10000110, 8'b10000101, 8'b10000100, 8'b10000011, 8'b10000010, 8'b10000001, 8'b10000000  // 7 .. 0
    };

    // the same matrix by rows: bit j of row i is set when position j is in row i
    localparam logic [8*72-1:0] HAMMING_ROWS = {
        72'hffffffffffffffffff,  // 7
        72'hff0000000000000000,  // 6
        72'h00ffffffff00000000,  // 5
        72'h00ffff0000ffff0000,  // 4
        72'h00ff00ff00ff00ff00,  // 3
        72'hf0f0f0f0f0f0f0f0f0,  // 2
        72'hcccccccccccccccccc,  // 1
        72'haaaaaaaaaaaaaaaaaa  // 0
    };

    // codeword position of each data bit
    localparam logic [64*8-1:0] HAMMING_DATA_POS = {
        8'd71, 8'd70, 8'd69, 8'd68, 8'd67, 8'd66, 8'd65, 8'd63, 8'd62, 8'd61, 8'd60, 8'd59, 8'd58, 8'd57, 8'd56, 8'd55,  // 63 .. 48
        8'd54, 8'd53, 8'd52, 8'd51, 8'd50, 8'd49, 8'd48, 8'd47, 8'd46, 8'd45, 8'd44, 8'd43, 8'd42, 8'd41, 8'd40, 8'd39,  // 47 .. 32
        8'd38, 8'd37, 8'd36, 8'd35, 8'd34, 8'd33, 8'd31, 8'd30, 8'd29, 8'd28, 8'd27, 8'd26, 8'd25, 8'd24, 8'd23, 8'd22,  // 31 .. 16
        8'd21, 8'd20, 8'd19, 8'd18, 8'd17, 8'd15, 8'd14, 8'd13, 8'd12, 8'd11, 8'd10, 8'd9, 8'd7, 8'd6, 8'd5, 8'd3  // 15 .. 0
    };

    // codeword position of each check bit
    localparam logic [8*8-1:0] HAMMING_CHECK_POS = {
        8'd64, 8'd32, 8'd16, 8'd8, 8'd4, 8'd2, 8'd1, 8'd0  // 7 .. 0
    };

    // data bits whose XOR is each check bit
    localparam logic [8*64-1:0] HAMMING_CHECK_MASKS = {
        64'hfe00000000000000,  // 7
        64'h01fffffffc000000,  // 6
        64'h01fffe0003fff800,  // 5
        64'h01fe01fe03fc07f0,  // 4
        64'hf1e1e1e1e3c3c78e,  // 3
        64'hcd9999999b33366d,  // 2
        64'hab55555556aaad5b,  // 1
        64'h972cd2d32da65cb7  // 0
    };
    // ---- HSIAO_DICT: Hsiao (72,64) in dictionary order, hsiao_72_64_dictionary() there

    localparam int HSIAO_DICT_N = 72;
    localparam int HSIAO_DICT_K = 64;
    localparam int HSIAO_DICT_R = 8;
    localparam int HSIAO_DICT_PARITY_ROW = -1;

    // column of the check matrix at each codeword position, row 0 in bit 0
    localparam logic [72*8-1:0] HSIAO_DICT_COLUMNS = {
        8'b10000000, 8'b01000000, 8'b00100000, 8'b00010000, 8'b00001000, 8'b00000100, 8'b00000010, 8'b00000001,  // 71 .. 64
        8'b01100111, 8'b10010111, 8'b01010111, 8'b00110111, 8'b10001111, 8'b01001111, 8'b00101111, 8'b00011111,  // 63 .. 56
        8'b11100000, 8'b11010000, 8'b10110000, 8'b01110000, 8'b11001000, 8'b10101000, 8'b01101000, 8'b10011000,  // 55 .. 48
        8'b01011000, 8'b00111000, 8'b11000100, 8'b10100100, 8'b01100100, 8'b10010100, 8'b01010100, 8'b00110100,  // 47 .. 40
        8'b10001100, 8'b01001100, 8'b00101100, 8'b00011100, 8'b11000010, 8'b10100010, 8'b01100010, 8'b10010010,  // 39 .. 32
        8'b01010010, 8'b00110010, 8'b10001010, 8'b01001010, 8'b00101010, 8'b00011010, 8'b10000110, 8'b01000110,  // 31 .. 24
        8'b00100110, 8'b00010110, 8'b00001110, 8'b11000001, 8'b10100001, 8'b01100001, 8'b10010001, 8'b01010001,  // 23 .. 16
        8'b00110001, 8'b10001001, 8'b01001001, 8'b00101001, 8'b00011001, 8'b10000101, 8'b01000101, 8'b00100101,  // 15 .. 8
        8'b00010101, 8'b00001101, 8'b10000011, 8'b01000011, 8'b00100011, 8'b00010011, 8'b00001011, 8'b00000111  // 7 .. 0
    };

    // the same matrix by rows: bit j of row i is set when position j is in row i
    localparam logic [8*72-1:0] HSIAO_DICT_ROWS = {
        72'h8048ed348d221a4420,  // 7
        72'h40a4daaa4a91152210,  // 6
        72'h2092b65926488c9108,  // 5
        72'h107171c711c4438884,  // 4
        72'h080f0fc0f03c207842,  // 3
        72'h04ff003ff003e007c1,  // 2
        72'h02ff00000fffe0003f,  // 1
        72'h01ff000000001fffff  // 0
    };

    // codeword position of each data bit
    localparam logic [64*8-1:0] HSIAO_DICT_DATA_POS = {
        8'd63, 8'd62, 8'd61, 8'd60, 8'd59, 8'd58, 8'd57, 8'd56, 8'd55, 8'd54, 8'd53, 8'd52, 8'd51, 8'd50, 8'd49, 8'd48,  // 63 .. 48
        8'd47, 8'd46, 8'd45, 8'd44, 8'd43, 8'd42, 8'd41, 8'd40, 8'd39, 8'd38, 8'd37, 8'd36, 8'd35, 8'd34, 8'd33, 8'd32,  // 47 .. 32
        8'd31, 8'd30, 8'd29, 8'd28, 8'd27, 8'd26, 8'd25, 8'd24, 8'd23, 8'd22, 8'd21, 8'd20, 8'd19, 8'd18, 8'd17, 8'd16,  // 31 .. 16
        8'd15, 8'd14, 8'd13, 8'd12, 8'd11, 8'd10, 8'd9, 8'd8, 8'd7, 8'd6, 8'd5, 8'd4, 8'd3, 8'd2, 8'd1, 8'd0  // 15 .. 0
    };

    // codeword position of each check bit
    localparam logic [8*8-1:0] HSIAO_DICT_CHECK_POS = {
        8'd71, 8'd70, 8'd69, 8'd68, 8'd67, 8'd66, 8'd65, 8'd64  // 7 .. 0
    };

    // data bits whose XOR is each check bit
    localparam logic [8*64-1:0] HSIAO_DICT_CHECK_MASKS = {
        64'h48ed348d221a4420,  // 7
        64'ha4daaa4a91152210,  // 6
        64'h92b65926488c9108,  // 5
        64'h7171c711c4438884,  // 4
        64'h0f0fc0f03c207842,  // 3
        64'hff003ff003e007c1,  // 2
        64'hff00000fffe0003f,  // 1
        64'hff000000001fffff  // 0
    };
    // ---- HSIAO_39_32: Hsiao (39,32), hsiao_39_32() there

    localparam int HSIAO_39_32_N = 39;
    localparam int HSIAO_39_32_K = 32;
    localparam int HSIAO_39_32_R = 7;
    localparam int HSIAO_39_32_PARITY_ROW = -1;

    // column of the check matrix at each codeword position, row 0 in bit 0
    localparam logic [72*8-1:0] HSIAO_39_32_COLUMNS = {
        8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000,  // 71 .. 64
        8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000,  // 63 .. 56
        8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000,  // 55 .. 48
        8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000, 8'b00000000,  // 47 .. 40
        8'b00000000, 8'b01000000, 8'b00100000, 8'b00010000, 8'b00001000, 8'b00000100, 8'b00000010, 8'b00000001,  // 39 .. 32
        8'b01110000, 8'b01101000, 8'b01011000, 8'b00111000, 8'b01100100, 8'b01010100, 8'b00110100, 8'b01001100,  // 31 .. 24
        8'b00101100, 8'b00011100, 8'b01010010, 8'b00110010, 8'b01001010, 8'b00101010, 8'b00011010, 8'b01000110,  // 23 .. 16
        8'b00100110, 8'b00010110, 8'b00001110, 8'b01100001, 8'b01010001, 8'b00110001, 8'b01001001, 8'b00101001,  // 15 .. 8
        8'b01000101, 8'b00100101, 8'b00010101, 8'b00001101, 8'b01000011, 8'b00100011, 8'b00010011, 8'b00001011  // 7 .. 0
    };

    // the same matrix by rows: bit j of row i is set when position j is in row i
    localparam logic [8*72-1:0] HSIAO_39_32_ROWS = {
        72'h000000000000000000,  // 7
        72'h0000000040ed291a88,  // 6
        72'h0000000020da949544,  // 5
        72'h0000000010b6724c22,  // 4
        72'h000000000871ce2311,  // 3
        72'h00000000040fc1e0f0,  // 2
        72'h0000000002003fe00f,  // 1
        72'h000000000100001fff  // 0
    };

    // codeword position of each data bit
    localparam logic [64*8-1:0] HSIAO_39_32_DATA_POS = {
        8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0,  // 63 .. 48
        8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0, 8'd0,  // 47 .. 32
        8'd31, 8'd30, 8'd29, 8'd28, 8'd27, 8'd26, 8'd25, 8'd24, 8'd23, 8'd22, 8'd21, 8'd20, 8'd19, 8'd18, 8'd17, 8'd16,  // 31 .. 16
        8'd15, 8'd14, 8'd13, 8'd12, 8'd11, 8'd10, 8'd9, 8'd8, 8'd7, 8'd6, 8'd5, 8'd4, 8'd3, 8'd2, 8'd1, 8'd0  // 15 .. 0
    };

    // codeword position of each check bit
    localparam logic [8*8-1:0] HSIAO_39_32_CHECK_POS = {
        8'd0, 8'd38, 8'd37, 8'd36, 8'd35, 8'd34, 8'd33, 8'd32  // 7 .. 0
    };

    // data bits whose XOR is each check bit
    localparam logic [8*64-1:0] HSIAO_39_32_CHECK_MASKS = {
        64'h0000000000000000,  // 7
        64'h00000000ed291a88,  // 6
        64'h00000000da949544,  // 5
        64'h00000000b6724c22,  // 4
        64'h0000000071ce2311,  // 3
        64'h000000000fc1e0f0,  // 2
        64'h00000000003fe00f,  // 1
        64'h0000000000001fff  // 0
    };
    // ---- one table, or one size, of the code a CODE number names

    function automatic int n_of(input int code);
        case (code)
            HSIAO: n_of = HSIAO_N;
            HAMMING: n_of = HAMMING_N;
            HSIAO_DICT: n_of = HSIAO_DICT_N;
            HSIAO_39_32: n_of = HSIAO_39_32_N;
            default: n_of = HSIAO_N;
        endcase
    endfunction

    function automatic int k_of(input int code);
        case (code)
            HSIAO: k_of = HSIAO_K;
            HAMMING: k_of = HAMMING_K;
            HSIAO_DICT: k_of = HSIAO_DICT_K;
            HSIAO_39_32: k_of = HSIAO_39_32_K;
            default: k_of = HSIAO_K;
        endcase
    endfunction

    function automatic int r_of(input int code);
        case (code)
            HSIAO: r_of = HSIAO_R;
            HAMMING: r_of = HAMMING_R;
            HSIAO_DICT: r_of = HSIAO_DICT_R;
            HSIAO_39_32: r_of = HSIAO_39_32_R;
            default: r_of = HSIAO_R;
        endcase
    endfunction

    function automatic int parity_row_of(input int code);
        case (code)
            HSIAO: parity_row_of = HSIAO_PARITY_ROW;
            HAMMING: parity_row_of = HAMMING_PARITY_ROW;
            HSIAO_DICT: parity_row_of = HSIAO_DICT_PARITY_ROW;
            HSIAO_39_32: parity_row_of = HSIAO_39_32_PARITY_ROW;
            default: parity_row_of = HSIAO_PARITY_ROW;
        endcase
    endfunction

    function automatic logic [MAX_N*MAX_R-1:0] columns_of(input int code);
        case (code)
            HSIAO: columns_of = HSIAO_COLUMNS;
            HAMMING: columns_of = HAMMING_COLUMNS;
            HSIAO_DICT: columns_of = HSIAO_DICT_COLUMNS;
            HSIAO_39_32: columns_of = HSIAO_39_32_COLUMNS;
            default: columns_of = HSIAO_COLUMNS;
        endcase
    endfunction

    function automatic logic [MAX_R*MAX_N-1:0] rows_of(input int code);
        case (code)
            HSIAO: rows_of = HSIAO_ROWS;
            HAMMING: rows_of = HAMMING_ROWS;
            HSIAO_DICT: rows_of = HSIAO_DICT_ROWS;
            HSIAO_39_32: rows_of = HSIAO_39_32_ROWS;
            default: rows_of = HSIAO_ROWS;
        endcase
    endfunction

    function automatic logic [MAX_K*8-1:0] data_pos_of(input int code);
        case (code)
            HSIAO: data_pos_of = HSIAO_DATA_POS;
            HAMMING: data_pos_of = HAMMING_DATA_POS;
            HSIAO_DICT: data_pos_of = HSIAO_DICT_DATA_POS;
            HSIAO_39_32: data_pos_of = HSIAO_39_32_DATA_POS;
            default: data_pos_of = HSIAO_DATA_POS;
        endcase
    endfunction

    function automatic logic [MAX_R*8-1:0] check_pos_of(input int code);
        case (code)
            HSIAO: check_pos_of = HSIAO_CHECK_POS;
            HAMMING: check_pos_of = HAMMING_CHECK_POS;
            HSIAO_DICT: check_pos_of = HSIAO_DICT_CHECK_POS;
            HSIAO_39_32: check_pos_of = HSIAO_39_32_CHECK_POS;
            default: check_pos_of = HSIAO_CHECK_POS;
        endcase
    endfunction

    function automatic logic [MAX_R*MAX_K-1:0] check_masks_of(input int code);
        case (code)
            HSIAO: check_masks_of = HSIAO_CHECK_MASKS;
            HAMMING: check_masks_of = HAMMING_CHECK_MASKS;
            HSIAO_DICT: check_masks_of = HSIAO_DICT_CHECK_MASKS;
            HSIAO_39_32: check_masks_of = HSIAO_39_32_CHECK_MASKS;
            default: check_masks_of = HSIAO_CHECK_MASKS;
        endcase
    endfunction

endpackage
