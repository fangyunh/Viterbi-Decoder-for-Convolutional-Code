// Viterbi Decoder for K=3, rate 1/2 convolutional code with soft decision decoding
// Parameters
parameter N = 32;           // Sequence length
parameter STATES = 4;       // Number of states (2^(K-1))
parameter Q_BITS = 3;       // Quantization bits
parameter BM_WIDTH = 4;     // Branch metric width (sum of two Q_BITS metrics)
parameter PM_WIDTH = 8;     // Path metric width

// Quantizer Module: Converts noisy 8-bit input to 3-bit soft decision value
module quantizer (
    input logic [7:0] noisy_in,  // 8-bit noisy input (e.g., fixed-point)
    output logic [Q_BITS-1:0] q_out  // 3-bit quantized output
);
    // Simple uniform quantization: map -256 to 255 range to 0-7
    // Assuming noisy_in is signed, -128 to 127 represents -1.5V to 1.5V
    always_comb begin
        if (noisy_in[7]) begin  // Negative values
            if (noisy_in >= 8'hC0)      q_out = 3'd0;  // < -0.75V
            else if (noisy_in >= 8'hE0) q_out = 3'd1;  // -0.75V to -0.5V
            else                        q_out = 3'd2;  // -0.5V to -0.25V
        end else begin          // Non-negative values
            if (noisy_in < 8'h20)       q_out = 3'd3;  // 0 to 0.25V
            else if (noisy_in < 8'h40)  q_out = 3'd4;  // 0.25V to 0.5V
            else if (noisy_in < 8'h60)  q_out = 3'd5;  // 0.5V to 0.75V
            else if (noisy_in < 8'h80)  q_out = 3'd6;  // 0.75V to 1.0V
            else                        q_out = 3'd7;  // > 1.0V
        end
    end
endmodule
//
//module tb_viterbi_decoder;
//
//    // Parameters
//    parameter N = 32;
//    parameter STATES = 4;
//    parameter Q_BITS = 3;
//    parameter BM_WIDTH = 4;
//    parameter PM_WIDTH = 8;
//
//    // Testbench signals
//    logic clk;
//    logic rst;
//    logic enable;
//    logic done;
//    logic [1:0] decisions [0:3];
//    logic [PM_WIDTH-1:0] new_pm [0:3];
//    logic [N-1:0] survivor_paths [0:3];
//    logic [PM_WIDTH-1:0] path_metrics [0:3];
//    logic [$clog2(N)-1:0] t;
//
//    // Intermediate signals for connecting modules
//    logic [7:0] noisy_in0, noisy_in1;              // Noisy input signals
//    logic [Q_BITS-1:0] q0, q1;                     // Quantized outputs
//    logic [BM_WIDTH-1:0] bm [0:7];                 // Branch metrics
//    logic [PM_WIDTH-1:0] pm [0:3];                 // Current path metrics
//
//    // Instantiate DUT (viterbi_decoder)
//    viterbi_decoder #(
//        .N(N),
//        .STATES(STATES),
//        .PM_WIDTH(PM_WIDTH)
//    ) dut (
//        .clk(clk),
//        .rst(rst),
//        .enable(enable),
//        .done(done),
//        .decisions(decisions),
//        .new_pm(new_pm),
//        .survivor_paths(survivor_paths),
//        .path_metrics(path_metrics),
//        .t(t)
//    );
//
//    // Instantiate supporting modules
//    quantizer quant0 (
//        .noisy_in(noisy_in0),
//        .q_out(q0)
//    );
//
//    quantizer quant1 (
//        .noisy_in(noisy_in1),
//        .q_out(q1)
//    );
//
//    bmc bmc_inst (
//        .q0(q0),
//        .q1(q1),
//        .bm(bm)
//    );
//
//    acs acs_inst (
//        .pm(path_metrics),
//        .bm(bm),
//        .new_pm(new_pm),
//        .decisions(decisions)
//    );
//
//    // Encoder function for K=3, rate 1/2, G0=[1,1,1], G1=[1,0,1]
//    function automatic void encode_sequence(input logic [31:0] seq, output logic [1:0] enc [0:31]);
//        logic [1:0] state = 2'b00;  // Initial state
//        for (int i = 0; i < 32; i++) begin
//            logic u = seq[31 - i];  // seq[31] is first bit (u0), seq[0] is last (u31)
//            logic p0 = u ^ state[1] ^ state[0];   // Parity bit 0
//            logic p1 = u ^ state[0];              // Parity bit 1
//            enc[i] = {p0, p1};                    // Store parity bits
//            state = {u, state[1]};                // Update state
//        end
//    endfunction
//
//    // Function to reverse bits of a 32-bit vector
//    function logic [31:0] reverse_bits(input logic [31:0] in);
//        for (int i = 0; i < 32; i++) begin
//            reverse_bits[i] = in[31 - i];
//        end
//    endfunction
//
//    // Clock generation (10ns period)
//    initial begin
//        clk = 0;
//        forever #5 clk = ~clk;
//    end
//
//    // Test sequences and control signals
//    logic [31:0] input_sequences [0:2];
//    logic [1:0] encoded_sequences [0:2][0:31];
//    logic [7:0] test_sequence [0:N-1][0:1];
//    logic [N-1:0] expected_path;
//    integer test_idx;
//
//    // Initialize test sequences
//    initial begin
//        // Define input sequences
//        input_sequences[0] = 32'hAAAAAAAA;  // Alternating 1s and 0s
//        input_sequences[1] = 32'hFFFF0000;  // First half 1s, second half 0s
//        input_sequences[2] = 32'h12345678;  // Pseudo-random pattern
//
//        // Encode all sequences
//        for (int i = 0; i < 3; i++) begin
//            encode_sequence(input_sequences[i], encoded_sequences[i]);
//        end
//
//        // Start with first test
//        test_idx = 0;
//        rst = 1;
//        enable = 0;
//        noisy_in0 = 8'h00;
//        noisy_in1 = 8'h00;
//        #10;
//        rst = 0;
//        #10;
//        run_test();
//    end
//
//    // Task to run a single test
//    task run_test();
//        for (int i = 0; i < N; i++) begin
//            logic [7:0] noise0;
//            logic [7:0] noise1;
//            noise0 = $random % 32 - 16;  // Noise range: -16 to 15
//            noise1 = $random % 32 - 16;
//            test_sequence[i][0] = (encoded_sequences[test_idx][i][0] ? 8'sh7F : 8'sh80) + noise0;
//            test_sequence[i][1] = (encoded_sequences[test_idx][i][1] ? 8'sh7F : 8'sh80) + noise1;
//        end
//        expected_path = input_sequences[test_idx];
//        enable = 1;
//    endtask
//
//    // Feed noisy inputs to decoder
//    always @(posedge clk) begin
//        if (enable && t < N) begin
//            noisy_in0 <= test_sequence[t][0];
//            noisy_in1 <= test_sequence[t][1];
//        end
//    end
//
//    // Check results and run next test
//    always @(posedge clk) begin
//        if (done) begin
//            logic [1:0] min_state;
//            logic [PM_WIDTH-1:0] min_pm;
//            logic [N-1:0] reversed_path;
//            min_state = 0;
//            min_pm = path_metrics[0];
//            for (int s = 1; s < STATES; s++) begin
//                if (path_metrics[s] < min_pm) begin
//                    min_pm = path_metrics[s];
//                    min_state = s;
//                end
//            end
//            reversed_path = reverse_bits(survivor_paths[min_state]);
//            $display("Test %0d:", test_idx);
//            if (reversed_path == expected_path) begin
//                $display("Time=%0t: Test Passed: Survivor path[%0d] = %h matches expected %h", 
//                         $time, min_state, reversed_path, expected_path);
//            end else begin
//                $display("Time=%0t: Test Failed: Survivor path[%0d] = %h, expected %h", 
//                         $time, min_state, reversed_path, expected_path);
//            end
//            // Dump all paths for debugging
//            for (int s = 0; s < STATES; s++) begin
//                logic [N-1:0] rev_path;
//                rev_path = reverse_bits(survivor_paths[s]);
//                $display("State %0d: Survivor Path = %h, Path Metric = %h", 
//                         s, rev_path, path_metrics[s]);
//            end
//
//            // Move to next test or finish
//            test_idx = test_idx + 1;
//            if (test_idx < 3) begin
//                enable = 0;
//                #20;
//                rst = 1;
//                #10;
//                rst = 0;
//                #10;
//                run_test();
//            end else begin
//                #10;
//                $finish;
//            end
//        end
//    end
//
//    // Simulation timeout
//    initial begin
//        #2000;
//        $display("Simulation timed out.");
//        $finish;
//    end
//
//endmodule


// Branch Metric Calculation (BMC) Module: Computes all branch metrics in parallel
module bmc (
    input logic [Q_BITS-1:0] q0,  // First quantized parity bit
    input logic [Q_BITS-1:0] q1,  // Second quantized parity bit
    output logic [BM_WIDTH-1:0] bm [0:7]  // 8 branch metrics
);
    // Bit metrics: for expected 0, metric = q; for expected 1, metric = 7 - q
    logic [BM_WIDTH-1:0] metric_0_q0, metric_1_q0, metric_0_q1, metric_1_q1;
    assign metric_0_q0 = {1'b0, q0};
    assign metric_1_q0 = 4'd7 - {1'b0, q0};
    assign metric_0_q1 = {1'b0, q1};
    assign metric_1_q1 = 4'd7 - {1'b0, q1};

    // Branch definitions based on standard K=3, G0=[1,1,1], G1=[1,0,1]
    // Transitions: state [s1,s0] with input u -> [u,s1]
    // Outputs: (u ⊕ s1 ⊕ s0, u ⊕ s0)
    // Branch 0: 00->00, u=0, out=00
    // Branch 1: 01->00, u=0, out=11
    // Branch 2: 00->10, u=1, out=11
    // Branch 3: 01->10, u=1, out=00
    // Branch 4: 10->01, u=0, out=10
    // Branch 5: 11->01, u=0, out=01
    // Branch 6: 10->11, u=1, out=01
    // Branch 7: 11->11, u=1, out=10
    assign bm[0] = metric_0_q0 + metric_0_q1;  // 00
    assign bm[1] = metric_1_q0 + metric_1_q1;  // 11
    assign bm[2] = metric_1_q0 + metric_1_q1;  // 11
    assign bm[3] = metric_0_q0 + metric_0_q1;  // 00
    assign bm[4] = metric_1_q0 + metric_0_q1;  // 10
    assign bm[5] = metric_0_q0 + metric_1_q1;  // 01
    assign bm[6] = metric_0_q0 + metric_1_q1;  // 01
    assign bm[7] = metric_1_q0 + metric_0_q1;  // 10
endmodule

// Add-Compare-Select (ACS) Module: Updates path metrics and selects predecessors in parallel
module acs (
    input logic [PM_WIDTH-1:0] pm [0:3],  // Current path metrics
    input logic [BM_WIDTH-1:0] bm [0:7],  // Branch metrics
    output logic [PM_WIDTH-1:0] new_pm [0:3],  // Updated path metrics
    output logic [1:0] decisions [0:3]  // Selected predecessor state index
);
    // State 00: from 00 (bm[0]), from 01 (bm[1])
    logic [PM_WIDTH-1:0] pm_00_to_00, pm_01_to_00;
    assign pm_00_to_00 = pm[0] + {4'b0, bm[0]};
    assign pm_01_to_00 = pm[1] + {4'b0, bm[1]};
    assign new_pm[0] = (pm_00_to_00 < pm_01_to_00) ? pm_00_to_00 : pm_01_to_00;
    assign decisions[0] = (pm_00_to_00 < pm_01_to_00) ? 2'd0 : 2'd1;

    // State 10: from 00 (bm[2]), from 01 (bm[3])
    logic [PM_WIDTH-1:0] pm_00_to_10, pm_01_to_10;
    assign pm_00_to_10 = pm[0] + {4'b0, bm[2]};
    assign pm_01_to_10 = pm[1] + {4'b0, bm[3]};
    assign new_pm[2] = (pm_00_to_10 < pm_01_to_10) ? pm_00_to_10 : pm_01_to_10;
    assign decisions[2] = (pm_00_to_10 < pm_01_to_10) ? 2'd0 : 2'd1;

    // State 01: from 10 (bm[4]), from 11 (bm[5])
    logic [PM_WIDTH-1:0] pm_10_to_01, pm_11_to_01;
    assign pm_10_to_01 = pm[2] + {4'b0, bm[4]};
    assign pm_11_to_01 = pm[3] + {4'b0, bm[5]};
    assign new_pm[1] = (pm_10_to_01 < pm_11_to_01) ? pm_10_to_01 : pm_11_to_01;
    assign decisions[1] = (pm_10_to_01 < pm_11_to_01) ? 2'd2 : 2'd3;

    // State 11: from 10 (bm[6]), from 11 (bm[7])
    logic [PM_WIDTH-1:0] pm_10_to_11, pm_11_to_11;
    assign pm_10_to_11 = pm[2] + {4'b0, bm[6]};
    assign pm_11_to_11 = pm[3] + {4'b0, bm[7]};
    assign new_pm[3] = (pm_10_to_11 < pm_11_to_11) ? pm_10_to_11 : pm_11_to_11;
    assign decisions[3] = (pm_10_to_11 < pm_11_to_11) ? 2'd2 : 2'd3;
endmodule

// Top Module: Viterbi Decoder
module viterbi_decoder (
    input logic clk,
    input logic rst,
    input logic enable,
    output logic done,
    // Other inputs/outputs omitted for brevity
    input logic [1:0] decisions [0:3], // Assuming STATES=4
    input logic [PM_WIDTH-1:0] new_pm [0:3],
    output logic [N-1:0] survivor_paths [0:3],
    output logic [PM_WIDTH-1:0] path_metrics [0:3],
    output logic [$clog2(N)-1:0] t
);
    parameter STATES = 4;
    parameter N = 32; // Example value, adjust as needed
    parameter PM_WIDTH = 8; // Example value, adjust as needed
	 
	 always_ff @(posedge clk or posedge rst) begin
		 if (rst) begin
			  for (int s = 0; s < STATES; s++) begin
					survivor_paths[s] <= 0;
					path_metrics[s] <= (s == 0) ? 0 : {PM_WIDTH{1'b1}};  // State 00 initial
			  end
			  t <= 0;
			  done <= 0;
		 end else if (enable && t < N) begin
			  path_metrics <= new_pm;
			  // Register exchange: Update survivor paths
			  for (int s = 0; s < STATES; s++) begin
					logic bit_to_append;  // Declare without initialization
					bit_to_append = (s >= 2) ? 1'b1 : 1'b0;  // Assign separately
					survivor_paths[s] <= {survivor_paths[decisions[s]][N-2:0], bit_to_append};
			  end
			  t <= t + 1;
			  if (t == N-1) done <= 1;
			  else done <= 0;
		 end else begin
			  done <= (t >= N);
		 end
	end

endmodule
