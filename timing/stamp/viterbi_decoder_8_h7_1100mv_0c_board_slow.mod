/*
 Copyright (C) 2020  Intel Corporation. All rights reserved.
 Your use of Intel Corporation's design tools, logic functions 
 and other software and tools, and any partner logic 
 functions, and any output files from any of the foregoing 
 (including device programming or simulation files), and any 
 associated documentation or information are expressly subject 
 to the terms and conditions of the Intel Program License 
 Subscription Agreement, the Intel Quartus Prime License Agreement,
 the Intel FPGA IP License Agreement, or other applicable license
 agreement, including, without limitation, that your use is for
 the sole purpose of programming logic devices manufactured by
 Intel and sold by Intel or its authorized distributors.  Please
 refer to the applicable agreement for further details, at
 https://fpgasoftware.intel.com/eula.
*/
MODEL
/*MODEL HEADER*/
/*
 This file contains Slow Corner delays for the design using part 5CGXFC7C7F23C8
 with speed grade 8_H7, core voltage 1.1V, and temperature 0 Celsius

*/
MODEL_VERSION "1.0";
DESIGN "viterbi_decoder";
DATE "04/17/2025 20:57:39";
PROGRAM "Quartus Prime";



INPUT clk;
INPUT rst;
INPUT enable;
INPUT decisions[3][0];
INPUT decisions[3][1];
INPUT decisions[2][1];
INPUT decisions[2][0];
INPUT decisions[0][0];
INPUT decisions[0][1];
INPUT decisions[1][1];
INPUT decisions[1][0];
INPUT new_pm[3][0];
INPUT new_pm[3][1];
INPUT new_pm[3][2];
INPUT new_pm[3][3];
INPUT new_pm[3][4];
INPUT new_pm[3][5];
INPUT new_pm[3][6];
INPUT new_pm[3][7];
INPUT new_pm[2][0];
INPUT new_pm[2][1];
INPUT new_pm[2][2];
INPUT new_pm[2][3];
INPUT new_pm[2][4];
INPUT new_pm[2][5];
INPUT new_pm[2][6];
INPUT new_pm[2][7];
INPUT new_pm[1][0];
INPUT new_pm[1][1];
INPUT new_pm[1][2];
INPUT new_pm[1][3];
INPUT new_pm[1][4];
INPUT new_pm[1][5];
INPUT new_pm[1][6];
INPUT new_pm[1][7];
INPUT new_pm[0][0];
INPUT new_pm[0][1];
INPUT new_pm[0][2];
INPUT new_pm[0][3];
INPUT new_pm[0][4];
INPUT new_pm[0][5];
INPUT new_pm[0][6];
INPUT new_pm[0][7];
OUTPUT done;
OUTPUT survivor_paths[3][0];
OUTPUT survivor_paths[3][1];
OUTPUT survivor_paths[3][2];
OUTPUT survivor_paths[3][3];
OUTPUT survivor_paths[3][4];
OUTPUT survivor_paths[3][5];
OUTPUT survivor_paths[3][6];
OUTPUT survivor_paths[3][7];
OUTPUT survivor_paths[3][8];
OUTPUT survivor_paths[3][9];
OUTPUT survivor_paths[3][10];
OUTPUT survivor_paths[3][11];
OUTPUT survivor_paths[3][12];
OUTPUT survivor_paths[3][13];
OUTPUT survivor_paths[3][14];
OUTPUT survivor_paths[3][15];
OUTPUT survivor_paths[3][16];
OUTPUT survivor_paths[3][17];
OUTPUT survivor_paths[3][18];
OUTPUT survivor_paths[3][19];
OUTPUT survivor_paths[3][20];
OUTPUT survivor_paths[3][21];
OUTPUT survivor_paths[3][22];
OUTPUT survivor_paths[3][23];
OUTPUT survivor_paths[3][24];
OUTPUT survivor_paths[3][25];
OUTPUT survivor_paths[3][26];
OUTPUT survivor_paths[3][27];
OUTPUT survivor_paths[3][28];
OUTPUT survivor_paths[3][29];
OUTPUT survivor_paths[3][30];
OUTPUT survivor_paths[3][31];
OUTPUT survivor_paths[2][0];
OUTPUT survivor_paths[2][1];
OUTPUT survivor_paths[2][2];
OUTPUT survivor_paths[2][3];
OUTPUT survivor_paths[2][4];
OUTPUT survivor_paths[2][5];
OUTPUT survivor_paths[2][6];
OUTPUT survivor_paths[2][7];
OUTPUT survivor_paths[2][8];
OUTPUT survivor_paths[2][9];
OUTPUT survivor_paths[2][10];
OUTPUT survivor_paths[2][11];
OUTPUT survivor_paths[2][12];
OUTPUT survivor_paths[2][13];
OUTPUT survivor_paths[2][14];
OUTPUT survivor_paths[2][15];
OUTPUT survivor_paths[2][16];
OUTPUT survivor_paths[2][17];
OUTPUT survivor_paths[2][18];
OUTPUT survivor_paths[2][19];
OUTPUT survivor_paths[2][20];
OUTPUT survivor_paths[2][21];
OUTPUT survivor_paths[2][22];
OUTPUT survivor_paths[2][23];
OUTPUT survivor_paths[2][24];
OUTPUT survivor_paths[2][25];
OUTPUT survivor_paths[2][26];
OUTPUT survivor_paths[2][27];
OUTPUT survivor_paths[2][28];
OUTPUT survivor_paths[2][29];
OUTPUT survivor_paths[2][30];
OUTPUT survivor_paths[2][31];
OUTPUT survivor_paths[1][0];
OUTPUT survivor_paths[1][1];
OUTPUT survivor_paths[1][2];
OUTPUT survivor_paths[1][3];
OUTPUT survivor_paths[1][4];
OUTPUT survivor_paths[1][5];
OUTPUT survivor_paths[1][6];
OUTPUT survivor_paths[1][7];
OUTPUT survivor_paths[1][8];
OUTPUT survivor_paths[1][9];
OUTPUT survivor_paths[1][10];
OUTPUT survivor_paths[1][11];
OUTPUT survivor_paths[1][12];
OUTPUT survivor_paths[1][13];
OUTPUT survivor_paths[1][14];
OUTPUT survivor_paths[1][15];
OUTPUT survivor_paths[1][16];
OUTPUT survivor_paths[1][17];
OUTPUT survivor_paths[1][18];
OUTPUT survivor_paths[1][19];
OUTPUT survivor_paths[1][20];
OUTPUT survivor_paths[1][21];
OUTPUT survivor_paths[1][22];
OUTPUT survivor_paths[1][23];
OUTPUT survivor_paths[1][24];
OUTPUT survivor_paths[1][25];
OUTPUT survivor_paths[1][26];
OUTPUT survivor_paths[1][27];
OUTPUT survivor_paths[1][28];
OUTPUT survivor_paths[1][29];
OUTPUT survivor_paths[1][30];
OUTPUT survivor_paths[1][31];
OUTPUT survivor_paths[0][0];
OUTPUT survivor_paths[0][1];
OUTPUT survivor_paths[0][2];
OUTPUT survivor_paths[0][3];
OUTPUT survivor_paths[0][4];
OUTPUT survivor_paths[0][5];
OUTPUT survivor_paths[0][6];
OUTPUT survivor_paths[0][7];
OUTPUT survivor_paths[0][8];
OUTPUT survivor_paths[0][9];
OUTPUT survivor_paths[0][10];
OUTPUT survivor_paths[0][11];
OUTPUT survivor_paths[0][12];
OUTPUT survivor_paths[0][13];
OUTPUT survivor_paths[0][14];
OUTPUT survivor_paths[0][15];
OUTPUT survivor_paths[0][16];
OUTPUT survivor_paths[0][17];
OUTPUT survivor_paths[0][18];
OUTPUT survivor_paths[0][19];
OUTPUT survivor_paths[0][20];
OUTPUT survivor_paths[0][21];
OUTPUT survivor_paths[0][22];
OUTPUT survivor_paths[0][23];
OUTPUT survivor_paths[0][24];
OUTPUT survivor_paths[0][25];
OUTPUT survivor_paths[0][26];
OUTPUT survivor_paths[0][27];
OUTPUT survivor_paths[0][28];
OUTPUT survivor_paths[0][29];
OUTPUT survivor_paths[0][30];
OUTPUT survivor_paths[0][31];
OUTPUT path_metrics[3][0];
OUTPUT path_metrics[3][1];
OUTPUT path_metrics[3][2];
OUTPUT path_metrics[3][3];
OUTPUT path_metrics[3][4];
OUTPUT path_metrics[3][5];
OUTPUT path_metrics[3][6];
OUTPUT path_metrics[3][7];
OUTPUT path_metrics[2][0];
OUTPUT path_metrics[2][1];
OUTPUT path_metrics[2][2];
OUTPUT path_metrics[2][3];
OUTPUT path_metrics[2][4];
OUTPUT path_metrics[2][5];
OUTPUT path_metrics[2][6];
OUTPUT path_metrics[2][7];
OUTPUT path_metrics[1][0];
OUTPUT path_metrics[1][1];
OUTPUT path_metrics[1][2];
OUTPUT path_metrics[1][3];
OUTPUT path_metrics[1][4];
OUTPUT path_metrics[1][5];
OUTPUT path_metrics[1][6];
OUTPUT path_metrics[1][7];
OUTPUT path_metrics[0][0];
OUTPUT path_metrics[0][1];
OUTPUT path_metrics[0][2];
OUTPUT path_metrics[0][3];
OUTPUT path_metrics[0][4];
OUTPUT path_metrics[0][5];
OUTPUT path_metrics[0][6];
OUTPUT path_metrics[0][7];
OUTPUT t[0];
OUTPUT t[1];
OUTPUT t[2];
OUTPUT t[3];
OUTPUT t[4];

/*Arc definitions start here*/
pos_decisions[0][0]__clk__setup:		SETUP (POSEDGE) decisions[0][0] clk ;
pos_decisions[0][1]__clk__setup:		SETUP (POSEDGE) decisions[0][1] clk ;
pos_decisions[1][0]__clk__setup:		SETUP (POSEDGE) decisions[1][0] clk ;
pos_decisions[1][1]__clk__setup:		SETUP (POSEDGE) decisions[1][1] clk ;
pos_decisions[2][0]__clk__setup:		SETUP (POSEDGE) decisions[2][0] clk ;
pos_decisions[2][1]__clk__setup:		SETUP (POSEDGE) decisions[2][1] clk ;
pos_decisions[3][0]__clk__setup:		SETUP (POSEDGE) decisions[3][0] clk ;
pos_decisions[3][1]__clk__setup:		SETUP (POSEDGE) decisions[3][1] clk ;
pos_enable__clk__setup:		SETUP (POSEDGE) enable clk ;
pos_new_pm[0][0]__clk__setup:		SETUP (POSEDGE) new_pm[0][0] clk ;
pos_new_pm[0][1]__clk__setup:		SETUP (POSEDGE) new_pm[0][1] clk ;
pos_new_pm[0][2]__clk__setup:		SETUP (POSEDGE) new_pm[0][2] clk ;
pos_new_pm[0][3]__clk__setup:		SETUP (POSEDGE) new_pm[0][3] clk ;
pos_new_pm[0][4]__clk__setup:		SETUP (POSEDGE) new_pm[0][4] clk ;
pos_new_pm[0][5]__clk__setup:		SETUP (POSEDGE) new_pm[0][5] clk ;
pos_new_pm[0][6]__clk__setup:		SETUP (POSEDGE) new_pm[0][6] clk ;
pos_new_pm[0][7]__clk__setup:		SETUP (POSEDGE) new_pm[0][7] clk ;
pos_new_pm[1][0]__clk__setup:		SETUP (POSEDGE) new_pm[1][0] clk ;
pos_new_pm[1][1]__clk__setup:		SETUP (POSEDGE) new_pm[1][1] clk ;
pos_new_pm[1][2]__clk__setup:		SETUP (POSEDGE) new_pm[1][2] clk ;
pos_new_pm[1][3]__clk__setup:		SETUP (POSEDGE) new_pm[1][3] clk ;
pos_new_pm[1][4]__clk__setup:		SETUP (POSEDGE) new_pm[1][4] clk ;
pos_new_pm[1][5]__clk__setup:		SETUP (POSEDGE) new_pm[1][5] clk ;
pos_new_pm[1][6]__clk__setup:		SETUP (POSEDGE) new_pm[1][6] clk ;
pos_new_pm[1][7]__clk__setup:		SETUP (POSEDGE) new_pm[1][7] clk ;
pos_new_pm[2][0]__clk__setup:		SETUP (POSEDGE) new_pm[2][0] clk ;
pos_new_pm[2][1]__clk__setup:		SETUP (POSEDGE) new_pm[2][1] clk ;
pos_new_pm[2][2]__clk__setup:		SETUP (POSEDGE) new_pm[2][2] clk ;
pos_new_pm[2][3]__clk__setup:		SETUP (POSEDGE) new_pm[2][3] clk ;
pos_new_pm[2][4]__clk__setup:		SETUP (POSEDGE) new_pm[2][4] clk ;
pos_new_pm[2][5]__clk__setup:		SETUP (POSEDGE) new_pm[2][5] clk ;
pos_new_pm[2][6]__clk__setup:		SETUP (POSEDGE) new_pm[2][6] clk ;
pos_new_pm[2][7]__clk__setup:		SETUP (POSEDGE) new_pm[2][7] clk ;
pos_new_pm[3][0]__clk__setup:		SETUP (POSEDGE) new_pm[3][0] clk ;
pos_new_pm[3][1]__clk__setup:		SETUP (POSEDGE) new_pm[3][1] clk ;
pos_new_pm[3][2]__clk__setup:		SETUP (POSEDGE) new_pm[3][2] clk ;
pos_new_pm[3][3]__clk__setup:		SETUP (POSEDGE) new_pm[3][3] clk ;
pos_new_pm[3][4]__clk__setup:		SETUP (POSEDGE) new_pm[3][4] clk ;
pos_new_pm[3][5]__clk__setup:		SETUP (POSEDGE) new_pm[3][5] clk ;
pos_new_pm[3][6]__clk__setup:		SETUP (POSEDGE) new_pm[3][6] clk ;
pos_new_pm[3][7]__clk__setup:		SETUP (POSEDGE) new_pm[3][7] clk ;
pos_decisions[0][0]__clk__hold:		HOLD (POSEDGE) decisions[0][0] clk ;
pos_decisions[0][1]__clk__hold:		HOLD (POSEDGE) decisions[0][1] clk ;
pos_decisions[1][0]__clk__hold:		HOLD (POSEDGE) decisions[1][0] clk ;
pos_decisions[1][1]__clk__hold:		HOLD (POSEDGE) decisions[1][1] clk ;
pos_decisions[2][0]__clk__hold:		HOLD (POSEDGE) decisions[2][0] clk ;
pos_decisions[2][1]__clk__hold:		HOLD (POSEDGE) decisions[2][1] clk ;
pos_decisions[3][0]__clk__hold:		HOLD (POSEDGE) decisions[3][0] clk ;
pos_decisions[3][1]__clk__hold:		HOLD (POSEDGE) decisions[3][1] clk ;
pos_enable__clk__hold:		HOLD (POSEDGE) enable clk ;
pos_new_pm[0][0]__clk__hold:		HOLD (POSEDGE) new_pm[0][0] clk ;
pos_new_pm[0][1]__clk__hold:		HOLD (POSEDGE) new_pm[0][1] clk ;
pos_new_pm[0][2]__clk__hold:		HOLD (POSEDGE) new_pm[0][2] clk ;
pos_new_pm[0][3]__clk__hold:		HOLD (POSEDGE) new_pm[0][3] clk ;
pos_new_pm[0][4]__clk__hold:		HOLD (POSEDGE) new_pm[0][4] clk ;
pos_new_pm[0][5]__clk__hold:		HOLD (POSEDGE) new_pm[0][5] clk ;
pos_new_pm[0][6]__clk__hold:		HOLD (POSEDGE) new_pm[0][6] clk ;
pos_new_pm[0][7]__clk__hold:		HOLD (POSEDGE) new_pm[0][7] clk ;
pos_new_pm[1][0]__clk__hold:		HOLD (POSEDGE) new_pm[1][0] clk ;
pos_new_pm[1][1]__clk__hold:		HOLD (POSEDGE) new_pm[1][1] clk ;
pos_new_pm[1][2]__clk__hold:		HOLD (POSEDGE) new_pm[1][2] clk ;
pos_new_pm[1][3]__clk__hold:		HOLD (POSEDGE) new_pm[1][3] clk ;
pos_new_pm[1][4]__clk__hold:		HOLD (POSEDGE) new_pm[1][4] clk ;
pos_new_pm[1][5]__clk__hold:		HOLD (POSEDGE) new_pm[1][5] clk ;
pos_new_pm[1][6]__clk__hold:		HOLD (POSEDGE) new_pm[1][6] clk ;
pos_new_pm[1][7]__clk__hold:		HOLD (POSEDGE) new_pm[1][7] clk ;
pos_new_pm[2][0]__clk__hold:		HOLD (POSEDGE) new_pm[2][0] clk ;
pos_new_pm[2][1]__clk__hold:		HOLD (POSEDGE) new_pm[2][1] clk ;
pos_new_pm[2][2]__clk__hold:		HOLD (POSEDGE) new_pm[2][2] clk ;
pos_new_pm[2][3]__clk__hold:		HOLD (POSEDGE) new_pm[2][3] clk ;
pos_new_pm[2][4]__clk__hold:		HOLD (POSEDGE) new_pm[2][4] clk ;
pos_new_pm[2][5]__clk__hold:		HOLD (POSEDGE) new_pm[2][5] clk ;
pos_new_pm[2][6]__clk__hold:		HOLD (POSEDGE) new_pm[2][6] clk ;
pos_new_pm[2][7]__clk__hold:		HOLD (POSEDGE) new_pm[2][7] clk ;
pos_new_pm[3][0]__clk__hold:		HOLD (POSEDGE) new_pm[3][0] clk ;
pos_new_pm[3][1]__clk__hold:		HOLD (POSEDGE) new_pm[3][1] clk ;
pos_new_pm[3][2]__clk__hold:		HOLD (POSEDGE) new_pm[3][2] clk ;
pos_new_pm[3][3]__clk__hold:		HOLD (POSEDGE) new_pm[3][3] clk ;
pos_new_pm[3][4]__clk__hold:		HOLD (POSEDGE) new_pm[3][4] clk ;
pos_new_pm[3][5]__clk__hold:		HOLD (POSEDGE) new_pm[3][5] clk ;
pos_new_pm[3][6]__clk__hold:		HOLD (POSEDGE) new_pm[3][6] clk ;
pos_new_pm[3][7]__clk__hold:		HOLD (POSEDGE) new_pm[3][7] clk ;
pos_clk__done__delay:		DELAY (POSEDGE) clk done ;
pos_clk__path_metrics[0][0]__delay:		DELAY (POSEDGE) clk path_metrics[0][0] ;
pos_clk__path_metrics[0][1]__delay:		DELAY (POSEDGE) clk path_metrics[0][1] ;
pos_clk__path_metrics[0][2]__delay:		DELAY (POSEDGE) clk path_metrics[0][2] ;
pos_clk__path_metrics[0][3]__delay:		DELAY (POSEDGE) clk path_metrics[0][3] ;
pos_clk__path_metrics[0][4]__delay:		DELAY (POSEDGE) clk path_metrics[0][4] ;
pos_clk__path_metrics[0][5]__delay:		DELAY (POSEDGE) clk path_metrics[0][5] ;
pos_clk__path_metrics[0][6]__delay:		DELAY (POSEDGE) clk path_metrics[0][6] ;
pos_clk__path_metrics[0][7]__delay:		DELAY (POSEDGE) clk path_metrics[0][7] ;
pos_clk__path_metrics[1][0]__delay:		DELAY (POSEDGE) clk path_metrics[1][0] ;
pos_clk__path_metrics[1][1]__delay:		DELAY (POSEDGE) clk path_metrics[1][1] ;
pos_clk__path_metrics[1][2]__delay:		DELAY (POSEDGE) clk path_metrics[1][2] ;
pos_clk__path_metrics[1][3]__delay:		DELAY (POSEDGE) clk path_metrics[1][3] ;
pos_clk__path_metrics[1][4]__delay:		DELAY (POSEDGE) clk path_metrics[1][4] ;
pos_clk__path_metrics[1][5]__delay:		DELAY (POSEDGE) clk path_metrics[1][5] ;
pos_clk__path_metrics[1][6]__delay:		DELAY (POSEDGE) clk path_metrics[1][6] ;
pos_clk__path_metrics[1][7]__delay:		DELAY (POSEDGE) clk path_metrics[1][7] ;
pos_clk__path_metrics[2][0]__delay:		DELAY (POSEDGE) clk path_metrics[2][0] ;
pos_clk__path_metrics[2][1]__delay:		DELAY (POSEDGE) clk path_metrics[2][1] ;
pos_clk__path_metrics[2][2]__delay:		DELAY (POSEDGE) clk path_metrics[2][2] ;
pos_clk__path_metrics[2][3]__delay:		DELAY (POSEDGE) clk path_metrics[2][3] ;
pos_clk__path_metrics[2][4]__delay:		DELAY (POSEDGE) clk path_metrics[2][4] ;
pos_clk__path_metrics[2][5]__delay:		DELAY (POSEDGE) clk path_metrics[2][5] ;
pos_clk__path_metrics[2][6]__delay:		DELAY (POSEDGE) clk path_metrics[2][6] ;
pos_clk__path_metrics[2][7]__delay:		DELAY (POSEDGE) clk path_metrics[2][7] ;
pos_clk__path_metrics[3][0]__delay:		DELAY (POSEDGE) clk path_metrics[3][0] ;
pos_clk__path_metrics[3][1]__delay:		DELAY (POSEDGE) clk path_metrics[3][1] ;
pos_clk__path_metrics[3][2]__delay:		DELAY (POSEDGE) clk path_metrics[3][2] ;
pos_clk__path_metrics[3][3]__delay:		DELAY (POSEDGE) clk path_metrics[3][3] ;
pos_clk__path_metrics[3][4]__delay:		DELAY (POSEDGE) clk path_metrics[3][4] ;
pos_clk__path_metrics[3][5]__delay:		DELAY (POSEDGE) clk path_metrics[3][5] ;
pos_clk__path_metrics[3][6]__delay:		DELAY (POSEDGE) clk path_metrics[3][6] ;
pos_clk__path_metrics[3][7]__delay:		DELAY (POSEDGE) clk path_metrics[3][7] ;
pos_clk__survivor_paths[0][0]__delay:		DELAY (POSEDGE) clk survivor_paths[0][0] ;
pos_clk__survivor_paths[0][1]__delay:		DELAY (POSEDGE) clk survivor_paths[0][1] ;
pos_clk__survivor_paths[0][2]__delay:		DELAY (POSEDGE) clk survivor_paths[0][2] ;
pos_clk__survivor_paths[0][3]__delay:		DELAY (POSEDGE) clk survivor_paths[0][3] ;
pos_clk__survivor_paths[0][4]__delay:		DELAY (POSEDGE) clk survivor_paths[0][4] ;
pos_clk__survivor_paths[0][5]__delay:		DELAY (POSEDGE) clk survivor_paths[0][5] ;
pos_clk__survivor_paths[0][6]__delay:		DELAY (POSEDGE) clk survivor_paths[0][6] ;
pos_clk__survivor_paths[0][7]__delay:		DELAY (POSEDGE) clk survivor_paths[0][7] ;
pos_clk__survivor_paths[0][8]__delay:		DELAY (POSEDGE) clk survivor_paths[0][8] ;
pos_clk__survivor_paths[0][9]__delay:		DELAY (POSEDGE) clk survivor_paths[0][9] ;
pos_clk__survivor_paths[0][10]__delay:		DELAY (POSEDGE) clk survivor_paths[0][10] ;
pos_clk__survivor_paths[0][11]__delay:		DELAY (POSEDGE) clk survivor_paths[0][11] ;
pos_clk__survivor_paths[0][12]__delay:		DELAY (POSEDGE) clk survivor_paths[0][12] ;
pos_clk__survivor_paths[0][13]__delay:		DELAY (POSEDGE) clk survivor_paths[0][13] ;
pos_clk__survivor_paths[0][14]__delay:		DELAY (POSEDGE) clk survivor_paths[0][14] ;
pos_clk__survivor_paths[0][15]__delay:		DELAY (POSEDGE) clk survivor_paths[0][15] ;
pos_clk__survivor_paths[0][16]__delay:		DELAY (POSEDGE) clk survivor_paths[0][16] ;
pos_clk__survivor_paths[0][17]__delay:		DELAY (POSEDGE) clk survivor_paths[0][17] ;
pos_clk__survivor_paths[0][18]__delay:		DELAY (POSEDGE) clk survivor_paths[0][18] ;
pos_clk__survivor_paths[0][19]__delay:		DELAY (POSEDGE) clk survivor_paths[0][19] ;
pos_clk__survivor_paths[0][20]__delay:		DELAY (POSEDGE) clk survivor_paths[0][20] ;
pos_clk__survivor_paths[0][21]__delay:		DELAY (POSEDGE) clk survivor_paths[0][21] ;
pos_clk__survivor_paths[0][22]__delay:		DELAY (POSEDGE) clk survivor_paths[0][22] ;
pos_clk__survivor_paths[0][23]__delay:		DELAY (POSEDGE) clk survivor_paths[0][23] ;
pos_clk__survivor_paths[0][24]__delay:		DELAY (POSEDGE) clk survivor_paths[0][24] ;
pos_clk__survivor_paths[0][25]__delay:		DELAY (POSEDGE) clk survivor_paths[0][25] ;
pos_clk__survivor_paths[0][26]__delay:		DELAY (POSEDGE) clk survivor_paths[0][26] ;
pos_clk__survivor_paths[0][27]__delay:		DELAY (POSEDGE) clk survivor_paths[0][27] ;
pos_clk__survivor_paths[0][28]__delay:		DELAY (POSEDGE) clk survivor_paths[0][28] ;
pos_clk__survivor_paths[0][29]__delay:		DELAY (POSEDGE) clk survivor_paths[0][29] ;
pos_clk__survivor_paths[0][30]__delay:		DELAY (POSEDGE) clk survivor_paths[0][30] ;
pos_clk__survivor_paths[0][31]__delay:		DELAY (POSEDGE) clk survivor_paths[0][31] ;
pos_clk__survivor_paths[1][0]__delay:		DELAY (POSEDGE) clk survivor_paths[1][0] ;
pos_clk__survivor_paths[1][1]__delay:		DELAY (POSEDGE) clk survivor_paths[1][1] ;
pos_clk__survivor_paths[1][2]__delay:		DELAY (POSEDGE) clk survivor_paths[1][2] ;
pos_clk__survivor_paths[1][3]__delay:		DELAY (POSEDGE) clk survivor_paths[1][3] ;
pos_clk__survivor_paths[1][4]__delay:		DELAY (POSEDGE) clk survivor_paths[1][4] ;
pos_clk__survivor_paths[1][5]__delay:		DELAY (POSEDGE) clk survivor_paths[1][5] ;
pos_clk__survivor_paths[1][6]__delay:		DELAY (POSEDGE) clk survivor_paths[1][6] ;
pos_clk__survivor_paths[1][7]__delay:		DELAY (POSEDGE) clk survivor_paths[1][7] ;
pos_clk__survivor_paths[1][8]__delay:		DELAY (POSEDGE) clk survivor_paths[1][8] ;
pos_clk__survivor_paths[1][9]__delay:		DELAY (POSEDGE) clk survivor_paths[1][9] ;
pos_clk__survivor_paths[1][10]__delay:		DELAY (POSEDGE) clk survivor_paths[1][10] ;
pos_clk__survivor_paths[1][11]__delay:		DELAY (POSEDGE) clk survivor_paths[1][11] ;
pos_clk__survivor_paths[1][12]__delay:		DELAY (POSEDGE) clk survivor_paths[1][12] ;
pos_clk__survivor_paths[1][13]__delay:		DELAY (POSEDGE) clk survivor_paths[1][13] ;
pos_clk__survivor_paths[1][14]__delay:		DELAY (POSEDGE) clk survivor_paths[1][14] ;
pos_clk__survivor_paths[1][15]__delay:		DELAY (POSEDGE) clk survivor_paths[1][15] ;
pos_clk__survivor_paths[1][16]__delay:		DELAY (POSEDGE) clk survivor_paths[1][16] ;
pos_clk__survivor_paths[1][17]__delay:		DELAY (POSEDGE) clk survivor_paths[1][17] ;
pos_clk__survivor_paths[1][18]__delay:		DELAY (POSEDGE) clk survivor_paths[1][18] ;
pos_clk__survivor_paths[1][19]__delay:		DELAY (POSEDGE) clk survivor_paths[1][19] ;
pos_clk__survivor_paths[1][20]__delay:		DELAY (POSEDGE) clk survivor_paths[1][20] ;
pos_clk__survivor_paths[1][21]__delay:		DELAY (POSEDGE) clk survivor_paths[1][21] ;
pos_clk__survivor_paths[1][22]__delay:		DELAY (POSEDGE) clk survivor_paths[1][22] ;
pos_clk__survivor_paths[1][23]__delay:		DELAY (POSEDGE) clk survivor_paths[1][23] ;
pos_clk__survivor_paths[1][24]__delay:		DELAY (POSEDGE) clk survivor_paths[1][24] ;
pos_clk__survivor_paths[1][25]__delay:		DELAY (POSEDGE) clk survivor_paths[1][25] ;
pos_clk__survivor_paths[1][26]__delay:		DELAY (POSEDGE) clk survivor_paths[1][26] ;
pos_clk__survivor_paths[1][27]__delay:		DELAY (POSEDGE) clk survivor_paths[1][27] ;
pos_clk__survivor_paths[1][28]__delay:		DELAY (POSEDGE) clk survivor_paths[1][28] ;
pos_clk__survivor_paths[1][29]__delay:		DELAY (POSEDGE) clk survivor_paths[1][29] ;
pos_clk__survivor_paths[1][30]__delay:		DELAY (POSEDGE) clk survivor_paths[1][30] ;
pos_clk__survivor_paths[1][31]__delay:		DELAY (POSEDGE) clk survivor_paths[1][31] ;
pos_clk__survivor_paths[2][0]__delay:		DELAY (POSEDGE) clk survivor_paths[2][0] ;
pos_clk__survivor_paths[2][1]__delay:		DELAY (POSEDGE) clk survivor_paths[2][1] ;
pos_clk__survivor_paths[2][2]__delay:		DELAY (POSEDGE) clk survivor_paths[2][2] ;
pos_clk__survivor_paths[2][3]__delay:		DELAY (POSEDGE) clk survivor_paths[2][3] ;
pos_clk__survivor_paths[2][4]__delay:		DELAY (POSEDGE) clk survivor_paths[2][4] ;
pos_clk__survivor_paths[2][5]__delay:		DELAY (POSEDGE) clk survivor_paths[2][5] ;
pos_clk__survivor_paths[2][6]__delay:		DELAY (POSEDGE) clk survivor_paths[2][6] ;
pos_clk__survivor_paths[2][7]__delay:		DELAY (POSEDGE) clk survivor_paths[2][7] ;
pos_clk__survivor_paths[2][8]__delay:		DELAY (POSEDGE) clk survivor_paths[2][8] ;
pos_clk__survivor_paths[2][9]__delay:		DELAY (POSEDGE) clk survivor_paths[2][9] ;
pos_clk__survivor_paths[2][10]__delay:		DELAY (POSEDGE) clk survivor_paths[2][10] ;
pos_clk__survivor_paths[2][11]__delay:		DELAY (POSEDGE) clk survivor_paths[2][11] ;
pos_clk__survivor_paths[2][12]__delay:		DELAY (POSEDGE) clk survivor_paths[2][12] ;
pos_clk__survivor_paths[2][13]__delay:		DELAY (POSEDGE) clk survivor_paths[2][13] ;
pos_clk__survivor_paths[2][14]__delay:		DELAY (POSEDGE) clk survivor_paths[2][14] ;
pos_clk__survivor_paths[2][15]__delay:		DELAY (POSEDGE) clk survivor_paths[2][15] ;
pos_clk__survivor_paths[2][16]__delay:		DELAY (POSEDGE) clk survivor_paths[2][16] ;
pos_clk__survivor_paths[2][17]__delay:		DELAY (POSEDGE) clk survivor_paths[2][17] ;
pos_clk__survivor_paths[2][18]__delay:		DELAY (POSEDGE) clk survivor_paths[2][18] ;
pos_clk__survivor_paths[2][19]__delay:		DELAY (POSEDGE) clk survivor_paths[2][19] ;
pos_clk__survivor_paths[2][20]__delay:		DELAY (POSEDGE) clk survivor_paths[2][20] ;
pos_clk__survivor_paths[2][21]__delay:		DELAY (POSEDGE) clk survivor_paths[2][21] ;
pos_clk__survivor_paths[2][22]__delay:		DELAY (POSEDGE) clk survivor_paths[2][22] ;
pos_clk__survivor_paths[2][23]__delay:		DELAY (POSEDGE) clk survivor_paths[2][23] ;
pos_clk__survivor_paths[2][24]__delay:		DELAY (POSEDGE) clk survivor_paths[2][24] ;
pos_clk__survivor_paths[2][25]__delay:		DELAY (POSEDGE) clk survivor_paths[2][25] ;
pos_clk__survivor_paths[2][26]__delay:		DELAY (POSEDGE) clk survivor_paths[2][26] ;
pos_clk__survivor_paths[2][27]__delay:		DELAY (POSEDGE) clk survivor_paths[2][27] ;
pos_clk__survivor_paths[2][28]__delay:		DELAY (POSEDGE) clk survivor_paths[2][28] ;
pos_clk__survivor_paths[2][29]__delay:		DELAY (POSEDGE) clk survivor_paths[2][29] ;
pos_clk__survivor_paths[2][30]__delay:		DELAY (POSEDGE) clk survivor_paths[2][30] ;
pos_clk__survivor_paths[2][31]__delay:		DELAY (POSEDGE) clk survivor_paths[2][31] ;
pos_clk__survivor_paths[3][0]__delay:		DELAY (POSEDGE) clk survivor_paths[3][0] ;
pos_clk__survivor_paths[3][1]__delay:		DELAY (POSEDGE) clk survivor_paths[3][1] ;
pos_clk__survivor_paths[3][2]__delay:		DELAY (POSEDGE) clk survivor_paths[3][2] ;
pos_clk__survivor_paths[3][3]__delay:		DELAY (POSEDGE) clk survivor_paths[3][3] ;
pos_clk__survivor_paths[3][4]__delay:		DELAY (POSEDGE) clk survivor_paths[3][4] ;
pos_clk__survivor_paths[3][5]__delay:		DELAY (POSEDGE) clk survivor_paths[3][5] ;
pos_clk__survivor_paths[3][6]__delay:		DELAY (POSEDGE) clk survivor_paths[3][6] ;
pos_clk__survivor_paths[3][7]__delay:		DELAY (POSEDGE) clk survivor_paths[3][7] ;
pos_clk__survivor_paths[3][8]__delay:		DELAY (POSEDGE) clk survivor_paths[3][8] ;
pos_clk__survivor_paths[3][9]__delay:		DELAY (POSEDGE) clk survivor_paths[3][9] ;
pos_clk__survivor_paths[3][10]__delay:		DELAY (POSEDGE) clk survivor_paths[3][10] ;
pos_clk__survivor_paths[3][11]__delay:		DELAY (POSEDGE) clk survivor_paths[3][11] ;
pos_clk__survivor_paths[3][12]__delay:		DELAY (POSEDGE) clk survivor_paths[3][12] ;
pos_clk__survivor_paths[3][13]__delay:		DELAY (POSEDGE) clk survivor_paths[3][13] ;
pos_clk__survivor_paths[3][14]__delay:		DELAY (POSEDGE) clk survivor_paths[3][14] ;
pos_clk__survivor_paths[3][15]__delay:		DELAY (POSEDGE) clk survivor_paths[3][15] ;
pos_clk__survivor_paths[3][16]__delay:		DELAY (POSEDGE) clk survivor_paths[3][16] ;
pos_clk__survivor_paths[3][17]__delay:		DELAY (POSEDGE) clk survivor_paths[3][17] ;
pos_clk__survivor_paths[3][18]__delay:		DELAY (POSEDGE) clk survivor_paths[3][18] ;
pos_clk__survivor_paths[3][19]__delay:		DELAY (POSEDGE) clk survivor_paths[3][19] ;
pos_clk__survivor_paths[3][20]__delay:		DELAY (POSEDGE) clk survivor_paths[3][20] ;
pos_clk__survivor_paths[3][21]__delay:		DELAY (POSEDGE) clk survivor_paths[3][21] ;
pos_clk__survivor_paths[3][22]__delay:		DELAY (POSEDGE) clk survivor_paths[3][22] ;
pos_clk__survivor_paths[3][23]__delay:		DELAY (POSEDGE) clk survivor_paths[3][23] ;
pos_clk__survivor_paths[3][24]__delay:		DELAY (POSEDGE) clk survivor_paths[3][24] ;
pos_clk__survivor_paths[3][25]__delay:		DELAY (POSEDGE) clk survivor_paths[3][25] ;
pos_clk__survivor_paths[3][26]__delay:		DELAY (POSEDGE) clk survivor_paths[3][26] ;
pos_clk__survivor_paths[3][27]__delay:		DELAY (POSEDGE) clk survivor_paths[3][27] ;
pos_clk__survivor_paths[3][28]__delay:		DELAY (POSEDGE) clk survivor_paths[3][28] ;
pos_clk__survivor_paths[3][29]__delay:		DELAY (POSEDGE) clk survivor_paths[3][29] ;
pos_clk__survivor_paths[3][30]__delay:		DELAY (POSEDGE) clk survivor_paths[3][30] ;
pos_clk__survivor_paths[3][31]__delay:		DELAY (POSEDGE) clk survivor_paths[3][31] ;
pos_clk__t[0]__delay:		DELAY (POSEDGE) clk t[0] ;
pos_clk__t[1]__delay:		DELAY (POSEDGE) clk t[1] ;
pos_clk__t[2]__delay:		DELAY (POSEDGE) clk t[2] ;
pos_clk__t[3]__delay:		DELAY (POSEDGE) clk t[3] ;
pos_clk__t[4]__delay:		DELAY (POSEDGE) clk t[4] ;

ENDMODEL
