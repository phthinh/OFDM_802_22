`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    14:02:46 12/04/2012 
// Design Name: 
// Module Name:    Pilots_Insert 
// Project Name: 
// Target Devices: 
// Tool versions: 
// Description: 
//
// Dependencies: 
//
// Revision: 
// Revision 0.01 - File Created
// Additional Comments: 
//
//////////////////////////////////////////////////////////////////////////////////
module Pilots_Insert(
 	input 			CLK_I, RST_I,
	input [31:0] 	DAT_I,
	input 			CYC_I, WE_I, STB_I, 
	output			ACK_O,
	
	output reg [31:0]	DAT_O,
	output reg		CYC_O, STB_O,
	output			WE_O,
	input				ACK_I	
    );
parameter P_P = 16'h7fff;	// +1 in Q1.15
parameter P_N = 16'h8001;	// -1 in Q1.15 
reg [1:0] alloc_vec 	 [0:8399];   // signed bit of real part of pilots,
initial $readmemh("./MY_SOURCES/Al_vec.txt", alloc_vec);
	 
reg 	[31:0]	idat;
wire [31:0]  odat;
reg			ival;	
wire 			out_halt, ena;
wire			datout_ack;


reg [10:0]	dat_cnt;
reg [12:0]	alloc_ptr;			// pointer of allocation vector
reg [6:0]	pilot_cnt;
reg			nul_insert_ena;		//inserting null symbol for guarding.
wire			pil_insert_ena;
wire			car_unactive;
wire [15:0] pil_Re;
wire [1:0]	cur_carrier;
wire 			sym_end;


assign 	out_halt   = (STB_O)&(CYC_O) & (~ACK_I);
assign 	datout_ack = STB_O & ACK_I;
assign 	ena 		= CYC_I & STB_I & WE_I;
assign 	ACK_O 	= ena & CYC_O & (~out_halt) & (~pil_insert_ena) & (~nul_insert_ena) & (~car_unactive);

always @(posedge CLK_I) begin
	if(RST_I) 			idat<= 2'b00;
	else if(ACK_O) 	idat <= DAT_I;
end
always @(posedge CLK_I) begin
	if(RST_I) 		ival <= 1'b0;
	else if(ena)	ival <= 1'b1;
	else				ival <= 1'b0;
end

always @(posedge CLK_I)
begin
	if(RST_I)	STB_O <= 1'b0;
	else if(ival|pil_insert_ena|nul_insert_ena|car_unactive)	STB_O <= 1'b1;
	else if(~ival) 														STB_O <= 1'b0;
end

reg [1:0] icyc;
always @(posedge CLK_I)
begin
	if(RST_I)		icyc <= 2'b00;		
	else				icyc <= {icyc[0],CYC_I};	
end

always @(posedge CLK_I)
begin
	if(RST_I)										CYC_O	<= 1'b0;			
	else if (icyc[1] & CYC_I & (~CYC_O))	CYC_O	<= 1'b1;
	else if (sym_end & (~CYC_I))				CYC_O	<= 1'b0;
end
assign odat = (nul_insert_ena|car_unactive)? 32'd0: (pil_insert_ena)? {16'd0, pil_Re} : DAT_I;
always @(posedge CLK_I)
begin
	if(RST_I)							DAT_O <= 32'b0;
	else if(ival & (~out_halt))	DAT_O <= odat;	
end
assign WE_O  = STB_O;	 

always@(posedge CLK_I)
begin
	if(RST_I)										dat_cnt	<= 11'd0;		
	else if(CYC_I & (~icyc[0]))				dat_cnt	<= 11'd0;
	else if(datout_ack)						   dat_cnt	<= dat_cnt + 1'b1;
end
assign sym_end = (dat_cnt == 11'd2047);

always@(posedge CLK_I)
begin
	if(RST_I)										pilot_cnt	<= 7'd0;			
	else if(CYC_I & (~icyc[0]))				pilot_cnt	<= 7'd0;
	else if(pil_insert_ena & (~out_halt))	pilot_cnt	<= pilot_cnt + 1'b1;
end

reg [14:0]	pil_seed;
wire 			cur_pil;		//current pilot


always@(posedge CLK_I)
begin
	if(RST_I)								pil_seed	<= 15'b011011100010101;			
	else if(CYC_I & (~icyc[0]))		pil_seed	<= 15'b011011100010101;
	else if(pil_insert_ena & CYC_O)	pil_seed	<= {pil_seed[13:0],cur_pil};
end
assign cur_pil = pil_seed[14]^pil_seed[13];
assign pil_Re = (cur_pil)? P_P : P_N;

always@(posedge CLK_I)
begin
	if(RST_I)											alloc_ptr	<= 13'd0;		
	else if(CYC_I & (~icyc[0]))					alloc_ptr	<= 13'd0;
	else if(datout_ack &(~nul_insert_ena))		alloc_ptr	<= alloc_ptr + 1'b1;
end

assign pil_insert_ena = (alloc_vec[alloc_ptr] == 2'b01)&(~nul_insert_ena);
assign car_unactive   = (alloc_vec[alloc_ptr] == 2'b00);


always@(posedge CLK_I)
begin
	if(RST_I)												nul_insert_ena  = 1'b0;		
	else if(CYC_I & (~icyc[0]))						nul_insert_ena  = 1'b0;	
	else if(icyc[0] & (~icyc[1]))						nul_insert_ena  = 1'b1;
	else if(icyc[1] & (~CYC_O))						nul_insert_ena  = 1'b0;
	else if(datout_ack & (dat_cnt == 11'd0))		nul_insert_ena  = 1'b0;
	else if(datout_ack & (dat_cnt == 11'd839))	nul_insert_ena  = 1'b1;	
	else if(datout_ack & (dat_cnt == 11'd1206))	nul_insert_ena  = 1'b0;	
	else if(icyc[0] & datout_ack & sym_end) 		nul_insert_ena  = 1'b1;	
end


endmodule 
