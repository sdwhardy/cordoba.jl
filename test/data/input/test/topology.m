%TNEP optimizor input file

function mpc = owpp_tnep_map
mpc.version = '2';
mpc.baseMVA = 100;

%% bus data
	%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
	mpc.bus = [
1 3 0 0 0 0 1 1 0 220 1 1.1 0.9;
2 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
3 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
4 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
5 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
6 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
7 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
8 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
9 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
10 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
11 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
12 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
13 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
14 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
15 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
16 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
17 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
18 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
19 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
20 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
21 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
22 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
23 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
24 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
25 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
26 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
27 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
28 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
29 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
30 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
31 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
32 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
33 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
34 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
35 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
36 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
37 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
38 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
39 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
40 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
41 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
42 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
43 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
44 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
45 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
46 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
47 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
48 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
49 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
50 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
51 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
52 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
53 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
54 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
55 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
56 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
57 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
58 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
59 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
60 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
61 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
62 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
63 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
];
%% generator data
	%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
	mpc.gen = [
9 0 0 0 0 1 100 1 60 0 0 0 0 0 0 0 0 0 0 0 1;
10 0 0 0 0 1 100 1 40 0 0 0 0 0 0 0 0 0 0 0 1;
11 0 0 0 0 1 100 1 160 0 0 0 0 0 0 0 0 0 0 0 1;
12 0 0 0 0 1 100 1 160 0 0 0 0 0 0 0 0 0 0 0 1;
13 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
14 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
15 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
16 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
17 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
18 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
19 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
20 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
21 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
22 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
23 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
24 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
25 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
26 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
27 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
28 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
29 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
30 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
31 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
32 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
33 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
34 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
35 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
36 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
37 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
38 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
39 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
40 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
41 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
42 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
43 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
44 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
45 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
46 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
47 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
48 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
49 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
50 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
51 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
52 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
53 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
54 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
55 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
56 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
57 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
58 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
59 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
60 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
61 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
62 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
63 0 0 0 0 1 100 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
];
%%Additional generator fields, type=0(wind generator), type=1(onshore market ie load and generator created)
	%column_names% 		type invest
	mpc.gen_type = [
0 210
0 210
0 210
0 210
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
];

	%% branch data
	%	fbus	tbus	r	       x	    b	  rateA	rateB	rateC	ratio	angle	status	angmin	angmax
	mpc.branch = [
1	2	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
1	3	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
1	4	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
1	5	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
1	6	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
1	7	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
1	8	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
1	9	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
1	10	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
1	11	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
1	12	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
2	3	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
2	4	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
2	5	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
2	6	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
2	7	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
2	8	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
2	9	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
2	10	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
2	11	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
2	12	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
3	4	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
3	5	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
3	6	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
3	7	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
3	8	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
3	9	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
3	10	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
3	11	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
3	12	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
4	5	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
4	6	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
4	7	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
4	8	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
4	9	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
4	10	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
4	11	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
4	12	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
5	6	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
5	7	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
5	8	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
5	9	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
5	10	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
5	11	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
5	12	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
6	7	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
6	8	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
6	9	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
6	10	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
6	11	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
6	12	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
7	8	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
7	9	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
7	10	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
7	11	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
7	12	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
8	9	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
8	10	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
8	11	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
8	12	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
9	10	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
9	11	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
9	12	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
10	11	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
10	12	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
11	12	0.010   0.100   0.00   0  0  0  0  0  1 -60  60;
];
%%Additional generator fields, type=0(wind generator), type=1(onshore market ie load and generator created)
	%column_names% 		cost
	mpc.branch_cost = [
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
1.14;
];

	%candidate branch data
	%column_names%	f_bus	t_bus	br_r	br_x	br_b	rate_a	rate_b	rate_c	tap	shift	br_status	angmin	angmax	mva	construction_cost
	mpc.ne_branch = [
1	2	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
1	3	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
1	4	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
1	5	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
1	6	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
1	7	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
1	8	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
1	9	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
1	10	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
1	11	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
1	12	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
2	3	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
2	4	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
2	5	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
2	6	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
2	7	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
2	8	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
2	9	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
2	10	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
2	11	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
2	12	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
3	4	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
3	5	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
3	6	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
3	7	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
3	8	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
3	9	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
3	10	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
3	11	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
3	12	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
4	5	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
4	6	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
4	7	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
4	8	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
4	9	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
4	10	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
4	11	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
4	12	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
5	6	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
5	7	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
5	8	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
5	9	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
5	10	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
5	11	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
5	12	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
6	7	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
6	8	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
6	9	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
6	10	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
6	11	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
6	12	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
7	8	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
7	9	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
7	10	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
7	11	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
7	12	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
8	9	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
8	10	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
8	11	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
8	12	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
9	10	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
9	11	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
9	12	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
10	11	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
10	12	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
11	12	0.004 0.04 0 1000 1000 1000 0 0 0 -60 60 100 10;
];
%%-----  OPF Data  -----%%
	%% generator cost data
	%	1	startup	shutdown	n	x1	y1	...	xn	yn
	%	2	startup	shutdown	n	c(n-1)	...	c0
	mpc.gencost = [
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
];
%% existing dc bus data
	%column_names%   busdc_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin Cdc
	mpc.busdc = [
1 1 0 1 300 1.1 0.9 0;
2 2 0 1 300 1.1 0.9 0;
3 3 0 1 300 1.1 0.9 0;
4 4 0 1 300 1.1 0.9 0;
5 5 0 1 300 1.1 0.9 0;
6 6 0 1 300 1.1 0.9 0;
7 7 0 1 300 1.1 0.9 0;
8 8 0 1 300 1.1 0.9 0;
9 9 0 1 300 1.1 0.9 0;
10 10 0 1 300 1.1 0.9 0;
11 11 0 1 300 1.1 0.9 0;
12 12 0 1 300 1.1 0.9 0;
13 13 0 1 300 1.1 0.9 0;
14 14 0 1 300 1.1 0.9 0;
15 15 0 1 300 1.1 0.9 0;
16 16 0 1 300 1.1 0.9 0;
17 17 0 1 300 1.1 0.9 0;
18 18 0 1 300 1.1 0.9 0;
19 19 0 1 300 1.1 0.9 0;
20 20 0 1 300 1.1 0.9 0;
21 21 0 1 300 1.1 0.9 0;
22 22 0 1 300 1.1 0.9 0;
23 23 0 1 300 1.1 0.9 0;
24 24 0 1 300 1.1 0.9 0;
25 25 0 1 300 1.1 0.9 0;
26 26 0 1 300 1.1 0.9 0;
27 27 0 1 300 1.1 0.9 0;
28 28 0 1 300 1.1 0.9 0;
29 29 0 1 300 1.1 0.9 0;
30 30 0 1 300 1.1 0.9 0;
31 31 0 1 300 1.1 0.9 0;
32 32 0 1 300 1.1 0.9 0;
33 33 0 1 300 1.1 0.9 0;
34 34 0 1 300 1.1 0.9 0;
35 35 0 1 300 1.1 0.9 0;
36 36 0 1 300 1.1 0.9 0;
37 37 0 1 300 1.1 0.9 0;
38 38 0 1 300 1.1 0.9 0;
39 39 0 1 300 1.1 0.9 0;
40 40 0 1 300 1.1 0.9 0;
41 41 0 1 300 1.1 0.9 0;
42 42 0 1 300 1.1 0.9 0;
43 43 0 1 300 1.1 0.9 0;
44 44 0 1 300 1.1 0.9 0;
45 45 0 1 300 1.1 0.9 0;
46 46 0 1 300 1.1 0.9 0;
47 47 0 1 300 1.1 0.9 0;
48 48 0 1 300 1.1 0.9 0;
49 49 0 1 300 1.1 0.9 0;
50 50 0 1 300 1.1 0.9 0;
51 51 0 1 300 1.1 0.9 0;
52 52 0 1 300 1.1 0.9 0;
53 53 0 1 300 1.1 0.9 0;
54 54 0 1 300 1.1 0.9 0;
55 55 0 1 300 1.1 0.9 0;
56 56 0 1 300 1.1 0.9 0;
57 57 0 1 300 1.1 0.9 0;
58 58 0 1 300 1.1 0.9 0;
59 59 0 1 300 1.1 0.9 0;
60 60 0 1 300 1.1 0.9 0;
61 61 0 1 300 1.1 0.9 0;
62 62 0 1 300 1.1 0.9 0;
63 63 0 1 300 1.1 0.9 0;
];
%% All Candidate equipment has equivalent existing infrastructure field as well ie  mpc.branchdc_ne <->  mpc.branchdc similar to AC
	%% Candidate DC buses here - refer to MATACDC 1.0 User’s Manual for description of fields
	%% candidate dc bus data
	%column_names%   busdc_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin  Cdc
	mpc.busdc_ne = [
	10000              4       0       1       300         1.1     0.9     0;
	];
%% existing dc branches
	%column_names%   fbusdc  tbusdc  r      l        c   rateA   rateB rateC cost status
	mpc.branchdc = [
1	2	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
1	3	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
1	4	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
1	5	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
1	6	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
1	7	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
1	8	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
1	9	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
1	10	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
1	11	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
1	12	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
2	3	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
2	4	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
2	5	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
2	6	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
2	7	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
2	8	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
2	9	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
2	10	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
2	11	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
2	12	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
3	4	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
3	5	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
3	6	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
3	7	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
3	8	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
3	9	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
3	10	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
3	11	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
3	12	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
4	5	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
4	6	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
4	7	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
4	8	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
4	9	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
4	10	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
4	11	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
4	12	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
5	6	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
5	7	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
5	8	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
5	9	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
5	10	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
5	11	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
5	12	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
6	7	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
6	8	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
6	9	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
6	10	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
6	11	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
6	12	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
7	8	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
7	9	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
7	10	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
7	11	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
7	12	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
8	9	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
8	10	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
8	11	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
8	12	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
9	10	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
9	11	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
9	12	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
10	11	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
10	12	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
11	12	0.0	0.00	 0.00  0	 		0	 0 0.125	 1;
];
%% Candidate Branches here - refer to MATACDC 1.0 Users Manual for description of fields
	 %% candidate branches
	 %column_names%   fbusdc  tbusdc  r      l        c   rateA   rateB   rateC status cost
	 mpc.branchdc_ne = [
1 2 0 0 0 10000 10000 10000 1 0;
1 3 0 0 0 10000 10000 10000 1 0;
1 4 0 0 0 10000 10000 10000 1 0;
1 5 0 0 0 10000 10000 10000 1 0;
1 6 0 0 0 10000 10000 10000 1 0;
1 7 0 0 0 10000 10000 10000 1 0;
1 8 0 0 0 10000 10000 10000 1 0;
1 9 0 0 0 10000 10000 10000 1 0;
1 10 0 0 0 10000 10000 10000 1 0;
1 11 0 0 0 10000 10000 10000 1 0;
1 12 0 0 0 10000 10000 10000 1 0;
2 3 0 0 0 10000 10000 10000 1 0;
2 4 0 0 0 10000 10000 10000 1 0;
2 5 0 0 0 10000 10000 10000 1 0;
2 6 0 0 0 10000 10000 10000 1 0;
2 7 0 0 0 10000 10000 10000 1 0;
2 8 0 0 0 10000 10000 10000 1 0;
2 9 0 0 0 10000 10000 10000 1 0;
2 10 0 0 0 10000 10000 10000 1 0;
2 11 0 0 0 10000 10000 10000 1 0;
2 12 0 0 0 10000 10000 10000 1 0;
3 4 0 0 0 10000 10000 10000 1 0;
3 5 0 0 0 10000 10000 10000 1 0;
3 6 0 0 0 10000 10000 10000 1 0;
3 7 0 0 0 10000 10000 10000 1 0;
3 8 0 0 0 10000 10000 10000 1 0;
3 9 0 0 0 10000 10000 10000 1 0;
3 10 0 0 0 10000 10000 10000 1 0;
3 11 0 0 0 10000 10000 10000 1 0;
3 12 0 0 0 10000 10000 10000 1 0;
4 5 0 0 0 10000 10000 10000 1 0;
4 6 0 0 0 10000 10000 10000 1 0;
4 7 0 0 0 10000 10000 10000 1 0;
4 8 0 0 0 10000 10000 10000 1 0;
4 9 0 0 0 10000 10000 10000 1 0;
4 10 0 0 0 10000 10000 10000 1 0;
4 11 0 0 0 10000 10000 10000 1 0;
4 12 0 0 0 10000 10000 10000 1 0;
5 6 0 0 0 10000 10000 10000 1 0;
5 7 0 0 0 10000 10000 10000 1 0;
5 8 0 0 0 10000 10000 10000 1 0;
5 9 0 0 0 10000 10000 10000 1 0;
5 10 0 0 0 10000 10000 10000 1 0;
5 11 0 0 0 10000 10000 10000 1 0;
5 12 0 0 0 10000 10000 10000 1 0;
6 7 0 0 0 10000 10000 10000 1 0;
6 8 0 0 0 10000 10000 10000 1 0;
6 9 0 0 0 10000 10000 10000 1 0;
6 10 0 0 0 10000 10000 10000 1 0;
6 11 0 0 0 10000 10000 10000 1 0;
6 12 0 0 0 10000 10000 10000 1 0;
7 8 0 0 0 10000 10000 10000 1 0;
7 9 0 0 0 10000 10000 10000 1 0;
7 10 0 0 0 10000 10000 10000 1 0;
7 11 0 0 0 10000 10000 10000 1 0;
7 12 0 0 0 10000 10000 10000 1 0;
8 9 0 0 0 10000 10000 10000 1 0;
8 10 0 0 0 10000 10000 10000 1 0;
8 11 0 0 0 10000 10000 10000 1 0;
8 12 0 0 0 10000 10000 10000 1 0;
9 10 0 0 0 10000 10000 10000 1 0;
9 11 0 0 0 10000 10000 10000 1 0;
9 12 0 0 0 10000 10000 10000 1 0;
10 11 0 0 0 10000 10000 10000 1 0;
10 12 0 0 0 10000 10000 10000 1 0;
11 12 0 0 0 10000 10000 10000 1 0;
];

	%% existing converters
	%column_names%   busdc_i busac_i type_dc type_ac P_g   Q_g  islcc  	Vtar 		rtf  xtf  transformer tm   bf  filter    rc      xc  reactor   basekVac Vmmax   Vmmin   Imax     status   LossA  LossB  LossCrec LossCinv   droop   Pdcset     Vdcset  dVdcset Pacmax Pacmin Qacmax   Qacmin cost
	mpc.convdc = [
1 1 2 3 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
2 2 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
3 3 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
4 4 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
5 5 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
6 6 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
7 7 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
8 8 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
9 9 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 51.5;
10 10 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 51.5;
11 11 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 51.5;
12 12 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 51.5;
13 13 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
14 14 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
15 15 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
16 16 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
17 17 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
18 18 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
19 19 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
20 20 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
21 21 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
22 22 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
23 23 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
24 24 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
25 25 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
26 26 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
27 27 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
28 28 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
29 29 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
30 30 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
31 31 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
32 32 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
33 33 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
34 34 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
35 35 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
36 36 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
37 37 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
38 38 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
39 39 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
40 40 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
41 41 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
42 42 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
43 43 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
44 44 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
45 45 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
46 46 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
47 47 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
48 48 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
49 49 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
50 50 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
51 51 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
52 52 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
53 53 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
54 54 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
55 55 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
56 56 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
57 57 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
58 58 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
59 59 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
60 60 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
61 61 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
62 62 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
63 63 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 40000 -40000 40000 -40000 19.25;
];

	%trans, filter, reactor, vmin vmax same as conv
	%% candidate converters
	%column_names%   busdc_i busac_i type_dc type_ac P_g   Q_g  islcc  Vtar    rtf   xtf  transformer tm   bf 	filter    rc     xc  reactor   basekVac Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv  droop     Pdcset    Vdcset  dVdcset Pacmax Pacmin Qacmax Qacmin cost
	mpc.convdc_ne = [
	                1000       2      1       1    400000    0   	0     1.0   0.001  0.1       0 			 1 	0.08 		0 		0.001   0.09 		0  				220    1.1     0.9    100000      1     	 0     0        0       0      0.0050    -52.7     1.0079     0  		1  -1   1    -1   1000000;
	];
 % hours
	 mpc.time_elapsed = 1.0

	 %% storage data
	 % storage_bus   ps   qs energy energy_rating charge_rating discharge_rating charge_efficiency discharge_efficiency thermal_rating    qmin   qmax    r    x p_loss q_loss status
	 mpc.storage = [
1 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
2 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
3 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
4 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
5 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
6 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
7 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
8 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
9 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
10 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
11 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
12 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
13 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
14 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
15 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
16 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
17 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
18 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
19 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
20 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
21 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
22 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
23 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
24 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
25 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
26 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
27 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
28 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
29 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
30 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
31 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
32 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
33 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
34 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
35 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
36 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
37 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
38 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
39 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
40 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
41 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
42 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
43 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
44 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
45 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
46 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
47 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
48 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
49 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
50 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
51 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
52 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
53 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
54 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
55 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
56 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
57 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
58 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
59 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
60 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
61 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
62 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
63 0 0 0 0 0 0 0.9 0.9 0 0 0 0 0 0 0 1;
];
%% cost 36.6
	 %% storage additional data
	 %column_names% max_energy_absorption stationary_energy_inflow stationary_energy_outflow self_discharge_rate cost
	 mpc.storage_extra = [
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
0 0 0 0.0001 39;
];
%% storage data
	 %column_names%   storage_bus ps 	qs 	energy  energy_rating charge_rating  discharge_rating  charge_efficiency  discharge_efficiency  thermal_rating  qmin  	qmax  	r  		x  p_loss  	q_loss  status eq_cost inst_cost co2_cost 	max_energy_absorption 	stationary_energy_inflow 	stationary_energy_outflow 	self_discharge_rate	 	cost_abs 	cost_inj on_off
	 mpc.ne_storage = [
	 1	 0.0 0.0 0.0 6400.0	1280.0 1600.0 0.9 0.9 16000	-3200.0	4480.0 0.1 0.0 0.0	0.0	 0 	100000 	100000	100000 	32000000 0 	0 1e-4 0 0 1;
	 											 											];
