#!/bin/bash
./pat_gen SamplePat.txt
./send.exe 1 SamplePat.txt_cfeb0_tmbc.pat  >& null
./send.exe 2 SamplePat.txt_cfeb1_tmbc.pat  >& null
./send.exe 3 SamplePat.txt_cfeb2_tmbc.pat  >& null
./send.exe 4 SamplePat.txt_cfeb3_tmbc.pat  >& null
./send.exe 5 SamplePat.txt_cfeb4_tmbc.pat  >& null
./send.exe 6 SamplePat.txt_cfeb5_tmbc.pat  >& null
./send.exe 7 SamplePat.txt_cfeb6_tmbc.pat  >& null
./send_d.exe
