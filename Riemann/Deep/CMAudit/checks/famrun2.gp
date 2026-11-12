read("Riemann/Deep/CMAudit/checks/famcm.gp");
OUT = "Riemann/Deep/CMAudit/checks/out_fam.txt";
famrun(8000, [64, 128, 256, 512, 1024, 2048, 4096], 5/6, OUT);
famrun(16000, [128, 256, 512, 1024, 2048, 4096, 8192], 5/6, OUT);
write(OUT, "DONE2");
