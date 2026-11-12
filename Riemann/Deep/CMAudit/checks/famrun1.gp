read("Riemann/Deep/CMAudit/checks/famcm.gp");
OUT = "Riemann/Deep/CMAudit/checks/out_fam.txt";
famrun(1000, [32, 64, 128, 256, 512], 5/6, OUT);
famrun(2000, [32, 64, 128, 256, 512, 1024], 5/6, OUT);
famrun(4000, [64, 128, 256, 512, 1024, 2048], 5/6, OUT);
write(OUT, "DONE1");
