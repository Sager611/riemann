\\ odlaudit: class numbers of quadratic fields, |D| <= 200000 (fundamental D); quadclassunit (GRH-conditional for D > 0 in PARI; test data only)
default(parisizemax, 120000000);
f = fileopen("/Users/adriansager/playground/riemann/.scratch/odlaudit/quadh.csv", "w");
for(D = -200000, 200000, if(D != 0 && D != 1 && isfundamental(D), filewrite(f, Str(D, ",", quadclassunit(D).no))));
fileclose(f);
quit;
