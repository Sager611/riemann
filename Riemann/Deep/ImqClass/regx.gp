v=readvec("runs/reg2000.extra"); bad=0; for(i=1,#v, if(qfbclassno(-v[i])<=2000, bad++)); print("extra entries ", #v, ", with h<=2000: ", bad);
