import sim_data
import osem
import plot_sim

for pf in ["tacs_AD.h5", "tacs_C.h5"]:
	sim_data.main(["--phantom_file",pf])
	
for pf in ["tacs_AD.npz", "tacs_C.npz"]:
	for sc in [0.1, 1.0, 10.0]:
		for frm in list(range(25)):
			osem.main(["--phantom_file",pf,"--num_seeds", str(20),"--scanner_sens", str(sc), "--frm",str(frm)])
			
flagfig = False

if flagfig:
	plot_sim.main()
