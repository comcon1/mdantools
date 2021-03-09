all: gdens.tsv
	
gdens.tsv: topol.tpr traj.xtc index.ndx
	../gdens_mda.py -s topol.tpr -i traj.xtc -b 5000:3000:45000 \
		-d 0.1:550:-10 -x index.ndx:POPC:MEMB:NH3 -r POPC:C2 -o gdens.tsv
