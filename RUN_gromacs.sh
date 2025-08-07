#!/bin/bash

gmx grompp -f em.mdp -c protein_solvate.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em 

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -pin on -nt 8 

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -pin on -nt 8

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_10ns.tpr
gmx mdrun -deffnm md_10ns -pin on -nt 8

for((s=20;s<=20000;s=s+10))
do
   rstpt=$((s - 10))
   gmx convert-tpr -s md_${rstpt}ns.tpr -o md_${s}ns.tpr -extend 10000
   gmx mdrun -cpi md_${rstpt}ns.cpt -deffnm md_${s}ns -noappend -pin on -nt 8 
done
