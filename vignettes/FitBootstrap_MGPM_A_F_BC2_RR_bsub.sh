for id in `seq 1 1 10`
do
mkdir -p ResultsBootstrap/MGPM_A_F_BC2_RR_BSID_$id
cd ResultsBootstrap/MGPM_A_F_BC2_RR_BSID_$id
if [ -f "Result_MGPM_A_F_BC2_RR_BSID_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
else
# The following would execute only if the bsub command is available on your system
# bsub -M 10000 -n 8 -W 3:59 -R ib R --vanilla --slave -f ../../FitBootstrap_MGPM_A_F_BC2_RR.R --args $id
# Should run on unix systems. Comment or remove this line if uncommenting the above.
R --vanilla --slave -f ../../FitBootstrap_MGPM_A_F_BC2_RR.R --args $id
fi
cd ../..
done
