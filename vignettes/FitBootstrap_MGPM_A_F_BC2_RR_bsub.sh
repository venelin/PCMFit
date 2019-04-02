for id in `seq 1 1 10`
do
mkdir -p ResultsBootstrap/MGPM_A_F_BC2_RR_BSID_$id
cd ResultsBootstrap/MGPM_A_F_BC2_RR_BSID_$id
if [ -f "Result_MGPM_A_F_BC2_RR_BSID_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
else
bsub -M 10000 -n 8 -W 3:59 -R ib sh R --vanilla --slave -f ../../FitBootstrap_MGPM_A_F_BC2_RR.R --args $id
fi
cd ../..
done
