for i in {2..4}
do 
cd /data/guang/SAM_evaluation_Paper/Trimming/Trim_Galore/Fake50bpWithAdapters_Stringency$i
bash TrimGaloreTrimmed_FakeSE50bpReads_Include0-25bpAdapter_Align_and_Evaluate.sh > TrimGaloreTrimmed_FakeSE50bpReads_Include0-25bpAdapter_Align_and_Evaluate\_$(date +"%m-%d-%Y").log 2>&1
done

