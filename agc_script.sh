mkdir agc
while read -r a b
do
  ln -s $a $b
done <(awk -F ',' 'NR>1 {print "data/currated_assemblies/"$1".fa.gz",agc/$5"_"$6".fa.gz"}' BPC_assemblies_20251020.csv)


REF=agc/Wagyu_haplotype1_v1.polished.fa.gz

agc create ${REF} -i <(ls agc *.fa.gz | grep -v "$REF" ) -t $THREADS -o BPC.agc
