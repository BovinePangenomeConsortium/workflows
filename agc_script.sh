mkdir -p agc

while read -r a b
do
  ln -s $a $b
done < <(awk -F ',' 'NR>1 {print "/cluster/work/pausch_bpc/bovine_pangenome/data/currated_assemblies/"$1".fa.gz","agc/"$5"_"$6".fa.gz"}' BPC_metadata.csv)

#agc create ${REF} -i <(ls *.fa.gz | grep -v "$REF" ) -t $THREADS -o BPC.agc
