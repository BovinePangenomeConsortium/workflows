import polars as pl
import sys

def pop_next_most_unique_sample(SNPs):
    singletons = (SNPs.with_columns(rowsum=pl.sum_horizontal(pl.all()))
                      .filter(pl.col('rowsum')==1)
                      .sum()
                      .drop('rowsum')
                      .unpivot(variable_name='sample',value_name='Unique SNPs')
                      .sample(shuffle=True,fraction=1) #to randomly handle ties in unique SNPs
                 )
    max_unique_SNPs = singletons.select(pl.col("Unique SNPs").arg_max()).item(0,0)
    max_unique_sample = singletons.select(pl.col("sample").get(max_unique_SNPs)).item(0,0)
    return {"Sample":max_unique_sample,"Unique SNPs":singletons.select(pl.col('Unique SNPs').get(max_unique_SNPs)).item(0,0)}, SNPs.drop(max_unique_sample)

def prioritise_samples(SNPs,N_samples):

    popped_samples = []
    while len(popped_samples) < N_samples:
        entry, SNPs = pop_next_most_unique_sample(SNPs)
        popped_samples.append(entry)

    return pl.DataFrame(popped_samples)

#hack to get number of samples, so we can override the schema
header_only = pl.read_csv(sys.argv[1],separator=' ',n_rows=1)

#assume all multiallelics have been split to multiple biallelics and are true/false input
SNPs = pl.read_csv(sys.argv[1],separator=' ',schema_overrides=[pl.Boolean]*len(header_only.columns))

N = len(SNPs.columns) if len(sys.argv) <= 3 else int(sys.argv[3])

prioritised_samples = prioritise_samples(SNPs,N)
prioritised_samples.write_csv(sys.argv[2])
