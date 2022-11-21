// impute.nf

process GETSAMPLEID {
  input:
    path input_vcfs

  output:
    path 'all_sampleIDs'

  script:
  """
  #!/bin/bash
  
  bcftools query -l $input_vcfs > all_sampleIDs
  """
}

process SPLIT {
  input:
    path input_vcfs
     
  
  output:
    path outDir+'chr$ch_$split_part.vcf.gz'


  script:
  """
  #!/bin/bash

  bcftools view -S $sample_name -Oz -o {output} $input 
  """
}
