#conda create -n flagcsnap python aria2 biopython cctyper dendropy diamond ete3 famsa fasttree infernal ncbi-datasets-cli padloc pandas=1.4.3 requests taxonkit numpy pyarrow==12.0.1 matplotlib defense-finder rich rich-click -c conda-forge -c bioconda -c padlocbio -c russel88
#micromamba activate flagcsnap

padloc --db-update
macsydata install -U -u --org mdmparis defense-finder-models
git clone https://github.com/pentamorfico/CRISPRCasTyper.git
cd CRISPRCasTyper
git checkout develop
pip install .
cd ..
yes | rm -r CRISPRCasTyper
mkdir -p $CONDA_PREFIX/db/flagcsnap
cp dbs/all.cm $CONDA_PREFIX/db/flagcsnap
cp dbs/type_strains_curated.txt $CONDA_PREFIX/db/flagcsnap
python -m pip install -e .
cd $CONDA_PREFIX/db/
mkdir -p taxonkit cm_models
cp all.cm cm_models/all.cm
cd taxonkit
curl -o "taxdump.tar.gz" "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz" 
tar -xvzf taxdump.tar.gz
