cp dbs/cm_models/all.cm $CONDA_PREFIX/db/fastgcsnap
cp dbs/type_strains_curated.txt $CONDA_PREFIX/db/
padloc --db-update
macsydata install -U -u --org mdmparis defense-finder-models
git clone https://github.com/pentamorfico/CRISPRCasTyper.git
cd CRISPRCasTyper
git checkout develop
pip install .
cd $CONDA_PREFIX/db/
mkdir -p taxonkit cm_models
cp all.cm cm_models/all.cm
cd taxonkit
curl -o "taxdump.tar.gz" "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz" 
tar -xvzf taxdump.tar.gz
