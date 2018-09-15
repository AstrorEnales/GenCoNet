# GenCoNet
Under development

## Pipeline

  1. Prerequisites
     - Most data downloads are automatic. However some databases like DrugBank need to be downloaded manually due to licensing or the need to register an account
       - Download DrugBank full release to `data/DrugBank/drugbank_all_full_database.xml.zip` and extract `full database.xml`
       - Execute the [DrugCentral queries](data/DrugCentral/README.md) and save them to the respective files
       - Download OMIM `genemap2.txt` to `data/OMIM/genemap2.txt`
     - Configure `data/config.json` to your Neo4j installation (bin path including admin tools)
  2. Execution
     1. Pre-processing
        - Run `mondo.py`
        - Run `drugbank.py`
        - Run `drugcentral.py`
        - Run `disgenet.py`
        - Run `gwas_catalog.py`
        - Run `hgnc.py`
        - Run `hpo.py` (WIP)
        - Run `med_rt.py` (WIP)
        - Run `omim.py`
        - Run `pubmed.py`
     2. Fusion
        - Run `fusion.py`
     3. Import database into Neo4j or directly use `[output-path]/graph.json`
