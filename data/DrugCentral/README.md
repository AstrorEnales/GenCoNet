Data extracted from [DrugCentral](http://drugcentral.org) published under [CC BY-SA](http://drugcentral.org/privacy).
Ursu, O., Holmes, J., Knockel, J., Bologa, C. G., Yang, J. J., Mathias, S. L., … Oprea, T. I. (2016). DrugCentral: online drug compendium. Nucleic Acids Research, 45(D1), D932–D939. https://doi.org/10.1093/nar/gkw993

Queries for data extraction:

    SELECT DISTINCT b.identifier AS drugbank_id, a.concept_name, a.snomed_conceptid, a.umls_cui
        FROM drugcentral.omop_relationship AS a
        LEFT JOIN identifier AS b ON a.struct_id=b.struct_id
        WHERE a.relationship_name = "indication" AND b.id_type = "DRUGBANK_ID" AND (a.umls_cui IS NOT NULL OR a.snomed_conceptid IS NOT NULL);
    
    SELECT DISTINCT b.identifier AS drugbank_id, a.concept_name, a.snomed_conceptid, a.umls_cui
        FROM drugcentral.omop_relationship AS a
        LEFT JOIN identifier AS b ON a.struct_id=b.struct_id
        WHERE a.relationship_name = "contraindication" AND b.id_type = "DRUGBANK_ID" AND (a.umls_cui IS NOT NULL OR a.snomed_conceptid IS NOT NULL);

    SELECT DISTINCT a.id as drugcentral_id, b.identifier as drugbank_id, c.identifier as rxnorm_id, a.name AS drug_name
    	FROM structures as a
        LEFT JOIN identifier AS b ON a.id=b.struct_id
        LEFT JOIN identifier AS c ON a.id=c.struct_id
        WHERE b.id_type = "DRUGBANK_ID" AND c.id_type = "RXNORM" ORDER BY b.identifier;